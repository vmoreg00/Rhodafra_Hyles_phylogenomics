#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=54
#SBATCH --time=168:00:00
#SBATCH --mem=50G
#SBATCH --job-name=hDNA_mapping_optim
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_optim-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_optim-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
REFGENOME=$WS/data/refgenome
DATA=$WS/data/fastqtrim
RESULTS=$WS/results/00_hDNA_mapper_opt
strategy=()
timing=()

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# create directories
if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS;
fi;

# Get reference mitogenome if it does not exist
## As it is the last chromosome, it can be extracted with grep
## -A every option ensures that every line after the hit will be printed
if [ ! -e $REFGENOME/DeiElpe2.1_refgenome_MITO.fna ]; then
	grep -A 1000 "^>OX457123.1" DeiElpe2.1_refgenome.fna \
	  > DeiElpe2.1_refgenome_MITO.fna;
fi;

# MIA mitochondrial assembly =================================================
if [ ! -d $RESULTS/00_MIA ]; then
	echo "Getting MIA assembly =========================================="
	mkdir -p $RESULTS/00_MIA;
	gunzip -c $DATA/11686_Rophe_h_NoAdapt_PEAR_Trim_Complex_U.fq.gz \
	  > $DATA/11686_Rophe_h_NoAdapt_PEAR_Trim_Complex_U.fastq;
	# subset
	lines=$(wc -l $DATA/11686_Rophe_h_NoAdapt_PEAR_Trim_Complex_U.fastq \
	          | awk '{printf "%.0f",($1 / 4)*0.1}' | awk '{print $0*4}');
	head -n $lines $DATA/11686_Rophe_h_NoAdapt_PEAR_Trim_Complex_U.fastq \
	  > $DATA/11686_Rophe_h_subset.fastq;
	mia -r $REFGENOME/DeiElpe2.1_refgenome_MITO.fna \
	    -f $DATA/11686_Rophe_h_subset.fastq \
	    -m $RESULTS/00_MIA/MIA_assembly \
	    -c;
fi;

# BWA-aln =====================================================================
echo "BWA-aln alignment =============================================="
if [ ! -e $REFGENOME/DeiElpe2.1_refgenome_MITO_index.bwt ]; then
        echo "Indexing the reference genome..."
        bwa index -p $REFGENOME/DeiElpe2.1_refgenome_MITO_index \
                  $REFGENOME/DeiElpe2.1_refgenome_MITO.fna;
fi;

for l in 20 1024; do
	for n in 0.01 0.001 0.0001; do
		echo "bwa aln -l $l -n $n -o 2 ..."
		start=$(date "+%s")
		bwa aln -t 104 \
			-l $l -n $n -o 2 \
			$REFGENOME/DeiElpe2.1_refgenome_MITO_index \
			$DATA/11686_Rophe_h_subset.fastq \
		  | bwa samse $REFGENOME/DeiElpe2.1_refgenome_MITO_index \
		              - \
		              $DATA/11686_Rophe_h_subset.fastq \
		  | samtools view -@ 20 -b -h \
		  | samtools sort -@ 20 -O bam -o $RESULTS/bwaaln_l$l\_n$n.bam;
		# index
		samtools index $RESULTS/bwaaln_l$l\_n$n.bam $RESULTS/bwaaln_l$l\_n$n.bam.bai;
		end=$(date "+%s")
		strategy+=("BWA-aln(l=$l;n=$n)")
		timing+=($((end-start)))
	done
done

# Index ref.genome with Bowtie2 ===============================================
echo "Bowtie2 alignment =============================================="
if [ ! -e $REFGENOME/DeiElpe2.1_refgenome_MITO_BT2index.1.bt2 ]; then
        echo "Indexing reference genome ..."
        bowtie2-build $REFGENOME/DeiElpe2.1_refgenome_MITO.fna \
                      $REFGENOME/DeiElpe2.1_refgenome_MITO_BT2index \
                      --threads 54;
fi;
for L in 20 22; do
	for mp in 4 5 6; do
		echo "bowtie2 --local -N 1 -D 15 -R 2 -L $L --mp $mp..."
		start=$(date "+%s")
                bowtie2 -x $REFGENOME/DeiElpe2.1_refgenome_MITO_BT2index \
                        -U $DATA/11686_Rophe_h_subset.fastq \
                        --local -N 1 -D 15 -R 2 \
                        -L $L --mp $mp \
                        --met-file $RESULTS/bowtie2_L$L\_mp$mp\_bowtie.log \
                        --threads 104 \
                  | samtools view -@ 54 -b -h \
                  | samtools sort -@ 54 -O bam -o $RESULTS/bowtie2_L$L\_mp$mp.bam;
		# index
                samtools index $RESULTS/bowtie2_L$L\_mp$mp.bam $RESULTS/bowtie2_L$L\_mp$mp.bam.bai
		end=$(date "+%s")
                strategy+=("Bowtie2(L=$L;mp=$mp)")
                timing+=($((end-start)))
	done
done

# Get summary =================================================================
for i in "${!timing[@]}"; do
	printf "%s\t%s\n" ${strategy[i]} ${timing[i]} >> $RESULTS/00_performance.tsv;
done;
python $WS/src/02_00.01_summary_optimization.py;

# Get sequences (and make phylogeny) ==========================================
echo ">MIA" > $RESULTS/consensus_sequences.fasta;
head -n 10 $RESULTS/00_MIA/MIA_assembly.4 | tail -n 1 | sed 's/SEQ //' \
  >> $RESULTS/consensus_sequences.fasta;

for BAM in $(ls $RESULTS/*.bam); do
	echo $BAM | cut -d "/" -f 8 | sed 's/^b/>b/' | sed 's/.bam//' \
	  >> $RESULTS/consensus_sequences.fasta
	samtools consensus --min-MQ 20 --ff UNMAP -f FASTA \
			   -T $REFGENOME/DeiElpe2.1_refgenome_MITO.fna \
			   $BAM \
	  | grep -v "^>" \
	  >> $RESULTS/consensus_sequences.fasta;
done

mafft $RESULTS/consensus_sequences.fasta > $RESULTS/consensus_sequences_aln.fasta
iqtree2 -s $RESULTS/consensus_sequences_aln.fasta \
	-T 20 --prefix $RESULTS/phylogeny
rm $DATA/11686_Rophe_h_NoAdapt_PEAR_Trim_Complex_U.fastq;

