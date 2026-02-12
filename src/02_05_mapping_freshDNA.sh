#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH --job-name=freshDNA_mapping
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/freshDNA_mapping-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/freshDNA_mapping-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
DATA=$WS/data/fastqtrim_fresh
REFGENOME=$WS/data/refgenome
RESULTS=$WS/results/02_freshDNA_alignment
SampleID=$(ls $DATA/*.gz -I "QC" | cut -d "/" -f 8 | sed 's/_[12].*//' | sort | uniq)
threads=52
# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# Index ref.genome with Bowtie2 ===============================================
if [ ! -e $REFGENOME/DeiElpe2.1_refgenome_BT2index.1.bt2 ]; then
	echo "Indexing reference genome ======================================"
	bowtie2-build $REFGENOME/DeiElpe2.1_refgenome.fna \
		      $REFGENOME/DeiElpe2.1_refgenome_BT2index \
		      --threads $threads;
fi;

# Mapping against the reference genome ========================================
if [ ! -d $RESULTS ]; then
	echo "Mapping the reads =============================================="
	mkdir -p $RESULTS;
	mkdir $RESULTS/QC;
fi;
for i in $SampleID; do
	if [ ! -e $RESULTS/$i.bam ]; then
		echo "Sample ${i}...";
		# map with bowtie2
		echo "    Mapping with bowtie..."
		bowtie2 -x $REFGENOME/DeiElpe2.1_refgenome_BT2index \
		        -1 $DATA/$i\_1.fq.gz -2 $DATA/$i\_2.fq.gz \
		        --local -N 0 \
		        --met-file $RESULTS/$i\bowtie.log \
		        --threads $threads \
		  | samtools view -F 4 -q 30 -@ $threads -h -b - \
		  | samtools sort -n -O BAM -@ $threads - \
		  | samtools fixmate -m -@ $threads - - \
		  | samtools sort -O BAM -@ $threads - \
		  | samtools markdup -r -@ $threads -O BAM - $RESULTS/$i.bam;
		# Index the bam file
		samtools index $RESULTS/$i.bam $RESULTS/$i.bam.bai;
		# computing statistics
                echo "    Retrieving summary statisitics..."
                samtools stats -@ 25 $RESULTS/$i.bam > $RESULTS/QC/$i.stats &
                samtools flagstats -@ 25 -O tsv $RESULTS/$i.bam > $RESULTS/QC/$i.flagstats.tsv &
                samtools coverage $RESULTS/$i.bam > $RESULTS/QC/$i.coverage.tsv &
		wait;

		# Getting mapping quality report
		echo "    Getting mapping quality report..."
                qualimap bamqc -bam $RESULTS/$i.bam -c -nr 10000 -nt 20 \
                               -outdir $RESULTS/QC/$i\_qualimap;
	fi;
done;

conda deactivate
