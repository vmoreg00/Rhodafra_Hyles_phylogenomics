#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=hDNA_mapping_merge
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_merge-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_merge-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
DATA=$WS/data/fastqtrim
RESULTS=$WS/results/01_hDNA_alignment
REFGENOME=$WS/data/refgenome
SampleID=$(ls $DATA/*.gz -I "QC" | cut -d "/" -f 8 | sed 's/_h.*/_h/' | sort | uniq)

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# Mapping against the reference genome ========================================
if [ ! -d $RESULTS/QC ]; then
	mkdir $RESULTS/QC;
	#mkdir $RESULTS/QC/Heterozygosity
fi;
echo -e "SampleID\tInitialReads\tCoverage\tMappedReads\n" > $RESULTS/QC/00_summary.stats;
for i in $SampleID; do
	if [ -e $RESULTS/$i\_merged.bam -a -e $RESULTS/$i\_unmerged.bam -a ! -e $RESULTS/$i\_end.bam ]; then
		echo "Sample ${i}...";
		# get some previous stats for each input bam
		for file in merged.bam unmerged.bam ; do
			STATS="${i}"
			## Initial reads
			STATS=$STATS"\t"$(\
			  grep $i $DATA/QC/03_trim.stats \
			  | head -n 1 | sed 's/,//g' | awk '{print $4}')
			## Coverage
	                STATS=$STATS"\t"$(\
			  samtools coverage $RESULTS/$i\_$file \
	                    | grep -v "^#" \
	                    | awk 'BEGIN{COV=0}; {if $1 != "CATKIA010000001.1" & $1 != "OX457123.1" {COV+=$6}}; END{print COV/30}')
	                ## number of mapped reads
	                STATS=$STATS$"\t"$(\
			  samtools stats -@ 16 $RESULTS/$i\_$file \
	                    | head -n 7 | tail -n 1)
			echo -e $STATS"\n" >> $RESULTS/QC/00_summary.stats;
		done
		# concatenate all bams
		echo "    ... concatenating bam"
		samtools cat $RESULTS/$i\_merged.bam \
			     $RESULTS/$i\_unmerged.bam \
                  | samtools view -F 4 -q 20 -@ $threads -h -b - \
		  | samtools sort -@ 32 -O bam -o $RESULTS/$i\_end.bam \
		  && rm $RESULTS/$i\_merged.bam $RESULTS/$i\_unmerged.bam;
                # QC (mapDamage2)
                echo "    ... estimating DNA damage"
                mapDamage -n 100000 -i $RESULTS/$i\_end.bam \
                          -r $REFGENOME/DeiElpe2.1_refgenome.fna \
                         --rescale --rescale-out $RESULTS/$i\_end_rescale.bam \
		  && rm $RESULTS/$i\_end.bam;
                mv results_$i\_end/ $RESULTS/QC/$i\_mapDamage/;

                # index
		echo "    ... indexing concatenated bam"
                samtools index -@ 32 $RESULTS/$i\_end_rescale.bam $RESULTS/$i\_end_rescale.bam.bai;

		# Get summary statistics (samtools)
		echo "    ... retrieving summary statisitics"
		samtools stats -@ 16 $RESULTS/$i\_end_rescale.bam > $RESULTS/QC/$i.stats;
                samtools flagstats -@ 16 -O tsv $RESULTS/$i\_end_rescale.bam > $RESULTS/QC/$i.flagstats.tsv;
                samtools coverage $RESULTS/$i\_end_rescale.bam > $RESULTS/QC/$i.coverage.tsv;

		# Qualimap
		echo "    ... Getting mapping quality report"
		qualimap bamqc -bam $RESULTS/$i\_end_rescale.bam -c -nr 10000 -nt 32 \
			       -outdir $RESULTS/QC/$i\_qualimap;

		# Estimate contamination ===
		echo "    ... estimating contamination"
		echo -ne "${i}\t" >> $RESULTS/QC/mitochondrion_heterozygosity.tsv;
		## get mitochondrial genome (will be used as reference)
                samtools consensus -r OX457123.1 -f FASTA $RESULTS/$i\_end_rescale.bam \
                  > $RESULTS/QC/Heterozygosity/$i\_mito.fna;
		## index mitochondrial genome
		bwa index -p $RESULTS/QC/Heterozygosity/$i\_mito \
			  $RESULTS/QC/Heterozygosity/$i\_mito.fna;
		## map reads against mitochondrial genome
		bwa aln -t 36 \
			$RESULTS/QC/Heterozygosity/$i\_mito \
			$DATA/$i\_NoAdapt_PEAR_Trim_Complex_U.fq.gz \
		  | bwa samse $RESULTS/QC/Heterozygosity/$i\_mito \
                              - \
                              $DATA/$i\_NoAdapt_PEAR_Trim_Complex_U.fq.gz \
                  | samtools view -@ 36 -b \
                  | samtools sort -@ 36 -o $RESULTS/QC/Heterozygosity/$i\_mito.bam;
		## variant calling and estimate heterozygosity
		## The grep command looks for "differences" against the mitochondrial genome
                bcftools mpileup -Ov --threads 36 \
                                 -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
                                 -f $RESULTS/QC/Heterozygosity/$i\_mito.fna \
				 $RESULTS/QC/Heterozygosity/$i\_mito.bam \
                  | bcftools call -mOv --threads 36 \
                                  --ploidy 2 \
                                  -V indels \
		  | bcftools view --no-header -Ov \
		  | cut -f 10 \
		  | grep -c "0/1\|1/0\|1/1" \
		  >> $RESULTS/QC/mitochondrion_heterozygosity.tsv;
		## remove intermediate files
		rm $RESULTS/QC/Heterozygosity/*
		echo "    DONE!"
	fi;
done;

conda deactivate
