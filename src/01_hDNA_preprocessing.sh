#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=25G
#SBATCH --job-name=hDNA_preprocessing
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_preprocessing-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_preprocessing-%j.err

# This script is based on nf-polish pipeline and adapted to Hyles samples
# In summary, for each sample, this script perform the QC, deduplicates,
# remove adapters, merge reads with pear, filters and trims by quality and
# remove low complexity reads.
#
# Despite all the intermediate files, only the three final files are kept:
#  - processed assembled fq.gz ($TRIM/$id\_NoAdapt_PEAR_Trim_Complex_U.fq.gz)
#  - processed unassembled fq.gz ($TRIM/$id\_NoAdapt_PEAR_Trim_Complex_R{1,2}.fq.gz)
#
# The assembled files should be treated as Single end reads, while the
# unassembled files should be treated as Paired End reads.


# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
RAW=$WS/data/fastqraw
TRIM=$WS/data/fastqtrim
threads=20

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

sampleID=$(ls $RAW/*.gz -I "QC" | cut -d "_" -f 1,2,3 | cut -d "/" -f 8 | sort | uniq)
rawFiles=$(ls $RAW/*.gz -I "QC")

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_pre

# First QC reports ============================================================
if [ ! -d $RAW/QC ]; then
	mkdir $RAW/QC;
fi;
echo "Performing Quality Control on Raw FASTQ files..."

# Do FastQC
fastqc -t $threads -o $RAW/QC -q $rawFiles;
# Summarize all reports with multqc
multiqc -f -o $RAW/QC -i "Raw reads Quality" -q $RAW/QC
# Get statistics with seqkit
seqkit stats --threads $threads $rawFiles > $RAW/QC/Raw-reads.stats;

echo "Done!"

# Trim, adapt and dedup =======================================================
# because of how SuperDeduper works, this step has to be run sequencially...
if [ ! -d $TRIM ]; then
	mkdir -p $TRIM;
fi;
if [ ! -d $TRIM/QC ]; then
        mkdir -p $TRIM/QC;
	mkdir -p $TRIM/fastp;
fi;
echo "Deduplicating and trimming adapters..."

for id in $sampleID; do
	echo $id
	# deduplicate (this will save PE_R1.fastq.gz and PE_R2.fastq.gz in wd)
	hts_SuperDeduper -1 $RAW/$id\_1.fq.gz -2 $RAW/$id\_2.fq.gz -f PE \
	                 --stats-file $TRIM/$id.stats;
	# Trim adapters and poly-g tails
	fastp -i PE_R1.fastq.gz -I PE_R2.fastq.gz \
	      -o $TRIM/$id\_NoAdapt_R1.fq.gz -O $TRIM/$id\_NoAdapt_R2.fq.gz \
	      --thread=$threads \
	      --trim_poly_g --detect_adapter_for_pe \
	      --html=report.html --json=report.json;
	# remove intermediate files
	rm PE_R1.fastq.gz PE_R2.fastq.gz report.html report.json;
done
# Get dedup_adapt statistics
trimFiles=$(ls $TRIM/*.gz)
seqkit stats --threads $threads $trimFiles > $TRIM/QC/01_dedup.stats;

echo "Done!"

# Merge =======================================================================
# At this step, I will keep all PEAR output files (assembled and unasembled)
echo "Merging reads with PEAR ..."

for id in $sampleID; do
	echo $id
	pear -f $TRIM/$id\_NoAdapt_R1.fq.gz \
	     -r $TRIM/$id\_NoAdapt_R2.fq.gz \
	     -o $TRIM/$id\_NoAdapt_Merged \
	     --p-value 0.0001 \
	     --min-overlap 20 \
	     -n 40 \
	     --threads $threads;
	# compress files
	gzip -c $TRIM/$id\_NoAdapt_Merged.assembled.fastq \
	  > $TRIM/$id\_NoAdapt_Merged.fq.gz &
	gzip -c $TRIM/$id\_NoAdapt_Merged.unassembled.forward.fastq \
	  > $TRIM/$id\_NoAdapt_Unass_fwd.fq.gz &
	gzip -c $TRIM/$id\_NoAdapt_Merged.unassembled.reverse.fastq \
	  > $TRIM/$id\_NoAdapt_Unass_rev.fq.gz &
	wait;
	# remove intermediate files
	rm $TRIM/$id\_NoAdapt_Merged.assembled.fastq \
	   $TRIM/$id\_NoAdapt_Merged.discarded.fastq \
	   $TRIM/$id\_NoAdapt_Merged.unassembled.forward.fastq \
	   $TRIM/$id\_NoAdapt_Merged.unassembled.reverse.fastq;
	# remove input files (no longer needed)
	rm $TRIM/$id\_NoAdapt_R1.fq.gz $TRIM/$id\_NoAdapt_R2.fq.gz;
done;
mergeFiles=$(ls $TRIM/*Merged*)
seqkit stats --threads $threads $mergeFiles > $TRIM/QC/02_merge.stats;

echo "Done!"

# Final trimming and filtering ================================================
# I'll trim both: assembled and unassembled fq.gz
# This part will
#   (1) remove the first and last 5 nucleotides of each read
#   (2) apply a 4-bp-sliding window from the tail of each read, cuting with Q<15
#   (3) remove any sequence with a complexity lower than 50%
#   (4) remove any sequence lower than 35 nt
echo "Trimming by quality ..."

for id in $sampleID; do
	# For unassembled pairs
	fastp -i $TRIM/$id\_NoAdapt_Unass_fwd.fq.gz -I $TRIM/$id\_NoAdapt_Unass_rev.fq.gz \
	      -o $TRIM/$id\_NoAdapt_PEAR_Trim_Complex_R1.fq.gz -O $TRIM/$id\_NoAdapt_PEAR_Trim_Complex_R2.fq.gz \
	      --disable_adapter_trimming --disable_trim_poly_g \
	      --trim_front1=5 --trim_front2=5 --trim_tail1=5 --trim_tail2=5 \
	      --cut_tail --cut_tail_window_size=4 --cut_tail_mean_quality=15 \
	      --low_complexity_filter --complexity_threshold=50 \
	      -l 35 --thread=4 \
	      --html=$TRIM/fastp/$id\_PE.html --json=$TRIM/fastp/$id\_PE.html --report_title="${id} PE" \
	  && rm $TRIM/$id\_NoAdapt_Unass_fwd.fq.gz $TRIM/$id\_NoAdapt_Unass_rev.fq.gz &
	# For assembled pairs
	fastp -i $TRIM/$id\_NoAdapt_Merged.fq.gz \
              -o $TRIM/$id\_NoAdapt_PEAR_Trim_Complex_U.fq.gz \
              --disable_adapter_trimming --disable_trim_poly_g \
              --trim_front1=5 --trim_tail1=5 \
              --cut_tail --cut_tail_window_size=4 --cut_tail_mean_quality=15 \
              --low_complexity_filter --complexity_threshold=50 \
              -l 35 --thread=3 \
              --html=$TRIM/fastp/$id\_SE.html --json=$TRIM/fastp/$id\_SE.html --report_title="${id} SE" \
	  && rm $TRIM/$id\_NoAdapt_Merged.fq.gz &
done;
wait;

trimFiles=$(ls $TRIM/*.fq.gz)
seqkit stats --threads $threads $trimFiles > $TRIM/QC/03_trim.stats;

echo "Done!"

# Final QC reports ============================================================
echo "Performing Quality Control on Trimmed FASTQ files..."

# Do FastQC
fastqc -t $threads -o $TRIM/QC -q $trimFiles;
# Summarize all reports with multqc
multiqc -f -o $TRIM/QC -i "Processed reads Quality" -q $TRIM/QC

echo "Done!"

conda deactivate
