#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=84
#SBATCH --time=3:00:00
#SBATCH --mem=90G
#SBATCH --job-name=freshDNA_preprocessing
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/freshDNA_preprocessing-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/freshDNA_preprocessing-%j.err

## Configuration hint:
##  - set 1 cpu per fq.gz file
##  - set at least 1GB mem per fq.gz.file

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
RAW=$WS/data/fastqraw_fresh
TRIM=$WS/data/fastqtrim_fresh
trimmomatic=$WS/src/Trimmomatic-0.39/trimmomatic-0.39.jar
adapters=$WS/src/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa
threads=84

sampleID=($(ls $RAW/*.gz | cut -d "_" -f 1,2,3 | cut -d "/" -f 8 | sort | uniq))

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/fastq_preprocessing

# First QC reports ============================================================
echo "Performing Quality Control on Raw FASTQ files..."
if [ ! -d $RAW/QC ]; then
	mkdir $RAW/QC;
	qc_raw=1
	rawFiles=($(ls $RAW/*.gz))
else
	# find processed files (so as to not process them again!)
	processedID=($(ls $RAW/QC/*_1_fastqc.zip | cut -d "_" -f 1,2,3 | cut -d "/" -f 9 | sort))
	joined_arrays=( "${sampleID[@]}" "${processedID[@]}" )
	not_processed=($(echo ${joined_arrays[@]} | tr " " "\n" | sort | uniq -c | grep "   1 " | sed 's/\ \ \+//' | cut -d " " -f 2))
        if [ ${#not_processed[@]} -le 1 -a ${not_processed[0]}=="" ]; then
		qc_raw=0
	else
		qc_raw=1
		rawFiles=()
		for i in ${not_processed[@]}; do
			rawFiles+=($(ls $RAW/$i*))
		done
	fi;
fi;

if [ $qc_raw -eq 1 ]; then
	echo ""	# Do FastQC
	fastqc -t $threads -o $RAW/QC -q ${rawFiles[@]};
	# Summarize all reports with multqc
	multiqc -f -o $RAW/QC -i "Raw reads Quality" -q $RAW/QC
	# Get statistics with seqkit
	seqkit stats --threads $threads ${rawFiles[@]} > $RAW/QC/Raw-reads.stats;

	echo "Done!"
else
	echo "All files were already processed. Nothing was done"
fi;

## trimming and filtering =====================================================
echo "Performing trimming and filtering with fastp..."

if [ ! -d $TRIM ]; then
	mkdir $TRIM;
	mkdir $TRIM/fastp
fi;

for i in ${sampleID[@]}; do
	if [ ! -e $TRIM/$i\_1.fq.gz ]; then
		echo "    Sample ${i}"
		fastp -i $RAW/$i\_1.fq.gz -I $RAW/$i\_2.fq.gz -o $TRIM/$i\_1.fq.gz -O $TRIM/$i\_2.fq.gz \
		      --trim_front1=10 --trim_front2=10 --trim_tail1=10 --trim_tail2=10 \
		      --cut_tail --cut_tail_window_size=10 --cut_tail_mean_quality=25 \
		      -q 25 -u 40 -e 25 \
		      -l 100 \
		      --thread 2 \
		      --html $TRIM/fastp/$i.html --json $TRIM/fastp/$i.json --report_title="${i}" &
	fi;
done;
echo "All executions have been launched. Waiting for them to finish..."
wait
echo "Done!"

## Final FASTQC ===============================================================
echo "Performing Quality Control on Trimmed FASTQ files..."

if [ ! -d $TRIM/QC ]; then
        mkdir $TRIM/QC;
        qc_trim=1
	trimFiles=($(ls $TRIM/*.gz -I "QC"))
else
        # find processed files (so as to not process them again!)
        processedID=($(ls $TRIM/QC/*_1_fastqc.zip | cut -d "_" -f 1,2,3 | cut -d "/" -f 9 | sort))
        joined_arrays=( "${sampleID[@]}" "${processedID[@]}" )
        not_processed=($(echo ${joined_arrays[@]} | tr " " "\n" | sort | uniq -c | grep "   1 " | sed 's/\ \ \+//' | cut -d " " -f 2))
        if [ ${#not_processed[@]} -le 1 -a ${not_processed[0]}=="" ]; then
                qc_trim=0
        else
                qc_trim=1
                trimFiles=()
                for i in ${not_processed[@]}; do
                        trimFiles+=($(ls $RAW/$i*))
                done
        fi;
fi;

if [ $qc_trim -eq 1 ]; then
	fastqc -t $threads -o $TRIM/QC -q $trimFiles;
	# Summarize all reports with multqc
	multiqc -f -o $TRIM/QC -i "Trimmed reads Quality" -q $TRIM/QC
	# Get statistics with seqkit
	seqkit stats --threads $threads $trimFiles > $TRIM/QC/Trim-reads.stats;
	echo "Done!"
else
        echo "All files were already processed. Nothing was done"
fi;

