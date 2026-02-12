#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --job-name=hDNA_mapping_unassembled
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_unassembled-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_unassembled-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
DATA=$WS/data/fastqtrim
REFGENOME=$WS/data/refgenome
RESULTS=$WS/results/01_hDNA_alignment
SampleID=$(ls $DATA/*.gz -I "QC" | cut -d "/" -f 8 | sed 's/_h.*/_h/' | sort | uniq)

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# Mapping against the reference genome ========================================
if [ ! -d $RESULTS ]; then
	echo "Mapping the reads =============================================="
	mkdir -p $RESULTS;
	mkdir $RESULTS/QC;
fi;
for i in $SampleID; do
	if [ ! -e $RESULTS/$i\_paired_unmapped_bwamem.bam -a ! -e $RESULTS/$i\_end.bam ]; then
		echo "Sample ${i}...";
		bowtie2 -x $REFGENOME/DeiElpe2.1_refgenome_BT2index \
                        -1 $DATA/$i\_NoAdapt_PEAR_Trim_Complex_R1.fq.gz \
			-2 $DATA/$i\_NoAdapt_PEAR_Trim_Complex_R2.fq.gz \
                        --local -N 1 -D 15 -R 2 \
                        -L 20 --mp 4 \
                        --met-file $RESULTS/bowtie2_L$L\_mp$mp\_bowtie.log \
                        --threads 32 \
                  | samtools view -@ 32 -b -h \
                  | samtools sort -@ 32 -O bam -o $RESULTS/$i\_unmerged.bam;
	fi;
done;

conda deactivate
