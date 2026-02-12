#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mem=8G
#SBATCH --account=p_hyleshawkmoths
#SBATCH --job-name=hDNA_mito_consensus
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mito_consensus-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mito_consensus-%j.err

module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
HDNA=$WS/results/01_hDNA_alignment
FRESH=$WS/results/02_freshDNA_alignment
REFGENOME=$WS/data/refgenome/DeiElpe2.1_refgenome.fna
RESULTS=$WS/results/05_MITO_phylogeny/00_consensus_fasta
hdna_files=$(ls $HDNA/*_rescale.bam)
fresh_files=$(ls -1 $FRESH/*.bam | grep -v "Deilephila")

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

if [ ! -d $RESULTS ]; then
        mkdir -p $RESULTS;
fi;

# Get consensus seqs
echo "Getting mitochondria consensus sequences from BAM alignments..."
for file in $hdna_files $fresh_files; do
        id=$(echo $file | sed 's/.*_alignment\///' | sed 's/_end_rescale.bam//' | sed 's/.bam//')
        if [ ! -e $RESULTS/$id.fa ]; then
                echo "Sample $id";
                samtools consensus -f FASTA \
                                   --region OX457123.1 \
                                   --reference $REFGENOME \
                                   --output $RESULTS/$id.fas \
                                   --min-BQ 15 --min-MQ 20 --min-depth 5 \
                                   -@ 10 \
                                   $file;
        fi;
done;
echo "DONE!"
echo "Output: $RESULTS"

# Concatenate files in different datasets
echo "Concatenating fasta files..."
if [ ! -e $RESULTS/mitogenomes_full.fa ]; then
	for fa in $(ls $RESULTS/*.fas | sed 's/.*00_consensus_fasta\///'); do
 		sed 's/OX457123.1/'"$fa"'/' $RESULTS/$fa \
		  | sed 's/.fas//' \
		  >> $RESULTS/mitogenomes_full.fa;
	done;
fi;
if [ ! -e $RESULTS/mitogenomes_1samp-x-sp.fa ]; then
	for fa in $(cut -d "/" -f 8 $WS/results/03_variant_calling/00_best_samples_per_species.txt \
                      | sed 's/.bam/.fas/' \
                      | sed 's/_end_rescale//' ); do
                sed 's/OX457123.1/'"$fa"'/' $RESULTS/$fa \
                  | sed 's/.fas//' \
                  >> $RESULTS/mitogenomes_1samp-x-sp.fa;
		if [[ ! "$fa" =~ _h\.fas$ ]]; then
			sed 's/OX457123.1/'"$fa"'/' $RESULTS/$fa \
	                  | sed 's/.fas//' \
	                  >> $RESULTS/mitogenomes_1samp-x-sp_noRhodafra.fa;
		fi;
	done
fi;
