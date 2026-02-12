#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mem=16G
#SBATCH --job-name=hDNA_phylo_consensus
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_consensus-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_consensus-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
DATA=$WS/results/03_variant_calling

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# Create SNP consensus for each dataset
## Deilephila is subsequently removed because their branches are too long and
## negatively impacts on phylogeny inference.
for dataset in full full_no508 1samp-x-sp 1samp-x-sp_noRhodafra strict; do
	if [ ! -e $DATA/variant_calling_$dataset\_1kb-thin_consensus.fa ]; then
		echo "Getting the consensus fasta from VCF.GZ file..."
		# consider also using --resolve-IUPAC to remove ambiguities
		python $WS/src/vcf2phylip.py -i $DATA/variant_calling_$dataset\_1kb-thin.vcf.gz \
					     --output-folder $DATA \
					     --output-prefix variant_calling_$dataset\_1kb-thin_consensus \
					     -f -n;
		# Shortening names and removing Deilephila
		## FASTA
		cat $DATA/variant_calling_$dataset\_1kb-thin_consensus.min4.fasta \
		  | awk -v FS="/" '{if($0 ~ /^>/) {gsub("_end_rescale.bam","", $8); print ">"$8} else print $0}' \
		  | awk -v FS="/" '{if($0 ~ /^>/) {gsub(".bam","", $0); print} else print $0}' \
                  | awk -v P=True '{ \
                      if($0 ~ /^>/){ \
                        if($0 ~ /Deilephila/){ \
                          P="False"; \
                        } else { \
                          P="True"; \
                        } \
                      } \
                      if(P=="True"){ \
                        print($0) \
                      } \
                    }' \
		  > $DATA/variant_calling_$dataset\_1kb-thin_consensus.fa \
		  && rm $DATA/variant_calling_$dataset\_1kb-thin_consensus.min4.fasta;
		## PHYLLIP
		cat $DATA/variant_calling_$dataset\_1kb-thin_consensus.min4.phy \
		  | awk '{gsub("_end_rescale.bam",""); print}' \
		  | awk '{gsub(".bam",""); print}' \
		  | awk '{gsub("/data/horse/ws/vimo762h-hDNAhyles/results/01_hDNA_alignment/",""); print}' \
		  | awk '{gsub("/data/horse/ws/vimo762h-hDNAhyles/results/02_freshDNA_alignment/",""); print}' \
		  | grep -v "Deilephila" \
		  | awk -v FS=" " '{if(NR==1){print $1-2" "$2}else{print $0}}' \
		  > $DATA/variant_calling_$dataset\_1kb-thin_consensus.phy \
		  && rm $DATA/variant_calling_$dataset\_1kb-thin_consensus.min4.phy;
		## NEXUX
		cat $DATA/variant_calling_$dataset\_1kb-thin_consensus.min4.nexus \
	          | awk '{gsub("_end_rescale.bam",""); print}' \
	          | awk '{gsub(".bam",""); print}' \
	          | awk '{gsub("/data/horse/ws/vimo762h-hDNAhyles/results/01_hDNA_alignment/",""); print}' \
	          | awk '{gsub("/data/horse/ws/vimo762h-hDNAhyles/results/02_freshDNA_alignment/",""); print}' \
		  | grep -v "Deilephila" \
		  | awk -F='[ =]' '{if(NR==3){print $1" "$2"="$3-2" "$4"="$5}else{print $0}}' \
	          > $DATA/variant_calling_$dataset\_1kb-thin_consensus.nexus \
	          && rm $DATA/variant_calling_$dataset\_1kb-thin_consensus.min4.nexus;
		echo "DONE!"
		echo "Output: $DATA/variant_calling_$dataset\_1kb-thin_consensus.fa"
	fi;
done;
