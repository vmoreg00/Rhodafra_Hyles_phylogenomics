#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=4:00:00
#SBATCH --job-name=hDNA_phylo_SVDquartets
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_SVDquartets-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_SVDquartets-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
MSA=$WS/results/03_variant_calling
RESULTS=$WS/results/04_SNP_phylogeny/03_SVDquartets

paup=$WS/src/paup4a168_centos64

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

if [ ! -d $RESULTS ]; then
        mkdir -p $RESULTS
fi;

# SVDquartets analysis (_no508 and _noRhodafra datasets not analyzed)
for dataset in full 1samp-x-sp strict; do
	echo "Analysing dataset $dataset..."
	# Simplify names to avoid problems
	if [ ! -e $RESULTS/snps_$dataset.nex ]; then
		echo "  - Getting input files..."
		cat $MSA/variant_calling_$dataset\_1kb-thin_consensus.nexus \
		  | sed 's/Hyles-/H/g' \
		  | sed 's/-nicaea//g' \
		  | sed 's/-catissima//g' \
		  | sed 's/-sheljuzkoi//g' \
		  | sed 's/Chaerocina-dohertyi/Chadoh/g' \
		  | sed 's/Euchloron-megaera/Eucmeg/g' \
		  | sed 's/Hippotion-celerio/Hipcel/g' \
		  | sed 's/Hippotion-echeclus/Hipech/g' \
		  | sed 's/Theretra-alecto/Theale/g' \
		  | sed 's/Theretra-japonica/Thejap/g' \
		  | sed 's/Theretra-silhetensis/Thesil/g' \
		  | sed 's/Xylophanes-falco/Xylfal/g' \
		  | sed 's/Xylophanes-porcus/Xylpor/g' \
		  | sed 's/Xylophanes-tersa/Xylter/g' \
		  > $RESULTS/snps_$dataset.nex
		# add the SVDquartets PAUP code to the nex file
		cat $WS/src/04_03_PAUP_code_$dataset.nex \
		  | sed 's@mrp.txt@'"$RESULTS"'\/'"$dataset"'_mrp.txt@g' \
		  | sed 's@quartet.qmc@'"$RESULTS"'\/'"$dataset"'_quartet.qmc@g' \
		  | sed 's@treefile.nex@'"$RESULTS"'\/'"$dataset"'_treefile.nex@g' \
		  | sed 's@treefile_consensus.nwk@'"$RESULTS"'\/'"$dataset"'_treefile_consensus.nwk@g' \
		  >> $RESULTS/snps_$dataset.nex;
	fi;
	# Fit SVD quartets
	if [ ! -e $RESULTS/$dataset\_treefile.nex ]; then
		echo "  - Fitting PAUP..."
		$paup -n $RESULTS/snps_$dataset.nex;
	fi;
        echo "DONE!"
        echo "Output: $RESULTS/$dataset\.*"
done;

echo -e "\n\n\n"
echo -e "------------------------------------------------------------\n\n"
echo "Session Info:"
echo "- PAUP"
$paup --version


