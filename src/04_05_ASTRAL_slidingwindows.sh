#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --time=4:00:00
#SBATCH --mem=16GB
#SBATCH --account=p_hyleshawkmoths
#SBATCH --job-name=hDNA_phylo_ASTRAL
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_ASTRAL-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_ASTRAL-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
MSA=$WS/results/03_variant_calling
RESULTS=$WS/results/04_SNP_phylogeny/05_ASTRAL

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/astral
export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general

if [ ! -d $RESULTS ]; then
        mkdir -p $RESULTS
fi;

# I only run this for the 1samp-x-sp and strict datasets. Technically, ASTRAL will not
# return a species tree if I provide several samples per species.
for dataset in 1samp-x-sp strict; do
	echo "==> Dataset: $dataset"
        # retrieving geno.gz file (change the headers for simplicity)
	# I also remove indels manually with grep because they are problematic in RAxML
        if [ ! -e $RESULTS/$dataset\_snps_subset_unphased.geno.gz ]; then
                echo "Getting the geno.gz file..."
                python $WS/src/genomics_general/VCF_processing/parseVCF.py \
		       --skipIndels \
                       -i $MSA/variant_calling_$dataset.vcf.gz \
                  | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
                  | sed 's/01_hDNA_alignment\///g' \
                  | sed 's/02_freshDNA_alignment\///g' \
                  | sed 's/\.bam//g' \
		  | awk '{if(NR==1){print $0} else {print toupper($0)}}' \
		  | grep -Ev "([ACGT]){2,}" \
                  | bgzip \
                  > $RESULTS/$dataset\_snps_subset_phased.geno.gz;
		# I want the geno file in 'diplo' format (not phased)
		# Remove Deilephila as it is a problematic outgroup (branches too long)
		python $WS/src/genomics_general/filterGenotypes.py \
		       -i $RESULTS/$dataset\_snps_subset_phased.geno.gz -if phased \
		       -o $RESULTS/$dataset\_snps_subset_unphased.geno.gz -of diplo \
		       --excludeSamples 1238_Deilephila-elpenor,502_Deilephila-porcellus \
		       --threads 20;
        fi;

        # sliding window phylogeny with RAxML
        if [ ! -d $RESULTS/$dataset\_raxml_sliding.trees ]; then
                echo "    - Getting 'gene' trees..."
		conda activate $WS/envirs/py2.7;
		export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general
                python2 $WS/src/genomics_general/phylo/raxml_sliding_windows.py \
                      -g $RESULTS/$dataset\_snps_subset_unphased.geno.gz \
                      --prefix $RESULTS/$dataset\_raxml_sliding \
                      --windType coordinate -w 50000 -S 50000 -M 100 -O 0 --model GTRCAT \
		      --raxml $WS/envirs/astral/bin/raxmlHPC \
		      --threads 104;
		conda deactivate;
		# remove NA trees to avoid failures in astral
                gunzip -c $RESULTS/$dataset\_raxml_sliding.trees.gz \
                  | awk '{if($0 != "NA"){print $0}}' \
                  > $RESULTS/$dataset\_raxml_sliding.trees;
	fi;
	# Run Astral
	if [ ! -e $RESULTS/$dataset\_astral.nwk ]; then
                echo "    - Inferring species tree..."
		# Run ASTRAL
		astral -i $RESULTS/$dataset\_raxml_sliding.trees \
		       -o $RESULTS/$dataset\_astral.nwk \
		       -t 2 --seed 666 \
		  2> $RESULTS/$dataset\_astral.log;
		# Run the same, but with polytomies test annotation
		astral -i $RESULTS/$dataset\_raxml_sliding.trees \
	               -o $RESULTS/$dataset\_astral_test_polytomies.nwk \
	               -t 2 --seed 666 \
	          2> $RESULTS/$dataset\_astral_test_polytomies.log;
	fi;
done;



echo -e "\n\n\n"
echo -e "------------------------------------------------------------\n\n"
echo "Session Info:"
echo "- parseVCF.py, filterGenotypes.py and raxml_sliding_windows.py retrieved from github.com/simonhmartin/genomics_general"
echo "- raxml-ng"
$WS/envirs/hDNA_map/bin/raxml-ng --version
echo "- ASTRAL version 5.7.8"
