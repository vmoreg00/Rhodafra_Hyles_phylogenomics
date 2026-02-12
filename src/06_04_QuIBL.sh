#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2:00:00
#SBATCH --account=p_hyleshawkmoths
#SBATCH --job-name=hDNA_introg_QuIBL
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_introg_QuIBL-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_introg_QuIBL-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3


# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
DATA=$WS/results/03_variant_calling
TREES=$WS/results/04_SNP_phylogeny/06_PhyloNet/1samp-x-sp_raxml_sliding_good.trees.gz
RESULTS=$WS/results/06_Introgression/04_QuIBL
QuIBL=$WS/src/QuIBL/QuIBL.py
QuIBL=$WS/src/QuIBL/cython_vers/QuIBL_cyth.py
summQuIBL=$WS/src/QuIBL/analysis/summaryFileAnalysis.R

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

if [ ! -d $RESULTS ]; then
        mkdir -p $RESULTS;
fi;

# Build 'gene' trees (raxml_sliding_windows)
# Subset the VCF
if [ ! -e $RESULTS/1samp-x-sp_snps_subset.vcf.gz ]; then
        bcftools view -h $DATA/variant_calling_1samp-x-sp.vcf.gz \
          | tail -n 1 \
          | sed 's/.*FORMAT\t//' \
          | tr '\t' '\n' \
          | grep "Hyles\|Rophe\|Rmars\|Xylophanes-tersa" \
          > $RESULTS/spp_list.txt;
        bcftools view -O z -S $RESULTS/spp_list.txt \
          $DATA/variant_calling_1samp-x-sp.vcf.gz \
          > $RESULTS/1samp-x-sp_snps_subset.vcf.gz;
fi;
# retrieving geno.gz file (change the headers for simplicity)
# I also remove indels manually with grep because they are problematic in RAxML
if [ ! -e $RESULTS/1samp-x-sp_snps_subset_unphased.geno.gz ]; then
        echo "Getting the geno.gz file..."
        python $WS/src/genomics_general/VCF_processing/parseVCF.py \
               --skipIndels \
               -i $RESULTS/1samp-x-sp_snps_subset.vcf.gz \
          | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
          | sed 's/01_hDNA_alignment\///g' \
          | sed 's/02_freshDNA_alignment\///g' \
          | sed 's/\.bam//g' \
          | awk '{if(NR==1){print $0} else {print toupper($0)}}' \
          | grep -Ev "([ACGT]){2,}" \
          | bgzip \
          > $RESULTS/1samp-x-sp_snps_subset_phased.geno.gz;
        # I want the geno file in 'diplo' format (not phased)
        # Remove Deilephila as it is a problematic outgroup (branches too long)
        conda activate $WS/envirs/py2.7;
        export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general
        python $WS/src/genomics_general/filterGenotypes.py \
               -i $RESULTS/1samp-x-sp_snps_subset_phased.geno.gz -if phased \
               -o $RESULTS/1samp-x-sp_snps_subset_unphased.geno.gz -of diplo \
               --threads 20;
        conda deactivate
fi;
# sliding window phylogeny with RAxML
# 10kb-sliding windows separated by 50kb to mitigate the effects of recombination
if [ ! -d $RESULTS/1samp-x-sp_raxml_sliding.trees.gz ]; then
        echo "    - Getting 'gene' trees..."
        conda activate $WS/envirs/py2.7;
        export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general
        python2 $WS/src/genomics_general/phylo/raxml_sliding_windows.py \
              -g $RESULTS/1samp-x-sp_snps_subset_unphased.geno.gz \
              --prefix $RESULTS/1samp-x-sp_raxml_sliding \
              --windType coordinate -w 10000 -S 60000 -M 10 -O 0 --model GTRCAT \
              --outgroup 6384_Xylophanes-tersa \
              --raxml $WS/envirs/astral/bin/raxmlHPC \
              --threads 40;
        conda deactivate;
        # remove NA trees to avoid failures in astral
        gunzip -c $RESULTS/1samp-x-sp_raxml_sliding.trees.gz \
          | awk '{if($0 != "NA"){print $0}}' \
          > $RESULTS/1samp-x-sp_raxml_sliding_good.trees;
fi;

# Prepare input file
if [ ! -e $RESULTS/QuIBL_input.txt ]; then
	# Download the templat
	wget https://raw.githubusercontent.com/miriammiyagi/QuIBL/refs/heads/master/Small_Test_Example/sampleInputFile.txt \
	  -O $RESULTS/QuIBL_input_template.txt;
	# Change the settings accordingly
	cat $RESULTS/QuIBL_input_template.txt \
	  | sed 's|treefile:.*|treefile: '$RESULTS'/1samp-x-sp_raxml_sliding_good.trees|g' \
	  | sed 's|totaloutgroup:.*|totaloutgroup: 6384_Xylophanes-tersa|g' \
	  | sed 's|maxcores:.*|maxcores: 40|g' \
	  | sed 's|OutputPath:.*|OutputPath: '$RESULTS'/QuIBL_results.csv|g' \
	  > $RESULTS/QuIBL_input.txt;
	# Remove the template
	rm $RESULTS/QuIBL_input_template.txt;
fi;

if [ ! -e $RESULTS/QuIBL_results.csv ]; then
	# Run QuIBL
	echo "Running QuIBL..."
	python $QuIBL $RESULTS/QuIBL_input.txt;
	Rscript --vanilla $summQuIBL -i $RESULTS/QuIBL_results.csv \
	                             -o $RESULTS/QuIBL_summary \
	                             -s -10 \
	                             -t 6384_Xylophanes-tersa;
fi;


echo "DONE!"
echo "Output: $RESULTS"
