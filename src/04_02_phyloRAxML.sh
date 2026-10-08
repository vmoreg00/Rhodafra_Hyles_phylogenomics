#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --time=16:00:00
#SBATCH --job-name=hDNA_phylo_RAXmL
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_RAXmL-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_RAXmL-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
MSA=$WS/results/03_variant_calling
RESULTS=$WS/results/04_SNP_phylogeny/02_RAxML

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS
fi;


# Under the model, I add the ASC_FELS option, as stated by Suissa et al. (2024)

for dataset in full full_no508 1samp-x-sp 1samp-x-sp_noRhodafra strict; do
        if [ ! -e $RESULTS/SNP_phylo_raxml_$dataset.raxml.log ]; then
		echo "Removing invariant sites..."
		$WS/envirs/hDNA_map/bin/python \
			 $WS/src/ascbias.py -p $MSA/variant_calling_$dataset\_1kb-thin_consensus.phy \
		                            -o $RESULTS/0_SNPs_$dataset\_noinv.phy;
		echo "Running RAxML-ng ($dataset dataset)..."
		raxml-ng --msa $RESULTS/0_SNPs_$dataset\_noinv.phy \
		         --data-type DNA \
			 --all \
			 --model GTR+G4+ASC_FELS{1000} \
			 --tree pars{10} \
			 --bs-trees 1000 \
			 --threads 52 \
			 --outgroup 6386_Theretra-alecto,6389_Theretra-japonica,6381_Theretra-silhetensis \
			 --prefix $RESULTS/SNP_phylo_raxml_$dataset;

		echo "Getting concordance factors (with IQtree)..."
                iqtree -t $RESULTS/SNP_phylo_raxml_$dataset.raxml.support \
                       --scf 100 -T 20 \
                       -s $RESULTS/0_SNPs_$dataset\_noinv.phy;
                echo "Estimating gene concordance factors (gCF) from 50 kb sliding-window trees..."
                ## get the geno file
                python $WS/src/genomics_general/VCF_processing/parseVCF.py \
                       --skipIndels \
                       -i $MSA/variant_calling_$dataset\_1kb-thin.vcf.gz \
                  | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
                  | sed 's/01_hDNA_alignment\///g' \
                  | sed 's/02_freshDNA_alignment\///g' \
                  | sed 's/\.bam//g' \
                  | awk '{if(NR==1){print $0} else {print toupper($0)}}' \
                  | grep -Ev "([ACGT]){2,}" \
                  | bgzip \
                  > $RESULTS/$dataset\_snps_subset_phased.geno.gz;
                conda activate $WS/envirs/py2.7;
                export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general
                python2 $WS/src/genomics_general/filterGenotypes.py \
                       -i $RESULTS/$dataset\_snps_subset_phased.geno.gz -if phased \
                       -o $RESULTS/$dataset\_snps_subset_unphased.geno.gz -of diplo \
                       --excludeSamples 1238_Deilephila-elpenor,502_Deilephila-porcellus \
                       --threads 20;
                ## get the sliding windows trees (non-overalping, 50kb)
                python2 $WS/src/genomics_general/phylo/raxml_sliding_windows.py \
                      -g $RESULTS/$dataset\_snps_subset_unphased.geno.gz \
                      --prefix $RESULTS/$dataset\_raxml_sliding \
                      --windType coordinate -w 50000 -S 50000 -M 10 -O 0 --model GTRCAT \
                      --raxml $WS/envirs/astral/bin/raxmlHPC \
                      --threads 104;
                conda deactivate;
                gunzip -c $RESULTS/$dataset\_raxml_sliding.trees.gz \
                  | awk '{if($0 != "NA"){print $0}}' \
                  | sed 's/_end_rescale//g' \
                  > $RESULTS/$dataset\_raxml_sliding.trees;
                ## get the gCF
                iqtree -t $RESULTS/SNP_phylo_raxml_$dataset.raxml.support \
                       --gcf $RESULTS/$dataset\_raxml_sliding.trees \
                       -s $RESULTS/0_SNPs_$dataset\_noinv.phy \
                       -T 20 \
                       --prefix $RESULTS/SNP_phylo_raxml_$dataset\_gCF;
                rm $RESULTS/$dataset\_snps_subset_* \
                   $RESULTS/$dataset\_raxml_sliding*;
		echo "DONE!"
		echo "Output: $RESULTS/SNP_phylo_raxml_$dataset.*"
	fi;
done;

echo -e "\n\n\n"
echo -e "------------------------------------------------------------\n\n"
echo "Session Info:"
echo "- ascbias.py retrieved from github.com/btmartin721/raxml_ascbias"
echo "- raxml-ng"
raxml-ng --version
echo "- iqtree"
iqtree --version

