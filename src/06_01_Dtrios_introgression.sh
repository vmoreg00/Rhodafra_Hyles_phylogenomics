#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=hDNA_Dtrios
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_Dtrios-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_Dtrios-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
VCF=$WS/results/03_variant_calling/variant_calling_full.vcf.gz
RESULTS=$WS/results/06_Introgression/01_Dtrios
Dsuite=/home/vimo762h/hDNAhyles/src/Dsuite/Build/Dsuite

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/introgression

if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS;
fi;


# create samples list. Only Rhodafra and Hyles (Xylophanes=Outgroup)
# Only use the full dataset
if [ ! -e $RESUTLS/Dtrios_list.txt ]; then
	# Get the file list used for Dtrios
	## NOTE: I have remove the sample 4349_Hyles-calida because I suspect
	## that it might be an hybrid between H. calida and H. perkinsi. In
	## runs in which I include that sample, the Fbranch is one for both
	## species, suggesting that they are identical.
	bcftools view -h $VCF \
	  | tail -n 1 \
	  | sed 's/.*FORMAT\t//' \
	  | tr '\t' '\n' \
          | grep 'Xylophanes\|_h\|Hyles' \
	  | grep -v '4349_Hyles-calida' \
          > $RESULTS/Dtrios_list.txt;
        # also make the pop file: sample_id\tgroup
	paste <(cat $RESULTS/Dtrios_list.txt) \
	      <(cat $RESULTS/Dtrios_list.txt \
                 | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
                 | sed 's/01_hDNA_alignment\///g' \
                 | sed 's/02_freshDNA_alignment\///g' \
                 | sed 's/\.bam//g' \
                 | awk -v FS="_" -v OFS="\t" '{print $2}' \
                 | sed 's/Xylophanes-falco\|Xylophanes-tersa\|Xylophanes-porcus/Outgroup/' \
                 | sed 's/Hyles-/H/g' \
                 | sed 's/-sheljuzkoi\|-catissima\|-nicaea//') \
	  > $RESULTS/Dtrios_popFile.txt
fi;

# Extract samples for this analyses
if [ ! -e $RESULTS/snps_subset.vcf.gz ]; then
        echo "Getting the VCF with the species of interest..."
        bcftools view -O z -S $RESULTS/Dtrios_list.txt \
                 $VCF > $RESULTS/snps_subset.vcf.gz;

fi;

###############################################################################
# Introgression Anlysis with Dtrios                                           #
###############################################################################

# The starting tree is based on SNPs phylogenies
# Only for representation purposes
# It does not affect the analyses
if [ ! -e $RESULTS/starting_tree.nwk ]; then
        echo '(Outgroup,((Rmars,Rophe),(Hlineata,((Hannei,Heuphorbiarum),(Hbiguttata,((Hgallii,Hnicaea),(Hvespertilio,(Hlivornicoides,(Hlivornica,(Hcalida,(Hperkinsi,Hwilsoni)))))))))));' \
          > $RESULTS/0_starting_tree_SNP.nwk;
        echo '(Outgroup,((Rophe,Rmars),(Hlineata,(Hlivornicoides,((Hannei,Heuphorbiarum),(Hbiguttata,((Hcalida,(Hperkinsi,Hwilsoni)),((Hnicaea,Hgallii),(Hvespertilio,Hlivornica)))))))));'\
          > $RESULTS/0_starting_tree_MITO.nwk;
fi;

# Run Dtrios and get plots
if [ ! -e $RESULTS/Dtrios_tree.txt ]; then
	for phylo in SNP MITO; do
		echo "Running Dtrios for Hyles and Rhodafra (basis: $phylo)..."
		$Dsuite Dtrios -c --ABBAclustering -o $RESULTS/$phylo\_Dtrios \
		               -t $RESULTS/0_starting_tree_$phylo.nwk \
		               $RESULTS/snps_subset.vcf.gz \
		               $RESULTS/Dtrios_popFile.txt;
		echo "Getting Fbrach summary ..."
		$Dsuite Fbranch -p 0.05 \
		                $RESULTS/0_starting_tree_$phylo.nwk \
		                $RESULTS/$phylo\_Dtrios_tree.txt \
		  > $RESULTS/$phylo\_Dtrios_Fbranch.txt;
                $Dsuite Fbranch -Z --Pb-matrix -p 0.05 \
                                $RESULTS/0_starting_tree_$phylo.nwk \
                                $RESULTS/$phylo\_Dtrios_tree.txt \
                  > $RESULTS/$phylo\_Dtrios_Fbranch_zscore.txt;
		echo "Getting the heatmaps ..."
		# Heatmap with phylogeny information
		python3 $WS/src/Dsuite/utils/dtools.py \
                        $RESULTS/$phylo\_Dtrios_Fbranch.txt \
                        $RESULTS/0_starting_tree_$phylo.nwk \
		        --color-cutoff 0.3 --dpi 333 \
		        -n $RESULTS/$phylo\_fbranch;
		# Heatmaps with no phylogeny information
		# I use the Dmin file in this case because this output file
		# does not rely on topology information
		## get the plot order
		cat $RESULTS/0_starting_tree_$phylo.nwk \
		  | sed 's/(\|;\|)//g' \
		  | tr "," "\n" \
		  | grep -v "Outgroup" \
		  > $RESULTS/$phylo\_plot_order.txt;
		## Substitute the p-value with the clustering_robust p-value (more conservative)
		## Ruby code is positional, so I can just change the order of columns
		cat $RESULTS/$phylo\_Dtrios_tree.txt \
		  | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$9,$7,$8,$6,$10,$11,$12}' \
		  > $RESULTS/$phylo\_Dtrios_tree_p_clust_robust.txt;
		## Make the plots
		ruby $WS/src/plot_d.rb $RESULTS/$phylo\_Dtrios_Dmin_p_clust_robust.txt \
		     $RESULTS/$phylo\_plot_order.txt 0.7 \
		      $RESULTS/$phylo\_Dtrios_D.svg;
		ruby $WS/src/plot_f4ratio.rb $RESULTS/$phylo\_Dtrios_Dmin_p_clust_robust.txt \
		     $RESULTS/$phylo\_plot_order.txt 0.7 \
		     $RESULTS/$phylo\_Dtrios_f4ratio.svg;
	done
fi;

echo "DONE!"
echo "Output: $RESULTS"
