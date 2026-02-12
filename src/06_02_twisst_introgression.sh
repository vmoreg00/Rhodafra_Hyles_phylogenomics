#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --time=48:00:00
#SBATCH --account=p_hyleshawkmoths
#SBATCH --job-name=hDNA_introgresion
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_introgresion-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_introgresion-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
VCF=$WS/results/03_variant_calling/variant_calling_full.vcf.gz
RESULTS=$WS/results/06_Introgression/02_twisst

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/introgression
export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general

if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS;
fi;

# create samples list for each scenario
if [ ! -e $RESUTLS/full_list.txt ]; then
	bcftools view -h $VCF \
	  | tail -n 1 \
	  | sed 's/.*FORMAT\t//' \
	  | tr '\t' '\n' \
	  > $RESULTS/full_list.txt;
fi;
if [ ! -e $RESUTLS/America_list.txt ]; then
	# Scenario one: Hyles species are able to hybridize in South America
	grep 'Xylophanes\|lineata\|annei\|3623_\|3722\|4379' \
	     $RESULTS/full_list.txt > $RESULTS/America_list.txt;
	# also make the pop file: sample_id\tgroup
	cat $RESULTS/America_list.txt \
	  | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
          | sed 's/01_hDNA_alignment\///g' \
          | sed 's/02_freshDNA_alignment\///g' \
          | sed 's/\.bam//g' \
	  | awk -v FS="_" -v OFS="\t" '{print $0"_A",$2"\n"$0"_B",$2}' \
	  | sed 's/euphorbiae/euphorbiarum/2' \
	  | sed 's/-falco\|-tersa\|-porcus//2' \
	  > $RESULTS/America_popFile.txt
fi;
if [ ! -e $RESUTLS/Africa_list.txt ]; then
	# Scenario two: Hyles and Rhodafra are able to hybridize in Africa
        grep 'Xylophanes\|Rophe_h\|Rmars_h\|biguttata\|livornica' \
             $RESULTS/full_list.txt > $RESULTS/Africa_list.txt;
	# also make the pop file: sample_id\tgroup
	cat $RESULTS/Africa_list.txt \
          | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
          | sed 's/01_hDNA_alignment\///g' \
          | sed 's/02_freshDNA_alignment\///g' \
          | sed 's/\.bam//g' \
          | awk -v FS="_" -v OFS="\t" '{print $0"_A",$2"\n"$0"_B",$2}' \
          | sed 's/-falco\|-tersa\|-porcus//2' \
          > $RESULTS/Africa_popFile.txt
fi;
if [ ! -e $RESUTLS/Hawaii_list.txt ]; then
	# Scenario three: Hyles species are able to hybridize in Hawaii
        grep 'Xylophanes\|lineata\|calida\|perkinsi\|wilsoni' \
             $RESULTS/full_list.txt > $RESULTS/Hawaii_list.txt;
	# also make the pop file: sample_id\tgroup
	cat $RESULTS/Hawaii_list.txt \
          | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
          | sed 's/01_hDNA_alignment\///g' \
          | sed 's/02_freshDNA_alignment\///g' \
          | sed 's/\.bam//g' \
          | awk -v FS="_" -v OFS="\t" '{print $0"_A",$2"\n"$0"_B",$2}' \
          | sed 's/-falco\|-tersa\|-porcus//2' \
          > $RESULTS/Hawaii_popFile.txt
fi;
if [ ! -e $RESUTLS/AncientILS_list.txt ]; then
        # Scenario four: There is an evident ILS between Rhodafra and H. lineata
	# but no signs of migration between Rhodafra and african Hyles species.
	# (H.gallii represent the rest of Hyles species; I take this species
	#  because its placement is very constant, and it is a derived group)
        grep 'Xylophanes\|Rophe_h\|Rmars_h\|lineata\|gallii' \
             $RESULTS/full_list.txt > $RESULTS/AncientILS_list.txt;
        # also make the pop file: sample_id\tgroup
        cat $RESULTS/AncientILS_list.txt \
          | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
          | sed 's/01_hDNA_alignment\///g' \
          | sed 's/02_freshDNA_alignment\///g' \
          | sed 's/\.bam//g' \
          | awk -v FS="_" -v OFS="\t" '{print $0"_A",$2"\n"$0"_B",$2}' \
          | sed 's/-falco\|-tersa\|-porcus//2' \
          > $RESULTS/AncientILS_popFile.txt
fi;

###############################################################################
# Introgression Analyses with TWISST ==========================================
###############################################################################
for scenario in America Africa Hawaii AncientILS; do
	echo "==============================================="
	echo "Testing introgression in species from $scenario"
	echo "==============================================="
	# define the groups
	if [ $scenario == "America" ]; then
		g1="Hyles-lineata";g2="Hyles-euphorbiarum";g3="Hyles-annei"
	elif [ $scenario == "Africa" ]; then
		g1="Rophe";g2="Rmars";g3="Hyles-biguttata";g4="Hyles-livornica"
	elif [ $scenario == "Hawaii" ]; then
		g1="Hyles-lineata";g2="Hyles-calida";g3="Hyles-perkinsi";g4="Hyles-wilsoni"
	else
		g1="Rophe";g2="Rmars";g3="Hyles-lineata";g4="Hyles-gallii"
	fi;
	# Create the directory
	if [ ! -d $RESULTS/$scenario ]; then
		mkdir -p $RESULTS/$scenario;
	fi;
	# Extract samples for this scenario
	if [ ! -e $RESULTS/$scenario/snps_subset.vcf.gz ]; then
		echo "Getting the VCF with the species from $scenario..."
		bcftools view -O z -S $RESULTS/$scenario\_list.txt \
			 $VCF > $RESULTS/$scenario/snps_subset.vcf.gz;

	fi;
	# retrieving geno.gz file (change the headers for simplicity)
	if [ ! -e $RESULTS/$scenario/snps_subset.geno.gz ]; then
		echo "Getting the geno.gz file..."
		python $WS/src/genomics_general/VCF_processing/parseVCF.py \
		       --skipIndels \
		       -i $RESULTS/$scenario/snps_subset.vcf.gz \
		  | sed 's/\/data\/horse\/ws\/vimo762h-hDNAhyles\/results\///g' \
		  | sed 's/01_hDNA_alignment\///g' \
		  | sed 's/02_freshDNA_alignment\///g' \
	          | sed 's/\.bam//g' \
                  | awk '{if(NR==1){print $0} else {print toupper($0)}}' \
                  | grep -Ev "([ACGT]){2,}" \
	          | bgzip \
	          > $RESULTS/$scenario/snps_subset.geno.gz;
	fi;

	# get the geno file (with genomics general)
	echo "ITERATING for different window sizes..."
	for ws in 20 50 100 500; do
		echo "==> Window size: $ws bp"
		if [ ! -d $RESULTS/$scenario/ws_$ws ]; then
			mkdir -p $RESULTS/$scenario/ws_$ws
			# sliding window phylogeny (genomics general)
			echo "    - Inferring tree..."
			conda activate $WS/envirs/py2.7
			export PYTHONPATH=$PYTHONPATH:$WS/src/genomics_general
		        python $WS/src/genomics_general/phylo/raxml_sliding_windows.py \
			      --windType sites -w $ws --minPerInd 4 \
		              -g $RESULTS/$scenario/snps_subset.geno.gz \
		              --prefix $RESULTS/$scenario/ws_$ws/01_raxml_sw_phylo \
			      --raxml $WS/envirs/astral/bin/raxmlHPC \
		              --model GTRCAT \
			      --threads 104;
			conda deactivate

			# twisst
			echo "Weighting..."
			if [ $scenario == "America" ]; then
				python $WS/src/twisst/twisst.py \
				       -t $RESULTS/$scenario/ws_$ws/01_raxml_sw_phylo.trees.gz \
			               -w $RESULTS/$scenario/ws_$ws/02_raxml_sw_phylo.weights.tsv.gz \
				       -g Xylophanes -g $g1 -g $g2 -g $g3 \
				       --groupsFile $RESULTS/$scenario\_popFile.txt \
			               --outgroup Xylophanes;
			else
			        python $WS/src/twisst/twisst.py \
                                       -t $RESULTS/$scenario/ws_$ws/01_raxml_sw_phylo.trees.gz \
                                       -w $RESULTS/$scenario/ws_$ws/02_raxml_sw_phylo.weights.tsv.gz \
                                       -g Xylophanes -g $g1 -g $g2 -g $g3 -g $g4 \
                                       --groupsFile $RESULTS/$scenario\_popFile.txt \
                                       --outgroup Xylophanes;
			fi;

			# plot
			echo "Plotting..."
			Rscript --vanilla $WS/src/06_02_plot_twisst.R $RESULTS/$scenario/ws_$ws;
		fi;
	done;
done;

echo "DONE!"
echo "Output: $RESULTS"
