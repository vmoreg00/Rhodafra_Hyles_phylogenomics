#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --time=168:00:00
#SBATCH --mem=250GB
#SBATCH --account=p_hyleshawkmoths
#SBATCH --job-name=hDNA_phylo_SNAPPER
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_SNAPPER-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_phylo_SNAPPER-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
MSA=$WS/results/03_variant_calling
RESULTS=$WS/results/04_SNP_phylogeny/04_SNAPPER
threads=$((SLURM_CPUS_PER_TASK*2))

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

if [ ! -d $RESULTS ]; then
        mkdir -p $RESULTS
fi;

# SNAPP only creates output files in the working directory. Go there!!
cd $RESULTS

# Create the constraint and divergence times file -----------------------------
if [ ! -e $RESULTS/0_constraints.tsv ]; then
	## 1-Root age (normal) : mean=10; SD=1
	## Took from Kawahara et al 2015 [10.1073/pnas.1416679112]
        echo -e \
          'normal(0,10,1)\tcrown\tXylophanes,Rmars,Rophe,Hannei,Hbiguttata,Hcalida,Heuphorbiarum,Hgallii,Hlineata,Hlivornica,Hlivornicoides,Hnicaea,Hperkinsi,Hvespertilio,Hwilsoni' \
          > $RESULTS/0_constraints.tsv;
	## 2-Hawaiian species age (normal) : mean=5.1; SD=1
	## Took from Hundsdoerfer et al 2017 [10.1111/zsc.12235]
        echo -e \
          'normal(0,3.14,1)\tcrown\tHwilsoni,Hcalida,Hperkinsi' \
          >> $RESULTS/0_constraints.tsv;
	## 3-calida-perkinsi age (normal) : mean=3.7; SD=0.612
	## Took from Hundsdoerfer et al 2017 [10.1111/zsc.12235]
        #echo -e \
        #  'normal(0,2.5,0.612)\tcrown\tHcalida,Hperkinsi' \
        #  >> $RESULTS/0_constraints.tsv;
fi;
# Get starting tree -----------------------------------------------------------
if [ ! -e $RESULTS/0_starting_tree.nwk ]; then
	# (taken from previous SNPs phylogenies of the present study)
	echo '(Xylophanes,((Rmars,Rophe),(Hlineata,((Hannei,Heuphorbiarum),(Hbiguttata,((Hgallii,Hnicaea),(Hvespertilio,(Hlivornicoides,(Hlivornica,(Hcalida,(Hperkinsi,Hwilsoni)))))))))));' \
	  > $RESULTS/0_starting_tree.nwk;
fi;
# I use SNAPPER because is more efficient
# I use the 1samp-x-sp dataset to (1) speed-up the process, (2) avoid
# missleading results due to missing data.
# In the same way, I use the strict dataset.
for dataset in 1samp-x-sp strict; do
	# get consensus in phylip format --------------------------------------
	if [ ! -e $RESULTS/$dataset\_good.phy ]; then
		# Use the 1kb thin dataset to avoid physical linkeage between loci.
		## get the autosomes
		regions="OX457094.1"
		for i in {096..122}; do
		        regions="$regions,OX457$i.1"
		done
		echo "Getting the phylip consensus..."
		## Subset vcf and get phy file
		bcftools view -O z -r $regions -o $RESULTS/$dataset\_1kb-thin.vcf.gz \
		              $MSA/variant_calling_$dataset\_1kb-thin.vcf.gz \
	          && python $WS/src/vcf2phylip.py -i $RESULTS/$dataset\_1kb-thin.vcf.gz \
	                                          --output-folder $RESULTS \
	                                          --output-prefix $RESULTS/$dataset \
		  && rm $RESULTS/$dataset\_1kb-thin.vcf.gz &&
		  grep -v 'Chaerocina\|Deilephila\|Hippotion\|Euchloron\|Theretra\|porcus\|tersa' \
	               $RESULTS/$dataset.min4.phy \
	            | cut -d "/" -f 8 \
	            | sed 's/.bam//g' \
	            | sed 's/_h_end_rescale//g' \
	            | sed 's/_/./g' \
	            | sed 's/Xylophanes-/X/g' \
	            | sed 's/Hyles-/H/g' \
	            | sed 's/-catissima//g' \
		    | sed 's/26/15/' \
	            > $RESULTS/$dataset\_good.phy &&
		  rm $RESULTS/$dataset.min4.phy;
	fi;
	# Prepare SNAPPER files -----------------------------------------------
	echo "Preparing SNAPPER files..."
	# Get the samples-species file
	if [ ! -e $RESULTS/samples_$dataset.tsv ]; then
	        echo -e "species\tsample" > $RESULTS/samples_$dataset.tsv;
	        paste \
	          <(cut -f 1 -d " " $RESULTS/$dataset\_good.phy \
	            | grep '.H\|.R\|.X' \
		    | sed 's/Xfalco/Xylophanes/g' \
	            | cut -d "." -f 2) \
	          <(cut -f 1 -d " " $RESULTS/$dataset\_good.phy \
	            | grep '.H\|.R\|.X') \
	          >> $RESULTS/samples_$dataset.tsv;
	fi;
        # Create the xml file for SNAPP
	## default chain length of 500000 is not enough to get estimates stability.
	## I use 1.000.000 instead.
        if [ ! -e $RESULTS/$dataset\_snapper.xml ]; then
                ruby $WS/src/snapp_prep.rb -p $RESULTS/$dataset\_good.phy \
                                           -t $RESULTS/samples_$dataset.tsv \
                                           -c $RESULTS/0_constraints.tsv \
                                           -l 1000000 \
                                           -s $RESULTS/0_starting_tree.nwk \
                                           -x $RESULTS/$dataset\_snapper.xml \
		                           -a SNAPPER \
                                           -o $RESULTS/$dataset\_snapper;
        fi;
        # Run beast/snapper ---------------------------------------------------
	echo "Running SNAPPER..."
        if [ ! -e $RESULTS/$dataset\_snapper.log ]; then
                # SNAPPER
                $WS/src/beast/bin/beast -threads $threads -seed 123 $RESULTS/$dataset\_snapper.xml;
        fi;

	# Get the final tree --------------------------------------------------
	if [ ! -e $RESULTS/$dataset\_snapper_final.tree ]; then
		$WS/src/beast/bin/treeannotator -burnin 10 -height mean \
		                                -file $RESULTS/$dataset\_snapper.trees \
		                                $RESULTS/$dataset\_snapper_final.tree;
	fi;

        # Run this if you want to check the statisitcs
        # tracer $RESULTS/$dataset\_snapper.log

	# Run this if you want to see the densitree
	# $WS/src/beast/bin/densitree $RESULTS/$dataset\_snapper.trees
done;
