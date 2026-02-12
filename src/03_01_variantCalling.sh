#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --job-name=hDNA_variant_calling
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_variant_calling-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_variant_calling-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
HDNA=$WS/results/01_hDNA_alignment
FRESH=$WS/results/02_freshDNA_alignment
REFGENOME=$WS/data/refgenome
RESULTS=$WS/results/03_variant_calling
hdna_files=$(ls $HDNA/*_rescale.bam)
fresh_files=$(ls $FRESH/*.bam)

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS;
	mkdir $RESULTS/QC;
fi;

# create the bam lists
if [ ! -e $RESULTS/hDNA_bamlist ]; then
	ls -1 $HDNA/*_rescale.bam > $RESULTS/hDNA_bamlist;
fi;
if [ ! -e $RESULTS/freshDNA_bamlist ]; then
        ls -1 $FRESH/*.bam > $RESULTS/freshDNA_bamlist;
fi;
if [ ! -e $RESULTS/full_bamlist ]; then
        cat $RESULTS/hDNA_bamlist $RESULTS/freshDNA_bamlist > $RESULTS/full_bamlist;
fi;
# Get the regions
# Only autosomes will be analyzed for simplicity.
# Deilephila elpernor chromosomes are OX457093 to OX457123, but
# OX457093, OX457095 and OX457123 correspond to Z, W and MT chromosomes, respectively.
regions="OX457094.1"
for i in {096..122}; do
	regions="$regions OX457$i.1"
done
echo "Running variant calling in chromosomes:"
for chr in $regions; do
	if [ ! -e $RESULTS/$chr.vcf ]; then
		echo "  - $chr ..."
		bcftools mpileup -Ov \
				 -b $RESULTS/full_bamlist \
				 -r $chr \
				 -f $REFGENOME/DeiElpe2.1_refgenome.fna \
		  | bcftools call -a GQ,GP -mOv --ploidy 2 \
		  		  -o $RESULTS/$chr.vcf &
	fi;
done
echo "Waiting for variant calling to be done..."
wait;
echo "DONE!"

# Concat, normalize and filter
# - variant sites within 5 nt from an indel will be removed
# - I only keep variants with coverage DP>120, quality QUAL>15,
#   missing ratio F_MISSING<0.2. Uncalled genotypes are excluded
if [ ! -e $RESULTS/variant_calling_full.vcf -a ! -e $RESULTS/variant_calling_full.vcf.gz ]; then
	echo "[All samples] Concatenating, sorting, normalizing, filtering..."
        mkdir $RESULTS/tmp;
        bcftools concat -O v $(ls $RESULTS/*1.vcf) \
	  | bcftools sort -m 75G -O v --temp-dir $RESULTS/tmp \
	  | bcftools norm -f $REFGENOME/DeiElpe2.1_refgenome.fna \
			  -m -any -O v \
          | bcftools filter -O v -g 5:indel,other \
          		    -i 'DP>120 && MAF>0.1 && QUAL>15 && F_MISSING<0.2 && GT!~"\."' \
                            -o $RESULTS/variant_calling_full.vcf;
	echo "[All samples] Thinning (1000bp)..."
	bcftools +prune -w 1000bp -n 1 -N maxAF -O v \
	  		-o $RESULTS/variant_calling_full_1kb-thin.vcf \
			$RESULTS/variant_calling_full.vcf;
        rm -r $RESULTS/tmp;
	echo "[All samples] Retrieving the per-sample missing rate";
        samps=($(bcftools view -h -O v $RESULTS/variant_calling_full.vcf \
                  | tail -n 1 | cut -f 10-))
        for i in $(seq 0 $((${#samps[@]} - 1))); do
                echo ${samps[$i]}
                printf "${samps[$i]}\t" >> $RESULTS/QC/missing_rate_per_sample.tsv;
                printf "$(echo ${samps[$i]} | cut -d '_' -f 4 | cut -d '-' -f 1-2 | sed 's/.bam//g')\t" \
                  >> $RESULTS/QC/missing_rate_per_sample.tsv
                bcftools view -H -O v $RESULTS/variant_calling_full.vcf \
                  | cut -f $(($i+11)) \
                  | cut -d ':' -f 1 \
                  | awk '/^\.\/\./ {NC++;} END{printf("%f\n",NC/(1.0*NR))}' \
                  >> $RESULTS/QC/missing_rate_per_sample.tsv;
        done;
	# get the best samples per species (those with the lowest missign rate)
        awk '{print $2,$3,$1}' $RESULTS/QC/missing_rate_per_sample.tsv \
          | sort \
          | awk -v prev="" '{if ($1 != prev) {print $3}; prev=$1}' \
          > $RESULTS/00_best_samples_per_species.txt
fi;

# Repeat the same process for the subsets =====================================
# 1. Full dataset without 508 sample (very low number of reads) - to check nb of SNPs
# 2. One sample per species (select only those with the lowest missing rate)
# 3. One sample per species removing Rhodafra (hDNA sampes) - to check effects on topology
if [ ! -e $RESULTS/variant_calling_full_no508.vcf -a ! -e $RESULTS/variant_calling_full_no508.vcf.gz ]; then
        echo "[Sample 508 (H. biguttata) removed] Concatenating, sorting, normalizing, filtering..."
        mkdir $RESULTS/tmp;
        bcftools concat -O v $(ls $RESULTS/*1.vcf) \
	  | bcftools view -s "^/data/horse/ws/vimo762h-hDNAhyles/results/02_freshDNA_alignment/508_Hyles-biguttata.bam" -O v \
          | bcftools sort -m 75G -O v --temp-dir $RESULTS/tmp \
          | bcftools norm -f $REFGENOME/DeiElpe2.1_refgenome.fna \
                          -m -any -O v \
          | bcftools filter -O v -g 5:indel,other \
                            -i 'DP>120 && MAF>0.1 && QUAL>15 && F_MISSING<0.2 && GT!~"\."' \
                            -o $RESULTS/variant_calling_full_no508.vcf;
        echo "[Sample 508 (H. biguttata) removed] Thinning (1000bp)..."
        bcftools +prune -w 1000bp -n 1 -N maxAF -O v \
                        -o $RESULTS/variant_calling_full_no508_1kb-thin.vcf \
                        $RESULTS/variant_calling_full_no508.vcf;
	rm -r $RESULTS/tmp
fi;
if [ ! -e $RESULTS/variant_calling_1samp-x-sp.vcf -a ! -e $RESULTS/variant_calling_1samp-x-sp.vcf.gz ]; then
        echo "[One sample per species] Concatenating, sorting, normalizing, filtering..."
        mkdir $RESULTS/tmp;
        bcftools concat -O v $(ls $RESULTS/*1.vcf) \
          | bcftools view -S $RESULTS/00_best_samples_per_species.txt -O v \
          | bcftools sort -m 75G -O v --temp-dir $RESULTS/tmp2 \
          | bcftools norm -f $REFGENOME/DeiElpe2.1_refgenome.fna \
                          -m -any -O v \
          | bcftools filter -O v -g 5:indel,other \
                            -i 'DP>120 && MAF>0.1 && QUAL>15 && F_MISSING<0.2 && GT!~"\."' \
                            -o $RESULTS/variant_calling_1samp-x-sp.vcf;
        echo "[One sample per species] Thinning (1000bp)..."
        bcftools +prune -w 1000bp -n 1 -N maxAF -O v \
                        -o $RESULTS/variant_calling_1samp-x-sp_1kb-thin.vcf \
                        $RESULTS/variant_calling_1samp-x-sp.vcf;
	rm -r $RESULTS/tmp
        echo "DONE!"
fi;
if [ ! -e $RESULTS/variant_calling_1samp-x-sp_noRhodafra.vcf -a ! -e $RESULTS/variant_calling_1samp-x-sp_noRhodafra.vcf.gz ]; then
        echo "[Rhodafra removed] Concatenating, sorting, normalizing, filtering..."
        mkdir $RESULTS/tmp;
        grep -v "_h" $RESULTS/00_best_samples_per_species.txt \
          > excludeRhodafra.txt
        bcftools concat -O v $(ls $RESULTS/*1.vcf) \
          | bcftools view -S excludeRhodafra.txt -O v \
          | bcftools sort -O v --temp-dir $RESULTS/tmp \
          | bcftools norm -f $REFGENOME/DeiElpe2.1_refgenome.fna \
                          -m -any -O v \
          | bcftools filter -O v -g 5:indel,other \
                            -i 'DP>120 && MAF>0.1 && QUAL>15 && F_MISSING<0.2 && GT!~"\."' \
                            -o $RESULTS/variant_calling_1samp-x-sp_noRhodafra.vcf;
        echo "[Rhodafra removed] Thinning (1000bp)..."
        bcftools +prune -w 1000bp -n 1 -N maxAF -O v \
                        -o $RESULTS/variant_calling_1samp-x-sp_noRhodafra_1kb-thin.vcf \
                        $RESULTS/variant_calling_1samp-x-sp_noRhodafra.vcf;
        rm -r $RESULTS/tmp excludeRhodafra.txt;
fi;
if [ ! -e $RESULTS/variant_calling_strict.vcf -a ! -e $RESULTS/variant_calling_strict.vcf.gz ]; then
	# As recommended by the Hawaiian team, stricter filtering rules are applied
	## - SNPs are removed if closed than 20nt to an indel
	## - DP is restricted between 60 (~half mean DP) and 350 (~3.5x mean DP) to avoid
	##   SNPs due to repetetive regions (missmapping)
	## - MAF is kept at 0.1
	## - QUAL threshold is increased to 30
	## - F_MISSING is decreased to 0.05
        echo "[Strict filtering] Concatenating, sorting, normalizing, filtering..."
	mkdir $RESULTS/tmp;
        bcftools concat -O v $(ls $RESULTS/*1.vcf) \
          | bcftools view -S $RESULTS/00_best_samples_per_species.txt -O v \
          | bcftools sort -m 75G -O v --temp-dir $RESULTS/tmp2 \
          | bcftools norm -f $REFGENOME/DeiElpe2.1_refgenome.fna \
                          -m -any -O v \
          | bcftools filter -O v -g 20:indel,other \
                            -i 'DP>60 && DP<=350 && MAF>0.1 && QUAL>30 && F_MISSING<0.05 && GT!~"\."' \
                            -o $RESULTS/variant_calling_strict.vcf;
        echo "[Strict removed] Thinning (1000bp)..."
        bcftools +prune -w 1000bp -n 1 -N maxAF -O v \
                        -o $RESULTS/variant_calling_strict_1kb-thin.vcf \
                        $RESULTS/variant_calling_strict.vcf;
        rm -r $RESULTS/tmp;
	echo "DONE!"
fi;
# Quality Control
if [ ! -e $RESULTS/QC/LDdecay.png ]; then
	$WS/src/PopLDdecay/bin/PopLDdecay -InVCF $RESULTS/variant_calling_full.vcf.gz \
	                                  -OutStat $RESULTS/QC/LDdecay;
	perl $WS/src/PopLDdecay/bin/Plot_OnePop.pl -inFile $RESULTS/QC/LDdecay.stat.gz \
	                                           -output $RESULTS/QC/LDdecay;
	# The LD stabilizes at ~0.1425. This asymptote is reached at ~50kb
fi;

echo "Getting summary statistics..."
for dataset in full full_no508 1samp-x-sp 1samp-x-sp_noRhodafra strict; do
	if [ ! -e $RESULTS/QC/variants_$dataset.stats ]; then
	        bcftools stats -F $REFGENOME/DeiElpe2.1_refgenome.fna \
	                       $RESULTS/variant_calling_$dataset.vcf \
	          > $RESULTS/QC/variants_$dataset.stats &
	fi;
done;
wait;
echo "DONE!"

echo "Summarizing statistics..."
for dataset in full full_no508 1samp-x-sp 1samp-x-sp_noRhodafra strict; do
	if [ ! -e $RESULTS/QC/QC_summary_$dataset.pdf ]; then
	        plot-vcfstats --prefix $RESULTS/tmp/ \
	                      --title Variant_calling_$dataset \
	                      --main-title "Variant Calling QC ($dataset)" \
	                      $RESULTS/QC/variants_$dataset.stats;
	        mv $RESULTS/tmp/summary.pdf $RESULTS/QC/QC_summary_$dataset.pdf;
	        rm -r $RESULTS/tmp;
	fi;
done
echo "DONE!"

# Compress to VCF.GZ format
echo "Converting into vcf.gz format..."
for dataset in full full_no508 1samp-x-sp 1samp-x-sp_noRhodafra strict; do
	if [ ! -e $RESULTS/variant_calling_$dataset.vcf.gz ]; then
	        bcftools view -O z -o $RESULTS/variant_calling_$dataset.vcf.gz \
	                      --threads 50 \
	                      $RESULTS/variant_calling_$dataset.vcf \
		  && rm $RESULTS/variant_calling_$dataset.vcf;
		bcftools view -O z -o $RESULTS/variant_calling_$dataset\_1kb-thin.vcf.gz \
	                      --threads 50 \
	                      $RESULTS/variant_calling_$dataset\_1kb-thin.vcf \
		  && rm $RESULTS/variant_calling_$dataset\_1kb-thin.vcf;
	fi;
done;
echo "DONE!"

# Index
echo "Indexing the VCF.GZ files..."
for dataset in full full_no508 1samp-x-sp 1samp-x-sp_noRhodafra strict; do
	if [ ! -e $RESULTS/variant_calling_$dataset.vcf.gz.csi ]; then
	        bcftools index -o $RESULTS/variant_calling_$dataset.vcf.gz.csi \
			       $RESULTS/variant_calling_$dataset.vcf.gz;
	        bcftools index -o $RESULTS/variant_calling_$dataset\_1kb-thin.vcf.gz.csi \
			       $RESULTS/variant_calling_$dataset\_1kb-thin.vcf.gz;
	fi;
done;
echo "DONE!"

conda deactivate


## NOTES ######################################################################
#                                                                             #
# Number of SNPs detected before filtering was very similar in all three      #
# datasets (all samples, no 508 sample, no Rhodafra samples).                 #
# Full dataset         2019598                                                #
# no 508 sample        2000321                                                #
# no Rhodafra samples  2023499                                                #
#                                                                             #
# Similar results were obtained after thining the variants                    #
# Full dataset         116008                                                 #
# no 508 sample        113994                                                 #
# no Rhodafra samples  113281                                                 #
#                                                                             #
# As these results are very similar, I will omit the "no508" dataset in the   #
# ongoing analyses. I will keep the "withoutRhodafra" dataset to check the    #
# effect on the topology of the phylogenetic analyses                         #
###############################################################################
