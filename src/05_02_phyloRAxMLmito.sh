#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH --account=p_hyleshawkmoths
#SBATCH --job-name=hDNA_mito_RAxML
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mito_RAxML-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mito_RAxML-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
DATA=$WS/results/05_MITO_phylogeny/00_consensus_fasta
MSA=$WS/results/05_MITO_phylogeny/00_alignments
RESULTS=$WS/results/05_MITO_phylogeny/02_RAxML

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS
fi;
if [ ! -d $MSA ]; then
	mkdir -p $MSA
fi;


if [ ! -e $MSA/aln_mafft_mitogenomes_full.fa \
     -o ! -e $MSA/aln_mafft_mitogenomes_1samp-x-sp.fa \
     -o ! -e $MSA/aln_mafft_mitogenomes_1samp-x-sp_noRhodafra.fa ]; then
	echo "Aligning mitogenomes..."
        mafft --thread 32 $DATA/mitogenomes_full.fa \
          > $MSA/aln_mafft_mitogenomes_full.fa &
        mafft --thread 32 $DATA/mitogenomes_1samp-x-sp.fa \
          > $MSA/aln_mafft_mitogenomes_1samp-x-sp.fa &
        mafft --thread 32 $DATA/mitogenomes_1samp-x-sp_noRhodafra.fa \
          > $MSA/aln_mafft_mitogenomes_1samp-x-sp_noRhodafra.fa &
	wait
	echo "DONE!"
	echo "Output: $RESULTS/aln_mafft_mitogenomes_*.fa"
fi;


for dataset in full 1samp-x-sp 1samp-x-sp_noRhodafra; do
	if [ ! -e $RESULTS/mito_raxml_$dataset.nwk ]; then
		echo "Running RAxML-ng for dataset $dataset ..."
		raxml-ng --msa $MSA/aln_mafft_mitogenomes_$dataset.fa \
		         --data-type DNA \
			 --all \
			 --model GTR+I+G4 \
			 --tree pars{10} \
			 --bs-trees 1000 \
			 --threads 30 \
			 --outgroup 6381_Theretra-silhetensis,6386_Theretra-alecto,6389_Theretra-japonica \
			 --prefix $RESULTS/mito_raxml_$dataset;

                echo "Getting concordance factors (with IQtree)..."
                iqtree -t $RESULTS/mito_raxml_$dataset.raxml.support \
                       --scf 100 -T 20 \
                       -s $MSA/aln_mafft_mitogenomes_$dataset.fa ;
                echo "Getting gene-concordance factores (with IQtree)..."
                ## Get genes from each sample's fasta file
                mkdir $RESULTS/mito_tmp/
                for fas in $(ls $DATA\/*.*fas); do
                        fas=$(basename "$fas")
                        seqkit subseq --bed $BED $DATA/$fas \
                                      -o $RESULTS/mito_tmp/$fas;
                        seqkit split --quiet --by-id $RESULTS/mito_tmp/$fas;
                        rm $RESULTS/mito_tmp/$fas
                done;
                for gene in $(awk '{print $3}' $BED); do
                        for fas in $(ls $RESULTS/mito_tmp/*.fas.split/*-$gene\_*_*.fas); do
                                sp=$(basename "${fas%.part_*}");
                                gene_i=$(basename "${fas#*part_}");
                                if [[ $dataset =~ 1samp-x-sp_noRhodafra && $sp =~ [0-9]*_R.* ]]; then continue; fi;
                                if [[ $dataset =~ 1samp-x-sp* ]] && ! grep -qF $sp $WS/results/03_variant_calling/00_best_samples_per_species.txt; then continue; fi;
                                cat $fas \
                                  | sed "s/^>OX.*/>$sp/g" \
                                  >> $RESULTS/mito_tmp/$gene_i;
                        done;
                done;
                ## Run RAxML for each gene
                for fas in $(ls $RESULTS/mito_tmp/*.fas); do
                        fas=$(basename $fas);
                        mafft --thread 5 --quiet $RESULTS/mito_tmp/$fas \
                          > $RESULTS/mito_tmp/ALN_$fas \
                          && raxml-ng --msa $RESULTS/mito_tmp/ALN_$fas \
                                      --data-type DNA --search \
                                      --model GTR+I+G4 --tree pars{10} \
                                      --threads 5 --log ERROR \
                                      --prefix $RESULTS/mito_tmp/TREE_$fas \
                          && cat $RESULTS/mito_tmp/TREE_$fas.raxml.bestTree \
                               >> $RESULTS/mito_raxml_$dataset\_geneTrees.nwk &
                done
                wait
                rm -r $RESULTS/mito_tmp/
                ## calcluate gCF
                iqtree -t $RESULTS/mito_raxml_$dataset.raxml.support \
                       --gcf $RESULTS/mito_raxml_$dataset\_geneTrees.nwk \
                       -s $MSA/aln_mafft_mitogenomes_$dataset.fa \
                       -T 20 \
                       --prefix $RESULTS/mito_raxml_$dataset\_gCF;

		echo "DONE!"
		echo "Output: $RESULTS/mito_raxml_$dataset.*"
	fi;
done;
