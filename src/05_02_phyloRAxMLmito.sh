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

		echo "DONE!"
		echo "Output: $RESULTS/mito_raxml_$dataset.*"
	fi;
done;
