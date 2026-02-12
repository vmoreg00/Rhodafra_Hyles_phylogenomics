#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0:20:00
#SBATCH --mem=10GB
#SBATCH --job-name=hDNA_HyDe
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_HyDe-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_HyDe-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
MSA=$WS/results/03_variant_calling/variant_calling_full_1kb-thin_consensus.phy
RESULTS=$WS/results/06_Introgression/05_HyDe

# environment loading
if [ ! -d $WS/envirs/py_venv_HyDe ]; then
        mkdir $WS/envirs
        python3 -m venv --system-site-package $WS/envirs/py_venv_HyDe
	source $WS/envirs/py_venv_HyDe/bin/activate
	pip install cython multiprocess numpy
	pip install phyde
else
	source $WS/envirs/py_venv_HyDe/bin/activate
fi;

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;
if [ ! -d $RESULTS ]; then
	mkdir -p $RESULTS;
fi;

# Create input files
## Sequence file
if [ ! -e $RESULTS/0_sequences.txt ]; then
	tail -n +2 $MSA \
	  | grep -v "Chaerocina\|Euchloron\|Theretra\|Hippotion" \
	  > $RESULTS/0_sequences.txt
fi;
## Map file
if [ ! -e $RESULTS/0_map.txt ]; then
	paste <(cut -f 1 -d " " $RESULTS/0_sequences.txt) \
	      <(cut -f 1 -d " " $RESULTS/0_sequences.txt \
	          | cut -f 2 -d "_" \
	          | sed 's/-catissima\|-nicaea\|-sheljuzkoi//2' \
	          | sed 's/Hyles-/H/' \
	          | sed 's/-falco\|-tersa\|-porcus//') \
	  > $RESULTS/0_map.txt
fi;

# Run HyDe ====================================================================
if [ ! -e $RESULTS/hyde-out.txt ]; then
	run_hyde_mp.py -i $RESULTS/0_sequences.txt -m $RESULTS/0_map.txt \
	               -o Xylophanes -n 46 -t 15 -s 108574 -j 20 --pvalue 0.001 \
	               --prefix $RESULTS/hyde;
fi;

echo "DONE!"
echo "Output: $RESULTS"
