#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:00
#SBATCH --mem=12G
#SBATCH --job-name=hDNA_mapping_index
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_index-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_mapping_index-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
REFGENOME=$WS/data/refgenome

# check workspace first
if [ ! -d $WS ]; then
	echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# Index the reference genome ==================================================
if [ ! -e $REFGENOME/DeiElpe2.1_refgenome_BT2index.1.bt2 ]; then
        echo "Indexing reference genome ..."
        bowtie2-build $REFGENOME/DeiElpe2.1_refgenome.fna \
                      $REFGENOME/DeiElpe2.1_refgenome_BT2index \
                      --threads 54;
fi;

conda deactivate
