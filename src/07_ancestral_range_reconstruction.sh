  GNU nano 5.6.1                                                                07_ancestral_range_reconstruction.sh                                                                          
#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=66:00:00
#SBATCH --mem=50G
#SBATCH --job-name=hDNA_biogeobears
#SBATCH --account=p_hyleshawkmoths
#SBATCH --output=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_biogeobears-%j.out
#SBATCH --error=/data/horse/ws/vimo762h-hDNAhyles/logs/hDNA_biogeobears-%j.err

# setup environment
module purge
module load release/24.10
module load GCCcore/13.3.0
module load Python/3.12.3
module load Anaconda3

# Variable definition
WS=/data/horse/ws/vimo762h-hDNAhyles
RESULTS=$WS/results/07_ancestral_range

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate $WS/envirs/hDNA_map

# check workspace first
if [ ! -d $WS ]; then
        echo "Error: Could not find the workspace ${WS}" && exit 1;
fi;

if [ ! -d $RESULTS ]; then
        mkdir -p $RESULTS
fi;


Rscript --vanilla $WS/src/07_biogeobears.R
