# Disentangling complex evolutionary histories with hDNA
## Case study with *Hyles* and *Rhodafra* (Lep.: Sphingidae)

This repository contains the code used for the analyses that were carried
out for the manuscript entitled **Disentangling complex evolutionary histories
with historical DNA. A case study on the phylogenomics of *Hyles* Hübner and
*Rhodafra* Rothschild & Jordan (Lepidoptera: Sphingidae)**

  * **Authors**: Moreno-González, V., Melichar, T., San Jose, M., Rubinoff, D.,
    Hundsdoerfer, A.K.
  * **Publication status**: in review

## Folder structure

```
./
 |
 +-- data/     # FASTQ files [not included; data can be donwladed from ENA]
 +-- envirs/   # conda and python environments.
 +-- logs/     # SLURM log files [not included; to be filled when the code is run]
 +-- results/  # [Not included; can be regenerated running the code]
 +-- src/      # bash scripts and other softwares used in the analyses
```

## Setup

### Data download
All sequence data was generated during the study and is publicaly available in the
ENA Archive (Project number: PRJEB103802).


| **Outgroup**                  | ***Hyles***                              | ***Rhodafra***         |
|:------------------------------|:-----------------------------------------|:-----------------------|
| *Chaerocina dohertyi* - 6374  | *H. annei* - 1332,6683,6947,6945         | *R. opheletes* - 11686 |
| *Euchloron megaera* - 6376    | *H. biguttata* - 508,510,4738            | *R. opheletes* - 12003 |
| *Hippotion echeclus* - 6382   | *H. calida* - 4349, ms852,Sa190la        | *R. opheletes* - 13597 |
| *Hippotion celerio* - 503     | *H. euphorbiarum* - 3623,3722,4379       | *R. opheletes* - 13598 |
| *Theretra alecto* - 6386      | *H. gallii* - 12608,3939,489             | *R. marshalli* - 13599 |
| *Theretra japonia* - 6389     | *H. lineata* - 4063,4340,6016            | *R. marshalli* - 13600 |
| *Theretra silhetensis* - 6381 | *H. livornica* - 13065,13097,13252,13270 |                        |
| *Xylophanes falco* - 6377     | *H. livornicoides* - 9001,8995,6408      |                        |
| *Xylophanes porcus* - 6391    | *H. nicaea* - 13121,13050,6046           |                        |
| *Xylophanes tersa* - 6384     | *H. perkinsi* - cg09,cg13                |                        |
|                               | *H. vespertilio* - 13068,13069,13215,628 |                        |
|                               | *H. wilsoni* - cg74,cg75                 |                        |


The following code allows to download the data and rename the files properly.

```
...
...
...

```

### Environments and software

In order to run the code in this repository, is it also necessary to download some
software and to create some conda and python environments. All packages installed
in each used environment is available in the folder `envirs`. The following code
install all these environments and download all necessary software to run the code


```
# Conda environments
## make sure conda is properly sourced
conda create -p envirs/astral --file envirs/astral_spec.txt
conda create -p envirs/fastq_preprocessing --file envirs/fastq_preprocessing_spec.txt
conda create -p envirs/hDNA_pre --file envirs/hDNA_pre_spec.txt
conda create -p envirs/hDNA_map --file envirs/hDNA_map_spec.txt
conda create -p envirs/introgression --file envirs/introgression_spec.txt
conda create -p envirs/py2.7 --file envirs/py2.7_spec.txt
# Python virtual environments
python -m venv envirs/py_venv_HyDe
source activate envirs/py_venv_HyDe/bin/activate
pip install -r envirs/py_venv_HyDe_req.txt
deactivate

# Other software
## Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -O src/Trimmomatic-0.39.zip
unzip src/Trimmomatic-0.39.zip
## genomics_general
git clone https://github.com/simonhmartin/genomics_general src/genomics_general
## PopLDdecay
git clone https://github.com/BGI-shenzhen/PopLDdecay src/PopLDdecay
## ascbias.py
wget https://github.com/btmartin721/raxml_ascbias/blob/master/ascbias.py -O src/ascbias.py
chmod +x src/ascbias.py
## Paup
wget https://phylosolutions.com/paup-test/paup4a168_centos64.gz -O src/paup4a168_centos64.gz
gunzip src/paup4a168_centos64.gz
chmod +x src/paup4a168_centos64
## For SNAPPER analyses
wget https://github.com/mmatschiner/snapp_prep/blob/master/snapp_prep.rb -O src/snapp_prep.rb
wget https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz -O src/BEAST.v2.7.7.Linux.x86.tgz 
tar fxz src/BEAST.v2.7.7.Linux.x86.tgz
## Dsuite and Dtrios scrips
git clone https://github.com/millanek/Dsuite.git src/Dsuite
wget https://github.com/mmatschiner/tutorials/blob/cbe81ecbb23536d58f2ab301c83898d2e74a8f73/analysis_of_introgression_with_snp_data/src/plot_d.rb -O src/plot_d.rb
wget https://github.com/mmatschiner/tutorials/blob/cbe81ecbb23536d58f2ab301c83898d2e74a8f73/analysis_of_introgression_with_snp_data/src/plot_f4ratio.rb -O src/plot_f4ratio.rb
wget https://github.com/mmatschiner/tutorials/blob/cbe81ecbb23536d58f2ab301c83898d2e74a8f73/analysis_of_introgression_with_snp_data/src/vcf2phylip.py -O src/vcf2phylip.py
## QuIBL
git clone https://github.com/miriammiyagi/QuIBL src/QuIBL
## Twisst
git clone https://github.com/simonhmartin/twisst src/twisst
```

### Important note
The analyses was run at the HPC at TUD and the scripts are configured
to work in that specific cluster. You must fix first the loaded modules
or remove them if you are not working on a Slurm machine. Besides,
you have to set the variable `WS` in every script accordingly.
All scrips are prepared to be run from whichever path.

## Executing the pipeline

The `.sh` files in `src` folder contain the code to execute all the analyses.
They are named with numbers as prefixes and are designed to be run in that
specific order:
 
Other `.py`, `.R` or `.nex` scripts are run within the `*.sh` files.
For example the script `06_02_plot_twisst.R` is run during the execution of
`06_02_twisst_introgression.sh`

  1. **Preprocessing and Quality Control**
      * `01_freshDNA_preprocessing.sh`
        Perform QC, filtering and trimming of fresh DNA raw reads
      * `01_hDNA_preprocessing.sh`
        Perform QC, filtering and trimming of hDNA raw reads
     
  1. **Mapping**
      * `02_00_hDNA_mapping_optimization.sh` & `02_00.01_summary_optimization.py`
        Test a set of parameters and mappers with hDNA samples to select the
        best option before mapping them.
      * `02_01_genome_index.sh`
        Index the reference genome
      * `02_02_mapping_mergedReads.sh`; `02_03_mapping_unassembledReads.sh`; `02_04_mergeResults.sh`
        Map the hDNA samples against the reference genome in two steps.
        First the paired & merged reads; then the paired & un-merged reads.
        Finnaly, it merges the assemblies into a single `bam` file
      * `02_05_mapping_freshDNA.sh`
        Map the reads of the fresh-DNA samples
        
  1. **Variant Calling**
      * `03_01_variantCalling.sh`
        Perform the variant calling with `bcftools` generating all three datasets:
        `full` (all samples), `1samp-x-sp` (only the sample with the highest coverage
        per species),, `1samp-x-sp_noRhodafra` (as 1samp-x-sp, but excluding
        *Rhodafra* [hDNA] samples) and `strict` (as 1samp-x-sp, but applying stricter
        filtering rules). Each dataset is also thinned, selecting only 1 SNP
        each 1000 bp.
        
  1. **Autosomes phylogenies**
      * `04_00_getConsensus.sh`
        Get the SNPs in formats `phyllip`, `fasta` and `nexus` to be used by different
        software.
      * `04_02_phyloRAxML.sh`
        Run the Maximum Likelihood phylogeny with RAxML for the concatenated matrix of
        autosome SNPs
      * `04_03_SVDquartets`; `04_03_PAUP_code_1samp-x-sp.nex`; `04_03_PAUP_code_full.nex`
        Run the SVDquartets phylogeny. The base PAUP code is in `.nex` scripts.
      * `04_04_SNAPP_bayesian.sh`
        Run the SNAPPER bayesian phylogeny and estimate node dates.
      * `04_05_ASTRAL_slidingwindows.sh`
        Get the "gene"-trees in sliding windows and run ASTRAL phylogeny
        
  1. **Mitochondrial phylogenies**
      * `05_00_getConsensusMito.sh`
        Get the consensus sequences of the mitochondria
      * `05_02_phyloRAxMLmito.sh`
        Align the mitchondrial genomes and run the Maximum Likelihood phylogeny with RAxML.
        
  1. **Introgression/ILS assesement**
      * `06_01_Dtrios_introgression.sh`
        Run the Dtrios analysis
      * `06_02_twisst_introgression.sh`; `06_02_plot_twisst.R`
        Run the TWISST analyisis for a reduced phylogeny with *Rhodafra* and two *Hyles* species
      * `06_04_QuIBL.sh`; `06_04.1_plotQuIBL.R`
        Run the QuIBL analysis.
      * `06_05_HyDe.sh`
        Run the HyDe analysis.

