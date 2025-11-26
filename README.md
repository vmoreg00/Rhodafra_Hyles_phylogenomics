# Shedding light on *Hyles* Hübner, 1819 and *Rhodafra* Rothschild & Jordan, 1903 (Lepidoptera: Sphingidae) evolutionary history using DNA from museum specimens

Scripts used for the phylogenomics análisis of the paper 
Moreno-González, V., Melichar, T., San Jose, M., Rubinoff, D. & Hundsdoerfer, A. K.
(in prep.) Shedding light on *Hyles* Hübner, 1819 and *Rhodafra* Rothschild & Jordan,
1903 (Lepidoptera: Sphingidae) evolutionary history using DNA from museum specimens

### Scripts (`src/`)
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
        per species), and `1samp-x-sp_noRhodafra` (as 1samp-x-sp, but excluding
        *Rhodafra* [hDNA] samples). Each dataset is also thinned, selecting only 1 SNP
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

### Results (`results/`)
