# script based on http://phylo.wikidot.com/biogeobears#script
# created 2026-05-18
# author: V. Moreno-Gonzalez

# Setup =======================================================================
library(ape)
library(optimx)
library(GenSA) 
library(rexpokit)
library(cladoRcpp)
library(snow)
library(parallel)
library(BioGeoBEARS)
library(ggplot2)
library(ggtree)

wd <- "/data/horse/ws/vimo762h-hDNAhyles"
resdir <- paste0(wd, "/results/07_ancestral_range")
tree_file <- paste0(resdir, "/0_time_calibrated_tree.nwk")
data_file <- paste0(wd, "/data/biogeography_data.phy")
dist_file <- paste0(resdir, "/0_distance_matrix.txt")
ncores <- 40

# Calibrate the tree ==========================================================
# Biogeobears requires a time-calibrated (ultrametric) tree. The bayesian
# phylogeny with time calibration was not performed with the full set of species.
# Hence, I'm using ape::chronos to get a relaxed-clock calibrated phylogeny.
# The calibration points are extracted from Couch et al. 2026
if ( ! file.exists(tree_file)) {
  in_tree <- read.tree(paste0(wd, "/results/04_SNP_phylogeny/02_RAxML/SNP_phylo_raxml_1samp-x-sp.raxml.support"))
  # # Check the nodes to apply the calibration points in
  # library(ggplot2)
  # ggplot() +
  #   geom_tree() +
  #   geom_nodelab(aes(label = 1:length(c(in_tree$node.label, in_tree$tip.label)))) +
  #   geom_tiplab(aes(label = label), parse = F, align = F)

  # Calibrate the tree running 10000 times
  set.seed(666)
  lengths <- matrix(data = NA, nrow = 1000, ncol = length(in_tree$edge.length))
  for (i in 1:1000) {
    tree_c <- 
      chronos(in_tree, model = "relaxed",
              calibration = data.frame(node = c(39, 33, 31, 30, 29, 25),
                                       age.min = c(2.3, 5.1, 6.4, 8.8, 8.3, 23.1),
                                       age.max = c(6.1, 11, 12.6, 16, 16.2, 28.1),
                                       soft.bounds = T),
              control = chronos.control(iter.max = 100),
              quiet = T)
    lengths[i,] <- tree_c$edge.length
    if (i %% 10 == 0) {
      message(i, "/1000 iterations")
    }
  }; rm(i, tree_c)
  # Get the mean edge length and 95% CI
  in_tree$edge.length <- colMeans(lengths, na.rm = T)
  in_tree$node.label <-
    as.list(data.frame(apply(lengths, 2, \(x){
      ci <- confint(lm(x ~ 1))
      range(x, na.rm = T)
    })))

  # Plot the tree
  ggplot(in_tree) +
    geom_tree() +
    geom_tiplab(aes(label = label), parse = F, align = F) +
    geom_range(range = "label", color='red', alpha=.6, size=2) +
    scale_x_continuous(name = "Million years before present",
                       breaks = seq(sum(in_tree$edge.length[1:2]), 0, by = -5), 
                       labels = seq(0, sum(in_tree$edge.length[1:2]), by = 5),
                       expand = c(.05, 0, .20, 0)) +
    theme_tree2() +
    theme(axis.title.x = element_text(face = "bold"))
  in_tree$node.label <- NULL
  write.tree(in_tree, tree_file)
}

# Get the distance matrix for the +x models ===================================
if ( ! file.exists(dist_file) ) {
  # Approximate distances (in km) between areas coasts
  ## Diagonal should have 0's
  ## Adjacent areas are set to 1
  distances <- matrix(c(
    #A(NA) B(SA)  C(HI)  D(MD)  E(AU)  F(AF)  G(PA)
    0,     1,     3500,  13000, 13000, 5500,  100,   # A (North America)
    1,     0,     8500,  9000,  13000, 3000,  6000,  # B (South America)
    3500,  8500,  0,     17000, 7000,  15000, 5000,  # C (Hawaii)
    13000, 9000,  17000, 0,     6500,  450,   3800,  # D (Madagascar)
    13000, 13000, 7000,  6500,  0,     7500,  100,   # E (Australia)
    5500,  3000,  15000, 450,   7500,  0,     1,     # F (Africa, South of Sahara)
    100,   6000,  5000,  3800,  100,   1,     0      # G (Palearctic, inc. N Africa)
  ), nrow = 7, byrow = TRUE)
  rownames(distances) <- colnames(distances) <- LETTERS[1:7]
  
  # rescale the matrix:
  # Divide by 1000; the distance matrix will represent distance in thousands km
  distances <- distances/1000
  write.table(distances, file = dist_file, quote = F,
              sep = "\t", row.names = F, col.names = T)
  write("\n\nEND", dist_file, append=T)
}
# Get tipranges ===============================================================
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=data_file)
max_range_size = 7

# Run DEC =====================================================================
if ( ! file.exists(paste0(resdir, "/11_DEC.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up DEC model
  # (nothing to do; defaults)

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resDEC = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resDEC, file = paste0(resdir, "/11_DEC.RDa"))
  pdf(paste0(resdir, "/11_DEC_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDEC)
  dev.off()
} else {
  load(paste0(resdir, "/11_DEC.RDa"))
}

# Run DEC+J ==================================================================
if ( ! file.exists(paste0(resdir, "/12_DECj.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

  # Check object
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run
  resDECj = bears_optim_run(BioGeoBEARS_run_object)
  # Save output
  save(resDECj, file = paste0(resdir, "/12_DECj.RDa"))
  pdf(paste0(resdir, "/12_DECj_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDECj)
  dev.off()
} else {
  load(paste0(resdir, "/12_DECj.RDa"))
}

# Run DEC+x ==================================================================
if ( ! file.exists(paste0(resdir, "/13_DECx.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$distsfn = dist_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

  # Set up DEC+x model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  # Add x as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0

  # Check object
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run
  resDECx = bears_optim_run(BioGeoBEARS_run_object)
  # Save output
  save(resDECx, file = paste0(resdir, "/13_DECx.RDa"))
  pdf(paste0(resdir, "/13_DECx_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDECx)
  dev.off()
} else {
  load(paste0(resdir, "/13_DECx.RDa"))
}

# Run DEC+x+j ================================================================
if ( ! file.exists(paste0(resdir, "/14_DECxj.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$distsfn = dist_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = FALSE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

  # Set up DEC+x model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = resDECj$output@params_table["j", "est"]
  xstart = resDECx$output@params_table["x", "est"]
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  # Add x as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart

  # Check object
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Run
  resDECxj = bears_optim_run(BioGeoBEARS_run_object)
  # Save output
  save(resDECxj, file = paste0(resdir, "/14_DECxj.RDa"))
  pdf(paste0(resdir, "/14_DECxj_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDECxj)
  dev.off()
} else {
  load(paste0(resdir, "/14_DECxj.RDa"))
}

# Run DIVALIKE ===============================================================
if ( ! file.exists(paste0(resdir, "/21_DIVALIKE.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE	# set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up DIVALIKE model
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resDIVALIKE = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resDIVALIKE, file = paste0(resdir, "/21_DIVALIKE.RDa"))
  pdf(paste0(resdir, "/21_DIVALIKE_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDIVALIKE)
  dev.off()
} else {
  load(paste0(resdir, "/21_DIVALIKE.RDa"))
}

# Run DIVALIKE+J ===============================================================
if ( ! file.exists(paste0(resdir, "/22_DIVALIKEj.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up DIVALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = 0.0001

  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resDIVALIKEj = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resDIVALIKEj, file = paste0(resdir, "/22_DIVALIKEj.RDa"))
  pdf(paste0(resdir, "/22_DIVALIKEj_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDIVALIKEj)
  dev.off()
} else {
  load(paste0(resdir, "/22_DIVALIKEj.RDa"))
}

# Run DIVALIKE+X ===============================================================
if ( ! file.exists(paste0(resdir, "/23_DIVALIKEx.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$distsfn = dist_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE	# set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up DIVALIKE+X model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]

  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

  # Add +x parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resDIVALIKEx = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resDIVALIKEx, file = paste0(resdir, "/23_DIVALIKEx.RDa"))
  pdf(paste0(resdir, "/23_DIVALIKEx_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDIVALIKEx)
  dev.off()
} else {
  load(paste0(resdir, "/23_DIVALIKEx.RDa"))
}


# Run DIVALIKE+X+J =============================================================
if ( ! file.exists(paste0(resdir, "/24_DIVALIKExj.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$distsfn = dist_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = FALSE                # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up DIVALIKE+X+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = resDIVALIKEj$outputs@params_table["j","est"]
  xstart = resDIVALIKEx$outputs@params_table["x","est"]
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

  # Add +x parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = jstart

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resDIVALIKExj = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resDIVALIKExj, file = paste0(resdir, "/24_DIVALIKExj.RDa"))
  pdf(paste0(resdir, "/24_DIVALIKExj_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resDIVALIKExj)
  dev.off()
} else {
  load(paste0(resdir, "/24_DIVALIKExj.RDa"))
}

# Run BAYAREALIKE =============================================================
if ( ! file.exists(paste0(resdir, "/31_BAYAREALIKE.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up BAYAREALIKE model
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resBAYAREALIKE = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resBAYAREALIKE, file = paste0(resdir, "/31_BAYAREALIKE.RDa"))
  pdf(paste0(resdir, "/31_BAYAREALIKE_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resBAYAREALIKE)
  dev.off()
} else {
  load(paste0(resdir, "/31_BAYAREALIKE.RDa"))
}

# Run BAYAREALIKE+J ===========================================================
if ( ! file.exists(paste0(resdir, "/32_BAYAREALIKEj.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE	# set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = 0.0001

  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resBAYAREALIKEj = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resBAYAREALIKEj, file = paste0(resdir, "/32_BAYAREALIKEj.RDa"))
  pdf(paste0(resdir, "/32_BAYAREALIKEj_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resBAYAREALIKEj)
  dev.off()
} else {
  load(paste0(resdir, "/32_BAYAREALIKEj.RDa"))
}

# Run BAYAREALIKE+X ===========================================================
if ( ! file.exists(paste0(resdir, "/33_BAYAREALIKEx.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$distsfn = dist_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE      # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE                 # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up BAYAREALIKE+X model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]

  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

  # Add +x parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0

  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resBAYAREALIKEx = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resBAYAREALIKEx, file = paste0(resdir, "/33_BAYAREALIKEx.RDa"))
  pdf(paste0(resdir, "/33_BAYAREALIKEx_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resBAYAREALIKEx)
  dev.off()
} else {
  load(paste0(resdir, "/33_BAYAREALIKEx.RDa"))
}

# Run BAYAREALIKE+X+J =========================================================
if ( ! file.exists(paste0(resdir, "/34_BAYAREALIKExj.RDa"))) {
  # Intitialize object
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree_file
  BioGeoBEARS_run_object$geogfn = data_file
  BioGeoBEARS_run_object$distsfn = dist_file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE	# set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Speed options
  BioGeoBEARS_run_object$on_NaN_error = -1e50           # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = FALSE                # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE              # if FALSE, use optim() instead of optimx();
  BioGeoBEARS_run_object$num_cores_to_use = ncores
  BioGeoBEARS_run_object$force_sparse = FALSE           # force_sparse=TRUE causes pathology & isn't much faster at this scale
  # Read and settings
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE           # get ancestral states from optim run

  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = resBAYAREALIKEj$outputs@params_table["j","est"]
  xstart = resBAYAREALIKEx$outputs@params_table["x","est"]

  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

  # Add +x parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0

  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.9999999

  # Run this to check inputs. Read the error messages if you get them!
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Run
  resBAYAREALIKExj = bears_optim_run(BioGeoBEARS_run_object)

  # Save output
  save(resBAYAREALIKExj, file = paste0(resdir, "/34_BAYAREALIKExj.RDa"))
  pdf(paste0(resdir, "/34_BAYAREALIKExj_tree.pdf"),
      width = 8, height = 8)
  plot_BioGeoBEARS_results(resBAYAREALIKExj)
  dev.off()
} else {
  load(paste0(resdir, "/34_BAYAREALIKExj.RDa"))
}


# CALCULATE SUMMARY STATISTICS TO COMPARE =====================================
# DEC,         DEC+J,         DEC+J,         DEC+J+X
# DIVALIKE,    DIVALIKE+J,    DIVALIKE+X,    DIVALIKE+J+X
# BAYAREALIKE, BAYAREALIKE+J, BAYAREALIKE+X, BAYAREALIKE+J+X

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

# Statistics -- DEC -----------------------------------------------------------
# Get LnL
LnL_2param     = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_3param_j   = get_LnL_from_BioGeoBEARS_results_object(resDECj)
LnL_3param_x   = get_LnL_from_BioGeoBEARS_results_object(resDECx)
LnL_4param_j_x = get_LnL_from_BioGeoBEARS_results_object(resDECxj)
# Perform LRT for pairs of nested models
lrt_null_vs_j = AICstats_2models(LnL_1 = LnL_3param_j, LnL_2 = LnL_2param,
                                 numparams1 = 3,       numparams2 = 2)
lrt_null_vs_x = AICstats_2models(LnL_1 = LnL_3param_x, LnL_2 = LnL_2param,
                                 numparams1 = 3,       numparams2 = 2)
lrt_j_vs_jx = AICstats_2models(LnL_1 = LnL_4param_j_x, LnL_2 = LnL_3param_j,
                               numparams1 = 4,         numparams2 = 3)
lrt_x_vs_jx = AICstats_2models(LnL_1 = LnL_4param_j_x, LnL_2 = LnL_3param_x,
                               numparams1 = 4,         numparams2 = 3)
# Get parameters
## Null model
par_n = extract_params_from_BioGeoBEARS_results_object(results_object = resDEC,
                                                       returnwhat = "table",
                                                       addl_params = c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + j
par_j = extract_params_from_BioGeoBEARS_results_object(results_object = resDECj,
                                                       returnwhat = "table",
                                                       addl_params=c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + x
par_x = extract_params_from_BioGeoBEARS_results_object(results_object = resDECx,
                                                       returnwhat = "table",
                                                       addl_params=c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + j + x
par_j_x = extract_params_from_BioGeoBEARS_results_object(results_object = resDECxj,
                                                         returnwhat = "table",
                                                         addl_params=c("j", "x"),
                                                         paramsstr_digits = 4)


# Merge results
restable = rbind(restable, "DEC" = par_n, "DEC+j" = par_j,
                 "DEC+x" = par_x, "DEC+j+x" = par_j_x)
teststable = rbind(teststable,
                   "DEC vs. DEC+j" = conditional_format_table(lrt_null_vs_j),
                   "DEC vs. DEC+x" = conditional_format_table(lrt_null_vs_x),
                   "DEC+j vs. DEC+j+x" = conditional_format_table(lrt_j_vs_jx),
                   "DEC+x vs. DEC+j+x" = conditional_format_table(lrt_x_vs_jx))

# Statistics -- DIVALIKE ------------------------------------------------------
# Get LnL
LnL_2param     = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_3param_j   = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
LnL_3param_x   = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEx)
LnL_4param_j_x = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKExj)
# Perform LRT for pairs of nested models
lrt_null_vs_j = AICstats_2models(LnL_1 = LnL_3param_j, LnL_2 = LnL_2param,
                                 numparams1 = 3,       numparams2 = 2)
lrt_null_vs_x = AICstats_2models(LnL_1 = LnL_3param_x, LnL_2 = LnL_2param,
                                 numparams1 = 3,       numparams2 = 2)
lrt_j_vs_jx = AICstats_2models(LnL_1 = LnL_4param_j_x, LnL_2 = LnL_3param_j,
                               numparams1 = 4,         numparams2 = 3)
lrt_x_vs_jx = AICstats_2models(LnL_1 = LnL_4param_j_x, LnL_2 = LnL_3param_x,
                               numparams1 = 4,         numparams2 = 3)
# Get parameters
## Null model
par_n = extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKE,
                                                       returnwhat = "table",
                                                       addl_params = c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + j
par_j = extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKEj,
                                                       returnwhat = "table",
                                                       addl_params=c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + x
par_x = extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKEx,
                                                       returnwhat = "table",
                                                       addl_params=c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + j + x
par_j_x = extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKExj,
                                                         returnwhat = "table",
                                                         addl_params=c("j", "x"),
                                                         paramsstr_digits = 4)


# Merge results
restable = rbind(restable, "DIVALIKE" = par_n, "DIVALIKE+j" = par_j,
                 "DIVALIKE+x" = par_x, "DIVALIKE+j+x" = par_j_x)
teststable = rbind(teststable,
                   "DIVALIKE vs. DIVALIKE+j" = conditional_format_table(lrt_null_vs_j),
                   "DIVALIKE vs. DIVALIKE+x" = conditional_format_table(lrt_null_vs_x),
                   "DIVALIKE+j vs. DIVALIKE+j+x" = conditional_format_table(lrt_j_vs_jx),
                   "DIVALIKE+x vs. DIVALIKE+j+x" = conditional_format_table(lrt_x_vs_jx))

# Statistics -- BAYAREALIKE ---------------------------------------------------
# Get LnL
LnL_2param     = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_3param_j   = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
LnL_3param_x   = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEx)
LnL_4param_j_x = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKExj)
# Perform LRT for pairs of nested models
lrt_null_vs_j = AICstats_2models(LnL_1 = LnL_3param_j, LnL_2 = LnL_2param,
                                 numparams1 = 3,       numparams2 = 2)
lrt_null_vs_x = AICstats_2models(LnL_1 = LnL_3param_x, LnL_2 = LnL_2param,
                                 numparams1 = 3,       numparams2 = 2)
lrt_j_vs_jx = AICstats_2models(LnL_1 = LnL_4param_j_x, LnL_2 = LnL_3param_j,
                               numparams1 = 4,         numparams2 = 3)
lrt_x_vs_jx = AICstats_2models(LnL_1 = LnL_4param_j_x, LnL_2 = LnL_3param_x,
                               numparams1 = 4,         numparams2 = 3)
# Get parameters
## Null model
par_n = extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKE,
                                                       returnwhat = "table",
                                                       addl_params = c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + j
par_j = extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKEj,
                                                       returnwhat = "table",
                                                       addl_params=c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + x
par_x = extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKEx,
                                                       returnwhat = "table",
                                                       addl_params=c("j", "x"),
                                                       paramsstr_digits = 4)
## Null model + j + x
par_j_x = extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKExj,
                                                         returnwhat = "table",
                                                         addl_params=c("j", "x"),
                                                         paramsstr_digits = 4)


# Merge results
restable = rbind(restable, "BAYAREALIKE" = par_n, "BAYAREALIKE+j" = par_j,
                 "BAYAREALIKE+x" = par_x, "BAYAREALIKE+j+x" = par_j_x)
teststable = rbind(teststable,
                   "BAYAREALIKE vs. BAYAREALIKE+j" = conditional_format_table(lrt_null_vs_j),
                   "BAYAREALIKE vs. BAYAREALIKE+x" = conditional_format_table(lrt_null_vs_x),
                   "BAYAREALIKE+j vs. BAYAREALIKE+j+x" = conditional_format_table(lrt_j_vs_jx),
                   "BAYAREALIKE+x vs. BAYAREALIKE+j+x" = conditional_format_table(lrt_x_vs_jx))

# Save ------------------------------------------------------------------------
teststable$alt <- sapply(strsplit(rownames(teststable), split = " vs. "), \(x) x[2])
teststable$null <- sapply(strsplit(rownames(teststable), split = " vs. "), \(x) x[1])
save(restable, file=paste0(resdir, "/restable_v1.RDa"))
save(teststable, file=paste0(resdir, "/teststable_v1.RDa"))
write.table(restable, file=paste0(resdir, "/restable.txt"), quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file=paste0(resdir, "/teststable.txt"),
            quote=FALSE, sep="\t", row.names = F)

# Final plot ------------------------------------------------------------------
library(dplyr)
library(tidyr)
options(device = "pdf")

# Read the tree
tr = read.tree(file=tree_file)
tr$tip.label <- gsub("^.*_Hyles-", "H. ", tr$tip.label)
tr$tip.label <- gsub("-catissima", "", tr$tip.label)
tr$tip.label <- gsub("^.*_Xylophanes-", "X. ", tr$tip.label)
tr$tip.label <- gsub("^.*_Euchloron-", "E. ", tr$tip.label)
tr$tip.label <- gsub("^.*_Theretra-", "T. ", tr$tip.label)
tr$tip.label <- gsub("^.*_Hippotion-", "Hi. ", tr$tip.label)
tr$tip.label <- gsub("^.*_Chaerocina-", "C. ", tr$tip.label)
tr$tip.label <- gsub("^.*Rmars_h", "R. marshali", tr$tip.label)
tr$tip.label <- gsub("^.*Rophe_h", "R. opheltes", tr$tip.label)

# Get data to plot
## Areas (A-F)
areas = getareas_from_tipranges_object(tipranges)
## states (=colnames of the matrix of states)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas,
                                                    maxareas=max_range_size,
                                                    include_null_range=TRUE)
## Probability states at nodes and corners
probs_nodes <-
  infprobs_to_probs_of_each_area(resDECx$ML_marginal_prob_each_state_at_branch_top_AT_node, 
                                 states_list = states_list_0based) %>%
  as.data.frame() %>%
  mutate(node = 1:47)
probs_corners <-
  infprobs_to_probs_of_each_area(resDECx$ML_marginal_prob_each_state_at_branch_bottom_below_node,
                                 states_list = states_list_0based) %>%
  as.data.frame()
colnames(probs_corners) <- paste0("C", 1:ncol(probs_corners))
probs_nodes <- cbind(probs_nodes, probs_corners)
# Plot 1 - Matrix of states
plot_matrix_states <-
  resDECx$ML_marginal_prob_each_state_at_branch_top_AT_node %>%
  data.frame() %>%
  cbind(., x = 1:47) %>%
  pivot_longer(cols = X1:X128) %>%
  mutate(name = as.numeric(gsub("X", "", name))) %>%
  ggplot() +
  geom_tile(aes(y = x, x = name, fill = value)) +
  geom_hline(yintercept=seq(.5, 47.5, by = 1), color = "gray") +
  geom_vline(xintercept=seq(.5,128.5, by = 1), color = "gray") +
  scale_x_continuous(name = "", breaks = 1:128, expand = c(0,0),
                     labels = sapply(states_list_0based, \(x) paste0(areas[x+1], collapse = "")),
                     sec.axis = sec_axis(transform = ~., name = "", breaks = 1:128,
                                         labels = sapply(states_list_0based, \(x) paste0(areas[x+1], collapse = "")))) +
  scale_y_continuous(name = "node", breaks = 1:47, expand = c(0,0),
                     sec.axis = sec_axis(transform = ~., name = "", breaks = 1:47)) +
  scale_fill_gradient(low = "white", high = "red3", limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust = .5, hjust = 1),
        axis.text.x.top = element_text(angle = 90, vjust = .5, hjust = 0)) +
  coord_equal() +
  #geom_hline(yintercept=31) +
  theme(legend.position = "top",
        legend.key.width = unit(3, "cm"),
        text = element_text(size = 8))
ggsave(paste0(resdir, "/plot_matrix_states.pdf"), plot_matrix_states,
       width = 30/2.54, height = 17/2.54)

plot_matrix_states_summary <-
  resDECx$ML_marginal_prob_each_state_at_branch_top_AT_node %>%
  data.frame() %>%
  cbind(., x = 1:47) %>%
  pivot_longer(cols = X1:X128) %>%
  mutate(name = as.numeric(gsub("X", "", name)),
         value = ifelse(value < 0.1, NA, value),
         name = factor(x = name, levels = 1:128, labels = sapply(states_list_0based, \(x) paste0(areas[x+1], collapse = "")))) %>%
  na.omit() %>%
  droplevels() %>%
  ggplot() +
  geom_tile(aes(y = x, x = name, fill = value)) +
  geom_hline(yintercept=seq(.5, 47.5, by = 1), color = "gray") +
  geom_vline(xintercept=seq(.5,24.5, by = 1), color = "gray") +
  scale_y_continuous(name = "node", breaks = 1:47, expand = c(0,0),
                     sec.axis = sec_axis(transform = ~., name = "", breaks = 1:47)) +
  scale_fill_gradient(low = "white", high = "red3", limits = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust = .5, hjust = 1),
        axis.text.x.top = element_text(angle = 90, vjust = .5, hjust = 0)) +
  coord_equal() +
  theme(legend.position	= "top", 
        legend.key.width = unit(1, "cm"), 
        text = element_text(size = 8))
ggsave(paste0(resdir, "/plot_matrix_states_summary.pdf"), plot_matrix_states_summary,
       width = 10/2.54, height = 17/2.54)

# Plot 2 - Tree with states
plot_tree <-
  ggplot(tr) %<+% probs_nodes +
  geom_tree() +
  geom_tiplab(offset = 1, fontface = 3) +
  # add the scales
  scale_x_continuous(name = "Million years before present",
                     breaks = seq(sum(tr$edge.length[1:2]), 0, by = -5),
                     labels = seq(0, sum(tr$edge.length[1:2]), by = 5),
                     expand = c(.05, 0, .20, 0)) +
  theme_tree2() +
  theme(axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank())
## Draw the nodes' matrices
node_matrices <- corner_matrices <- list(nrow(probs_nodes))
for (i in 1:nrow(probs_nodes)) {
  # Node
  node_matrices[[probs_nodes$node[i]]] <-
    probs_nodes[i,] %>%
    select(V1:node) %>%
    pivot_longer(V1:V7) %>%
    mutate(x1 = case_when(name %in% c("V1", "V2") ~ -2,
                          name == "V6" ~ -1,
                          name == "V7" ~ -2/3,
                          name == "V4" ~ 0,
                          name == "V3" ~ 2/3,
                          name == "V5" ~ 1),
           x2 = case_when(name %in% c("V1", "V7", "V3") ~ x1 + (4/3),
                          name %in% c("V2", "V6", "V4", "V5") ~ x1 + 1),
           y1 = case_when(name %in% c("V1", "V7", "V3") ~ 0,
                          name %in% c("V2", "V6", "V4", "V5") ~ -1),
           y2 = y1 + 1) %>%
    ggplot() +
    geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2,
                  group = name, fill = value), color = "black") +
    scale_fill_gradient(name = "", low = "white", high = "firebrick") +
    theme_void() +
    theme(legend.position = "none") +
    coord_equal()
  # Corner
  corner_matrices[[probs_nodes$node[i]]] <-
    probs_nodes[i,] %>%
    select(node:C7) %>%
    pivot_longer(C1:C7) %>%
    mutate(x1 = case_when(name %in% c("C1", "C2") ~ -2,
                          name == "C6" ~ -1,
                          name == "C7" ~ -2/3,
                          name == "C4" ~ 0,
                          name == "C3" ~ 2/3,
                          name == "C5" ~ 1),
           x2 = case_when(name %in% c("C1", "C7", "C3") ~ x1 + (4/3),
                          name %in% c("C2", "C6", "C4", "C5") ~ x1 + 1),
           y1 = case_when(name %in% c("C1", "C7", "C3") ~ 0,
                          name %in% c("C2", "C6", "C4", "C5") ~ -1),
           y2 = y1 + 1) %>%
    ggplot() +
    geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2,
                  group = name, fill = value), color = "black") +
    scale_fill_gradient(name = "", low = "white", high = "firebrick") +
    theme_void() +
    theme(legend.position = "none") +
    coord_equal()
  # Add to the tree
  parent <- plot_tree@data$parent[plot_tree@data$node == probs_nodes$node[i]]
  plot_tree <- 
    plot_tree +
    annotation_custom(ggplotGrob(node_matrices[[probs_nodes$node[i]]]),
                      xmin = plot_tree@data$x[plot_tree@data$node == probs_nodes$node[i]]-.5,
                      xmax = plot_tree@data$x[plot_tree@data$node == probs_nodes$node[i]]+.5,
                      ymin = plot_tree@data$y[plot_tree@data$node == probs_nodes$node[i]]-.25, 
                      ymax = plot_tree@data$y[plot_tree@data$node == probs_nodes$node[i]]+.25) +
    annotation_custom(ggplotGrob(corner_matrices[[probs_nodes$node[i]]]),
                      xmin = plot_tree@data$x[plot_tree@data$node == parent]-.3,
                      xmax = plot_tree@data$x[plot_tree@data$node == parent]+.3,
                      ymin = plot_tree@data$y[plot_tree@data$node == probs_nodes$node[i]]-.2, 
                      ymax = plot_tree@data$y[plot_tree@data$node == probs_nodes$node[i]]+.2)
}; rm(i, parent)
ggsave(paste0(resdir, "/plot_tree.pdf"), plot_tree,
       width = 50, height = 28, units = "cm", dpi = 444)

plot_tree_w_nodes <- plot_tree +
  # add node number
  geom_text(data = ~ subset(., isTip),
            aes(x = x+.75, y = y, label = node)) +
  geom_text(data = ~ subset(., !isTip),
            aes(x = x-.6, y = y+.35, label = node)) 
ggsave(paste0(resdir, "/plot_tree_w_nodes.pdf"), plot_tree_w_nodes,
       width = 50, height = 28, units = "cm", dpi = 444)
