####################################################################
# Compare performance of CV-TMLE vs. AIPW vs. IPW-GLM in realistic sample size
# PLASMODE SIMULATION
#
# Created: 2020 9 March
# Modified: 2020 9 June
#
# Updated with fixed plasmode simulation function (PlasmodeContNew)
#
####################################################################

library(tidyverse)
library(ggplot2)
library(haven)
library(readxl)
library(table1)
# library(ggbeeswarm)
# library(ltmle)
# library(Plasmode)
library(SuperLearner)
#library(ctmle)
library(sl3)
library(tmle3)
library(tictoc)
library(doParallel)
library(doRNG)
library(foreach)
library(randomForest)
rm(list = ls())

cluster <- T;  # Whether we want to run on cluster or the local machine
is.timing <- T # Whether we want to record the run time
parallel.for.DC <- F # This is variable is controlling whether we use parallelism in DC_single fitting
N_sims <- 100# this should <= plas_sim_N
if (cluster == T){
  args = commandArgs(TRUE)
  print(args)
  est.method <- args[[1]];
  non.par <- ifelse(args[[2]]=="non-par",T,F)
  N_sims <- as.numeric(args[[3]])
  is.timing <- !is.na(args[[4]])
  doIPW = 0; doLASSO=0;
  doAIPW=0; doDCAIPW=0
  doManuTMLE=0; doDCTMLE=0; doGComp=0
  num_cf=5
  doShortTMLE = 0
  control=SuperLearner.CV.control(V=2)
  if (est.method=="non-DC"){
    doIPW = 1; doLASSO=1;
    doAIPW=1; doManuTMLE=1;
  }else if (est.method=="DCTMLE"){
    doDCTMLE = 1
  }else if (est.method=="DCAIPW"){
    doDCAIPW = 1
  }else if (est.method=="TMLE"){
    doManuTMLE=1
  }else if (est.method=="AIPW"){
    doAIPW=1
  }else if (est.method=="IPW"){
    doIPW=1
  }else if (est.method=="GComp"){
    doGComp=1
  }
  
  
  if (is.timing==T){
    no_cores <- as.numeric(args[[5]])
  }else{
    no_cores <- detectCores(all.tests = T) - 2
  }
}else{
  non.par <- F
  est.method <- ''
  no_cores <- detectCores(all.tests = T) - 2
  setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no4correct")
}

set.seed(42782)
options(tibble.print_max = 40, tibble.print_min = 30)

registerDoParallel(cores=no_cores)
print(no_cores)

# Set simulation parameters
{
  sims.ver <- "plas"
  # sims.ver <- "5var"
  # sims.ver <- "5var.then.plas"
  # sims.ver <- "5var.hard"
  
  Effect_Size <- 6.6
  
  plas.seed <- 1111
  
  ########
  # parameters for plasmode simulation
  #######
  # p=331 for the whole set 
  data.ver <- "FULL"
  # data.ver <- "13"
  size <- 1178
  # size <-200
  plas_sim_N <- 500; use.subset <- F
  generateA <- T
  estimateWithMore <- F; p.sim = 100;p.est <- 50
  randVar = F # Do we permute variable order?
  isHandPick = T; idx.handpick <- c(1,2,5,18, 217, 2:40)
  interact=T; 
  interact.w.exp = F
  est.interact <- T # Estimate with first order interaction?
  
  
  ########
  # parameters for the simulation in Section 3
  #######
  # Nsets <- 10000
  # Nsamp <- 3000
  Nsets <- 500
  Nsamp <- 600
}

source("20200803-Sims-Function.R")
cat(sims.ver)
sims.obj <- general.sim(sims.ver)

if (sims.ver == "plas"|sims.ver == "5var.then.plas"){
  plas <- sims.obj$plas
  plas_sims <- sims.obj$plas_sims
  vars <- sims.obj$vars
}else{
  sim_boots <- sims.obj
}

# Plot the regression coefficient in Plasmode simulation
plot(plas_sims$TrueOutBeta)
# Output the OR and PS model used for Plasmode
print("Simulation: ")
print(paste0("OR form: ",plas_sims$outForm))
print(paste0("PS form: ",plas_sims$expForm))

##################################
## PARALLELIZE ANALYSES
##################################
# Set param
{
  # DC implementation
  source("20200705-DCDR-Functions.R")
  # getRES function
  source("20200720-Algos-code.R")
  
  
  
  # specify which set of learners for SL
  if (non.par == T){
    ### NON-SMOOTH
    short_tmle_lib <- SL_param_list
    tmle_lib <- lrnr_SL
    aipw_lib <- SL.lib
  }
  else{
    # SMOOTH
    short_tmle_lib <- SL_list
    tmle_lib <- lrnr_SL_param
    aipw_lib <- SL.param
  }
  
  N_sims<- 100
  
  # errorhandling="stop"
  errorhandling="remove"
  
  if (cluster != T){
    doIPW = 0; doLASSO=0;
    doAIPW=0; doDCAIPW=0
    doManuTMLE=1; doDCTMLE=0
    num_cf=5; doGComp=0
    doShortTMLE = 0
    control=SuperLearner.CV.control(V=2)
  }
}




#########
# run
##########
ptm0 <- proc.time()
source("20200904-run-sim-code.R")
ptm1 <- proc.time()

print(paste("It takes",round((ptm1 - ptm0)[3],2),"s"))
print(paste("Cores used:",no_cores))
print(paste("estmate method:",est.method))
print(paste("plas.seed=",plas.seed, "generateA=",generateA, "sims.ver=",sims.ver))
print(paste0("is.timing==",is.timing))
print(paste("./RDataFiles/result-",est.method,"-",args[[2]],"-timing.RData",sep=""))
print(length(boot1))

if (is.timing == F){
  save(boot1, file=paste("./RDataFiles/result-",est.method,"-",args[[2]],".RData",sep="")) 
}else{
  save(boot1, file=paste("./RDataFiles/result-",est.method,"-",args[[2]],"-timing.RData",sep=""))
  tm <- round((ptm1 - ptm0)[3],2)
  tm <- data.frame(tm)
  save(tm, file=paste("./RDataFiles/time-",est.method,"-",no_cores,"cores.RData",sep=""))
}
 

