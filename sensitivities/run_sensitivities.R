#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Calculation of energy fluxes in 48 lake food webs applying sensitivty
#
#=========================================================================================================================


# This code runs the sensitivity analyses for the different ecosystem types. 

rm(list=ls())

## To run this code, a local working directory must be set where all accompanying data and source code are lodged ##
setwd("~/BEF-in-Food-Webs/")



source("Food_web_functions.r")
library(fluxweb); library(cheddar); library(igraph); library(RColorBrewer); 
library(colorspace); library(sp); library(dplyr); library(tibble)

library(future.apply)

options(stringsAsFactors=FALSE)

#### Omnivory function ####
Omnivory.species = function(i, fw, TL){
  # computes omnivory of species i
  # TL: vector of all species' TLs
  if (TL[i] == 1) {
    omn = 0
    return(omn)
  }
  prey = fw[,i] != 0
  omn = 1/sum(fw[,i]) * sum(fw[prey,i] * (TL[prey] - (TL[i] - 1))^2)
  return(omn)
}

### source the functions for the different ecosystem types ###
source("sensitivities/1. FuSED Lakes_sensitivity.r")
source("sensitivities/2. FuSED Marine_sensitivity.R")
source("sensitivities/3. FuSED Soils_sensitivity.R")
source("sensitivities/4. FuSED Streams_sensitivity.R")

#### global variables to be used ####

boltz <- 0.00008617343
T0 <- 273.15 + 20   
perday <- 60*60*24


### set number of replicates to run ###
n = 1500

### generate set of n parameter combination ###
params = data.frame(
  n = 1:n,
  X.exp = rnorm(n,0.71,0.02),
  X.temp = rnorm(n,0.69,0.02),
  eff.exp.inv = rnorm(n,2.266,0.139),
  eff.exp.prod = rnorm(n,0.179,0.161),
  eff.exp.det = rnorm(n,-1.67,0.223),
  eff.temp.prod = rnorm(n,0.164,0.05),
  eff.temp.anim = NA,
  eff.temp.det = NA
)
### temperature effects for intercepts are the same, correct:
params$eff.temp.anim = params$eff.temp.prod
params$eff.temp.det = params$eff.temp.prod

### wrapper to introduce tryCatch
#  for some parameter combination flux calculation is not doable:

wrapper.fun = function(pars, fun){
  tryCatch(fun(as.numeric(pars)), 
           error = function(e){return(NA)})
}

# results.lakes = apply(params, 1, flux.lakes)
# results.marine = apply(params, 1, wrapper.fun, flux.marine)
# results.soils = apply(params, 1, flux.soils)
# results.streams = apply(params, 1, flux.streams)

# # parallelised version (check compatibility with windows systems)
nb.proc = availableCores() - 2 # see number of processors locally
plan(multisession, workers = nb.proc)
results.lakes = future_apply(params, 1, wrapper.fun, flux.lakes)
results.marine = future_apply(params, 1, wrapper.fun, flux.marine)
results.soils = future_apply(params, 1, wrapper.fun, flux.soils)
results.streams = future_apply(params, 1, wrapper.fun, flux.streams)
# 
# results.lakes = 
#   results.marine = 
#   results.soils = 
#   results.streams = 

plan(sequential)
# save in a R object:
save(results.lakes, results.marine, results.soils, results.streams, file = "sensitivities/sensitivity.Rdata")
