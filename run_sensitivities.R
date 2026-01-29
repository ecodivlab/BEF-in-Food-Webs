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
n = 5

### generate set of n parameter combination ###
params = data.frame(
  n = 1:n,
  X.exp = rnorm(n,0,1)
)


#### run the flux calculation for each parameter combination ####
results.lakes = apply(params, 1, flux.lakes)
results.marine = apply(params, 1, flux.marine)
results.soils = apply(params, 1, flux.soils)
results.streams = apply(params, 1, flux.streams)

# # parallelised version (check compatibility with windows systems)
# results.lakes = future_apply(params, flux.lakes)
# results$marine = future_apply(params, flux.marine)
# results$soils = future_apply(params, flux.soils)
# results$streams = future_apply(params, flux.streams)
