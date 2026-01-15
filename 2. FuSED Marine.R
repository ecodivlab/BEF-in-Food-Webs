#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Calculation of energy fluxes in 131 marine food webs
#
#=========================================================================================================================

# This code calculates energy fluxes based on organism metabolic rates (estimated from body mass and temperature), 
# trophic assimilation efficiency (based on resource type and temperature), and food web structure, following methods 
# described in Gauzens et al. (2019) DOI: 10.1111/2041-210X.13109
# The code is run on two datasets: intertidal rockpools and Gulf of Riga Baltic Sea.  
# Code was developed on R version 4.5.0

rm(list=ls())

## To run this code, a local working directory must be set where all accompanying data and source code are lodged ##
#setwd()

source("Food_web_functions.r")
library(fluxweb); library(igraph); library(dplyr); library(cheddar); library(colorspace)
options(stringsAsFactors=FALSE)


#### Intertidal Rockpools ####

boltz <- 0.00008617343
T0 <- 273.15 + 20   
meta.RP <- read.csv("Intertidalrockpools/Intertidalrockypools_metadata.csv")
meta.RP$study_ID <- rep('Intertidal rockpools', nrow(meta.RP))
perday <- 60*60*24
web <- unique(meta.RP$FW_name)

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


for(i in 1:length(web)){
  
  temp.K = meta.RP$temperature_C[meta.RP$FW_name == web[i]] + 273.15
  temp.arr = (temp.K-T0)/(boltz*(temp.K)*T0)
  matrix <-as.matrix(read.csv(paste("Intertidalrockpools/Intertidalrockypools_matrix_",web[i],".csv",sep=""),sep=",",header=T))
  row.names(matrix) <- colnames(matrix)
  diag(matrix) <- 0   # remove cannibalistic links
  igraph <- graph_from_adjacency_matrix(matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard')) #calculates trophic Dissimilarity

  attributes <- read.csv(paste("Intertidalrockpools/Intertidalrockypools_spAttributes_",web[i],".csv",sep=""))
  attributes$bodymass <- attributes$bodymass*4 # *4 convert dry to fresh weight
  attributes$abundance <- attributes$biomass/attributes$bodymass
  attributes$metabolic_type[attributes$metabolic_type=="invertebrate"] <- "ectotherm invertebrate"
  
  animals <- attributes$metabolic_type != "primary producer"  
  plants <- attributes$metabolic_type == "primary producer"

  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                                - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
  attributes$efficiencies[animals] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[plants] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(matrix)
  attributes$omnivory = sapply(1:nrow(matrix), Omnivory.species, matrix, attributes$TL, simplify = TRUE)
  attributes$omnivory[attributes$TL == 1] <- NA
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
  #meta.RP$NAbiomass[i] <- sum(is.na(attributes$biomass))
  ##Fill missing consumer biomass to equal biomass of 1 individual (only 15 out of 3837 nodes, or 0.39% of all nodes)
  attributes$biomass[is.na(attributes$biomass)] <- attributes$bodymass[is.na(attributes$biomass)]

  eat.plants = colSums(matrix[plants, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.plants & eat.animals
  
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes 
  # Total flux
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  # Per capita flux
  PC.flux <- sweep(flux, 1, attributes$abundance, FUN = "/")
  
  ## Food web stability 
  # create vector that encodes for species types (animal, detritus or plant)
  met.types = rep('animal', nrow(matrix))
  met.types[plants] = 'plant'
  # compute stability (leading eigenvalue of the Jacobian)
  stab = stability.value(val.mat = flux,
                         biomasses = attributes$biomass,
                         efficiencies = attributes$efficiencies,
                         metabolic.types = met.types,
                         ef.level = "prey"
  )
  
  ## save stability value:
  meta.RP$stability[i] = stab

  ## Food web fluxes
  meta.RP$flux[i] <- sum(flux)
  meta.RP$second.consumption[i] <- sum(flux[animals])
  meta.RP$prim.consumption[i] <- sum(flux[plants,])  
  # Per capita flux
  meta.RP$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)

  ## Food web metrics
  meta.RP$S[i] <- Number.of.species(matrix)
  meta.RP$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
  meta.RP$MaxTL[i] <- max(attributes$TL)
  meta.RP$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.RP$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.RP$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
  
  meta.RP$L[i] <- Number.of.links(matrix)
  meta.RP$LD[i] <- Link.density(matrix)
  meta.RP$C[i] <- Connectance(matrix)
  meta.RP$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
}




#### Fluxweb analysis Baltic Sea ####

meta.BS <- read.csv("baltic_sea/baltic_metadata.csv")
meta.BS$study_ID <- rep('Baltic Sea', nrow(meta.BS))
web <- unique(meta.BS$FW_name)

for(i in 1:length(web)){
  
  temp.K = meta.BS$temperature_C[meta.BS$FW_name == web[i]] + 273.15
  temp.arr = ((temp.K)-T0)/(boltz*(temp.K)*T0)
  matrix <-as.matrix(read.csv(paste("baltic_sea/baltic_matrix_",web[i],".csv",sep=""),sep=",",header=T, row.names = 1))
  diag(matrix) <- 0   # remove cannibalistic links
  igraph <- graph_from_adjacency_matrix(matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  attributes <- read.csv(paste("baltic_sea/baltic_spAttributes_",web[i],".csv",sep=""))
  attributes <- attributes[attributes$species %in% row.names(matrix), ] 
  attributes$bodymass <- attributes$bodymass
  attributes$abundance <- attributes$biomass/attributes$bodymass
  attributes$metabolic_type[attributes$metabolic_type=="invertebrate"] <- "ectotherm invertebrate"
  
  animals <- attributes$metabolic_type == "ectotherm invertebrate" | attributes$metabolic_type == "ectotherm vertebrate"
  plants <- attributes$metabolic_type == "primary producer"
  detritus <- attributes$metabolic_type == "detritus"
  
  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                              - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
  
  attributes$efficiencies[attributes$metabolic_type=="ectotherm invertebrate"|attributes$metabolic_type=="ectotherm vertebrate"] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="primary producer"] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="detritus"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(matrix)
  attributes$omnivory = sapply(1:nrow(matrix), Omnivory.species, matrix, attributes$TL, simplify = TRUE)
  attributes$omnivory[attributes$TL == 1] <- NA
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
   
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes 
  # Total flux
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  # Per capita flux
  PC.flux <- sweep(flux, 1, attributes$abundance, FUN = "/")
  
  ## Food web stability 
  # create vector that encodes for species types (animal, detritus or plant)
  met.types = rep('animal', nrow(matrix))
  met.types[plants] = 'plant'
  met.types[detritus] = 'detritus'
  # compute stability (leading eigenvalue of the Jacobian)
  stab = stability.value(val.mat = flux,
                         biomasses = attributes$biomass,
                         efficiencies = attributes$efficiencies,
                         metabolic.types = met.types,
                         ef.level = "prey"
  )
  
  # save stability value:
  meta.BS$stability[i] = stab

  ## Food web fluxes  
  meta.BS$flux[i] <- sum(flux)
  meta.BS$second.consumption[i] <- sum(flux[eat.animals])
  meta.BS$prim.consumption[i] <- sum(flux[basals,])  
  # Per capita flux
  meta.BS$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)

  ## Food web metrics
  meta.BS$S[i] <- Number.of.species(matrix)
  meta.BS$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
  meta.BS$MaxTL[i] <- max(attributes$TL)
  meta.BS$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.BS$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.BS$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
  
  meta.BS$L[i] <- Number.of.links(matrix)
  meta.BS$LD[i] <- Link.density(matrix)
  meta.BS$C[i] <- Connectance(matrix)
  meta.BS$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
}


commcols <- intersect(names(meta.RP), names(meta.BS))
meta.Marine <- bind_rows(select(meta.RP, all_of(commcols)),
                          select(meta.BS, all_of(commcols)))

write.csv(meta.Marine, file = 'meta.Marine.csv', row.names = F)
