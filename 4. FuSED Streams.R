rm(list=ls())
setwd("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\BEF-in-Food-Webs")
source("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\BEF-in-Food-Webs\\Food_web_functions.r")
library(fluxweb); library(igraph); library(dplyr); library(cheddar); library(colorspace)
options(stringsAsFactors=FALSE)


#### Fluxweb analysis Icelandic Streams ####

meta.IS <- read.csv("IcelandicStreams/IcelandicStreams_metadata.csv")
meta.IS$study_ID <- rep('Hengill valley', nrow(meta.IS))
perday <- 60*60*24
boltz <- 0.00008617343
T0 <- 273.15 + 20   
web <- unique(meta.IS$FW_name)

for(i in 1:length(web)){
  
  matrix <- bin.matrix <- as.matrix(read.csv(paste("IcelandicStreams/IcelandicStreams_matrix_",web[i],".csv",sep=""),header=T))
  diag(matrix) <- diag(bin.matrix) <- 0   ## Make sure there are no cannibalistic links
  bin.matrix[bin.matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(bin.matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  attributes <- read.csv(paste("IcelandicStreams/IcelandicStreams_spAttributes_",web[i],".csv",sep=""))
  attributes$losses[is.na(attributes$losses)] <- 0 # Metabolic losses already calculated as ln(I)=ln(i_o)+0.71*ln(M)-E(1/kT)
  attributes$TL <- TL(matrix)
  
  animals <- attributes$species_type == "animal"   
  plants <- attributes$species_type == "plant"
  detritus <- attributes$species_type == "detritus" 
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
  
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  
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
  
  ## save stability value:
  meta.IS$stability[i] = stab
 
  ## food web fluxes
  meta.IS$flux[i] <- sum(flux)
  meta.IS$second.consumption[i] <- sum(flux[attributes$species_type=="animal",])
  meta.IS$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
  
  ## Food web metrics
  meta.IS$S[i] <- Number.of.species(bin.matrix)
  meta.IS$MaxTL[i] <- max(attributes$TL)
  meta.IS$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.IS$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.IS$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
}


#### Fluxweb analysis UK Streams ####
meta.UK <- read.csv("UKstreams/UKstreams_metadata.csv")
meta.UK$study_ID <- rep('UK chalk streams', nrow(meta.UK))
web <- unique(meta.UK$FW_name)

for(i in 1:length(web)){
  
  temp.K = meta.UK$temperature[meta.UK$FW_name == web[i]] + 273.15
  temp.arr = ((temp.K)-T0)/(boltz*(temp.K)*T0)
  
  matrix <- bin.matrix <- as.matrix(read.csv(paste("UKstreams/UKstreams_matrix_",web[i],".csv",sep=""),header=T))
  diag(matrix) <- diag(bin.matrix) <- 0   # remove cannibalistic links
  bin.matrix[bin.matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(bin.matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  attributes <- read.csv(paste("UKstreams/UKstreams_spAttributes_",web[i],".csv",sep=""))
  attributes$TL <- TL(matrix)
  attributes$bodymass <- (attributes$bodymass*4) # Convert body mass dry to fresh mass 
  
  animals <- attributes$species_type == "animal"   
  plants <- attributes$species_type == "plant"
  detritus <- attributes$species_type == "detritus" 
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
  
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  primary.cons.and.omnivores = omnivores | prim.cons
  
  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                              - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
  attributes$losses[is.na(attributes$losses)] <- 0
  
  attributes$efficiencies[attributes$metabolic_type=="ectotherm invertebrate"|attributes$metabolic_type=="ectotherm vertebrate"] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="plant"] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="detritus"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  
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
  
  ## save stability value:
  meta.UK$stability[i] = stab
  
  ## food web fluxes
  meta.UK$flux[i] <- sum(flux)
  meta.UK$second.consumption[i] <- sum(flux[attributes$species_type=="animal",])
  meta.UK$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])

  ## Food web metrics
  meta.UK$S[i] <- Number.of.species(bin.matrix)
  meta.UK$MaxTL[i] <- max(attributes$TL)
  meta.UK$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.UK$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.UK$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
}


#### River data from Brauns, Hall, and Rene-Mor et al. #### 
meta.Brauns <- read.csv("Rivers_Mario/rivers_metadata.csv")
meta.Brauns$study_ID <- rep(c('Elbe','Hall et al','Montsant'),times=c(3,4,4))
web <- unique(meta.Brauns$FW_name)

for(i in 1:length(web)){
  
  temp.K = meta.Brauns$temperature + 273.15
  temp.arr = ((temp.K)-T0)/(boltz*(temp.K)*T0)
  matrix <-as.matrix(read.csv(paste("Rivers_Mario/rivers_matrix_",web[i],".csv",sep=""),sep=",",header=T, row.names = 1, check.names = F))
  matrix[is.na(matrix)] = 0
  diag(matrix) <- 0   # remove cannibalistic links
  bin.matrix <- matrix
  bin.matrix[matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(bin.matrix, mode="directed", )
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  attributes <- read.csv(paste("Rivers_Mario/rivers_spAttributes_",web[i],".csv",sep=""))
  attributes <- attributes[attributes$species %in% colnames(bin.matrix), ] 
  attributes %>% arrange(factor(species, levels = row.names(matrix)))
  attributes$bodymass <- (attributes$bodymass*4)/1000 #Convert body mass in mg to g and dry to fresh mass conversion
  attributes$abundance <- attributes$abundance
  
  animals <- attributes$metabolic_type == "ectotherm invertebrate"
  plants <- attributes$metabolic_type == "plant"| attributes$metabolic_type == "phytoplankton"
  detritus <- attributes$metabolic_type == "detritus"| attributes$metabolic_type == "fungi"
  
  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                              - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
  
  attributes$efficiencies[attributes$metabolic_type=="ectotherm invertebrate" | attributes$metabolic_type=="ectotherm vertebrate"] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="plant"| attributes$metabolic_type == "phytoplankton"] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="detritus" | attributes$metabolic_type=="fungi" | attributes$metabolic_type=="heterotrophic bacteria"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(matrix)
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
  
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, biomasses = attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  
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
  
  ## save stability value:
  meta.Brauns$stability[i] = stab

  ## food web fluxes
  meta.Brauns$flux[i] <- sum(flux)
  meta.Brauns$second.consumption[i] <- sum(flux[eat.animals])
  meta.Brauns$prim.consumption[i] <- sum(flux[basals,])

  ## Food web metrics
  meta.Brauns$S[i] <- Number.of.species(bin.matrix)
  meta.Brauns$MaxTL[i] <- max(attributes$TL)
  meta.Brauns$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.Brauns$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.Brauns$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
}




#### Fluxweb analysis Brazilian streams ####
meta.Brazil <- read.csv("Brazilian streams_Saito/CANANEIA_metadata.csv")
meta.Brazil$study_ID <- rep('Cananeia', nrow(meta.Brazil))
web <- unique(meta.Brazil$FW_name)

for(i in 1:length(web)){
  
  temp.K = meta.Brazil$temperature[meta.Brazil$FW_name == web[i]] + 273.15
  temp.arr = ((temp.K)-T0)/(boltz*(temp.K)*T0)
  
  matrix <- bin.matrix <- as.matrix(read.csv(paste("Brazilian streams_Saito/CANANEIA_matrix_",web[i],".csv",sep=""),header=T, row.names=1))
  diag(matrix) <- diag(bin.matrix) <- 0   ## Make sure there are no cannibalistic links
  bin.matrix[bin.matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(bin.matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  attributes <- read.csv(paste("Brazilian streams_Saito/CANANEIA_spAttributes_",web[i],".csv",sep=""))
  attributes$TL <- TL(matrix)
  attributes$bodymass <- (attributes$bodymass*4)/1000 #Convert body mass in mg to g and dry to fresh mass conversion
  
  animals <- attributes$metabolic_type == "ectotherm invertebrate" | attributes$metabolic_type == "ectotherm vertebrate"   
  plants <- attributes$metabolic_type == "plant"
  detritus <- attributes$metabolic_type == "detritus"
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
  
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  primary.cons.and.omnivores = omnivores | prim.cons
  
  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                               - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
  attributes$losses[is.na(attributes$losses)] <- 0
  
  attributes$efficiencies[attributes$metabolic_type=="ectotherm invertebrate"|attributes$metabolic_type=="ectotherm vertebrate"] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="plant"] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="detritus"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  
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
  
  ## save stability value:
  meta.Brazil$stability[i] = stab

  ## food web fluxes
  meta.Brazil$flux[i] <- sum(flux)
  meta.Brazil$second.consumption[i] <- sum(flux[eat.animals])
  meta.Brazil$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
 
  ## Food web metrics
  meta.Brazil$S[i] <- Number.of.species(bin.matrix)
  meta.Brazil$MaxTL[i] <- max(attributes$TL)
  meta.Brazil$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.Brazil$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.Brazil$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
}


commcols1 <- intersect(names(meta.IS), names(meta.UK))
commcols2 <- intersect(commcols1, names(meta.Brauns))
commcols <- intersect(commcols2, names(meta.Brazil))
meta.Streams <- bind_rows(select(meta.IS, all_of(commcols)),
                          select(meta.UK, all_of(commcols)),
                          select(meta.Brauns, all_of(commcols)),
                          select(meta.Brazil, all_of(commcols)))


# #### Calculate sample coverage using iNEXT ####
# library(iNEXT)
# library(ggplot2)
# 
# # Get the full set of unique species across all files
# files.IS <- list.files("IcelandicStreams", pattern = "^IcelandicStreams_spAttributes_", full.names = TRUE)
# files.UK <- list.files("UKstreams", pattern = "^UKstreams_spAttributes_", full.names = TRUE)
# files.Brauns <- list.files("Rivers_Mario", pattern = "^rivers_spAttributes_", full.names = TRUE)
# files.Brazil <- list.files("Brazilian streams_Saito", pattern = "^CANANEIA_spAttributes_", full.names = TRUE)
# all_files <- c(files.IS, files.UK, files.Brauns, files.Brazil)
# folder_labels <- c(
#   rep("IS", length(files.IS)),
#   rep("UK", length(files.UK)),
#   rep("Brauns", length(files.Brauns)),
#   rep("Brazil", length(files.Brazil))
# )
# 
# all_taxa <- unique(unlist(lapply(all_files, function(f) {
#   dat <- read.csv(f, stringsAsFactors = FALSE)
#   unique(dat$species)
# })))
# 
# all_taxa <- sort(all_taxa)
# 
# # Create abundance matrix
# abundance_matrix <- matrix(0, nrow = length(all_files), ncol = length(all_taxa))
# colnames(abundance_matrix) <- all_taxa
# base_names <- basename(all_files)
# rownames(abundance_matrix) <- sub(
#   "^(IcelandicStreams_spAttributes_|UKstreams_spAttributes_|rivers_spAttributes_|CANANEIA_spAttributes_)(.*)\\.csv$", "\\2", 
#   base_names)
# 
# # Loop over files to build rows
# for (i in seq_along(all_files)) {
#   community <- read.csv(all_files[i], stringsAsFactors = FALSE)
#   community <- community[!is.na(community$abundance), ]
#   sp <- as.character(community$species)
#   abund <- community$abundance
#   abundance_matrix[i, sp] <- abund
# }
# abundance_matrix[is.na(abundance_matrix)] <- 0
# species_counts <- rowSums(abundance_matrix > 0)
# 
# 
# # create list for iNEXT
# abund_list <- apply(abundance_matrix, 1, ceiling)
# abund_list <- as.list(data.frame(abund_list))
# names(abund_list) <- rownames(abundance_matrix)
# 
# #Calculate sample coverage
# coverage <- DataInfo(abund_list, datatype = "abundance")  #iNEXT
# colnames(coverage)[colnames(coverage) == 'Assemblage']  <- "FW_name"
# 
# meta.Streams$SC <- coverage$SC[match(meta.Streams$FW_name, coverage$FW_name)]
write.csv(meta.Streams, file = 'meta.Streams.csv', row.names = F)
