rm(list=ls())
setwd("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\Data")
source("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\code\\Food_web_functions.r")
library(fluxweb); library(igraph); library(dplyr); library(cheddar); library(colorspace)
options(stringsAsFactors=FALSE)

## Put the food webs into FUSED format
boltz <- 0.00008617343
T0 <- 273.15 + 20   ## Centre temperatures around 20 degrees K across all datasets
meta.RP <- read.csv("Intertidalrockpools/Intertidalrockypools_metadata.csv")
#meta.RP <- meta.RP[19:26,] #Only run for plotting Mozambique food webs

meta.RP$study_ID <- rep('Intertidal rockpools', nrow(meta.RP))
perday <- 60*60*24
web <- unique(meta.RP$FW_name)


metadata = read.csv('Intertidalrockpools/Intertidalrockypools_metadata.csv')
#metadata <- metadata[19:26,] #Only run for plotting Mozambique food webs

for(i in 1:length(web)){
  
  temperature = metadata$temperature_C[metadata$FW_name == web[i]]
  temp.kT = temperature + 273.15
  temp.arr = ((273.15+temperature)-T0)/(boltz*(273.15+temperature)*T0)
  matrix <-as.matrix(read.csv(paste("Intertidalrockpools/Intertidalrockypools_matrix_",web[i],".csv",sep=""),sep=",",header=T))
  row.names(matrix) <- colnames(matrix)
  # diag(matrix) <- 0.01*diag(matrix)   ## Down-weight preference of cannibalistic links
  diag(matrix) <- 0   ## Make sure there are no cannibalistic links
  # matrix[matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard')) #calculates trophic Dissimilarity

  
  attributes <- read.csv(paste("Intertidalrockpools/Intertidalrockypools_spAttributes_",web[i],".csv",sep=""))
  attributes$bodymass <- attributes$bodymass
  attributes$abundance <- attributes$biomass/attributes$bodymass*4 # *4 converts dry weight to fresh weight
  attributes$metabolic_type[attributes$metabolic_type=="invertebrate"] <- "ectotherm invertebrate"
  
  animals <- attributes$metabolic_type != "primary producer"   ## "attributes" is just whatever you called your species attributes data frame
  plants <- attributes$metabolic_type == "primary producer"

  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.kT))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$efficiencies[animals] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[plants] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  #attributes$efficiencies[attributes$species_type=="detritus"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.kT))
  
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(matrix)
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass
  attributes$biomass[is.na(attributes$biomass)] <- attributes$bodymass[is.na(attributes$biomass)] ##Fill missing consumer biomass to = body mass (i.e. 1 individual)

  #SP1 doesn't run with raw biomass (magnitudes too large?)
  #This is fixed with log-transforming and rescaling to all positive values (min of 1)
  attributes$biomass[web[i]=='SP1'] <- (log(attributes$biomass) + abs(min(log(attributes$biomass)))) +1
  
  eat.plants = colSums(matrix[plants, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.plants & eat.animals
  
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  
  ############# stability ###############
  # I first need a vector that encodes for species types (animal, detritus or plant)
  # here I assume that the vectors created L97-99 contain this information
  met.types = rep('animal', nrow(matrix))
  met.types[plants] = 'plant'
  #met.types[detritus] = 'detritus'
  # Then I can compute the stability (leading eigenvalue of the Jacobian)
  stab = stability.value(val.mat = flux,
                         biomasses = attributes$biomass,
                         efficiencies = attributes$efficiencies,
                         metabolic.types = met.types,
                         ef.level = "prey"
  )
  
  # save stability value:
  meta.RP$stability[i] = stab
  ########################################

  meta.RP$flux[i] <- sum(flux)
  meta.RP$second.consumption[i] <- sum(flux[animals])
  meta.RP$prim.consumption[i] <- sum(flux[plants,])

  meta.RP$biomass[i] <- sum(attributes$biomass, na.rm=T)
  meta.RP$animal.biomass[i] <- sum(attributes$biomass[animals], na.rm=T)
  meta.RP$plant.biomass[i] <- sum(attributes$biomass[plants], na.rm=T)

  ## Food web metrics
  meta.RP$S[i] <- Number.of.species(matrix)
  meta.RP$L[i] <- Number.of.links(matrix)
  meta.RP$LD[i] <- Link.density(matrix)
  meta.RP$C[i] <- Connectance(matrix)
  #meta.RP$IG.predation[i] <- Connectance(matrix[animals,])
  meta.RP$generality[i] <- Gen.sd(matrix)
  meta.RP$vulnerability[i] <- Vul.sd(matrix[!basals, !basals])
  meta.RP$MeanTL[i] <- mean(attributes$TL)
  meta.RP$MaxTL[i] <- max(attributes$TL)
  #meta.RP$modularity[i] <- cluster_spinglass(igraph)$modularity
  meta.RP$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.RP$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.RP$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
  meta.RP$prim.cons.richness[i] = sum(primary.cons.and.omnivores) # primary consumer richness
  meta.RP$sec.cons.richness[i] = sum(sec.cons)  # secondary consumer richness
}




#### Fluxweb analysis Baltic Sea ####
meta.BS <- read.csv("baltic_sea/baltic_metadata.csv")
meta.BS$study_ID <- rep('Baltic Sea', nrow(meta.BS))
web <- unique(meta.BS$FW_name)

for(i in 1:length(web)){
  
  temperature = meta.BS$temperature_C[meta.BS$FW_name == web[i]]
  temp.kT = temperature + 273.15
  temp.arr = ((273.15+temperature)-T0)/(boltz*(273.15+temperature)*T0)
  matrix <-as.matrix(read.csv(paste("baltic_sea/baltic_matrix_",web[i],".csv",sep=""),sep=",",header=T, row.names = 1))
  # diag(matrix) <- 0.01*diag(matrix)   ## Down-weight preference of cannibalistic links
  diag(matrix) <- 0   ## Make sure there are no cannibalistic links
  # matrix[matrix>0] <- 1
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
                                                                                - 0.69/(boltz*(temp.kT))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                              - 0.69/(boltz*(temp.kT))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
  
  attributes$efficiencies[attributes$metabolic_type=="ectotherm invertebrate"|attributes$metabolic_type=="ectotherm vertebrate"] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="primary producer"] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[attributes$metabolic_type=="detritus"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
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
  flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  bioms.losses=F, bioms.prefs = T) * perday
  
  ############# stability ###############
  # I first need a vector that encodes for species types (animal, detritus or plant)
  # here I assume that the vectors created L97-99 contain this information
  met.types = rep('animal', nrow(matrix))
  met.types[plants] = 'plant'
  met.types[detritus] = 'detritus'
  # Then I can compute the stability (leading eigenvalue of the Jacobian)
  stab = stability.value(val.mat = flux,
                         biomasses = attributes$biomass,
                         efficiencies = attributes$efficiencies,
                         metabolic.types = met.types,
                         ef.level = "prey"
  )
  
  # save stability value:
  meta.BS$stability[i] = stab
  ########################################
  
  meta.BS$flux[i] <- sum(flux)
  meta.BS$second.consumption[i] <- sum(flux[eat.animals])
  meta.BS$prim.consumption[i] <- sum(flux[basals,])

  meta.BS$biomass[i] <- sum(attributes$biomass, na.rm=T)
  meta.BS$animal.biomass[i] <- sum(attributes$biomass[animals], na.rm=T)
  meta.BS$plant.biomass[i] <- sum(attributes$biomass[plants], na.rm=T)

  ## Food web metrics
  meta.BS$S[i] <- Number.of.species(matrix)
  meta.BS$L[i] <- Number.of.links(matrix)
  meta.BS$LD[i] <- Link.density(matrix)
  meta.BS$C[i] <- Connectance(matrix)
  #meta.BS$IG.predation[i] <- Connectance(matrix[animals,])
  meta.BS$generality[i] <- Gen.sd(matrix)
  meta.BS$vulnerability[i] <- Vul.sd(matrix[!basals, !basals])
  meta.BS$MeanTL[i] <- mean(attributes$TL)
  meta.BS$MaxTL[i] <- max(attributes$TL)
  #meta.BS$modularity[i] <- cluster_spinglass(igraph)$modularity
  meta.BS$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.BS$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.BS$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
  meta.BS$prim.cons.richness[i] = sum(primary.cons.and.omnivores) # primary consumer richness
  meta.BS$sec.cons.richness[i] = sum(sec.cons)  # secondary consumer richness
}


commcols <- intersect(names(meta.RP), names(meta.BS))
meta.Marine <- bind_rows(select(meta.RP, all_of(commcols)),
                          select(meta.BS, all_of(commcols)))


#### Calculate sample coverage ####
library(ggplot2)
# Get the full set of unique species across all files
files.RP <- list.files("Intertidalrockpools", pattern = "^Intertidalrockypools_spAttributes_", full.names = TRUE)
files.BS <- list.files("baltic_sea", pattern = "^baltic_spAttributes_", full.names = TRUE)
all_files <- c(files.RP, files.BS)

all_taxa <- unique(unlist(lapply(all_files, function(f) {
  dat <- read.csv(f, stringsAsFactors = FALSE)
  unique(dat$species)
})))

all_taxa <- sort(all_taxa)

# Create abundance matrix
abundance_matrix <- matrix(0, nrow = length(all_files), ncol = length(all_taxa))
colnames(abundance_matrix) <- all_taxa
base_names <- basename(all_files)
rownames(abundance_matrix) <- sub("^(Intertidalrockypools_spAttributes_|baltic_spAttributes_)(.*)\\.csv$", "\\2", base_names) 

# Loop over files to build rows
for (i in seq_along(all_files)) {
  community <- read.csv(all_files[i], stringsAsFactors = FALSE)
  sp <- as.character(community$species)
  abund <- community$biomass/community$bodymass
  abundance_matrix[i, sp] <- abund
}
abundance_matrix[is.na(abundance_matrix)] <- 0
species_counts <- rowSums(abundance_matrix > 0)

# create list for iNEXT
abund_list <- apply(abundance_matrix, 1, ceiling)
abund_list <- as.list(data.frame(abund_list))
names(abund_list) <- rownames(abundance_matrix)

coverage <- DataInfo(abund_list, datatype = "abundance")  #iNEXT
colnames(coverage)[colnames(coverage) == 'Assemblage']  <- "FW_name"

meta.Marine$SC <- coverage$SC[match(meta.Marine$FW_name, coverage$FW_name)]
write.csv(meta.Marine, file = 'meta.Marine.csv', row.names = F)
