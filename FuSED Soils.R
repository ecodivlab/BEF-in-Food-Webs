rm(list=ls())
setwd("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\Data")
source("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\code\\Food_web_functions.r")
library(fluxweb); library(igraph); library(dplyr); library(cheddar); library(colorspace)
options(stringsAsFactors=FALSE)

## Put the food webs into FUSED format
boltz <- 0.00008617343
T0 <- 273.15 + 20   ## Use 20 degrees as standardised temperature across all dataset
meta.BioExp <- read.csv("Biodiversity Exploratories/explo_metadata_copy.csv")
#meta.BioExp <- meta.BioExp[1:29,] #Only run for plotting first 29 food webs
meta.BioExp$study_ID <- rep('Biodiversity exploratories', nrow(meta.BioExp))
perday <- 60*60*24
web <- unique(meta.BioExp$FW_name)

#### Biodiversity Exploratories
metadata = read.csv('Biodiversity Exploratories/explo_metadata_copy.csv')

for(i in 1:length(web)){
  
  temperature = metadata$temperature[metadata$FW_name == web[i]]
  temp.kT = temperature + 273.15
  temp.arr = ((273.15+temperature)-T0)/(boltz*(273.15+temperature)*T0)
  matrix <-as.matrix(read.csv(paste("Biodiversity Exploratories/explo_matrix_",web[i],".csv",sep=""),sep=",",header=T))
  row.names(matrix) <- colnames(matrix)
  matrix[,c('bacteria', 'mykorrhiza', 'saphrophytes')] = 0
  row.names(matrix) <- colnames(matrix)
  # diag(matrix) <- 0.01*diag(matrix)   ## Down-weight preference of cannibalistic links
  diag(matrix) <- 0   ## Make sure there are no cannibalistic links
  # matrix[matrix>0] <- 1

  attributes <- read.csv(paste("Biodiversity Exploratories/explo_spAttributes_",web[i],".csv",sep=""))
  attributes$species <- make.names(attributes$species)
  attributes$bodymass <- attributes$bodymass/1000 #Convert body mass in mg to g 
  attributes$abundance <- attributes$biomass/attributes$bodymass
  attributes$metabolic_type[attributes$metabolic_type=="invertebrate"] <- "ectotherm invertebrate"
  attributes$biomass[attributes$metabolic_type=="primary producer"|attributes$metabolic_type=="dead organic material"|
                     attributes$metabolic_type=="heterotrophic bacteria"|attributes$metabolic_type=="heterotrophic fungi"
                     ] <- mean(attributes$biomass, na.rm=TRUE) #Set biomass values for basal resources to mean of other groups
  
  to.remove <- attributes$species[is.na(attributes$biomass)]
  attributes <- attributes[!(attributes$species %in% to.remove), ] #Remove consumers with no biomass
  matrix <- matrix[!(row.names(matrix) %in% to.remove), !(colnames(matrix) %in% to.remove) ] #Remove consumers with no biomass
  
  igraph <- graph_from_adjacency_matrix(matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  #Now need to remove this list from the matrix and attributes! 
  
  animals <- attributes$metabolic_type == "ectotherm invertebrate"   
  plants <- attributes$metabolic_type == "primary producer"
  detritus <- attributes$metabolic_type == "heterotrophic fungi" | attributes$metabolic_type == "heterotrophic bacteria" | 
              attributes$metabolic_type == "dead organic material"
  
  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.kT))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$efficiencies[animals] <- exp(2.266)*exp(0.164*temp.arr) / (1 + ((exp(2.266)*exp(0.164*temp.arr))))
  attributes$efficiencies[plants] <-  exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[detritus] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(matrix)
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  #eat.plants = colSums(matrix[plants, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, losses=attributes$losses, efficiencies=attributes$efficiencies, bioms.losses=F, bioms.prefs = T,
                  biomass = attributes$biomass) * perday
  
  ############# stability ###############
  # I first need a vector that encodes for species types (animal, detritus or plant)
  # here I assume that the vectors created L97-99 contain this information
  met.types = rep('animal', nrow(matrix))
  met.types[plants] = 'plant'
  met.types[detritus] = 'detritus'
  # Then I can compute the stability
  stab = stability.value(val.mat = flux,
                         biomasses = attributes$biomass,
                         efficiencies = attributes$efficiencies,
                         metabolic.types = met.types,
                         ef.level = "prey"
  )
  
  # save stability value:
  meta.BioExp$stability[i] = stab
  ########################################
  
  meta.BioExp$flux[i] <- sum(flux)
  meta.BioExp$second.consumption[i] <- sum(flux[animals])
  meta.BioExp$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])

  meta.BioExp$biomass[i] <- sum(attributes$biomass, na.rm=T)
  meta.BioExp$animal.biomass[i] <- sum(attributes$biomass[animals], na.rm=T)
  meta.BioExp$plant.biomass[i] <- sum(attributes$biomass[plants], na.rm=T)

  ## Food web metrics
  meta.BioExp$S[i] <- Number.of.species(matrix)
  meta.BioExp$L[i] <- Number.of.links(matrix)
  meta.BioExp$LD[i] <- Link.density(matrix)
  meta.BioExp$C[i] <- Connectance(matrix)
  meta.BioExp$generality[i] <- Gen.sd(matrix)
  meta.BioExp$vulnerability[i] <- Vul.sd(matrix)
  meta.BioExp$MeanTL[i] <- mean(attributes$TL)
  meta.BioExp$MaxTL[i] <- max(attributes$TL)
  #meta.BioExp$modularity[i] <- cluster_spinglass(igraph)$modularity
  meta.BioExp$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.BioExp$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.BioExp$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
  meta.BioExp$prim.cons.richness[i] = sum(primary.cons.and.omnivores) # primary consumer richness
  meta.BioExp$sec.cons.richness[i] = sum(sec.cons)  # secondary consumer richness
}




#### Russian forest soils
meta.Russian = read.csv('Russian soils/Russianforests_metadata.csv')
meta.Russian$study_ID <- rep('Russian soils', nrow(meta.Russian))
web <- unique(meta.Russian$FW_name)

for(i in 1:length(web)){
  
  temperature = meta.Russian$temperature[meta.Russian$FW_name == web[i]]
  temp.kT = temperature + 273.15
  temp.arr = ((273.15+temperature)-T0)/(boltz*(273.15+temperature)*T0)
  matrix <-as.matrix(read.csv(paste("Russian soils/Russianforests_matrix_",web[i],".csv",sep=""),sep=",",header=T))
  row.names(matrix) <- colnames(matrix)
  # diag(matrix) <- 0.01*diag(matrix)   ## Down-weight preference of cannibalistic links
  matrix <- t(matrix)
  diag(matrix) <- 0   ## Make sure there are no cannibalistic links
  matrix[matrix>0] <- 1
  matrix[matrix<1] <- 0
  
  attributes <- read.csv(paste("Russian soils/Russianforests_spAttributes_",web[i],".csv",sep=""))
  attributes$species <- make.names(attributes$taxon)
  attributes$bodymass <- (attributes$bodymass*4) #Convert body mass dry to fresh mass 
  attributes$abundance <- attributes$abundance
  attributes$metabolic_type[attributes$metabolic_type=="invertebrate"] <- "ectotherm invertebrate"
  attributes$biomass[attributes$metabolic_type=="primary producer"|attributes$metabolic_type=="dead organic material"|
                       attributes$metabolic_type=="heterotrophic bacteria"|attributes$metabolic_type=="heterotrophic fungi"
  ] <- mean(attributes$biomass, na.rm=TRUE) #Set biomass values for basal resources to mean of other groups
  
  to.remove <- attributes$species[is.na(attributes$biomass)]
  attributes <- attributes[!(attributes$species %in% to.remove), ] #Remove consumers with no biomass
  matrix <- matrix[!(row.names(matrix) %in% to.remove), !(colnames(matrix) %in% to.remove) ] #Remove consumers with no biomass
  
  igraph <- graph_from_adjacency_matrix(matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  
  animals <- attributes$metabolic_type == "ectotherm invertebrate"   
  plants <- attributes$metabolic_type == "primary producer"
  detritus <- attributes$metabolic_type == "heterotrophic fungi" | attributes$metabolic_type == "heterotrophic bacteria" | 
    attributes$metabolic_type == "dead organic material"
  
  attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((0.71 * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                - 0.69/(boltz*(temp.kT))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
  attributes$efficiencies[animals] <- exp(2.266)*exp(0.164*temp.arr) / (1 + ((exp(2.266)*exp(0.164*temp.arr))))
  attributes$efficiencies[plants] <-  exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
  attributes$efficiencies[detritus] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
  
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(matrix)
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  #eat.plants = colSums(matrix[plants, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  
  primary.cons.and.omnivores = omnivores | prim.cons
  
  ## Calculating fluxes
  flux <- fluxing(mat=matrix, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                  biomass = attributes$biomass, bioms.losses=F, bioms.prefs = T) * perday
  
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
  meta.Russian$stability[i] = stab
  ########################################
  
  meta.Russian$flux[i] <- sum(flux)
  meta.Russian$second.consumption[i] <- sum(flux[animals])
  meta.Russian$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
  
  meta.Russian$biomass[i] <- sum(attributes$biomass, na.rm=T)
  meta.Russian$animal.biomass[i] <- sum(attributes$biomass[animals], na.rm=T)
  meta.Russian$plant.biomass[i] <- sum(attributes$biomass[plants], na.rm=T)

  ## Food web metrics
  meta.Russian$S[i] <- Number.of.species(matrix)
  meta.Russian$L[i] <- Number.of.links(matrix)
  meta.Russian$LD[i] <- Link.density(matrix)
  meta.Russian$C[i] <- Connectance(matrix)
  meta.Russian$generality[i] <- Gen.sd(matrix)
  meta.Russian$vulnerability[i] <- Vul.sd(matrix)
  meta.Russian$MeanTL[i] <- mean(attributes$TL)
  meta.Russian$MaxTL[i] <- max(attributes$TL)
  #meta.Russian$modularity[i] <- cluster_spinglass(igraph)$modularity
  meta.Russian$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta.Russian$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta.Russian$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
  meta.Russian$prim.cons.richness[i] = sum(primary.cons.and.omnivores) # primary consumer richness
  meta.Russian$sec.cons.richness[i] = sum(sec.cons)  # secondary consumer richness
}

commcols <- intersect(names(meta.BioExp), names(meta.Russian))
meta.Soils <- bind_rows(select(meta.BioExp, all_of(commcols)),
                         select(meta.Russian, all_of(commcols)))


#### Calculate sample coverage using iNEXT ####
library(iNEXT)
library(ggplot2)

# Get the full set of unique species across all files
files.BioExp <- list.files("Biodiversity Exploratories", pattern = "^explo_spAttributes_", full.names = TRUE)
files.Russian <- list.files("Russian soils", pattern = "^Russianforests_spAttributes_", full.names = TRUE)
all_files <- c(files.BioExp, files.Russian)
folder_labels <- c(
  rep("BioExp", length(files.BioExp)),
  rep("Russian", length(files.Russian))
)

all_taxa <- unique(unlist(mapply(function(f, label) {
  dat <- read.csv(f, stringsAsFactors = FALSE)
  if (label == "BioExp") {
    unique(dat$species)
  } else if (label == "Russian") {
    unique(dat$taxon)
  }
},
f = all_files,
label = folder_labels,
SIMPLIFY = FALSE)
))

all_taxa <- sort(all_taxa)

# Create abundance matrix
abundance_matrix <- matrix(0, nrow = length(all_files), ncol = length(all_taxa))
colnames(abundance_matrix) <- all_taxa
base_names <- basename(all_files)
rownames(abundance_matrix) <- sub("^(explo_spAttributes_|Russianforests_spAttributes_)(.*)\\.csv$", "\\2", base_names) #need to fix naming of baltic sea sites...

# Loop over files to build rows
for (i in seq_along(all_files)) {
  community <- read.csv(all_files[i], stringsAsFactors = FALSE)
  community <- community[!is.na(community$biomass) & !is.na(community$bodymass) & community$bodymass > 0, ]
  if ("species" %in% names(community)) {sp <- as.character(community$species)} 
  else if ("taxon" %in% names(community)) {sp <- as.character(community$taxon)} 
  else {stop(paste("No 'species' or 'taxon' column in file:", all_files[i]))}
  abund <- community$biomass/community$bodymass
  abundance_matrix[i, sp] <- abund
}
abundance_matrix[is.na(abundance_matrix)] <- 0
species_counts <- rowSums(abundance_matrix > 0)


# create list for iNEXT
abund_list <- apply(abundance_matrix, 1, ceiling)
abund_list <- as.list(data.frame(abund_list))
names(abund_list) <- rownames(abundance_matrix)

#Calculate sample coverage
coverage <- DataInfo(abund_list, datatype = "abundance")  #iNEXT
colnames(coverage)[colnames(coverage) == 'Assemblage']  <- "FW_name"

meta.Soils$SC <- coverage$SC[match(meta.Soils$FW_name, coverage$FW_name)]
write.csv(meta.Soils, file = 'meta.Soils.csv', row.names = F)
