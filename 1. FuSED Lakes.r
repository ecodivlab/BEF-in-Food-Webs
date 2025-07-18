rm(list=ls())
setwd("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\BEF-in-Food-Webs")
source("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\BEF-in-Food-Webs\\Food_web_functions.r")
library(fluxweb); library(cheddar); library(igraph); library(RColorBrewer); library(colorspace); library(sp); library(dplyr); library(tibble)
options(stringsAsFactors=FALSE)

#### Adirondack and Maggiore lakes ####

boltz <- 0.00008617343
T0 <- 273.15 + 20   
perday <- 60*60*24
meta <- read.csv("Lakes/Lakes_metadata.csv")
lakes <- unique(meta$FW_name)


for(i in 1:length(lakes)) {
  
  # imports
    lake <- NULL; nodes <- NULL; matrix<-NULL
    lake <- lakes[i]
    nodes <- read.csv(paste0("Lakes/Lakes_spAttributes_",lake,".csv"))
     
    # corrections to match database
    if (meta$study_ID[i]=='Adirondack lakes'){
    nodes$taxonomy[nodes$taxonomy=="Salmo trutta"] <- "Salmo rutta" # should be trutta but match to incorrect edgelist
    nodes$taxonomy[nodes$taxonomy=="Anabaena flos-aquae"] <- "Anabaena flos.aquae"
    }
      
    matrix <- as.matrix(read.csv(paste0("Lakes/Lakes_matrix_",lake,".csv"), row.names = 1))

  # check order
    nodes <- nodes[order(nodes$taxonomy),]
    matrix <- matrix[order(rownames(matrix)),order(colnames(matrix))]
    
  # remove taxa that do not have matches, and with 0 density
    unmatched.rows <- NULL; fishy <- NULL; detr <- NULL
    unmatched.rows <- which(is.na(nodes$Count.of.Lake_species)) # which are missing? Could not be matched with density from report
    eggs <- which(grepl("fish eggs", nodes$taxonomy)) # but keep fish eggs
    detr <- which(grepl("detritus", nodes$taxonomy))  # but keep detritus
    nauplii <- which(grepl("nauplii", nodes$taxonomy)) 
    unmatched.rows <- unmatched.rows[!(unmatched.rows %in% c(detr))] ## <--- add 'eggs' here to keep eggs

    if(length(unmatched.rows)>0) {
         nodes <- nodes[-unmatched.rows,]
         matrix <- matrix[-unmatched.rows,-unmatched.rows]
    }
    
  # temperature
    temp.K <- meta$temperature_C[meta$FW_name==lake] + 273.15
    temp.arr <- ((temp.K)-T0)/(boltz*(temp.K)*T0)
    
  # losses
    nodes$losses <- 0
    nodes$losses[nodes$metabolic.type=="invertebrate"] <- exp((0.71 * log(nodes$mass.mean.g.[nodes$metabolic.type=="invertebrate"]) + 17.17) - 0.69/(boltz*(temp.K))) * nodes$density.comb[nodes$metabolic.type=="invertebrate"]
    nodes$losses[nodes$metabolic.type=="ectotherm vertebrate"] <- exp((0.71 * log(nodes$mass.mean.g.[nodes$metabolic.type=="ectotherm vertebrate"]) + 18.47) - 0.69/(boltz*(temp.K))) * nodes$density.comb[nodes$metabolic.type=="ectotherm vertebrate"]
    
        ## Losses for fish eggs?? - unknown density - currently these are omitted anyway
        nodes$losses[is.na(nodes$losses)] <- 0
  
  # efficiencies
    nodes$efficiencies[nodes$metabolic.type=="invertebrate" | nodes$metabolic.type=="ectotherm vertebrate"] <- exp(2.266)*exp(0.164*temp.arr) / (1 + exp(2.266)*exp(0.164*temp.arr)) #animal
    nodes$efficiencies[nodes$metabolic.type=="primary producer"] <- exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr)) # plant
    nodes$efficiencies[nodes$metabolic.type=="detritus"] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr)) #detritus

    # preferences (biomass matrix)
    B.matrix <- matrix * nodes$biomass
    # missing detritus node given sum of phytoplankton biomass for each predator (i.e. equal weight to green/brown origin)
    missing.detritus <- nodes$metabolic.type=="detritus"
    B.matrix[missing.detritus,] <- rep(colSums(B.matrix[nodes$metabolic.type=="primary producer",], na.rm=T) / sum(missing.detritus), each=sum(missing.detritus))
    B.matrix[matrix!=0 & B.matrix==0] <- 1 # for when detritus is the only food source


  ## Fluxweb analysis
  diag(matrix) <- 0   # remove cannibalistic links
  bin.matrix <- matrix
  bin.matrix[bin.matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(bin.matrix, mode="directed")
  sim.mat = 1-(similarity(igraph, mode = 'in', method = 'jaccard'))
  attributes <- nodes
  attributes$losses[is.na(attributes$losses)] <- 0
  attributes$TL <- TL(bin.matrix)[,1]
  
  animals <- attributes$metabolic.type == "invertebrate" | attributes$metabolic.type == "ectotherm vertebrate" 
  plants <- attributes$metabolic.type == "primary producer"
  detritus <- attributes$metabolic.type == "detritus" 
  
  prim.cons = attributes$TL == 2 # primary consumers
  basals = attributes$TL == 1 # basal species
  sec.cons = attributes$TL > 2  # secondary consumers
  
  attributes$biomass[attributes$TL==1] <- mean(attributes$biomass[attributes$TL>1], na.rm=T) ## Set basal resource biomass to equal mean consumer biomass

  eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
  eat.animals = colSums(matrix[animals,]) > 0
  omnivores = eat.basals & eat.animals
  primary.cons.and.omnivores = omnivores | prim.cons

  ## Calculating fluxes
  if(meta$study_ID[i]=='Lake Maggiore'){
  flux <- fluxing(mat=bin.matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies,
                  bioms.losses=F, bioms.prefs = T) * perday}
  else if (meta$study_ID[i]=='Adirondack lakes'){
  flux <- fluxing(mat=bin.matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies,
                  bioms.losses=F, bioms.prefs = T) * perday * 1000000 #*1000000 to convert from ml to cubic meters for Adirondack Lakes
  }
  
  ## Food web stability 
  # create vector that encodes for species types (animal, detritus or plant)
  met.types = rep('animal', nrow(bin.matrix))
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
  meta$stability[i] = stab

  ## Food web fluxes
  meta$flux[i] <- sum(flux)
  meta$second.consumption[i] <- sum(flux[animals])
  meta$second.consumption[meta$second.consumption==0] <- NA
  meta$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])

  ## Food web metrics
  meta$S[i] <- Number.of.species(bin.matrix)
  meta$MaxTL[i] <- max(TL(bin.matrix))
  meta$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
  meta$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
  meta$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
}
#view(meta[,-c(2:32, 34)])

meta.Lakes <- meta[,-c(2:32, 34)]

meta.Lakes <- filter(meta.Lakes, second.consumption > 0) #remove food webs with no secondary consumers (Brooke Trout and South Lake. 


#### Calculate sample coverage using iNEXT #### 
library(iNEXT)
library(ggplot2)

# # Get the full set of unique species across all files
# files <- list.files("Lakes", pattern = "^Lakes_spAttributes_", full.names = TRUE)
# 
# all_taxa <- unique(unlist(lapply(files, function(f) {
#   dat <- read.csv(f, stringsAsFactors = FALSE)
#   unique(dat$taxonomy)
# })))
# 
# all_taxa <- sort(all_taxa)
# 
# # Create abundance matrix
# abundance_matrix <- matrix(0, nrow = length(files), ncol = length(all_taxa))
# colnames(abundance_matrix) <- all_taxa
# rownames(abundance_matrix) <- sub("^Lakes_spAttributes_(.*)\\.csv$", "\\1", basename(files))
# 
# # Loop over files to build rows
# for (i in seq_along(files)) {
#   community <- read.csv(files[i], stringsAsFactors = FALSE)
#   sp <- as.character(community$taxonomy)
#   abund <- community$density.comb
#   abundance_matrix[i, sp] <- abund
# }
# abundance_matrix <- abundance_matrix[!row.names(abundance_matrix) %in% c("Brook trout lake", "South Lake"),]
# abundance_matrix[is.na(abundance_matrix)] <- 0
# species_counts <- rowSums(abundance_matrix > 0)
# 
# # create list for iNEXT
# abund_list <- apply(abundance_matrix, 1, ceiling)
# abund_list <- as.list(data.frame(abund_list))
# names(abund_list) <- rownames(abundance_matrix)
# 
# coverage <- DataInfo(abund_list, datatype = "abundance")  #iNEXT
# colnames(coverage)[colnames(coverage) == 'Assemblage']  <- "FW_name"
# 
# meta.Lakes$SC <- coverage$SC[match(meta.Lakes$FW_name, coverage$FW_name)]
write.csv(meta.Lakes, file = 'meta.Lakes.csv', row.names = F)


