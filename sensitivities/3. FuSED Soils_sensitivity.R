#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Calculation of energy fluxes in 65 soil food webs
#
#=========================================================================================================================


flux.soils = function(params){
  
  # access specific parameter combination
  repl = params[1]
  X.exp = params[2]
  X.temp = params[3]
  eff.exp.inv = params[4]
  eff.exp.prod = params[5]
  eff.exp.det = params[6]
  eff.temp.inv = params[7]
  eff.temp.prod = params[8]
  eff.temp.det = params[9]
  
  #### Biodiversity Exploratories ####
  meta.BioExp <- read.csv("Biodiversity Exploratories/explo_metadata_copy.csv")
  meta.BioExp$study_ID <- rep('Biodiversity exploratories', nrow(meta.BioExp))
  web <- unique(meta.BioExp$FW_name)
  
  
  for(i in 1:length(web)){
    
    temp.K = meta.BioExp$temperature[meta.BioExp$FW_name == web[i]] + 273.15
    temp.arr = ((temp.K)-T0)/(boltz*(temp.K)*T0)
    matrix <-as.matrix(read.csv(paste("Biodiversity Exploratories/explo_matrix_",web[i],".csv",sep=""),sep=",",header=T))
    row.names(matrix) <- colnames(matrix)
    matrix[,c('bacteria', 'mykorrhiza', 'saphrophytes')] = 0
    row.names(matrix) <- colnames(matrix)
    diag(matrix) <- 0   # remove cannibalistic links
    
    attributes <- read.csv(paste("Biodiversity Exploratories/explo_spAttributes_",web[i],".csv",sep=""))
    attributes$species <- make.names(attributes$species)
    attributes$bodymass <- attributes$bodymass/1000 # Convert body mass in mg to g 
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
    
    animals <- attributes$metabolic_type == "ectotherm invertebrate"   
    plants <- attributes$metabolic_type == "primary producer"
    detritus <- attributes$metabolic_type == "heterotrophic fungi" | attributes$metabolic_type == "heterotrophic bacteria" | 
      attributes$metabolic_type == "dead organic material"
    
    attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((X.exp * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                  - X.temp/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
    attributes$efficiencies[animals] <- exp(eff.exp.inv)*exp(eff.temp.inv*temp.arr) / (1 + exp(eff.exp.inv)*exp(eff.temp.inv*temp.arr)) #animal
    attributes$efficiencies[plants] <- exp(eff.exp.prod)*exp(eff.temp.inv*temp.arr) / (1 + exp(eff.exp.prod)*exp(eff.temp.inv*temp.arr)) # plant
    attributes$efficiencies[detritus] <- exp(eff.exp.det)*exp(eff.temp.inv*temp.arr) / (1 + exp(eff.exp.det)*exp(eff.temp.inv*temp.arr)) #detritus
    
    
    attributes$losses[is.na(attributes$losses)] <- 0
    attributes$TL <- TL(matrix)
    attributes$omnivory = sapply(1:nrow(matrix), Omnivory.species, matrix, attributes$TL, simplify = TRUE)
    attributes$omnivory[attributes$TL == 1] <- NA
    
    prim.cons = attributes$TL == 2 # primary consumers
    basals = attributes$TL == 1 # basal species
    sec.cons = attributes$TL > 2  # secondary consumers
    
    eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
    eat.animals = colSums(matrix[animals,]) > 0
    omnivores = eat.basals & eat.animals
    
    primary.cons.and.omnivores = omnivores | prim.cons
    
    ## Calculating fluxes
    # Total flux
    flux <- fluxing(mat=matrix, losses=attributes$losses, efficiencies=attributes$efficiencies, bioms.losses=F, bioms.prefs = T,
                    biomass = attributes$biomass) * perday
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
    
    ## save stability value:
    meta.BioExp$stability[i] = stab
    
    ## Food web fluxes
    meta.BioExp$flux[i] <- sum(flux)
    meta.BioExp$second.consumption[i] <- sum(flux[animals])
    meta.BioExp$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
    # Per capita flux
    meta.BioExp$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta.BioExp$S[i] <- Number.of.species(matrix)
    meta.BioExp$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
    meta.BioExp$MaxTL[i] <- max(attributes$TL)
    meta.BioExp$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta.BioExp$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta.BioExp$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web  
    
    meta.BioExp$L[i] <- Number.of.links(matrix)
    meta.BioExp$LD[i] <- Link.density(matrix)
    meta.BioExp$C[i] <- Connectance(matrix)
    meta.BioExp$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta.BioExp$modularity[i] = modularity(igraph, membership(community))
  }
  
  
  
  
  #### Russian forest soils ####
  meta.Russian = read.csv('Russian soils/Russianforests_metadata.csv')
  meta.Russian$study_ID <- rep('Russian soils', nrow(meta.Russian))
  web <- unique(meta.Russian$FW_name)
  
  for(i in 1:length(web)){
    
    temp.K = meta.Russian$temperature[meta.Russian$FW_name == web[i]] + 273.15
    temp.arr = ((temp.K)-T0)/(boltz*(temp.K)*T0)
    matrix <-as.matrix(read.csv(paste("Russian soils/Russianforests_matrix_",web[i],".csv",sep=""),sep=",",header=T))
    row.names(matrix) <- colnames(matrix)
    matrix <- t(matrix)
    diag(matrix) <- 0   # remove cannibalistic links
    
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
                                                                                  - 0.69/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
    attributes$efficiencies[animals] <- exp(2.266)*exp(0.164*temp.arr) / (1 + ((exp(2.266)*exp(0.164*temp.arr))))
    attributes$efficiencies[plants] <-  exp(0.179)*exp(0.164*temp.arr) / (1 + exp(0.179)*exp(0.164*temp.arr))
    attributes$efficiencies[detritus] <- exp(-1.670)*exp(0.164*temp.arr) / (1 + exp(-1.670)*exp(0.164*temp.arr))
    
    attributes$losses[is.na(attributes$losses)] <- 0
    attributes$TL <- TL(matrix)
    attributes$omnivory = sapply(1:nrow(matrix), Omnivory.species, matrix, attributes$TL, simplify = TRUE)
    attributes$omnivory[attributes$TL == 1] <- NA
    
    prim.cons = attributes$TL == 2 # primary consumers
    basals = attributes$TL == 1 # basal species
    sec.cons = attributes$TL > 2  # secondary consumers
    
    eat.basals = colSums(matrix[basals, , drop = FALSE])  > 0
    eat.animals = colSums(matrix[animals,]) > 0
    omnivores = eat.basals & eat.animals
    
    primary.cons.and.omnivores = omnivores | prim.cons
    
    ## Calculating fluxes
    # Total flux
    flux <- fluxing(mat=matrix, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                    biomass = attributes$biomass, bioms.losses=F, bioms.prefs = T) * perday
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
    
    ## save stability value:
    meta.Russian$stability[i] = stab
    
    ## food web fluxes
    meta.Russian$flux[i] <- sum(flux)
    meta.Russian$second.consumption[i] <- sum(flux[animals])
    meta.Russian$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
    # Per capita flux
    meta.Russian$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta.Russian$S[i] <- Number.of.species(matrix)
    meta.Russian$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
    meta.Russian$MaxTL[i] <- max(attributes$TL)
    meta.Russian$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta.Russian$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta.Russian$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
    
    meta.Russian$L[i] <- Number.of.links(matrix)
    meta.Russian$LD[i] <- Link.density(matrix)
    meta.Russian$C[i] <- Connectance(matrix)
    meta.Russian$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta.Russian$modularity[i] = modularity(igraph, membership(community))
  }
  
  commcols <- intersect(names(meta.BioExp), names(meta.Russian))
  meta.Soils <- bind_rows(select(meta.BioExp, all_of(commcols)),
                          select(meta.Russian, all_of(commcols)))
  
  meta.Soils$replicate = repl
  
  return(meta.Soils)
}
