#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Calculation of energy fluxes in 74 stream food webs
#
#=========================================================================================================================


flux.streams = function(params){
  
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
  
  
  #### Fluxweb analysis Icelandic Streams ####
  meta.IS <- read.csv("IcelandicStreams/IcelandicStreams_metadata.csv")
  meta.IS$study_ID <- rep('Hengill valley', nrow(meta.IS))
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
    attributes$omnivory = sapply(1:nrow(bin.matrix), Omnivory.species, bin.matrix, attributes$TL, simplify = TRUE)
    attributes$omnivory[attributes$TL == 1] <- NA
    
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
    # Total fluxes
    flux <- fluxing(mat=matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
                    bioms.losses=F, bioms.prefs = T) * perday
    # Per capita flux
    PC.flux <- sweep(flux, 1, (attributes$biomass/attributes$bodymass), FUN = "/")
    
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
    meta.IS$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta.IS$S[i] <- Number.of.species(bin.matrix)
    meta.IS$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
    meta.IS$MaxTL[i] <- max(attributes$TL)
    meta.IS$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta.IS$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta.IS$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
    
    meta.IS$L[i] <- Number.of.links(bin.matrix)
    meta.IS$LD[i] <- Link.density(bin.matrix)
    meta.IS$C[i] <- Connectance(bin.matrix)
    meta.IS$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta.IS$modularity[i] = modularity(igraph, membership(community))
    
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
    attributes$omnivory = sapply(1:nrow(bin.matrix), Omnivory.species, bin.matrix, attributes$TL, simplify = TRUE)
    attributes$omnivory[attributes$TL == 1] <- NA
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
    
    attributes$losses[attributes$metabolic_type=="ectotherm invertebrate"] <- exp((X.exp * log(attributes$bodymass[attributes$metabolic_type=="ectotherm invertebrate"]) + 17.17) 
                                                                                  - X.temp/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm invertebrate"]
    attributes$losses[attributes$metabolic_type=="ectotherm vertebrate"] <- exp((X.exp * log(attributes$bodymass[attributes$metabolic_type=="ectotherm vertebrate"]) + 18.47) 
                                                                                - X.temp/(boltz*(temp.K))) * attributes$abundance[attributes$metabolic_type=="ectotherm vertebrate"]
    attributes$losses[is.na(attributes$losses)] <- 0
    
    attributes$efficiencies[animals] <- exp(eff.exp.inv)*exp(eff.temp.inv*temp.arr) / (1 + exp(eff.exp.inv)*exp(eff.temp.inv*temp.arr)) #animal
    attributes$efficiencies[plants] <- exp(eff.exp.prod)*exp(eff.temp.inv*temp.arr) / (1 + exp(eff.exp.prod)*exp(eff.temp.inv*temp.arr)) # plant
    attributes$efficiencies[detritus] <- exp(eff.exp.det)*exp(eff.temp.inv*temp.arr) / (1 + exp(eff.exp.det)*exp(eff.temp.inv*temp.arr)) #detritus
    
    
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
    
    ## save stability value:
    meta.UK$stability[i] = stab
    
    ## food web fluxes
    meta.UK$flux[i] <- sum(flux)
    meta.UK$second.consumption[i] <- sum(flux[animals,])
    meta.UK$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
    # Per capita flux
    meta.UK$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta.UK$S[i] <- Number.of.species(bin.matrix)
    meta.UK$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
    meta.UK$MaxTL[i] <- max(attributes$TL)
    meta.UK$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta.UK$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta.UK$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
    
    meta.UK$L[i] <- Number.of.links(bin.matrix)
    meta.UK$LD[i] <- Link.density(bin.matrix)
    meta.UK$C[i] <- Connectance(bin.matrix)
    meta.UK$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta.UK$modularity[i] = modularity(igraph, membership(community))
  }
  
  
  #### River data from Brauns, Hall, and Rene-Mor et al. #### 
  meta.Brauns <- read.csv("Rivers_Mario/rivers_metadata.csv")
  meta.Brauns$study_ID <- rep(c('Elbe','Hall et al','Montsant'),times=c(3,4,4))
  web <- unique(meta.Brauns$FW_name)
  temperatures.K = meta.Brauns$temperature + 273.15
  for(i in 1:length(web)){
    
    temp.K = temperatures.K[i]
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
    attributes$omnivory = sapply(1:nrow(bin.matrix), Omnivory.species, bin.matrix, attributes$TL, simplify = TRUE)
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
    flux <- fluxing(mat=matrix, biomasses = attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies, 
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
    
    ## save stability value:
    meta.Brauns$stability[i] = stab
    
    
    ## food web fluxes
    meta.Brauns$flux[i] <- sum(flux)
    meta.Brauns$second.consumption[i] <- sum(flux[eat.animals])
    meta.Brauns$prim.consumption[i] <- sum(flux[basals,])
    # Per capita flux
    meta.Brauns$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta.Brauns$S[i] <- Number.of.species(bin.matrix)
    meta.Brauns$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
    meta.Brauns$MaxTL[i] <- max(attributes$TL)
    meta.Brauns$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta.Brauns$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta.Brauns$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
    
    meta.Brauns$L[i] <- Number.of.links(bin.matrix)
    meta.Brauns$LD[i] <- Link.density(bin.matrix)
    meta.Brauns$C[i] <- Connectance(bin.matrix)
    meta.Brauns$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta.Brauns$modularity[i] = modularity(igraph, membership(community))
  }
  meta.Brauns$C <- as.numeric(meta.Brauns$C)
  
  
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
    attributes$omnivory = sapply(1:nrow(bin.matrix), Omnivory.species, bin.matrix, attributes$TL, simplify = TRUE)
    attributes$omnivory[attributes$TL == 1] <- NA
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
    
    ## save stability value:
    meta.Brazil$stability[i] = stab
    
    ## food web fluxes
    meta.Brazil$flux[i] <- sum(flux)
    meta.Brazil$second.consumption[i] <- sum(flux[animals,])
    meta.Brazil$prim.consumption[i] <- sum(flux[plants,], flux[detritus,])
    # Per capita flux
    meta.Brazil$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta.Brazil$S[i] <- Number.of.species(bin.matrix)
    meta.Brazil$density[i] <- sum(attributes$abundance[animals], na.rm = TRUE)
    meta.Brazil$MaxTL[i] <- max(attributes$TL)
    meta.Brazil$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta.Brazil$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta.Brazil$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
    
    meta.Brazil$L[i] <- Number.of.links(bin.matrix)
    meta.Brazil$LD[i] <- Link.density(bin.matrix)
    meta.Brazil$C[i] <- Connectance(bin.matrix)
    meta.Brazil$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta.Brazil$modularity[i] = modularity(igraph, membership(community))
  }
  
  
  commcols1 <- intersect(names(meta.IS), names(meta.UK))
  commcols2 <- intersect(commcols1, names(meta.Brauns))
  commcols <- intersect(commcols2, names(meta.Brazil))
  meta.Streams <- bind_rows(select(meta.IS, all_of(commcols)),
                            select(meta.UK, all_of(commcols)),
                            select(meta.Brauns, all_of(commcols)),
                            select(meta.Brazil, all_of(commcols)))
  
  meta.Streams$replicate = repl
  
  return(meta.Streams)
}
