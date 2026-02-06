#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Calculation of energy fluxes in 48 lake food webs applying sensitivty
#
#=========================================================================================================================


#### Adirondack and Maggiore lakes ####


flux.lakes = function(params){
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
  
  meta <- read.csv("Lakes/Lakes_metadata.csv")
  lakes <- unique(meta$FW_name)
  
  #### Calculate food web properties and energy fluxes ####
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
    consumers <- (bin.matrix)[,1]
    attributes$omnivory = sapply(1:nrow(bin.matrix), Omnivory.species, bin.matrix, attributes$TL, simplify = TRUE)
    attributes$omnivory[attributes$TL == 1] <- NA
    
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
    # Total flux
    if(meta$study_ID[i]=='Lake Maggiore'){
      flux <- fluxing(mat=bin.matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies,
                      bioms.losses=F, bioms.prefs = T) * perday}
    else if (meta$study_ID[i]=='Adirondack lakes'){
      flux <- fluxing(mat=bin.matrix, biomasses=attributes$biomass, losses=attributes$losses, efficiencies=attributes$efficiencies,
                      bioms.losses=F, bioms.prefs = T) * perday * 1000000 #*1000000 to convert from ml to cubic meters for Adirondack Lakes
    }
    # Per capita flux
    PC.flux <- sweep(flux, 1, nodes$density.comb, FUN = "/")
    
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
    # Per capita flux
    meta$PC.predation[i] <- mean(PC.flux[animals], na.rm = TRUE)
    
    ## Food web metrics
    meta$S[i] <- Number.of.species(bin.matrix)
    meta$density[i] <- sum(attributes$density.comb[animals], na.rm = TRUE)
    meta$MaxTL[i] <- max(TL(bin.matrix))
    meta$sim.sec.cons[i] = mean(sim.mat[sec.cons, sec.cons]) # trophic similarity for secondary consumers
    meta$sim.prim.cons[i] = mean(sim.mat[primary.cons.and.omnivores, primary.cons.and.omnivores])
    meta$sim.total[i] = mean(sim.mat[!basals, !basals]) # trophic similarity for whole food web
    # cat(Number.of.links(bin.matrix), '\n')
    meta$L[i] <- Number.of.links(bin.matrix)
    meta$LD[i] <- Link.density(bin.matrix)
    meta$C[i] <- Connectance(bin.matrix)
    meta$omnivory[i] <- mean(attributes$omnivory, na.rm=T)
    
    # modularity
    community = cluster_walktrap(igraph)
    meta$modularity[i] = modularity(igraph, membership(community))
  }
  #view(meta[,-c(2:32, 34)])
  meta$L = as.numeric(meta$L) #not clear why converted to string just here
  meta.Lakes <- meta[,-c(2:32, 34)]
  
  meta.Lakes <- filter(meta.Lakes, second.consumption > 0) #remove food webs with no secondary consumers (Brooke Trout and South Lake). 
  
  meta.Lakes$replicate = repl
  
  return(meta.Lakes)
  
}


