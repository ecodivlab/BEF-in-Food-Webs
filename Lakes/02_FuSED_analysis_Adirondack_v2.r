

rm(list=ls())
setwd("C:/Users/darren/OneDrive - University of Canberra/Workflows/FuSED/Adirondack flux analysis")
source("C:/Users/darren/OneDrive - University of Canberra/Workflows/FuSED/Adirondack flux analysis/Food_web_functions.r")
library(fluxweb); library(cheddar); library(igraph); library(RColorBrewer)
options(stringsAsFactors=FALSE)

## Put the food webs into FUSED format
boltz <- 0.00008617343
T0 <- 273.15 + 20   ## Use 20 degrees as standardised temperature across all dataset
perday <- 60*60*24
meta <- read.csv("metadata/Adirondack_metadata.csv")
meta$FW_name <- meta$FW_name_database
lakes <- unique(meta$FW_name)
lakes <- lakes[-50] # remove Wolf Lake - missing pages in report, no density data


for(i in 1:length(lakes)) {
  
  # imports
    lake <- NULL; nodes <- NULL; matrix<-NULL
    lake <- lakes[i]
    nodes <- read.csv(paste0("sp_att/Adirondack_spAttributes_",lake,".csv"))
      # corrections to match database
      nodes$taxonomy[nodes$taxonomy=="Salmo trutta"] <- "Salmo rutta" # should be trutta but match to incorrect edgelist
      nodes$taxonomy[nodes$taxonomy=="Anabaena flos-aquae"] <- "Anabaena flos.aquae"
      
    matrix <- as.matrix(read.csv(paste0("adj_mat/Adirondack_matrix_",lake,".csv"), row.names = 1))

  # check order
    nodes <- nodes[order(nodes$taxonomy),]
    matrix <- matrix[order(colnames(matrix)),order(colnames(matrix))]
    print(i)
    print(all(colnames(matrix) == gsub(" ", ".", nodes$taxonomy)))
    
  # remove taxa that do not have matches, and with 0 density
    unmatched.rows <- NULL; fishy <- NULL; detr <- NULL
    unmatched.rows <- which(is.na(nodes$Count.of.Lake_species)) # which are missing? Could not be matched with density from report
    eggs <- which(grepl("fish eggs", nodes$Lake_species)) # but keep fish eggs
    detr <- which(grepl("detritus", nodes$Lake_species))  # but keep detritus
    nauplii <- which(grepl("nauplii", nodes$Lake_species)) 
    unmatched.rows <- unmatched.rows[!(unmatched.rows %in% c(detr))] ## <--- add 'eggs' here to keep eggs
    
    phyto0 <- which(is.na(nodes$density) & nodes$species_type == "phytoplankton") # phytoplankton with 0 density reported
    unmatched.rows <- c(unmatched.rows, phyto0)
    
    if(length(unmatched.rows)>0) {
         nodes <- nodes[-unmatched.rows,]
         matrix <- matrix[-unmatched.rows,-unmatched.rows]
    }
    print(all(colnames(matrix) == gsub(" ", ".", nodes$taxonomy)))
    
  # temperature
    temp <- meta$mean_temp_C[meta$FW_name_database==lake]
    temp.kT <- ((273.15+temp)-T0)/(boltz*(273.15+temp)*T0)
    
  # losses
    nodes$losses <- 0
    nodes$losses[nodes$metabolic.type=="invertebrate"] <- exp((0.71 * log(nodes$mass.mean.g.[nodes$metabolic.type=="invertebrate"]) + 17.17) - 0.69/(boltz*(273.15+temp))) * nodes$density.comb[nodes$metabolic.type=="invertebrate"]
    nodes$losses[nodes$metabolic.type=="ectotherm vertebrate"] <- exp((0.71 * log(nodes$mass.mean.g.[nodes$metabolic.type=="ectotherm vertebrate"]) + 18.47) - 0.69/(boltz*(273.15+temp))) * nodes$density.comb[nodes$metabolic.type=="ectotherm vertebrate"]
    
        ## Loses for fish eggs?? - unknown density - currently these are omitted anyway
        nodes$losses[is.na(nodes$losses)] <- 0
  
  # efficiencies
    nodes$efficiencies[nodes$species_type=="zooplankton" | nodes$species_type=="fish"] <- exp(2.266)*exp(0.164*temp.kT) / (1 + exp(2.266)*exp(0.164*temp.kT)) #animal
    nodes$efficiencies[nodes$species_type=="phytoplankton"] <- exp(0.179)*exp(0.164*temp.kT) / (1 + exp(0.179)*exp(0.164*temp.kT)) # plant
    nodes$efficiencies[nodes$species_type=="detritus"] <- exp(-1.670)*exp(0.164*temp.kT) / (1 + exp(-1.670)*exp(0.164*temp.kT)) #detritus

  # preferences (biomass matrix)
    B.matrix <- matrix * nodes$biomass
   
    # missing detritus node given sum of phytoplankton biomass for each predator (i.e. equal weight to green/brown origin)
    missing.detritus <- nodes$species_type=="detritus"
    B.matrix[missing.detritus,] <- rep(colSums(B.matrix[nodes$species_type=="phytoplankton",], na.rm=T) / sum(missing.detritus), each=sum(missing.detritus))
    B.matrix[matrix!=0 & B.matrix==0] <- 1 # for when detritus is the only food source
      
    # missing fish eggs
        #missing.fish <- nodes$taxonomy=="fish eggs"
        #B.matrix[missing.fish,] <- 
    

  ## Fluxweb analysis
  diag(B.matrix) <- 0   ## Make sure there are no cannibalistic links
  bin.matrix <- B.matrix
  bin.matrix[bin.matrix>0] <- 1
  igraph <- graph_from_adjacency_matrix(bin.matrix, mode="directed")
  attributes <- nodes
  attributes$losses[is.na(attributes$losses)] <- 0

  ## Calculating fluxes
  flux <- fluxing(mat=B.matrix, losses=attributes$losses, efficiencies=attributes$efficiencies, bioms.losses=F, bioms.prefs = F) * perday
  meta$flux[i] <- sum(flux)
  meta$animal.flux[i] <- sum(flux[attributes$species_type=="zooplankton" | attributes$species_type=="fish",])
  meta$plant.flux[i] <- sum(flux[attributes$species_type=="phytoplankton",])
  meta$detritus.flux[i] <- sum(flux[attributes$species_type=="detritus",])
  meta$ratio.flux[i] <- meta$animal.flux[i] / sum(meta$plant.flux[i] + meta$detritus.flux[i])
  
  ## Biomass distributions
  meta$biomass[i] <- sum(attributes$biomass, na.rm=T)
  meta$animal.biomass[i] <- sum(attributes$biomass[attributes$species_type=="zooplankton" | attributes$species_type=="fish"], na.rm=T)
  meta$plant.biomass[i] <- sum(attributes$biomass[attributes$species_type=="phytoplankton"], na.rm=T)
  #meta$detritus.biomass[i] <- sum(attributes$biomass[attributes$species_type=="detritus"], na.rm=T)
  meta$ratio.biomass[i] <- meta$animal.biomass[i] / meta$plant.biomass[i]
  meta$MNslope[i] <- coef(lm(log10(attributes$density.comb[attributes$taxonomy!="detritus"]) ~ log10(attributes$mass.mean.g.[attributes$taxonomy!="detritus"])))[2]
  
  ## Food web metrics
  meta$S[i] <- Number.of.species(bin.matrix)
  meta$L[i] <- Number.of.links(bin.matrix)
  meta$LD[i] <- Link.density(bin.matrix)
  meta$C[i] <- Connectance(bin.matrix)
  meta$B[i] <- Bottom.Intermediate.Top(bin.matrix)$Proportions.of.each[1]
  meta$I[i] <- Bottom.Intermediate.Top(bin.matrix)$Proportions.of.each[2]
  meta$T[i] <- Bottom.Intermediate.Top(bin.matrix)$Proportions.of.each[3]
  meta$generality[i] <- Gen.sd(bin.matrix)
  meta$vulnerability[i] <- Vul.sd(bin.matrix)
  meta$MeanTL[i] <- mean(TL(bin.matrix))
  meta$MaxTL[i] <- max(TL(bin.matrix))
  #meta$modularity[i] <- cluster_spinglass(igraph)$modularity
  
  # plot fluxes
    # community object for prey-averaged trophic level
    flux.igraph <- graph_from_adjacency_matrix(flux, mode="directed", weighted=TRUE)
  
    nodes.ched <- data.frame("node"=colnames(matrix))
    edge.list <- data.frame(as_edgelist(igraph))
    colnames(edge.list) = c("resource", "consumer")
    ched <- Community(nodes=nodes.ched, properties=list(title=i), trophic.links=edge.list)
    patl <-   PreyAveragedTrophicLevel(ched)
    
    V(flux.igraph)$size <- (nodes$biomass - min(nodes$biomass,na.rm=T)) / (max(nodes$biomass,na.rm=T) - min(nodes$biomass,na.rm=T))*20
    V(flux.igraph)$size[is.na(V(flux.igraph)$size)] <- 10
    #V(flux.igraph)$size <- 10
    
    cols <- brewer.pal(12, "Paired")
    V(flux.igraph)$color[nodes$species_type=="detritus"] <- cols[11]
    V(flux.igraph)$color[nodes$species_type=="phytoplankton"] <- cols[3]
    V(flux.igraph)$color[nodes$species_type=="zooplankton"] <- cols[7]
    V(flux.igraph)$color[nodes$species_type=="fish"] <- cols[5]
    
    V(flux.igraph)$label.color[nodes$species_type=="detritus"] <- cols[12]
    V(flux.igraph)$label.color[nodes$species_type=="phytoplankton"] <- cols[4]
    V(flux.igraph)$label.color[nodes$species_type=="zooplankton"] <- cols[8]
    V(flux.igraph)$label.color[nodes$species_type=="fish"] <- cols[6]
    
    E(flux.igraph)$weight <- E(flux.igraph)$weight*10000
    
    ys <- patl
    ys[patl==1] <- jitter(ys[patl==1],8)
    ys[patl==2] <- jitter(ys[patl==2],8)
    xs <- rep(0,nrow(nodes))
    xs[patl==1] <- seq(0.2,1,length.out=sum(patl==1))
    xs[patl>1] <- seq(0,1,length.out=sum(patl>1))
    xs[patl>2] <- seq(0,1,length.out=sum(patl>2))
    xs[nodes$taxonomy=="benthic detritus"] <- 0
    ys[degree(flux.igraph)==0] <- 0.5
    
    png(paste0("C:/Users/darren/OneDrive - University of Canberra/Workflows/FuSED/Adirondack flux analysis/fluxplots/",lake,".png"), width=15, height=15, units="cm", res=400)
    
    plot(flux.igraph, layout=cbind(xs,ys), edge.arrow.size=0.25, edge.width=E(flux.igraph)$weight,
         vertex.color = V(flux.igraph)$color, vertex.label.cex=0.5)
  
    dev.off()
  
}


# ## Relationships between flux and temperature
plot(log10(meta$flux) ~ meta$mean_temp_C)
plot(meta$mean_temp_C ~ meta$max.depth..m.) # need to think about how temperature/depth is affecting the rates

# plot(log10(meta$ratio.flux) ~ meta$mean_temp_C)
# plot(log10(meta$animal.flux) ~ meta$mean_temp_C)
# plot(log10(meta$plant.flux) ~ meta$mean_temp_C)
# 
# ## Relationships between flux and biomass
plot(log10(meta$flux*1000000) ~ meta$biomass)
# plot(log10(meta$ratio.flux) ~ meta$biomass)
# plot(log10(meta$animal.flux) ~ meta$biomass)
# plot(log10(meta$plant.flux) ~ meta$biomass)
# 
# ## Relationships between flux and node richness
plot(log10(meta$flux) ~ meta$S)
# plot(log10(meta$ratio.flux) ~ meta$S)
# plot(log10(meta$animal.flux) ~ meta$S)
# plot(log10(meta$plant.flux) ~ meta$S)
# 
# ## Relationships between flux and connectance
# plot(log10(meta$flux) ~ meta$C)
# plot(log10(meta$ratio.flux) ~ meta$C)
# plot(log10(meta$animal.flux) ~ meta$C)
# plot(log10(meta$plant.flux) ~ meta$C)
# 
# ## Relationships between flux and mean prey-averaged trophic level
plot(log10(meta$flux) ~ meta$MeanTL)
# plot(log10(meta$ratio.flux) ~ meta$MeanTL)
# plot(log10(meta$animal.flux) ~ meta$MeanTL)
# plot(log10(meta$plant.flux) ~ meta$MeanTL)
# 
# ## Relationships between flux and maximum prey-averaged trophic level
plot(log10(meta$flux) ~ meta$MaxTL)
# plot(log10(meta$ratio.flux) ~ meta$MaxTL)
# plot(log10(meta$animal.flux) ~ meta$MaxTL)
# plot(log10(meta$plant.flux) ~ meta$MaxTL)
# 
# ## Relationships between flux and modularity
# plot(log10(meta$flux) ~ meta$modularity)
# plot(log10(meta$ratio.flux) ~ meta$modularity)
# plot(log10(meta$animal.flux) ~ meta$modularity)
# plot(log10(meta$plant.flux) ~ meta$modularity)
# 
# ## Relationships between flux and LD
plot(log10(meta$flux) ~ meta$LD)
# plot(log10(meta$ratio.flux) ~ meta$modularity)
# plot(log10(meta$animal.flux) ~ meta$modularity)
# plot(log10(meta$plant.flux) ~ meta$modularity)
