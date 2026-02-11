## Load packages ##
library(tidyverse); library(ggeffects); library(gridExtra); library(piecewiseSEM);
library(patchwork); library(nlme); library(grid); library(car); library(rempsyc); library(ggpattern);
library(ggh4x); library(scales); library(ggtext); library(ggplot2)


rm(list = ls())

########## read and shape metadata for the different ecosystems ####################
NPP.proxy <- read.csv("NDVI and Chlorophyll-a/data/proxy-npp.csv") ## import NDVI & Chl-a data


############################################################


# retrieve flux results 
load("sensitivities/sensitivity.Rdata")

for (i in 1: length(results.lakes)){
  
  meta.Marine <- results.marine[[i]]
  meta.Marine <- meta.Marine %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
    rename(NPP.proxy = avg)
  meta.Marine$NPP.proxy <- (meta.Marine$NPP.proxy - min(meta.Marine$NPP.proxy, na.rm = TRUE)) /
    (max(meta.Marine$NPP.proxy, na.rm = TRUE) - min(meta.Marine$NPP.proxy, na.rm = TRUE))
  meta.Marine <- meta.Marine %>% mutate(NPP.scale = logit(NPP.proxy))
  meta.Marine <- meta.Marine %>% mutate(NPP.scale2 = NPP.scale^2)
  
  meta.Soils <- results.soils[[i]]
  meta.Soils <- meta.Soils %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
    rename(temperature_C = temperature, NPP.proxy = avg)
  meta.Soils <- meta.Soils %>% mutate(NPP.scale = logit(NPP.proxy))
  meta.Soils <- meta.Soils %>% mutate(NPP.scale2 = NPP.scale^2)
  
  meta.Streams <- results.streams[[i]]
  meta.Streams <- meta.Streams %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
    rename(temperature_C = temperature, NPP.raw = avg)
  meta.Streams <- meta.Streams %>% mutate(NPP.proxy = ifelse(NPP.raw < 0, 0, NPP.raw), # replace negative value at Cananeia SP6 with 0
                                          NPP.scale = logit(NPP.proxy))
  meta.Streams <- meta.Streams %>% mutate(NPP.scale2 = NPP.scale^2)
  
  meta.Lakes <- results.lakes[[i]]
  meta.Lakes <- meta.Lakes %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
    rename(NPP.proxy = avg)
  meta.Lakes <- meta.Lakes %>% mutate(NPP.scale = logit(NPP.proxy))
  meta.Lakes <- meta.Lakes %>% mutate(NPP.scale2 = NPP.scale^2)
  
  ## Compile data sets for cross-ecosystem analysis ##
  commcols <- intersect(names(meta.Marine), names(meta.Soils))
  commcols <- intersect(commcols, names(meta.Streams))
  commcols <- intersect(commcols, names(meta.Lakes))
  
  all_data <- bind_rows(select(meta.Marine, all_of(commcols)),
                        select(meta.Soils, all_of(commcols)),
                        select(meta.Streams, all_of(commcols)),
                        select(meta.Lakes, all_of(commcols)))
  
  ## Predation
  S.predation_Global <- lme(log10(second.consumption) ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data)
  S.predation_Globala=update(S.predation_Global, random = ~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID))
  S.predation_Globalb=update(S.predation_Global, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S))) #Best model
  S.predation_Global <- S.predation_Globalb
  
  ## Primary consumption
  S.prim.cons_Global <- lme(log10(prim.consumption) ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data)
  S.prim.cons_Globala=update(S.prim.cons_Global, weights=varIdent(form=~1|study_ID)) #Best model
  S.prim.cons_Globalb=update(S.prim.cons_Global, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)))
  S.prim.cons_Global = S.prim.cons_Globala
  
  all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                    levels = c("Marine", "Soils", "Streams", "Lakes"))
  
  ### Maximal model
  SEM.all <- psem(
    lme(log(S) ~ NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
    lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data),
    lme(logit(sim.prim.cons) ~ log(MaxTL) + log(S), weights=varIdent(form=~1|study_ID), random=~1|ecosystem.type/study_ID, data=all_data),
    lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S), random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
    lme(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), random=~1|ecosystem.type/study_ID, 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
    lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), random=~1|ecosystem.type/study_ID, data=all_data),
    
    logit(sim.prim.cons) %~~% logit(sim.sec.cons),
    log(prim.consumption) %~~% logit(sim.sec.cons),
    log(second.consumption) %~~% logit(sim.prim.cons),
    log(second.consumption) %~~% log(prim.consumption)
  )
  
  xx = summary(SEM.all)
  
  
  
  
}
