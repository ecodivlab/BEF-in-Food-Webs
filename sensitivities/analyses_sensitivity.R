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
  
  all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                    levels = c("Marine", "Soils", "Streams", "Lakes"))
  
  
  ## Cross-ecosystem analysis
  SEM.all2 <- psem(
    lme(log(S) ~ NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
    lme(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
    lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), data=all_data),
    lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S), random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
    lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data),
    lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale + logit(sim.sec.cons), random=~1|ecosystem.type/study_ID, data=all_data),
    
    logit(sim.prim.cons) %~~% logit(sim.sec.cons),
    log(prim.consumption) %~~% logit(sim.sec.cons),
    log(second.consumption) %~~% logit(sim.prim.cons),
    log(second.consumption) %~~% log(prim.consumption)
  )
  
  #Extracts r-squared for primary & secondary consumption, global P-value, and p-values for paths to functions included in minimum adequate SEM (excluding NPP)
  xx = summary(SEM.all2)
  output.all = c(rsquared(SEM.all2)[5,6], rsquared(SEM.all2)[6,6], xx$Cstat[,3],
                 xx$coefficients[10,7], xx$coefficients[12,7], xx$coefficients[13,7], xx$coefficients[1,7]) 
  names(output.all) = c('r2.prim.consumtion', 'r2.predation', 'SEM.P_value',
                        'P_taxa_prim.consumption', 'P_MaxTL_predation', 'P_taxa_predation', 'P_sim.sec.cons_predation')
  
  
  
  ## Marine
  SEM.Marine3 <- psem(
    
    gls(log(S) ~ NPP.scale + NPP.scale2, data=meta.Marine),
    gls(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, data=meta.Marine),
    gls(logit(sim.prim.cons) ~ log(S), data=meta.Marine),
    gls(logit(sim.sec.cons) ~ log(S) + log(MaxTL) + NPP.scale, data=meta.Marine),
    gls(log(prim.consumption) ~ log(S) + NPP.scale + NPP.scale2, weights=varPower(form=~NPP.scale), data=meta.Marine),
    gls(log(second.consumption) ~ log(MaxTL) + NPP.scale + NPP.scale2, weights=varPower(form=~S), data=meta.Marine),
    
    logit(sim.prim.cons) %~~% logit(sim.sec.cons),
    log(prim.consumption) %~~% logit(sim.sec.cons),
    log(second.consumption) %~~% logit(sim.prim.cons),
    log(second.consumption) %~~% log(prim.consumption)
  )
  
  #Extracts r-squared for primary & secondary consumption, global P-value, and p-values for paths to functions included in minimum adequate SEM (excluding NPP)
  xx = summary(SEM.Marine3)
  output.Marine = c(rsquared(SEM.Marine3)[5,5], rsquared(SEM.Marine3)[6,5], xx$Cstat[,3],
                    xx$coefficients[10,7], xx$coefficients[13,7]) 
  names(output.Marine) = c('r2.prim.consumtion', 'r2.predation', 'SEM.P_value',
                           'P_taxa_prim.consumption', 'P_MaxTL_predation')
  
  
  
  ## Soils
  SEM.Soils3 <- psem(
    
    gls(log(S) ~ NPP.scale + NPP.scale2, data=meta.Soils),
    gls(log(MaxTL) ~ log(S), data=meta.Soils),
    gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Soils),
    gls(logit(sim.sec.cons) ~ NPP.scale + log(S), data=meta.Soils),
    gls(log(prim.consumption) ~ log(S) + NPP.scale, weights=varIdent(form=~1|study_ID), data=meta.Soils),
    gls(log(second.consumption) ~ log(S) + log(MaxTL), weights=varIdent(form=~1|study_ID), data=meta.Soils),
    
    logit(sim.prim.cons) %~~% logit(sim.sec.cons),
    log(prim.consumption) %~~% logit(sim.sec.cons),
    log(second.consumption) %~~% logit(sim.prim.cons),
    log(second.consumption) %~~% log(prim.consumption)
  )
  
  #Extracts r-squared for primary & secondary consumption, global P-value, and p-values for paths to functions included in minimum adequate SEM (excluding NPP)
  xx = summary(SEM.Soils3)
  output.Soils = c(rsquared(SEM.Soils3)[5,5], rsquared(SEM.Soils3)[6,5], xx$Cstat[,3],
                   xx$coefficients[8,7], xx$coefficients[10,7], xx$coefficients[11,7]) 
  names(output.Soils) = c('r2.prim.consumtion', 'r2.predation', 'SEM.P_value',
                          'P_taxa_prim.consumption', 'P_taxa_predation', 'P_MaxTL_predation')
  
  
  
  ## Streams
  SEM.Streams3 <- psem(
    
    lme(log(S) ~ NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
    lme(log(MaxTL) ~ log(S), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
    lme(logit(sim.prim.cons) ~ log(MaxTL), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
    lme(logit(sim.sec.cons) ~ log(S), random=~1|study_ID, weights=varExp(), data=meta.Streams),
    lme(log(prim.consumption) ~ log(S) , random=~1|study_ID, data=meta.Streams),
    lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), random=~1|study_ID, data=meta.Streams),
    
    logit(sim.prim.cons) %~~% logit(sim.sec.cons),
    log(prim.consumption) %~~% logit(sim.sec.cons),
    log(second.consumption) %~~% logit(sim.prim.cons),
    log(second.consumption) %~~% log(prim.consumption)
  )
  
  #Extracts r-squared for primary & secondary consumption, global P-value, and p-values for paths to functions included in minimum adequate SEM (excluding NPP)
  xx = summary(SEM.Streams3)
  output.Streams = c(rsquared(SEM.Streams3)[5,6], rsquared(SEM.Streams3)[6,6], xx$Cstat[,3],
                     xx$coefficients[5,7], xx$coefficients[6,7], xx$coefficients[7,7]) 
  names(output.Streams) = c('r2.prim.consumtion', 'r2.predation', 'SEM.P_value',
                            'P_taxa_prim.consumption', 'P_MaxTL_predation', 'P_sim.sec.cons_predation')
  
  
  
  ## Lakes
  SEM.Lakes3 <- psem(
    
    gls(log(S) ~ NPP.scale, data=meta.Lakes),
    gls(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, weights=varPower(form=~S), data=meta.Lakes),
    gls(logit(sim.prim.cons) ~ log(S), weights=varExp(form=~S), data=meta.Lakes),
    gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varExp(form=~log(MaxTL)), data=meta.Lakes),
    gls(log(prim.consumption) ~ log(MaxTL) + log(S), weights=varPower(), data=meta.Lakes),
    gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S), data=meta.Lakes),
    
    logit(sim.prim.cons) %~~% logit(sim.sec.cons),
    log(prim.consumption) %~~% logit(sim.sec.cons),
    log(second.consumption) %~~% logit(sim.prim.cons),
    log(second.consumption) %~~% log(prim.consumption)
  )
  
  #Extracts r-squared for primary & secondary consumption, global P-value, and p-values for paths to functions included in minimum adequate SEM (excluding NPP)
  xx = summary(SEM.Lakes3)
  output.Lakes = c(rsquared(SEM.Lakes3)[5,5], rsquared(SEM.Lakes3)[6,5], xx$Cstat[,3],
                   xx$coefficients[7,7], xx$coefficients[8,7], xx$coefficients[9,7], xx$coefficients[10,7]) 
  names(output.Lakes) = c('r2.prim.consumtion', 'r2.predation', 'SEM.P_value',
                          'P_MaxTL_prim.consumption', 'P_taxa_prim.consumption', 'P_sim.sec.cons_predation', 'P_taxa_predation')
  
  
}
