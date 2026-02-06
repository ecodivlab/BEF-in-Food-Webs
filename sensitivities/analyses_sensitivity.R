## Load packages ##
library(tidyverse); library(ggeffects); library(gridExtra); library(piecewiseSEM);
library(patchwork); library(nlme); library(grid); library(car); library(rempsyc); library(ggpattern);
library(ggh4x); library(scales); library(ggtext); library(ggplot2)


rm(list = ls())

########## read and shape metadata for the different ecosystems ####################

NPP.proxy <- read.csv("NDVI and Chlorophyll-a/data/proxy-npp.csv") ## import NDVI & Chl-a data
meta.Marine <- read.csv('meta.Marine.csv')  
meta.Marine <- meta.Marine %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(NPP.proxy = avg)
meta.Marine$NPP.proxy <- (meta.Marine$NPP.proxy - min(meta.Marine$NPP.proxy, na.rm = TRUE)) /
  (max(meta.Marine$NPP.proxy, na.rm = TRUE) - min(meta.Marine$NPP.proxy, na.rm = TRUE))
meta.Marine <- meta.Marine %>% mutate(NPP.scale = logit(NPP.proxy))
meta.Marine <- meta.Marine %>% mutate(NPP.scale2 = NPP.scale^2)

meta.Soils <- read.csv('meta.Soils.csv')
meta.Soils <- meta.Soils %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(temperature_C = temperature, NPP.proxy = avg)
meta.Soils <- meta.Soils %>% mutate(NPP.scale = logit(NPP.proxy))
meta.Soils <- meta.Soils %>% mutate(NPP.scale2 = NPP.scale^2)

meta.Streams <- read.csv('meta.Streams.csv')
meta.Streams <- meta.Streams %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(temperature_C = temperature, NPP.raw = avg)
meta.Streams <- meta.Streams %>% mutate(NPP.proxy = ifelse(NPP.raw < 0, 0, NPP.raw), # replace negative value at Cananeia SP6 with 0
                                        NPP.scale = logit(NPP.proxy))
meta.Streams <- meta.Streams %>% mutate(NPP.scale2 = NPP.scale^2)

meta.Lakes <- read.csv('meta.Lakes.csv')
meta.Lakes <- meta.Lakes %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(NPP.proxy = avg)
meta.Lakes <- meta.Lakes %>% mutate(NPP.scale = logit(NPP.proxy))
meta.Lakes <- meta.Lakes %>% mutate(NPP.scale2 = NPP.scale^2)



############################################################


# retrieve flux results 
load("sensitivities/sensitivity.Rdata")

for (i in 1: length(results.lakes)){
  
  ## Compile data sets for cross-ecosystem analysis ##
  commcols <- intersect(names(results.marine[[i]]), names(results.soils[[i]]))
  commcols <- intersect(commcols, names(results.streams[[i]]))
  commcols <- intersect(commcols, names(results.lakes[[i]]))
  
  
  all_data <- bind_rows(select(results.marine[[i]], all_of(commcols)),
                        select(results.soils[[i]], all_of(commcols)),
                        select(results.streams[[i]], all_of(commcols)),
                        select(results.lakes[[i]], all_of(commcols)))
  
  
  ## Predation
  S.predation_Global <- lme(log10(second.consumption) ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data)
  plot(S.predation_Global, which=1)
  qqnorm(S.predation_Global)
  summary(S.predation_Global)
  
  S.predation_Globala=update(S.predation_Global, random = ~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID))
  anova(S.predation_Global,S.predation_Globala)
  plot(S.predation_Globala)
  qqnorm(S.predation_Globala)
  summary(S.predation_Globala)
  
  S.predation_Globalb=update(S.predation_Global, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S))) #Best model
  anova(S.predation_Global,S.predation_Globala,S.predation_Globalb)
  plot(S.predation_Globalb)
  qqnorm(S.predation_Globalb)
  summary(S.predation_Globalb)
  
  S.predation_Global <- S.predation_Globalb
  
  ## Primary consumption
  S.prim.cons_Global <- lme(log10(prim.consumption) ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data)
  plot(S.prim.cons_Global, which=1)
  qqnorm(S.prim.cons_Global)
  summary(S.prim.cons_Global)
  
  S.prim.cons_Globala=update(S.prim.cons_Global, weights=varIdent(form=~1|study_ID)) #Best model
  anova(S.prim.cons_Global,S.prim.cons_Globala)
  plot(S.prim.cons_Globala)
  qqnorm(S.prim.cons_Globala)
  summary(S.prim.cons_Globala)
  
  S.prim.cons_Globalb=update(S.prim.cons_Global, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)))
  # anova(S.prim.cons_Global,S.prim.cons_Globala,S.prim.cons_Globalb)
  # plot(S.prim.cons_Globalb)
  # qqnorm(S.prim.cons_Globalb)
  # summary(S.prim.cons_Globalb)
  
  S.prim.cons_Global = S.prim.cons_Globala
  
  all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                    levels = c("Marine", "Soils", "Streams", "Lakes"))
  
  
  
  
  
  
}