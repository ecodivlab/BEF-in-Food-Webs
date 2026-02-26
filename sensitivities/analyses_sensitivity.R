## Load packages ##
library(tidyverse); library(ggeffects); library(gridExtra); library(piecewiseSEM);
library(patchwork); library(nlme); library(grid); library(car); library(rempsyc); library(ggpattern);
library(ggh4x); library(scales); library(ggtext); library(ggplot2); 
library(future.apply)


rm(list = ls())
options(digits = 5)


########## read and shape metadata for the different ecosystems ####################
NPP.proxy <- read.csv("NDVI and Chlorophyll-a/data/proxy-npp.csv") ## import NDVI & Chl-a data


############################################################


# retrieve flux results 
load("sensitivities/sensitivity.Rdata")


# run SEM for all flux replicates
SEM.sensitivity = function(i, results.lakes, results.marine, results.soils, results.streams){
  results = list()
  meta.Marine <- results.marine[[i]]
  meta.Marine$xxx = NA
  # print(meta.Marine)
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
    lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
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
  
  results$sensitivity.all = output.all
  meta.Marine[1,1]
  # weird interaction between gls and apply... have to make it global if in a function???
  meta.Marine2 <<- meta.Marine
 
  ## Marine
  SEM.Marine3 <- psem(
    
    gls(log(S) ~ NPP.scale + NPP.scale2, data=meta.Marine2),
    gls(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, data=meta.Marine2),
    gls(logit(sim.prim.cons) ~ log(S), data=meta.Marine2),
    gls(logit(sim.sec.cons) ~ log(S) + log(MaxTL) + NPP.scale, data=meta.Marine2),
    gls(log(prim.consumption) ~ log(S) + NPP.scale + NPP.scale2, weights=varPower(form=~NPP.scale), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Marine2),
    gls(log(second.consumption) ~ log(MaxTL) + NPP.scale + NPP.scale2, weights=varPower(form=~S), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Marine2),
    
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
  results$sensitivity.Marine = output.Marine

  
  ## Soils
  
  meta.Soils2 <<-meta.Soils
  SEM.Soils3 <- psem(
    
    gls(log(S) ~ NPP.scale + NPP.scale2, data=meta.Soils2),
    gls(log(MaxTL) ~ log(S), data=meta.Soils2),
    gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Soils2),
    gls(logit(sim.sec.cons) ~ NPP.scale + log(S), data=meta.Soils2),
    gls(log(prim.consumption) ~ log(S) + NPP.scale, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Soils2),
    gls(log(second.consumption) ~ log(S) + log(MaxTL), weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Soils2),
    
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

  results$sensitivity.Soils = output.Soils
  
  ## Streams
  SEM.Streams3 <- psem(
    
    lme(log(S) ~ NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Streams),
    lme(log(MaxTL) ~ log(S), random=~1|study_ID, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Streams),
    lme(logit(sim.prim.cons) ~ log(MaxTL), random=~1|study_ID, weights=varIdent(form=~1|study_ID), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Streams),
    lme(logit(sim.sec.cons) ~ log(S), random=~1|study_ID, weights=varExp(), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Streams),
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
  results$sensitivity.Streams = output.Streams
  
  
  ## Lakes
  meta.Lakes2 <<- meta.Lakes
  SEM.Lakes3 <- psem(
    
    gls(log(S) ~ NPP.scale, data=meta.Lakes2),
    gls(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, weights=varPower(form=~S), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Lakes2),
    gls(logit(sim.prim.cons) ~ log(S), weights=varExp(form=~S), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Lakes2),
    gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varExp(form=~log(MaxTL)), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Lakes2),
    gls(log(prim.consumption) ~ log(MaxTL) + log(S), weights=varPower(), 
        control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=meta.Lakes2),
    gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S), data=meta.Lakes2),
    
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
  results$sensitivity.Lakes = output.Lakes
  return(results)
}

# because of the replicative approach, some SEM might fail. using the wrapper function 
# with a trycatch will ensure that NA are returned in that case (and not errors breaking everything)
SEM.wrapper = function(i, results.lakes, results.marine, results.soils, results.streams){
  tryCatch(SEM.sensitivity(i, results.lakes, results.marine, results.soils, results.streams),
           error = function(e){
             xx = list()
             xx$sensitivity.all = NA
             xx$sensitivity.Marine = NA
             xx$sensitivity.Soils = NA
             xx$sensitivity.Streams = NA
             xx$sensitivity.Lakes = NA
             return(xx)
             })
}

##########################################################################
###### end of functions now different code options to run them ############


stack.vectors = function(name, combined.results) {
  as.data.frame(do.call(rbind, future_lapply(combined.results, `[[`, name)))
}

##### basic loop
n = 10
n = min(length(results.lakes), length(results.marine), length(results.soils), length(results.streams))
n = 1450

start_time <- Sys.time()
combined.results = list()
for (i in 1:n){
  combined.results[[i]] = SEM.wrapper(i, results.lakes, results.marine, results.soils, results.streams)
  rm(meta.Marine2)
  rm(meta.Soils2)
  rm(meta.Lakes2)
}
end.time = Sys.time()
timing = end.time - start_time
#check average duration of one iteration
timing
timing / n

combined.results <- combined.results[!sapply(combined.results, function(el) any(is.na(el$sensitivity.all)))]

#average time of succesful iterations:
timing / length(combined.results)

#proportion that returned NAs:
length(combined.results) / n

# xx.all = stack.vectors("sensitivity.all", xx)

sensitivity.all <- stack.vectors("sensitivity.all", combined.results)[1:1000,]
sensitivity.Marine <- stack.vectors("sensitivity.Marine", combined.results)[1:1000,]
sensitivity.Soils <- stack.vectors("sensitivity.Soils", combined.results)[1:1000,]
sensitivity.Streams <- stack.vectors("sensitivity.Streams", combined.results)[1:1000,]
sensitivity.Lakes <- stack.vectors("sensitivity.Lakes", combined.results)[1:1000,]

#save sensitivity data frames as RData files
save(sensitivity.all, sensitivity.Marine, sensitivity.Soils, sensitivity.Streams, sensitivity.Lakes, file = "sensitivities/results_sensitivity.Rdata")


#### Plot results ####
library(ggdist)
theme_set(theme_ggdist())

#All data
sensitivity.all.Pvalue <- sensitivity.all %>%
  pivot_longer(cols = starts_with("SEM.P") | starts_with("P_"), names_to = "parameter.output",  values_to = "parameter.value")

All.P = sensitivity.all.Pvalue %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  geom_vline(xintercept = 0.05, linewidth = .75, linetype = "dashed", colour = "red") + 
  labs(x = "P value", y = "") + 
  scale_y_discrete(labels = c("MaxTL->\npredation", "trophic complementarity->\npredation", "taxon richness->\npredation",
                              "taxon richness->\nprimary consumption", "SEM global P"))


sensitivity.all.R2 <- sensitivity.all %>%
  pivot_longer(cols = starts_with("r2."), names_to = "parameter.output",  values_to = "parameter.value")

All.R2 = sensitivity.all.R2 %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  labs(x = "R²", y = "") + 
  scale_y_discrete(labels = c("predation", "primary consumption"))


#Marine
sensitivity.Marine.Pvalue <- sensitivity.Marine %>%
  pivot_longer(cols = starts_with("SEM.P") | starts_with("P_"), names_to = "parameter.output",  values_to = "parameter.value")

Marine.P = sensitivity.Marine.Pvalue %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  geom_vline(xintercept = 0.05, linewidth = .75, linetype = "dashed", colour = "red") + 
  labs(x = "P value", y = "") + 
  scale_y_discrete(labels = c("MaxTL->\npredation", "taxon richness->\nprimary consumption", "SEM global P"))


sensitivity.Marine.R2 <- sensitivity.Marine %>%
  pivot_longer(cols = starts_with("r2."), names_to = "parameter.output",  values_to = "parameter.value")

Marine.R2 = sensitivity.Marine.R2 %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  labs(x = "R²", y = "") + 
  scale_y_discrete(labels = c("predation", "primary consumption"))


#Soils
sensitivity.Soils.Pvalue <- sensitivity.Soils %>%
  pivot_longer(cols = starts_with("SEM.P") | starts_with("P_"), names_to = "parameter.output",  values_to = "parameter.value")

Soils.P = sensitivity.Soils.Pvalue %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top")+
  geom_vline(xintercept = 0.05, linewidth = .75, linetype = "dashed", colour = "red") + 
  labs(x = "P value", y = "") + 
  scale_y_discrete(labels = c("MaxTL->\npredation", "taxon richness->\npredation", 
                              "taxon richness->\nprimary consumption", "SEM global P"))


sensitivity.Soils.R2 <- sensitivity.Soils %>%
  pivot_longer(cols = starts_with("r2."), names_to = "parameter.output",  values_to = "parameter.value")

Soils.R2 = sensitivity.Soils.R2 %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  labs(x = "R²", y = "") + 
  scale_y_discrete(labels = c("predation", "primary consumption"))


#Streams
sensitivity.Streams.Pvalue <- sensitivity.Streams %>%
  pivot_longer(cols = starts_with("SEM.P") | starts_with("P_"), names_to = "parameter.output",  values_to = "parameter.value")

Streams.P = sensitivity.Streams.Pvalue %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  geom_vline(xintercept = 0.05, linewidth = .75, linetype = "dashed", colour = "red") + 
  labs(x = "P value", y = "") + 
  scale_y_discrete(labels = c("MaxTL->\npredation", "trophic complementarity->\npredation", 
                              "taxon richness->\nprimary consumption", "SEM global P"))


sensitivity.Streams.R2 <- sensitivity.Streams %>%
  pivot_longer(cols = starts_with("r2."), names_to = "parameter.output",  values_to = "parameter.value")

Streams.R2 = sensitivity.Streams.R2 %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  labs(x = "R²", y = "") + 
  scale_y_discrete(labels = c("predation", "primary consumption"))


##Lakes
sensitivity.Lakes.Pvalue <- sensitivity.Lakes %>%
  pivot_longer(cols = starts_with("SEM.P") | starts_with("P_"), names_to = "parameter.output",  values_to = "parameter.value")

Lakes.P = sensitivity.Lakes.Pvalue %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  geom_vline(xintercept = 0.05, linewidth = .75, linetype = "dashed", colour = "red") + 
  labs(x = "P value", y = "") + 
  scale_y_discrete(labels = c("MaxTL->\nprimary consumption", "trophic complementarity->\npredation", 
                              "taxon richness->\npredation", "taxon richness->\nprimary consumption", "SEM global P"))


sensitivity.Lakes.R2 <- sensitivity.Lakes %>%
  pivot_longer(cols = starts_with("r2."), names_to = "parameter.output",  values_to = "parameter.value")

Lakes.R2 = sensitivity.Lakes.R2 %>% ggplot(aes(x = parameter.value, y = parameter.output)) + 
  stat_slabinterval(side = "top") +
  labs(x = "R²", y = "") + 
  scale_y_discrete(labels = c("predation", "primary consumption"))


library(gridExtra)

grid.arrange(
  All.P, All.R2, Marine.P, Marine.R2, Soils.P, Soils.R2, Streams.P, Streams.R2, Lakes.P, Lakes.R2,
  ncol = 2
)

