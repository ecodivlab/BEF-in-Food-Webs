#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Supplementary Information and Extended Data (figures and tables)
#
#=========================================================================================================================
 
# All models presented in Barnes et al. are run according to the described methods in the manuscript.
# The data accompanying this script are available on the Zenodo repository DOI: 
# Code was developed on R version 4.5.0

## Load packages ##
library(tidyverse); library(ggeffects); library(gridExtra); library(piecewiseSEM);
library(patchwork); library(nlme); library(grid); library(car); library(rempsyc); library(ggpattern);
library(ggh4x); library(scales); library(ggtext); library(ggrain); library(glmmTMB); 

## Clear environment and read ecosystem-specific food web data sets ##
rm(list=ls())
options(scipen = 999)

## To run this code, a local working directory must be set where all accompanying data and source code are lodged ##
setwd()


NPP.proxy <- read.csv("NDVI and Chlorophyll-a/data/proxy-npp.csv") ## import NDVI & Chl-a data
meta.Marine <- read.csv('meta.Marine.csv')  
meta.Marine <- meta.Marine %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(NPP.proxy = avg)
meta.Marine$NPP.proxy <- (meta.Marine$NPP.proxy - min(meta.Marine$NPP.proxy, na.rm = TRUE)) /
  (max(meta.Marine$NPP.proxy, na.rm = TRUE) - min(meta.Marine$NPP.proxy, na.rm = TRUE))
meta.Marine <- meta.Marine %>% mutate(NPP.scale = logit(NPP.proxy))
meta.Marine <- meta.Marine %>% mutate(NPP.scale2 = NPP.scale^2)
meta.Marine <- meta.Marine %>% mutate(stability.inv = 1-log10(stability))

meta.Soils <- read.csv('meta.Soils.csv')
meta.Soils <- meta.Soils %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(temperature_C = temperature, NPP.proxy = avg)
meta.Soils <- meta.Soils %>% mutate(NPP.scale = logit(NPP.proxy))
meta.Soils <- meta.Soils %>% mutate(NPP.scale2 = NPP.scale^2)
meta.Soils <- meta.Soils %>% mutate(stability.inv = 1-log10(stability))

meta.Streams <- read.csv('meta.Streams.csv')
meta.Streams <- meta.Streams %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(temperature_C = temperature, NPP.proxy = avg)
colnames(meta.Streams)[c(5,22)] <- c("temperature_C", "NPP.raw")
meta.Streams <- meta.Streams %>% mutate(NPP.proxy = ifelse(NPP.raw < 0, 0, NPP.raw), # replace negative value at Cananeia SP6 with 0
                                        NPP.scale = logit(NPP.proxy))
meta.Streams <- meta.Streams %>% mutate(NPP.scale2 = NPP.scale^2)
meta.Streams <- meta.Streams %>% mutate(stability.inv = 1-log10(stability))

meta.Lakes <- read.csv('meta.Lakes.csv')
meta.Lakes <- meta.Lakes %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
  rename(NPP.proxy = avg)
meta.Lakes <- meta.Lakes %>% mutate(NPP.scale = logit(NPP.proxy))
meta.Lakes <- meta.Lakes %>% mutate(NPP.scale2 = NPP.scale^2)
meta.Lakes <- meta.Lakes %>% mutate(stability.inv = 1-log10(stability))


## Compile data sets for cross-ecosystem analysis ##
commcols <- intersect(names(meta.Marine), names(meta.Soils))
commcols <- intersect(commcols, names(meta.Streams))
commcols <- intersect(commcols, names(meta.Lakes))


all_data <- bind_rows(select(meta.Marine, all_of(commcols)),
                      select(meta.Soils, all_of(commcols)),
                      select(meta.Streams, all_of(commcols)),
                      select(meta.Lakes, all_of(commcols)))

## graphics settings ##
set_theme(base=theme_classic(base_size = 10))
all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                  levels = c("Marine", "Soils", "Streams", "Lakes"))

#### Bivariate relationship total flux ~ taxon richness (Extended Data Fig. 1a) #####

## Total flux
S.flux_Global <- lme(log10(flux) ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data, na.action=na.omit)
plot(S.flux_Global, which=1)
qqnorm(S.flux_Global)
summary(S.flux_Global)

S.flux_Globala=update(S.flux_Global, random = ~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
                      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200)) #Best model
anova(S.flux_Global,S.flux_Globala)
plot(S.flux_Globala)
qqnorm(S.flux_Globala)
summary(S.flux_Globala)

S.flux_Globalb=update(S.flux_Global, weights=varPower(form=~S))
anova(S.flux_Global,S.flux_Globala,S.flux_Globalb)
plot(S.flux_Globalb)
qqnorm(S.flux_Globalb)
summary(S.flux_Globalb)

S.flux_Global = S.flux_Globalb


## Graph BEF for total flux 
Total_flux_Global=ggpredict(S.flux_Global, terms = "S")
S.total.sjp_global <- ggplot(Total_flux_Global, aes(x, predicted)) + 
  ylab(bquote('Total energy flux ('~J~day^-1*')')) + xlab("Multitrophic taxon richness") + 
  geom_line(linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = all_data, aes(x = S, y = flux, colour=ecosystem.type), size=2, alpha=0.4) +
  labs(colour = ("Ecosystem type")) +
  scale_y_continuous(breaks=breaks_log(n = 6), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="top", plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"), text=element_text(size=14)) 

## Produce Supplementary Figure 1 scatterplot ##
#ggsave("Supplementary fig 1 scatterplot.svg", S.total.sjp_global, width = 16, height = 16, units = "cm")







#### BEF relationships for p.c. predation and stability (Extended Data Fig. 2a) ####
#### per capita predation

PC.predation_Global = lme(log10(PC.predation) ~ log10(S), random=~1|study_ID, data=all_data, method='ML') 
  plot(PC.predation_Global, which=1)
  qqnorm(PC.predation_Global)
  summary(PC.predation_Global)
PC.predation_Globala = lme(log10(PC.predation) ~ log10(S), weights=varIdent(form=~1|study_ID), random=~1|ecosystem.type/study_ID, data=all_data, method='ML') #Best model
  plot(PC.predation_Globala, which=1)
  qqnorm(PC.predation_Globala)
  summary(PC.predation_Globala)
anova(PC.predation_Global, PC.predation_Globala)

PC.predationrate_Global=ggpredict(PC.predation_Globala, terms = "S")
PC.predation_Global <- ggplot(PC.predationrate_Global, aes(x, predicted)) + 
  ylab(expression("Predation" ~ "(" * italic("per capita") * ")")) + 
  geom_line(aes(linetype=group, color='#489FA759'), linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = all_data, aes(x = S, y = PC.predation, color='#489FA759'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", axis.text.x = element_blank(),
        plot.title = element_text(size = 14),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) +
  scale_colour_identity()

#### Stability
S.stability_Global <- lme(stability.inv ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data, na.action=na.omit)
plot(S.stability_Global, which=1)
qqnorm(S.stability_Global)
summary(S.stability_Global)

S.stability_Globala=update(S.stability_Global, random = ~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
                           control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200)) 
anova(S.stability_Global,S.stability_Globala)
plot(S.stability_Globala)
qqnorm(S.stability_Globala)
summary(S.stability_Globala)

S.stability_Globalb=update(S.stability_Global, weights=varComb(varIdent(form=~1|study_ID), varExp())) #Best model
anova(S.stability_Global,S.stability_Globala,S.stability_Globalb)
plot(S.stability_Globalb)
qqnorm(S.stability_Globalb)
summary(S.stability_Globalb)

S.stability_Globalc = glmmTMB(stability.inv ~ log10(S) + (1|ecosystem.type/study_ID), family = t_family(), data = all_data)
AIC(S.stability_Global,S.stability_Globala,S.stability_Globalb,S.stability_Globalc)
summary(S.stability_Globalc)

S.stability_Global = S.stability_Globalc
## Graph BEF for total flux 
Stability_Global=ggpredict(S.stability_Global, terms = "S")
S.StabilityTotal_Global <- ggplot(Stability_Global, aes(x, predicted)) + 
  ylab(expression(Stability~(log[10]~-lambda[max]))) + xlab("Taxon richness") + 
  geom_line(aes(linetype=group, color='#F5C2699E'), linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = all_data, aes(x = S, y = stability.inv, color='#F5C2699E'), size=2, alpha=0.5) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(x = "log10")+
  theme(legend.position="none", plot.margin = unit(c(4, 5.5, 5.5, 5.5), "pt")) +
  scale_colour_identity()

## Produce Supplementary Figure 3 ##
Richness_main_SuppInfo <- grid.arrange(patchworkGrob(PC.predation_Global / S.StabilityTotal_Global)) 

ggsave("Supplementary Figure 3 scatterplots.png", Richness_main_SuppInfo, width = 8, height = 10.5, units = "cm")





#### Supplementary Fig 2. Raincloud plot of ecosystem and food web properties ####
all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                  levels = c("Marine", "Soils", "Streams", "Lakes"))

NPP <- ggplot(all_data, aes(x = ecosystem.type, y = NPP.scale, fill = ecosystem.type)) +
  geom_rain(alpha = .3, cov='ecosystem.type', rain.side = 'r',
            boxplot.args.pos = list(width = .15, position = position_nudge(x = -.2), linewidth=0.25),
            violin.args.pos =  list(width = .8, position = position_nudge(x = .15), linewidth=0.25)) + 
  ylab('Net primary productivity') + xlab('') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) 

Temperature <- ggplot(all_data, aes(x = ecosystem.type, y = temperature_C, fill = ecosystem.type)) +
  geom_rain(alpha = .3, cov='ecosystem.type', rain.side = 'r',
            boxplot.args.pos = list(width = .15, position = position_nudge(x = -.2), linewidth=0.25),
            violin.args.pos =  list(width = .8, position = position_nudge(x = .15), linewidth=0.25)) + 
  ylab('Temperature (°C)') + xlab('') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"))

Diversity <- ggplot(all_data, aes(x = ecosystem.type, y = S, fill = ecosystem.type)) +
  geom_rain(alpha = .3, cov='ecosystem.type', rain.side = 'r',
            boxplot.args.pos = list(width = .15, position = position_nudge(x = -.2), linewidth=0.25),
            violin.args.pos =  list(width = .8, position = position_nudge(x = .15), linewidth=0.25)) + 
  ylab('Taxon richness') + xlab('') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"))

max.TL <- ggplot(all_data, aes(x = ecosystem.type, y = MaxTL, fill = 	ecosystem.type)) +
  geom_rain(alpha = .3, cov='ecosystem.type', rain.side = 'r',
            boxplot.args.pos = list(width = .15, position = position_nudge(x = -.2), linewidth=0.25),
            violin.args.pos =  list(width = .8, position = position_nudge(x = .15), linewidth=0.25)) + 
  ylab('Maximum trophic level') + xlab('Ecosystem type') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) 

prim.similarity <- ggplot(all_data, aes(x = ecosystem.type, y = sim.prim.cons, fill = 	ecosystem.type)) +
  geom_rain(alpha = .3, cov='ecosystem.type', rain.side = 'r',
            boxplot.args.pos = list(width = .15, position = position_nudge(x = -.2), linewidth=0.25),
            violin.args.pos =  list(width = .8, position = position_nudge(x = .15), linewidth=0.25)) +  
  ylab('1\u00B0 consumer dissimilarity') + xlab('Ecosystem type') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"))

pred.similarity <- ggplot(all_data, aes(x = ecosystem.type, y = sim.sec.cons, fill = 	ecosystem.type)) +
  geom_rain(alpha = .3, cov='ecosystem.type', rain.side = 'r',
            boxplot.args.pos = list(width = .15, position = position_nudge(x = -.2), linewidth=0.25),
            violin.args.pos =  list(width = .8, position = position_nudge(x = .15), linewidth=0.25)) + 
  ylab('Predator dissimilarity') + xlab('Ecosystem type') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"))

FW.Properties <- Richness_main <- grid.arrange(patchworkGrob((NPP | Temperature | Diversity)/
                                                               (max.TL | prim.similarity | pred.similarity))) 

ggsave(filename = "FW.Properties.png", FW.Properties, width = 200, height = 130, unit = "mm", dpi = 300)




#### Predation/prim. consumption relationships with environmental temperature and NPP ####

marine.a = ggplot(meta.Marine, aes(x = NPP.scale, y = log10(second.consumption))) +
  ylab('Predation') + xlab("NPP") + ggtitle("Marine") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") +
  theme(plot.title = element_text(face = "bold"))
marine.b = ggplot(meta.Marine, aes(x = temperature_C, y = log10(second.consumption))) +
  ylab('Predation') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") 
marine.c = ggplot(meta.Marine, aes(x = NPP.scale, y = log10(prim.consumption))) +
  ylab('Primary consumption') + xlab("NPP") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 
marine.d = ggplot(meta.Marine, aes(x = temperature_C, y = log10(prim.consumption))) +
  ylab('Primary consumption') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 

soils.a = ggplot(meta.Soils, aes(x = NPP.scale, y = log10(second.consumption))) +
  ylab('') + xlab("NPP") + ggtitle("Soils") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") +
  theme(plot.title = element_text(face = "bold"))
soils.b = ggplot(meta.Soils, aes(x = temperature_C, y = log10(second.consumption))) +
  ylab('') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") 
soils.c = ggplot(meta.Soils, aes(x = NPP.scale, y = log10(prim.consumption))) +
  ylab('') + xlab("NPP") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 
soils.d = ggplot(meta.Soils, aes(x = temperature_C, y = log10(prim.consumption))) +
  ylab('') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 

streams.a = ggplot(meta.Streams, aes(x = NPP.scale, y = log10(second.consumption))) +
  ylab('') + xlab("NPP") + ggtitle("Streams") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") +
  theme(plot.title = element_text(face = "bold")) 
streams.b = ggplot(meta.Streams, aes(x = temperature_C, y = log10(second.consumption))) +
  ylab('') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") 
streams.c = ggplot(meta.Streams, aes(x = NPP.scale, y = log10(prim.consumption))) +
  ylab('') + xlab("NPP") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 
streams.d = ggplot(meta.Streams, aes(x = temperature_C, y = log10(prim.consumption))) +
  ylab('') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 

lakes.a = ggplot(meta.Lakes, aes(x = NPP.scale, y = log10(second.consumption))) +
  ylab('') + xlab("NPP") + ggtitle("Lakes") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") +
  theme(plot.title = element_text(face = "bold")) 
lakes.b = ggplot(meta.Lakes, aes(x = temperature_C, y = log10(second.consumption))) +
  ylab('') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#C257579E") + geom_smooth(method = "gam", color="#C257579E") 
lakes.c = ggplot(meta.Lakes, aes(x = NPP.scale, y = log10(prim.consumption))) +
  ylab('') + xlab("NPP") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 
lakes.d = ggplot(meta.Lakes, aes(x = temperature_C, y = log10(prim.consumption))) +
  ylab('') + xlab("Temperature (°C)") +
  geom_point(alpha = 0.6, color="#3A67AE9E") + geom_smooth(method = "gam", color="#3A67AE9E") 

ecosystem_combined <- grid.arrange(patchworkGrob(
  (marine.a | soils.a | streams.a | lakes.a) /
  (marine.b | soils.b | streams.b | lakes.b) /
  (marine.c | soils.c | streams.c | lakes.c) /
  (marine.d | soils.d | streams.d | lakes.d)))


#### Reanalyse bivariate relationships with consumer density covariate ####

## Predation
S.predation_Global <- lme(log10(second.consumption) ~ log10(S) + log(density), random = ~1|ecosystem.type/study_ID, data=all_data)
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
S.prim.cons_Global <- lme(log10(prim.consumption) ~ log10(S) + log(density), random = ~1|ecosystem.type/study_ID, data=all_data)
plot(S.prim.cons_Global, which=1)
qqnorm(S.prim.cons_Global)
summary(S.prim.cons_Global)

S.prim.cons_Globala=update(S.prim.cons_Global, weights=varIdent(form=~1|study_ID))
anova(S.prim.cons_Global,S.prim.cons_Globala)
plot(S.prim.cons_Globala)
qqnorm(S.prim.cons_Globala)
summary(S.prim.cons_Globala)

S.prim.cons_Globalb=update(S.prim.cons_Global, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), 
                           control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200)) #Best model
anova(S.prim.cons_Global, S.prim.cons_Globala, S.prim.cons_Globalb)
plot(S.prim.cons_Globalb)
qqnorm(S.prim.cons_Globalb)
summary(S.prim.cons_Globalb)

S.prim.cons_Global = S.prim.cons_Globalb


#### Ecosystem-specific analyses 

## MARINE ##
S.predation_MAR <- gls(log10(second.consumption) ~ log10(S) + log(density), data=meta.Marine)
plot(S.predation_MAR, which=1)
qqnorm(S.predation_MAR)
summary(S.predation_MAR)

S.predation_MARa=update(S.predation_MAR, weights=varIdent(form=~1|study_ID))
anova(S.predation_MAR,S.predation_MARa)
plot(S.predation_MARa)
qqnorm(S.predation_MARa)
summary(S.predation_MARa)

S.predation_MARb=update(S.predation_MAR, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S))) #Best model
anova(S.predation_MAR,S.predation_MARa,S.predation_MARb)
plot(S.predation_MARb)
qqnorm(S.predation_MARb)
summary(S.predation_MARb)

S.predation_MAR <- S.predation_MARb

S.prim.cons_MAR <- gls(log10(prim.consumption) ~ log10(S) + log(density), data=meta.Marine) #Best model
plot(S.prim.cons_MAR, which=1)
qqnorm(S.prim.cons_MAR)
summary(S.prim.cons_MAR)

# S.prim.cons_MARa=update(S.prim.cons_MAR, weights=varIdent(form=~1|study_ID))
# anova(S.prim.cons_MAR,S.prim.cons_MARa)
# plot(S.prim.cons_MARa)
# qqnorm(S.prim.cons_MARa)
# summary(S.prim.cons_MARa)
# 
# S.prim.cons_MARb=update(S.prim.cons_MAR, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)))
# anova(S.prim.cons_MAR,S.prim.cons_MARa,S.prim.cons_MARb)
# plot(S.prim.cons_MARb)
# qqnorm(S.prim.cons_MARb)
# summary(S.prim.cons_MARb)




## SOIL ##
S.predation_SOIL <- gls(log10(second.consumption) ~ log10(S) + log(density), data=meta.Soils)
plot(S.predation_SOIL, which=1)
qqnorm(S.predation_SOIL)
summary(S.predation_SOIL)

S.predation_SOILa=update(S.predation_SOIL, weights=varIdent(form=~1|study_ID)) #Best model
anova(S.predation_SOIL,S.predation_SOILa)
plot(S.predation_SOILa)
qqnorm(S.predation_SOILa)
summary(S.predation_SOILa)

S.predation_SOIL <- S.predation_SOILa


S.prim.cons_SOIL <- gls(log10(prim.consumption) ~ log10(S) + log(density), data=meta.Soils)
plot(S.prim.cons_SOIL, which=1)
qqnorm(S.prim.cons_SOIL)
summary(S.prim.cons_SOIL)

S.prim.cons_SOILa=update(S.prim.cons_SOIL, weights=varIdent(form=~1|study_ID)) #Best model
anova(S.prim.cons_SOIL,S.prim.cons_SOILa)
plot(S.prim.cons_SOILa)
qqnorm(S.prim.cons_SOILa)
summary(S.prim.cons_SOILa)

S.prim.cons_SOIL <- S.prim.cons_SOILa



## LAKES ## 
S.predation_LAKE <- gls(log10(second.consumption) ~ log10(S) + log(density), data=meta.Lakes)
plot(S.predation_LAKE, which=1)
qqnorm(S.predation_LAKE)
summary(S.predation_LAKE)

S.predation_LAKEa=update(S.predation_LAKE, weights=varPower(form=~S))
anova(S.predation_LAKE,S.predation_LAKEa)
plot(S.predation_LAKEa)
qqnorm(S.predation_LAKEa)
summary(S.predation_LAKEa)

S.predation_LAKE <- S.predation_LAKEa


S.prim.cons_LAKE <- gls(log10(prim.consumption) ~ log10(S) + log(density), data=meta.Lakes) #Best model
plot(S.prim.cons_LAKE, which=1)
qqnorm(S.prim.cons_LAKE)
summary(S.prim.cons_LAKE)

S.prim.cons_LAKEa=update(S.prim.cons_LAKE, weights=varPower(form=~S))
anova(S.prim.cons_LAKE,S.prim.cons_LAKEa)
plot(S.prim.cons_LAKEa)
qqnorm(S.prim.cons_LAKEa)
summary(S.prim.cons_LAKEa)

S.prim.cons_LAKEb=update(S.prim.cons_LAKE, weights=varIdent(form=~1|study_ID)) #Best model
anova(S.prim.cons_LAKE,S.prim.cons_LAKEa, S.prim.cons_LAKEb)
plot(S.prim.cons_LAKEb)
qqnorm(S.prim.cons_LAKEb)
summary(S.prim.cons_LAKEb)


## STREAMS ##
# S.predation_STREAM <- gls(log10(second.consumption) ~ log10(S), data=meta.Streams)
# plot(S.predation_STREAM, which=1)
# qqnorm(S.predation_STREAM)
# summary(S.predation_STREAM)

S.predation_STREAM.mm <- lme(log10(second.consumption) ~ log10(S) + log(density), random=~1|study_ID, data=meta.Streams) #Best model
plot(S.predation_STREAM.mm, which=1)
qqnorm(S.predation_STREAM.mm)
summary(S.predation_STREAM.mm)

# anova(S.predation_STREAM.mm, S.predation_STREAM)


# S.prim.cons_STREAM <- gls(log10(prim.consumption) ~ log10(S), data=meta.Streams)
# plot(S.prim.cons_STREAM, which=1)
# qqnorm(S.prim.cons_STREAM)
# summary(S.prim.cons_STREAM)

S.prim.cons_STREAM.mm <- lme(log10(prim.consumption) ~ log10(S) + log(density), random=~1|study_ID, data=meta.Streams) #Best model
plot(S.prim.cons_STREAM.mm, which=1)
qqnorm(S.predation_STREAM.mm)
summary(S.prim.cons_STREAM.mm)




#####################################
#### Structural Equation Models #####
#####################################


## Create vectors for effect size plotting ##
flux <- c(rep(c(rep("stability", 3), rep("PC.predation", 3)),2))
FW_prop <- rep(c("Taxon richness","Max TL","Trophic dissim."), 4)
effect.type <- factor(c(rep("direct", 6),rep("indirect", 6)), levels=c('indirect', 'direct'))

#### Cross-ecosystem analysis with total flux (Extended Data Fig. 1b) ####


mod2 = lme(logit(sim.total) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|study_ID, data=all_data)
  plot(mod2, which=1)
  qqnorm(mod2)
  summary(mod2)


mod4 = lme(log(flux) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|study_ID, data=all_data, method = 'ML') #Best model
  plot(mod4, which=1)
  qqnorm(mod4)
  summary(mod4)
mod4a = lme(log(flux) ~ log(MaxTL) + log(S) + NPP.scale, weights=varComb(varIdent(form=~1|study_ID)), random=~1|study_ID, data=all_data, method = 'ML') #Best model
  plot(mod4a, which=1)
  qqnorm(mod4a)
  summary(mod4a)
mod4b = lme(log(flux) ~ log(MaxTL) + log(S) + NPP.scale, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~MaxTL)), random=~1|study_ID, data=all_data, method = 'ML') #Best model
  plot(mod4b, which=1)
  qqnorm(mod4b)
  summary(mod4b)
mod4c = lme(log(flux) ~ log(MaxTL) + log(S) + poly(NPP.scale,2,raw=TRUE), weights=varComb(varIdent(form=~1|study_ID), varExp(form=~MaxTL)), random=~1|study_ID, data=all_data, method = 'ML') #Best model
  plot(mod4c, which=1)
  qqnorm(mod4c)
  summary(mod4c)

anova(mod4, mod4a, mod4b, mod4c)


### Maximal model
SEM.all.total <- psem(
  
  lme(log(S) ~ NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(logit(sim.total) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), random=~1|study_ID, data=all_data),
  lme(log(flux) ~ logit(sim.total) + log(MaxTL), weights=varComb(varIdent(form=~1|study_ID)), random=~1|study_ID, data=all_data)
)

summary(SEM.all.total)

### Min adequate model
SEM.all.total2 <- psem(
  
  lme(log(S) ~ NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(logit(sim.total) ~ log(MaxTL), random=~1|study_ID, data=all_data),
  lme(log(flux) ~ log(MaxTL) + log(S) + NPP.scale + NPP.scale2, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~MaxTL)), random=~1|study_ID, data=all_data)
)

summary(SEM.all.total2)
results.all.totalSEM <- summary(SEM.all.total2)$coefficients[,c(1:5, 8, 7)]
names(results.all.totalSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
fun <- function(x) {
  formatC(x, format = "f", digits = 3)
}
all.totalSEM_table <- nice_table(results.all.totalSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(all.totalSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/all.totalSEM_table.docx")




####calculate std. effect size for quadratic variables sensu Henseler et al. 2012 (https://doi.org/10.1057/ejis.2011.36)
#Refit SEM.all2 without NPP.scale quadratic variable
SEM.all.total2_no.NPP_R2 <- rsquared(psem(
  lme(log(S) ~ 1, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(MaxTL) ~ log(S), random=~1|ecosystem.type/study_ID, data=all_data),
  lme(logit(sim.total) ~ log(MaxTL), random=~1|study_ID, data=all_data),
  lme(log(flux) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varExp(form=~MaxTL)), random=~1|study_ID, data=all_data)
))
SEM.all.total2_R2 <- rsquared(SEM.all.total2)

#Calculate Cohen's f² = (R²full - R²reduced) / (1 - R²full)
SEM.all.total2_Q.Std.Est <- data.frame(SEM.all.total2_R2[,1],
                                 (SEM.all.total2_R2[,6] - SEM.all.total2_no.NPP_R2[,6]) / (1 - SEM.all.total2_R2[,6]),
                                 row.names=NULL)
colnames(SEM.all.total2_Q.Std.Est) <- c("response", "F^2")






#### Test stream SEM with density ####
SEM.Streams.density <- psem(
  
  lme(log(S) ~ NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(log(density) ~ NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(log(MaxTL) ~ log(S), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(logit(sim.prim.cons) ~ log(MaxTL), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(logit(sim.sec.cons) ~ log(S) + log(density), random=~1|study_ID, weights=varExp(), data=meta.Streams),
  lme(log(prim.consumption) ~ log(density) + logit(sim.prim.cons), random=~1|study_ID, data=meta.Streams),
  lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons) + log(density), random=~1|study_ID, data=meta.Streams),
  
  log(S) %~~% log(density),
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Streams.density)

results.streamsSEM <- summary(SEM.Streams.density)$coefficients[c(1:11, 16),c(1:5, 8, 7)]
names(results.streamsSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
streamsSEM_table <- nice_table(results.streamsSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(streamsSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/streamsSEM_table with density.docx")


## Create results dataframe for summary boxplots ##
std.effect <- c(
  SR.direct.pred <- 0,
  MTL.direct.pred <- results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(second.consumption)' & results.streamsSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.pred <- results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(second.consumption)' & results.streamsSEM$Predictor == 'logit(sim.sec.cons)'],
  SR.direct.prim <- 0,
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  
  SR.indirect.pred <- 
    results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'logit(sim.sec.cons)' & results.streamsSEM$Predictor == 'log(S)'] *
    results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(second.consumption)' & results.streamsSEM$Predictor == 'logit(sim.sec.cons)'],
  MTL.indirect.pred <- 0,
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- 0,
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0
)

flux <- c(rep(c(rep("predation", 3), rep("primary consumption",3)),2))
FW_prop <- rep(c("Taxon richness","Max TL","Trophic dissim."), 4)
effect.type <- factor(c(rep("direct", 6),rep("indirect", 6)), levels=c('indirect', 'direct'))

eff.table_streams <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot
streams.effects.density <- ggplot(eff.table_streams, 
                          aes(x = flux, y = std.effect, fill = flux, pattern = effect.type)) + 
  ylab('Effect size') + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5) +
  geom_bar_pattern(stat = "identity", position = "stack", pattern_color = 'White', pattern_fill = "white", width = .98,
                   pattern_alpha = 0.7, pattern_angle = 55, pattern_density = 0.75, pattern_spacing = 0.1, pattern_size=0) +
  scale_fill_manual(values = c("#C257579E", "#3A67AE9E")) + 
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  coord_flip(ylim=c(-0.5, 0.5)) + scale_y_continuous(labels = c(-0.5, '',0.0, '',0.5)) +
  facet_grid(fct_relevel(FW_prop,'Taxon richness', 'Max TL', 'Trophic dissim.')~., scales = "free_y", 
             labeller = label_wrap_gen(width=10), switch = "y") +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), 
        strip.clip = "off", strip.background = element_blank(), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

ggsave("Streams effects_density.png", streams.effects.density, width = 6, height = 9, units = "cm", bg='transparent')




#### Cross-ecosystem SEM with stability and PC.predation (Extended Data Fig. 1b and c) ####

mod.PC.predation = lme(log(PC.predation) ~ logit(sim.sec.cons) + log(MaxTL), random=~1|study_ID, data=all_data, method='ML') 
  plot(mod.PC.predation, which=1)
  qqnorm(mod.PC.predation)
  summary(mod.PC.predation)
mod.PC.predationa = lme(log(PC.predation) ~ logit(sim.sec.cons) + log(MaxTL), weights=varIdent(form=~1|study_ID), random=~1|ecosystem.type/study_ID, data=all_data, method='ML') #Best model
  plot(mod.PC.predationa, which=1)
  qqnorm(mod.PC.predationa)
  summary(mod.PC.predationa)
anova(mod.PC.predation, mod.PC.predationa)

modStability = lme(stability.inv ~ log(PC.predation) + log(prim.consumption) + logit(sim.prim.cons), random=~1|study_ID, data=all_data, method='ML') 
  plot(modStability, which=1)
  qqnorm(modStability)
  summary(modStability)
modStabilitya = lme(stability.inv ~ log(PC.predation) + log(prim.consumption) + logit(sim.prim.cons), weights=varComb(varPower(form=~S), varExp()), random=~1|ecosystem.type/study_ID, data=all_data, method='ML') # Best model
  plot(modStabilitya, which=1)
  qqnorm(modStabilitya)
  summary(modStabilitya)
modStabilityb = lme(stability.inv ~ log(PC.predation) + log(prim.consumption) + logit(sim.prim.cons), weights=varIdent(form=~1|ecosystem.type), random=~1|ecosystem.type/study_ID, data=all_data, method='ML') 
  plot(modStabilityb, which=1)
  qqnorm(modStabilityb)
  summary(modStabilityb)
modStabilityc = lme(stability.inv ~ log(PC.predation) + log(prim.consumption) + logit(sim.prim.cons), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), random=~1|ecosystem.type/study_ID, data=all_data, method='ML') 
  plot(modStabilityc, which=1)
  qqnorm(modStabilityc)
  hist(modStabilityc$residuals)
  summary(modStabilityc)
modStabilityd = glmmTMB(stability.inv ~ log(PC.predation) + log(prim.consumption) + logit(sim.prim.cons) + (1|ecosystem.type/study_ID), family = t_family(), data = all_data, REML = FALSE)
  summary(modStabilityd)
AIC(modStability, modStabilitya, modStabilityb, modStabilityc,modStabilityd)


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
  lme(log(PC.predation) ~ logit(sim.sec.cons) + log(MaxTL), weights=varPower(form=~S), random=~1|ecosystem.type/study_ID, data=all_data),
  glmmTMB(stability.inv ~ log(PC.predation) + log(prim.consumption) + logit(sim.prim.cons) + (1|ecosystem.type/study_ID), family = t_family(), data = all_data),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption),
  log(PC.predation) %~~% log(prim.consumption),
  log(PC.predation) %~~% log(second.consumption)
)

summary(SEM.all)


### Min adequate model
SEM.all2 <- psem(
  lme(log(S) ~ NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), data=all_data),
  lme(logit(sim.sec.cons) ~ NPP.scale + log(MaxTL) + log(S), random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, 
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale + logit(sim.sec.cons), random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(PC.predation) ~ logit(sim.sec.cons) + log(MaxTL) + NPP.scale, weights=varIdent(form=~1|study_ID), random=~1|ecosystem.type/study_ID,
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  glmmTMB(stability.inv ~ log(S) + logit(sim.prim.cons) + log(PC.predation) + (1|ecosystem.type/study_ID), family = t_family(),
          start = list(psi = log(4)), map = list(psi = factor(NA)), data = all_data),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption),
  log(PC.predation) %~~% log(prim.consumption),
  log(PC.predation) %~~% log(second.consumption)
)

summary(SEM.all2)


####Calculate standardised path coefficients for student-t regressions following (Grace 2006: standardized coef = beta * SDx/SDy)
stability.mod <- glmmTMB(stability.inv ~ log(S) + logit(sim.prim.cons) + log(PC.predation) + (1|ecosystem.type/study_ID), family = t_family(),
        start = list(psi = log(4)), map = list(psi = factor(NA)), data = all_data)

betas_S <- fixef(stability.mod)$cond["log(S)"]       
betas_sim.prim.cons <- fixef(stability.mod)$cond["logit(sim.prim.cons)"] 
betas_PC.predation <- fixef(stability.mod)$cond["log(PC.predation)"] 

sd_x_S <- sd(log(all_data$S)) # SD of predictors 
sd_x_sim.prim.cons <- sd(logit(all_data$sim.prim.cons)) 
sd_x_PC.predation <- sd(log(all_data$PC.predation))  

lp <- predict(stability.mod, type = "link") 
sd_y <- sd(lp)  # SD of response on link scale (Lefcheck 2016)
beta_std_S <- betas_S * (sd_x_S / sd_y)  
betas_std_sim.prim.cons <- betas_sim.prim.cons * (sd_x_sim.prim.cons / sd_y)
betas_std_PC.predation <- betas_PC.predation * (sd_x_PC.predation / sd_y)

results.allSEM <- summary(SEM.all2)$coefficients[c(1:22),c(1:5, 8, 7)]
names(results.allSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
fun <- function(x) {
  formatC(x, format = "f", digits = 3)
}
allSEM_table_stability <- nice_table(results.allSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

flextable::save_as_docx(allSEM_table_stability, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/allSEM_table_stability.docx")


####calculate std. effect size for quadratic variables sensu Henseler et al. 2012 (https://doi.org/10.1057/ejis.2011.36)
#Refit SEM.all2 without NPP.scale quadratic variable
SEM.all.total2_no.NPP_R2 <- rsquared(psem(
  lme(log(S) ~ 1, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(MaxTL) ~ log(S), random=~1|ecosystem.type/study_ID, data=all_data),
  lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), data=all_data),
  lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S), random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, 
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale + logit(sim.sec.cons), random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(PC.predation) ~ logit(sim.sec.cons) + log(MaxTL), weights=varPower(form=~S), random=~1|ecosystem.type/study_ID,
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  glmmTMB(stability.inv ~ log(S) + logit(sim.prim.cons) + log(PC.predation) + (1|ecosystem.type/study_ID), family = t_family(),
          start = list(psi = log(4)), map = list(psi = factor(NA)), data = all_data),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption),
  log(PC.predation) %~~% log(prim.consumption),
  log(PC.predation) %~~% log(second.consumption)
))
SEM.all.total2_R2 <- rsquared(SEM.all2)

#Calculate Cohen's f² = (R²full - R²reduced) / (1 - R²full)
SEM.all.total2_Q.Std.Est <- data.frame(SEM.all.total2_R2[,1],
                                       (SEM.all.total2_R2[,6] - SEM.all.total2_no.NPP_R2[,6]) / (1 - SEM.all.total2_R2[,6]),
                                       row.names=NULL)
colnames(SEM.all.total2_Q.Std.Est) <- c("response", "F^2")



#### Create results dataframe for summary boxplots 
results.allSEM$`Std. Estimate` <- as.numeric(results.allSEM$`Std. Estimate`)
std.effect <- c(
  SR.direct.stab <- beta_std_S,
  MTL.direct.stab <- 0,
  SIM.direct.stab <- betas_std_sim.prim.cons,
  SR.direct.PCpred <- 0,
  MTL.direct.PCpred <- results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.PCpred <- results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'],
  
  SR.indirect.stab <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(MaxTL)' & results.allSEM$Predictor == 'log(S)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'log(MaxTL)'] *
    betas_std_PC.predation +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(MaxTL)' & results.allSEM$Predictor == 'log(S)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(MaxTL)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'] *
    betas_std_PC.predation +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(S)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'] *
    betas_std_PC.predation +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.prim.cons)' & results.allSEM$Predictor == 'log(S)'] *
    betas_std_sim.prim.cons,
  MTL.indirect.stab <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'log(MaxTL)'] *
    betas_std_PC.predation +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(MaxTL)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'] *
    betas_std_PC.predation,
  SIM.indirect.stab <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'] *
    betas_std_PC.predation,
  SR.indirect.PCpred <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(MaxTL)' & results.allSEM$Predictor == 'log(S)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'log(MaxTL)'] +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(MaxTL)' & results.allSEM$Predictor == 'log(S)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(MaxTL)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'] +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(S)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'],
  MTL.indirect.PCpred <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(MaxTL)'] *
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(PC.predation)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'],
  SIM.indirect.PCpred <- 0
)
eff.table_all <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot for Fig. 3 ##
global.effects <- ggplot(eff.table_all, 
                         aes(x = flux, y = std.effect, fill = flux, pattern = effect.type)) + 
  ylab('Effect size') + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5) +
  geom_bar_pattern(stat = "identity", position = "stack", pattern_color = 'White', pattern_fill = "white", width = .98,
                   pattern_alpha = 0.7, pattern_angle = 55, pattern_density = 0.75, pattern_spacing = 0.1, pattern_size=0) +
  scale_fill_manual(values = c("#489FA759", "#F5C2699E")) + 
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  coord_flip(ylim=c(-0.5, 0.5)) + scale_y_continuous(labels = c(-0.5, '',0.0, '',0.5)) +
  facet_grid(fct_relevel(FW_prop,'Taxon richness', 'Max TL', 'Trophic dissim.')~., scales = "free_y", labeller = label_wrap_gen(width=10)) +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), 
        strip.clip = "off", strip.background = element_blank(), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))


#ggsave("Global effects Fig S3.png", global.effects, width = 6, height = 9, units = "cm", bg='transparent')






