#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Main text analyses and figures
#
#=========================================================================================================================

# This code generates bivariate and SEM models testing for the effect of taxa richness, maximum trophic level, and
# trophic dissimilarity on total primary consumption and predation rates (derived from estimated energy fluxes) in
# 318 food webs, spanning marine, lake, stream, and soil food webs. 
# 
# All models presented in Barnes et al. are run according to the described methods in the manuscript.
# The data accompanying this script are available on the Zenodo repository DOI: 
# Code was developed on R version 4.5.0

## Load packages ##
library(tidyverse); library(ggeffects); library(gridExtra); library(piecewiseSEM);
library(patchwork); library(nlme); library(grid); library(car); library(rempsyc); library(ggpattern);
library(ggh4x); library(scales); library(ggtext); library(ggrain)

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
  
meta.Soils <- read.csv('meta.Soils.csv')
  meta.Soils <- meta.Soils %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
    rename(temperature_C = temperature, NPP.proxy = avg)
  meta.Soils <- meta.Soils %>% mutate(NPP.scale = logit(NPP.proxy))
  meta.Soils <- meta.Soils %>% mutate(NPP.scale2 = NPP.scale^2)
                                      
meta.Streams <- read.csv('meta.Streams.csv')
  meta.Streams <- meta.Streams %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name") %>%
    rename(temperature_C = temperature, NPP.proxy = avg)
  colnames(meta.Streams)[c(5,22)] <- c("temperature_C", "NPP.raw")
  meta.Streams <- meta.Streams %>% mutate(NPP.proxy = ifelse(NPP.raw < 0, 0, NPP.raw), # replace negative value at Cananeia SP6 with 0
                  NPP.scale = logit(NPP.proxy))
  meta.Streams <- meta.Streams %>% mutate(NPP.scale2 = NPP.scale^2)
  
meta.Lakes <- read.csv('meta.Lakes.csv')
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


#write.csv(all_data, file="Master dataset_FuSED.csv", row.names = F)

## Set ggplot2 theme ##
set_theme(base=theme_classic(base_size = 10))


############################################################################
#### Construct linear models and graphs for bivariate BEF relationships ####
############################################################################


#### Cross-ecosystem analyses #####

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

# S.prim.cons_Globalb=update(S.prim.cons_Global, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)))
# anova(S.prim.cons_Global,S.prim.cons_Globala,S.prim.cons_Globalb)
# plot(S.prim.cons_Globalb)
# qqnorm(S.prim.cons_Globalb)
# summary(S.prim.cons_Globalb)

S.prim.cons_Global = S.prim.cons_Globala

all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                  levels = c("Marine", "Soils", "Streams", "Lakes"))


## Graph BEF for predation and primary consumption 
predation_global=ggpredict(S.predation_Global, terms = "S")
# calculate predicted increase in fluxes with S
abs.increase = round(max(predation_global$predicted) - min(predation_global$predicted)) #predation
percent.increase = round((max(predation_global$predicted) - min(predation_global$predicted)) / 
                           min(predation_global$predicted) * 100) #predation

S.predation.sjp_global <- ggplot(predation_global, aes(x, predicted)) + 
  ylab('Predation') + 
  geom_line(aes(linetype=group, color='#C257579E'), linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = all_data, aes(x = S, y = second.consumption, color='#C257579E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", axis.text.x = element_blank(),
        plot.title = element_text(size = 14),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) +
  scale_colour_identity()

prim.cons_global=ggpredict(S.prim.cons_Global, terms = "S")
# calculate predicted increase in fluxes with S
abs.increase = round(max(prim.cons_global$predicted) - min(prim.cons_global$predicted)) #primary consumption
percent.increase = round((max(prim.cons_global$predicted) - min(prim.cons_global$predicted)) / 
                           min(prim.cons_global$predicted) * 100) #primary consumption

S.prim.cons.sjp_global <- ggplot(prim.cons_global, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color='#3A67AE9E'), linewidth=1.5) + 
  ylab('Primary consumption') + xlab("Taxon richness") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = all_data, aes(x = S, y = prim.consumption, color='#3A67AE9E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", plot.margin = unit(c(4, 5.5, 5.5, 5.5), "pt")) +
  scale_colour_identity()


#(max(prim.cons_global$predicted)-min(prim.cons_global$predicted))/min(prim.cons_global$predicted)*100
#(max(predation_global$predicted)-min(predation_global$predicted))/min(predation_global$predicted)*100


## Produce Figure 3 ##
Richness_main <- grid.arrange(patchworkGrob(S.predation.sjp_global / S.prim.cons.sjp_global)) 

ggsave("Figure 3 scatterplots.svg", Richness_main, width = 8, height = 10.5, units = "cm")




#### Ecosystem-specific analyses #####

## MARINE ##
S.predation_MAR <- gls(log10(second.consumption) ~ log10(S), data=meta.Marine)
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

# S.predation_MAR.mm <- lme(log10(second.consumption) ~ log10(S), random=~1|study_ID, data=meta.Marine)
# plot(S.predation_MAR.mm, which=1)
# qqnorm(S.predation_MAR.mm)
# summary(S.predation_MAR.mm)
# 
# anova(S.predation_MAR.mm, S.predation_MAR, S.predation_MARa, S.predation_MARb)

S.prim.cons_MAR <- gls(log10(prim.consumption) ~ log10(S), data=meta.Marine) #Best model
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
# 
# S.prim.cons_MAR.mm <- lme(log10(prim.consumption) ~ log10(S), random=~1|study_ID, data=meta.Marine)
# plot(S.prim.cons_MAR.mm, which=1)
# qqnorm(S.prim.cons_MAR)
# summary(S.prim.cons_MAR.mm)
# 
# anova(S.prim.cons_MAR.mm, S.prim.cons_MAR, S.prim.cons_MARa, S.prim.cons_MARb)


predation_MAR=ggpredict(S.predation_MAR, terms = "S")
S.predation.sjp_MAR <- ggplot(predation_MAR, aes(x, predicted)) + 
  ylab('Predation') + ggtitle("Marine (n = 131)") +
  geom_line(aes(linetype=group, color='#C257579E'), linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Marine, aes(x = S, y = second.consumption, color='#C257579E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", axis.text.x = element_blank(),
        plot.title = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) +
  scale_colour_identity()


prim.cons_MAR=ggpredict(S.prim.cons_MAR, terms = "S")
S.prim.cons.sjp_MAR <- ggplot(prim.cons_MAR, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color='#3A67AE9E'), linewidth=1.5) + 
  ylab('Primary<br>consumption') + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Marine, aes(x = S, y = prim.consumption, color='#3A67AE9E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", plot.margin = unit(c(4, 5.5, 0, 5.5), "pt"), 
        axis.title.x = element_blank(), axis.title.y = element_markdown()) +
  scale_colour_identity()


## SOIL ##
S.predation_SOIL <- gls(log10(second.consumption) ~ log10(S), data=meta.Soils)
plot(S.predation_SOIL, which=1)
qqnorm(S.predation_SOIL)
summary(S.predation_SOIL)

S.predation_SOILa=update(S.predation_SOIL, weights=varIdent(form=~1|study_ID)) #Best model
anova(S.predation_SOIL,S.predation_SOILa)
plot(S.predation_SOILa)
qqnorm(S.predation_SOILa)
summary(S.predation_SOILa)

S.predation_SOIL <- S.predation_SOILa

# S.predation_SOILb=update(S.predation_SOILa, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S))) #Best model
# anova(S.predation_SOIL,S.predation_SOILa, S.predation_SOILb)
# plot(S.predation_SOILb)
# qqnorm(S.predation_SOILb)
# summary(S.predation_SOILb)

# S.predation_SOIL.mm <- lme(log10(second.consumption) ~ log10(S), random=~1|study_ID, data=meta.Soils)
# plot(S.predation_SOIL.mm, which=1)
# qqnorm(S.predation_SOIL.mm)
# summary(S.predation_SOIL.mm)
# 
# anova(S.predation_SOIL.mm, S.predation_SOIL, S.predation_SOILa)

S.prim.cons_SOIL <- gls(log10(prim.consumption) ~ log10(S), data=meta.Soils)
plot(S.prim.cons_SOIL, which=1)
qqnorm(S.prim.cons_SOIL)
summary(S.prim.cons_SOIL)

S.prim.cons_SOILa=update(S.prim.cons_SOIL, weights=varIdent(form=~1|study_ID)) #Best model
anova(S.prim.cons_SOIL,S.prim.cons_SOILa)
plot(S.prim.cons_SOILa)
qqnorm(S.prim.cons_SOILa)
summary(S.prim.cons_SOILa)

S.prim.cons_SOIL <- S.prim.cons_SOILa

# S.prim.cons_SOILb=update(S.prim.cons_SOILa, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S))) #Best model
# anova(S.prim.cons_SOIL,S.prim.cons_SOILa, S.prim.cons_SOILb)
# plot(S.prim.cons_SOILb)
# qqnorm(S.prim.cons_SOILb)
# summary(S.prim.cons_SOILb)

# S.prim.cons_SOIL.mm <- lme(log10(prim.consumption) ~ log10(S), random=~1|study_ID, data=meta.Soils)
# plot(S.prim.cons_SOIL.mm, which=1)
# qqnorm(S.prim.cons_SOIL.mm)
# summary(S.prim.cons_SOIL.mm)
# 
# anova(S.prim.cons_SOIL.mm, S.prim.cons_SOIL, S.prim.cons_SOILa)

predation_SOIL=ggpredict(S.predation_SOIL, terms = "S")
S.predation.sjp_SOIL <- ggplot(predation_SOIL, aes(x, predicted)) +  ggtitle("Soils (n = 65)") +
  geom_line(aes(linetype=group, color='#C257579E'), linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Soils, aes(x = S, y = second.consumption, color='#C257579E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=9), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", axis.text.x = element_blank(),
        plot.title = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) +
  scale_colour_identity()

prim.cons_SOIL=ggpredict(S.prim.cons_SOIL, terms = "S")
S.prim.cons.sjp_SOIL <- ggplot(prim.cons_SOIL, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color='#3A67AE9E'), linewidth=1.5) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Soils, aes(x = S, y = prim.consumption, color='#3A67AE9E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=9), expand = expansion(mult = c(0.05, 0.1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", plot.margin = unit(c(4, 5.5, 0, 5.5),"pt"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_colour_identity()


## LAKES ## 
S.predation_LAKE <- gls(log10(second.consumption) ~ log10(S), data=meta.Lakes)
plot(S.predation_LAKE, which=1)
qqnorm(S.predation_LAKE)
summary(S.predation_LAKE)

S.predation_LAKEa=update(S.predation_LAKE, weights=varPower(form=~S))
anova(S.predation_LAKE,S.predation_LAKEa)
plot(S.predation_LAKEa)
qqnorm(S.predation_LAKEa)
summary(S.predation_LAKEa)

S.predation_LAKEb=update(S.predation_LAKE, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S))) #Best model
anova(S.predation_LAKE,S.predation_LAKEa, S.predation_LAKEb)
plot(S.predation_LAKEb)
qqnorm(S.predation_LAKEb)
summary(S.predation_LAKEb)

S.predation_LAKE <- S.predation_LAKEb

S.prim.cons_LAKE <- gls(log10(prim.consumption) ~ log10(S), data=meta.Lakes) #Best model
plot(S.prim.cons_LAKE, which=1)
qqnorm(S.prim.cons_LAKE)
summary(S.prim.cons_LAKE)

# S.prim.cons_LAKEa=update(S.prim.cons_LAKE, weights=varPower(form=~S))
# anova(S.prim.cons_LAKE,S.prim.cons_LAKEa)
# plot(S.prim.cons_LAKEa)
# qqnorm(S.prim.cons_LAKEa)
# summary(S.prim.cons_LAKEa)

predation_LAKE=ggpredict(S.predation_LAKE, terms = "S")
S.predation.sjp_LAKE <- ggplot(predation_LAKE, aes(x, predicted)) +  ggtitle("Lakes (n = 48)") +
  geom_line(aes(linetype=group, color='#C257579E'), linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Lakes, aes(x = S, y = second.consumption, color='#C257579E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = math_format(format = log10)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", axis.text.x = element_blank(),
        plot.title = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) +
  scale_colour_identity()

prim.cons_LAKE=ggpredict(S.prim.cons_LAKE, terms = "S")
S.prim.cons.sjp_LAKE <- ggplot(prim.cons_LAKE, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color='#3A67AE9E'), linewidth=1.5) +
  scale_linetype_manual(values = c("dashed")) +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, fill = '#3A67AE9E', alpha = 0.15) +
  geom_point(data = meta.Lakes, aes(x = S, y = prim.consumption, color='#3A67AE9E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), expand = expansion(mult = c(0.5, .1)), labels = math_format(format = log10)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", plot.margin = unit(c(4, 5.5, 0, 5.5), "pt"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_colour_identity()


## STREAMS ##
# S.predation_STREAM <- gls(log10(second.consumption) ~ log10(S), data=meta.Streams)
# plot(S.predation_STREAM, which=1)
# qqnorm(S.predation_STREAM)
# summary(S.predation_STREAM)

S.predation_STREAM.mm <- lme(log10(second.consumption) ~ log10(S), random=~1|study_ID, data=meta.Streams) #Best model
plot(S.predation_STREAM.mm, which=1)
qqnorm(S.predation_STREAM.mm)
summary(S.predation_STREAM.mm)

# anova(S.predation_STREAM.mm, S.predation_STREAM)


# S.prim.cons_STREAM <- gls(log10(prim.consumption) ~ log10(S), data=meta.Streams)
# plot(S.prim.cons_STREAM, which=1)
# qqnorm(S.prim.cons_STREAM)
# summary(S.prim.cons_STREAM)

S.prim.cons_STREAM.mm <- lme(log10(prim.consumption) ~ log10(S), random=~1|study_ID, data=meta.Streams) #Best model
plot(S.prim.cons_STREAM.mm, which=1)
qqnorm(S.predation_STREAM.mm)
summary(S.prim.cons_STREAM.mm)

#anova(S.prim.cons_STREAM, S.prim.cons_STREAM.mm)

predation_STREAM=ggpredict(S.predation_STREAM.mm, terms = "S", interval = "confidence")
S.predation.sjp_STREAM <- ggplot(predation_STREAM, aes(x, predicted)) +  ggtitle("Streams (n = 74)") +
  geom_line(aes(linetype="group", color='#C257579E'), linewidth=1.5) + 
  #scale_linetype_manual(values = c("dashed")) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Streams, aes(x = S, y = second.consumption, color='#C257579E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", axis.text.x = element_blank(),
        plot.title = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) +
  scale_colour_identity()

prim.cons_STREAM=ggpredict(S.prim.cons_STREAM.mm, terms = "S", interval = "confidence")
S.prim.cons.sjp_STREAM <- ggplot(prim.cons_STREAM, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color='#3A67AE9E'), linewidth=1.5) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = meta.Streams, aes(x = S, y = prim.consumption, color='#3A67AE9E'), size=2, alpha=0.5) +
  scale_y_continuous(breaks=breaks_log(n = 5), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(y = "log10", x = "log10")+
  theme(legend.position="none", plot.margin = unit(c(4, 5.5, 0, 5.5), "pt"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_colour_identity()




#### Produce BEF graphs for Fig. 4 ####

Richness_ECOSYSTEM <- grid.arrange(patchworkGrob(S.predation.sjp_MAR / S.prim.cons.sjp_MAR |
  S.predation.sjp_SOIL / S.prim.cons.sjp_SOIL |
  S.predation.sjp_STREAM / S.prim.cons.sjp_STREAM |
  S.predation.sjp_LAKE / S.prim.cons.sjp_LAKE), bottom=textGrob("Multitrophic taxon richness", gp=gpar(fontsize=10))) 

ggsave(filename = "Figure 4 scatterplots.svg", Richness_ECOSYSTEM, width = 200, height = 80, unit = "mm", dpi = 300)




# ########################################################################################################
# #### Structural Equation Models ########################################################################
# ########################################################################################################


## Create vectors for effect size plotting ##
flux <- c(rep(c(rep("predation", 3), rep("primary consumption",3)),2))
FW_prop <- rep(c("Taxon richness","Max TL","Trophic dissim."), 4)
effect.type <- factor(c(rep("direct", 6),rep("indirect", 6)), levels=c('indirect', 'direct'))



#### Cross-ecosystem SEM ####
##Explore simple linearity of simple NPP relationships
ggplot(all_data, aes(x = NPP.scale, y = log(S))) +
  geom_point(data = all_data, aes(alpha = 0.6, colour=ecosystem.type)) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(all_data, aes(x = NPP.scale, y = log(MaxTL))) +
  geom_point(data = all_data, aes(alpha = 0.6, colour=ecosystem.type)) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(all_data, aes(x = NPP.scale, y = logit(sim.prim.cons))) +
  geom_point(data = all_data, aes(alpha = 0.6, colour=ecosystem.type)) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(all_data, aes(x = NPP.scale, y = logit(sim.sec.cons))) +
  geom_point(data = all_data, aes(alpha = 0.6, colour=ecosystem.type)) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(all_data, aes(x = NPP.scale, y = log(prim.consumption))) +
  geom_point(data = all_data, aes(alpha = 0.6, colour=ecosystem.type)) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(all_data, aes(x = NPP.scale, y = log(second.consumption))) +
  geom_point(data = all_data, aes(alpha = 0.6, colour=ecosystem.type)) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) 

modS = lme(log(S) ~ NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data, method="ML")
  plot(modS, which=1)
  qqnorm(modS)
modSa = lme(log(S) ~ poly(NPP.scale,2,raw=TRUE), random=~1|ecosystem.type/study_ID, data=all_data, method="ML") #best model
  plot(modSa, which=1)
  qqnorm(modSa)
anova(modS, modSa)
summary(modSa)


mod1 = lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data, method="ML")
  plot(mod1, which=1)
  qqnorm(mod1)
mod1a = lme(log(MaxTL) ~ log(S) + poly(NPP.scale,2,raw=TRUE), random=~1|ecosystem.type/study_ID, data=all_data, method="ML") #best model
  plot(mod1, which=1)
  qqnorm(mod1)
anova(mod1, mod1a)
summary(mod1a)


mod2 = lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data, method="ML")
  plot(mod2)
  qqnorm(mod2)
mod2a = lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights = varIdent(form=~1|study_ID), data=all_data, method="ML") #best model
  plot(mod2a)
  qqnorm(mod2a)
mod2b = lme(logit(sim.prim.cons) ~ log(S) + poly(NPP.scale,2,raw=TRUE), weights=varIdent(form=~1|study_ID), 
            control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), random=~1|ecosystem.type/study_ID, data=all_data, method="ML")
  plot(mod2b)
  qqnorm(mod2b)
anova(mod2, mod2a, mod2b)
summary(mod2a)

mod3 = lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data, method="ML") 
  plot(mod3)
  qqnorm(mod3)
mod3a = lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, weights=varIdent(form=~1|study_ID), random=~1|ecosystem.type/study_ID,  # Best model
            control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data, method="ML") 
  plot(mod3a)
  qqnorm(mod3a)
anova(mod3, mod3a)
summary(mod3a)

mod4 = lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|study_ID, data=all_data, method="ML") # Best model
  plot(mod4, which=1)
  qqnorm(mod4)
mod4a = lme(log(prim.consumption) ~ log(S) + poly(NPP.scale,2,raw=TRUE), random=~1|ecosystem.type/study_ID, data=all_data, method="ML")
  plot(mod4a, which=1)
  qqnorm(mod4a)
anova(mod4, mod4a)
summary(mod4)

mod5 = lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale + logit(sim.sec.cons), random=~1|study_ID, data=all_data, method="ML")  # Best model
  plot(mod5, which=1)
  qqnorm(mod5)
mod5a = lme(log(second.consumption) ~ log(MaxTL) + log(S) + poly(NPP.scale,2,raw=TRUE) + logit(sim.sec.cons), random=~1|ecosystem.type/study_ID, data=all_data, method="ML")
  plot(mod5a, which=1)
  qqnorm(mod5a)
anova(mod5, mod5a)
summary(mod5)



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

summary(SEM.all)


### Min adequate model
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

summary(SEM.all2)

results.allSEM <- summary(SEM.all2)$coefficients[c(1:15,19),c(1:5, 8, 7)]
names(results.allSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
fun <- function(x) {
  formatC(x, format = "f", digits = 3)
}
allSEM_table <- nice_table(results.allSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")


#### calculate std. effect size for quadratic variables sensu Henseler et al. 2012 (https://doi.org/10.1057/ejis.2011.36)
#Refit SEM.all2 without NPP.scale quadratic variable
SEM.all_no.NPP_R2 <- rsquared(psem(
  lme(log(S) ~ 1, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(MaxTL) ~ log(S), random=~1|ecosystem.type/study_ID, data=all_data),
  lme(logit(sim.prim.cons) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), data=all_data),
  lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S), random=~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID), 
      control=nlmeControl(opt = "nlminb",maxIter = 200,msMaxIter=200), data=all_data),
  lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, data=all_data),
  lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale + logit(sim.sec.cons), random=~1|ecosystem.type/study_ID, data=all_data),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
))
SEM.all2_R2 <- rsquared(SEM.all2)

#Calculate Cohen's f² = (R²full - R²reduced) / (1 - R²full)
SEM.all2_Q.Std.Est <- data.frame(SEM.all2_R2[,1],
  (SEM.all2_R2[,6] - SEM.all_no.NPP_R2[,6]) / (1 - SEM.all2_R2[,6]),
  row.names=NULL)
colnames(SEM.all2_Q.Std.Est) <- c("response", "F^2")

#### Create results dataframe for summary boxplots 
std.effect <- c(
  SR.direct.pred <- results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'log(S)'],
  MTL.direct.pred <- results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.pred <- results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'],
  SR.direct.prim <- results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(prim.consumption)' & results.allSEM$Predictor == 'log(S)'],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  
  SR.indirect.pred <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(MaxTL)' & results.allSEM$Predictor == 'log(S)']*
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'log(MaxTL)'] +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(MaxTL)' & results.allSEM$Predictor == 'log(S)']*
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(MaxTL)']*
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'] +
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(S)']*
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'],
  MTL.indirect.pred <- 
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'logit(sim.sec.cons)' & results.allSEM$Predictor == 'log(MaxTL)']*
    results.allSEM$`Std. Estimate`[results.allSEM$Response == 'log(second.consumption)' & results.allSEM$Predictor == 'logit(sim.sec.cons)'],
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- 0,
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0
)
eff.table_all <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot for Fig. 3 ##
global.effects <- ggplot(eff.table_all, 
                         aes(x = flux, y = std.effect, fill = flux, pattern = effect.type)) + 
  ylab('Effect size') + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5) +
  geom_bar_pattern(stat = "identity", position = "stack", pattern_color = 'White', pattern_fill = "white", width = .98,
                   pattern_alpha = 0.7, pattern_angle = 55, pattern_density = 0.75, pattern_spacing = 0.1, pattern_size=0) +
  scale_fill_manual(values = c("#C257579E", "#3A67AE9E")) + 
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  coord_flip(ylim=c(-0.5, 0.5)) + scale_y_continuous(labels = c(-0.5, '',0.0, '',0.5)) +
  facet_grid(fct_relevel(FW_prop,'Taxon richness', 'Max TL', 'Trophic dissim.')~., scales = "free_y", labeller = label_wrap_gen(width=10)) +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), 
        strip.clip = "off", strip.background = element_blank(), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))


ggsave("Global effects.svg", global.effects, width = 6, height = 9, units = "cm", bg='transparent')



#### Construct ecosystem-specific SEMs as in Fig. 4 ####


#### MARINE ####
##Explore simple linearity of simple NPP relationships
ggplot(meta.Marine, aes(x = NPP.scale, y = log(S))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Marine, aes(x = NPP.scale, y = log(MaxTL))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Marine, aes(x = NPP.scale, y = logit(sim.prim.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Marine, aes(x = NPP.scale, y = logit(sim.sec.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Marine, aes(x = NPP.scale, y = log(prim.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Marine, aes(x = NPP.scale, y = log(second.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) 


mod0 = gls(log(S) ~ NPP.scale, data=meta.Marine, method='ML')
  plot(mod0, which=1)
  qqnorm(mod0)
  summary(mod0)
mod0a = gls(log(S) ~ poly(NPP.scale, degree=2, raw=TRUE), data=meta.Marine, method='ML') # Best model
  plot(mod0a, which=1)
  qqnorm(mod0a)
  summary(mod0a)
anova(mod0, mod0a)


mod1 = gls(log(MaxTL) ~ log(S) + NPP.scale, data=meta.Marine, method='ML')
  plot(mod1, which=1)
  qqnorm(mod1)
  summary(mod1)
mod1a = gls(log(MaxTL) ~ log(S) + poly(NPP.scale,2,raw=TRUE), data=meta.Marine, method='ML') # Best model
  plot(mod1a, which=1)
  qqnorm(mod1a)
  summary(mod1a)
anova(mod1, mod1a)


mod2 = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Marine, method='ML') # Best model
  plot(mod2, which=1)
  qqnorm(mod2)
  summary(mod2)
mod2a = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + poly(NPP.scale,2,raw=TRUE), data=meta.Marine, method='ML')
  plot(mod2a, which=1)
  qqnorm(mod2a)
  summary(mod2a)
anova(mod2, mod2a)


mod3 = gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Marine, method='ML') # Best model
  plot(mod3, which=1)
  qqnorm(mod3)
  summary(mod3)
mod3a = gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + poly(NPP.scale,2,raw=TRUE), data=meta.Marine, method='ML')
  plot(mod3a, which=1)
  qqnorm(mod3a)
  summary(mod3a)
anova(mod3, mod3a)


mod4 = gls(log(prim.consumption) ~ log(S) + NPP.scale, data=meta.Marine, method='ML')
  plot(mod4, which=1)
  qqnorm(mod4)
  summary(mod4)
mod4a = gls(log(prim.consumption) ~ log(S) + NPP.scale, weights=varPower(form=~NPP.scale), data=meta.Marine, method='ML')  # Best model
  plot(mod4a, which=1)
  qqnorm(mod4a)
  summary(mod4a)
mod4b = gls(log(prim.consumption) ~ log(S) + poly(NPP.scale,2,raw=TRUE), weights=varPower(form=~NPP.scale), data=meta.Marine, method='ML')
  plot(mod4b, which=1)
  qqnorm(mod4b)
  summary(mod4b)
anova(mod4, mod4a, mod4b)


mod5 = gls(log(second.consumption) ~ log(MaxTL) + NPP.scale, data=meta.Marine, method='ML')
  plot(mod5, which=1)
  qqnorm(mod5)
  summary(mod5)
mod5a = gls(log(second.consumption) ~ log(MaxTL) + NPP.scale, weights=varPower(form=~S), data=meta.Marine, method='ML')
  plot(mod5a, which=1)
  qqnorm(mod5a)
  summary(mod5a)
mod5b = gls(log(second.consumption) ~ log(MaxTL) + poly(NPP.scale,2,raw=TRUE), weights=varPower(form=~S), data=meta.Marine, method='ML') # Best model
  plot(mod5b, which=1)
  qqnorm(mod5b)
  summary(mod5b)
anova(mod5, mod5a, mod5b)


### Maximal model
SEM.Marine2 <- psem(
  
  gls(log(S) ~ NPP.scale + NPP.scale2, data=meta.Marine),
  gls(log(MaxTL) ~ log(S) + NPP.scale + NPP.scale2, data=meta.Marine),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale + NPP.scale2, data=meta.Marine),
  gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale + NPP.scale2, data=meta.Marine),
  gls(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), weights=varPower(form=~NPP.scale), data=meta.Marine),
  gls(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), weights=varPower(form=~S), data=meta.Marine),

  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Marine2)

### Min adequate model
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

summary(SEM.Marine3)
results.marineSEM <- summary(SEM.Marine3)$coefficients[c(1:15,19),c(1:5, 8, 7)]
names(results.marineSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
fun <- function(x) {
  formatC(x, format = "f", digits = 3)
}
marineSEM_table <- nice_table(results.marineSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")


####calculate std. effect size for quadratic variables sensu Henseler et al. 2012 (https://doi.org/10.1057/ejis.2011.36)
#Refit SEM without NPP.scale quadratic variable
SEM.Marine3_no.NPP_R2 <- rsquared(psem(
  
  gls(log(S) ~ 1, data=meta.Marine),
  gls(log(MaxTL) ~ log(S), data=meta.Marine),
  gls(logit(sim.prim.cons) ~ log(S), data=meta.Marine),
  gls(logit(sim.sec.cons) ~ log(S) + log(MaxTL) + NPP.scale, data=meta.Marine),
  gls(log(prim.consumption) ~ log(S), weights=varPower(form=~NPP.scale), data=meta.Marine),
  gls(log(second.consumption) ~ log(MaxTL), weights=varPower(form=~S), data=meta.Marine),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
))
SEM.Marine3_R2 <- rsquared(SEM.Marine3)

#Calculate Cohen's f² = (R²full - R²reduced) / (1 - R²full)
SEM.Marine3_Q.Std.Est <- data.frame(SEM.Marine3_R2[,1],
                                 (SEM.Marine3_R2[,5] - SEM.Marine3_no.NPP_R2[,5]) / (1 - SEM.Marine3_R2[,5]),
                                 row.names=NULL)
colnames(SEM.Marine3_Q.Std.Est) <- c("response", "F^2")

#### Create results dataframe for summary boxplots 
std.effect <- c(
  SR.direct.pred <- 0,
  MTL.direct.pred <- results.marineSEM$`Std. Estimate`[results.marineSEM$Response == 'log(second.consumption)' & results.marineSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.pred <- 0,
  SR.direct.prim <- results.marineSEM$`Std. Estimate`[results.marineSEM$Response == 'log(prim.consumption)' & results.marineSEM$Predictor == 'log(S)'],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  
  SR.indirect.pred <- 
    results.marineSEM$`Std. Estimate`[results.marineSEM$Response == 'log(MaxTL)' & results.marineSEM$Predictor == 'log(S)']*
    results.marineSEM$`Std. Estimate`[results.marineSEM$Response == 'log(second.consumption)' & results.marineSEM$Predictor == 'log(MaxTL)'],
  MTL.indirect.pred <- 0,
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- 0,
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0
)
eff.table_marine <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot
marine.effects <- ggplot(eff.table_marine, 
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

ggsave("Marine effects.svg", marine.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### SOIL ####
##Explore simple linearity of simple NPP relationships
ggplot(meta.Soils, aes(x = NPP.scale, y = log(S))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Soils, aes(x = NPP.scale, y = log(MaxTL))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Soils, aes(x = NPP.scale, y = logit(sim.prim.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Soils, aes(x = NPP.scale, y = logit(sim.sec.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Soils, aes(x = NPP.scale, y = log(prim.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Soils, aes(x = NPP.scale, y = log(second.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2))

mod0 = gls(log(S) ~ NPP.scale, data=meta.Soils, method='ML')
  plot(mod0, which=1)
  qqnorm(mod0)
  summary(mod0)
mod0a = gls(log(S) ~ poly(NPP.scale,2,raw=TRUE), data=meta.Soils, method='ML') # Best model
  plot(mod0a, which=1)
  qqnorm(mod0a)
  summary(mod0a)
anova(mod0, mod0a)


mod1 = gls(log(MaxTL) ~ log(S), data=meta.Soils, method='ML') # Best model
  plot(mod1, which=1)
  qqnorm(mod1)
  summary(mod1)
mod1a = gls(log(MaxTL) ~ log(S), weights=varIdent(form=~1|study_ID), data=meta.Soils, method='ML')
  plot(mod1a, which=1)
  qqnorm(mod1a)
  summary(mod1a)
anova(mod1, mod1a)


mod2 = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Soils, method='ML')
  plot(mod2, which=1)
  qqnorm(mod2)
  summary(mod2)
mod2a = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), weights=varIdent(form=~1|study_ID), data=meta.Soils, method='ML') # Best model
  plot(mod2a, which=1)
  qqnorm(mod2a)
  summary(mod2a)
anova(mod2, mod2a)


mod3 = gls(logit(sim.sec.cons) ~ NPP.scale, data=meta.Soils, method='ML') # Best model
  plot(mod3, which=1)
  qqnorm(mod3)
  summary(mod3)
mod3a = gls(logit(sim.sec.cons) ~ poly(NPP.scale,2,raw=TRUE), data=meta.Soils, method='ML')
  plot(mod3a, which=1)
  qqnorm(mod3a)
  summary(mod3a)
anova(mod3, mod3a)


mod4 = gls(log(prim.consumption) ~ log(S) + NPP.scale, data=meta.Soils, method='ML')
  plot(mod4, which=1)
  qqnorm(mod4)
  summary(mod4)
mod4a = gls(log(prim.consumption) ~ log(S) + NPP.scale, weights=varIdent(form=~1|study_ID), data=meta.Soils, method='ML') # Best model
  plot(mod4a, which=1)
  qqnorm(mod4a)
  summary(mod4a)
mod4b = gls(log(prim.consumption) ~ log(S) + poly(NPP.scale,2,raw=TRUE), weights=varIdent(form=~1|study_ID), data=meta.Soils, method='ML')
  plot(mod4b, which=1)
  qqnorm(mod4b)
  summary(mod4b)
anova(mod4, mod4a, mod4b)


mod5 = gls(log(second.consumption) ~ log(S) + log(MaxTL), data=meta.Soils, method='ML')
  plot(mod5, which=1)
  qqnorm(mod5)
  summary(mod5)
mod5a = gls(log(second.consumption) ~ log(S) + log(MaxTL), weights=varIdent(form=~1|study_ID), data=meta.Soils, method='ML')  # Best model
  plot(mod5a, which=1)
  qqnorm(mod5a)
  summary(mod5a)
anova(mod5, mod5a)


### Maximal model
SEM.Soils2 <- psem(

  gls(log(S) ~ NPP.scale + NPP.scale2, data=meta.Soils),
  gls(log(MaxTL) ~ log(S) + NPP.scale, data=meta.Soils),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Soils),
  gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Soils),
  gls(log(prim.consumption) ~ logit(sim.prim.cons) + log(MaxTL), weights=varIdent(form=~1|study_ID), data=meta.Soils),
  gls(log(second.consumption) ~ logit(sim.sec.cons) + log(MaxTL), weights=varIdent(form=~1|study_ID), data=meta.Soils),

  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Soils2)


### Min adequate model
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

summary(SEM.Soils3)
results.soilsSEM <- summary(SEM.Soils3)$coefficients[c(1:11, 15),c(1:5, 8, 7)]
names(results.soilsSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
soilsSEM_table <- nice_table(results.soilsSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")


####calculate std. effect size for quadratic variables sensu Henseler et al. 2012 (https://doi.org/10.1057/ejis.2011.36)
#Refit SEM without NPP.scale quadratic variable
SEM.Soils3_no.NPP_R2 <- rsquared(psem(
  
  gls(log(S) ~ 1, data=meta.Soils),
  gls(log(MaxTL) ~ log(S), data=meta.Soils),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Soils),
  gls(logit(sim.sec.cons) ~ NPP.scale + log(S), data=meta.Soils),
  gls(log(prim.consumption) ~ log(S) + NPP.scale, weights=varIdent(form=~1|study_ID), data=meta.Soils),
  gls(log(second.consumption) ~ log(S) + log(MaxTL), weights=varIdent(form=~1|study_ID), data=meta.Soils),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
))
SEM.Soils3_R2 <- rsquared(SEM.Soils3)

#Calculate Cohen's f² = (R²full - R²reduced) / (1 - R²full)
SEM.Soils3_Q.Std.Est <- data.frame(SEM.Soils3_R2[,1],
                                    (SEM.Soils3_R2[,5] - SEM.Soils3_no.NPP_R2[,5]) / (1 - SEM.Soils3_R2[,5]),
                                    row.names=NULL)
colnames(SEM.Soils3_Q.Std.Est) <- c("response", "F^2")


#### Create results dataframe for summary boxplots 
std.effect <- c(
  SR.direct.pred <- results.soilsSEM$`Std. Estimate`[results.soilsSEM$Response == 'log(second.consumption)' & results.soilsSEM$Predictor == 'log(S)'],
  MTL.direct.pred <- results.soilsSEM$`Std. Estimate`[results.soilsSEM$Response == 'log(second.consumption)' & results.soilsSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.pred <- 0,
  SR.direct.prim <- results.soilsSEM$`Std. Estimate`[results.soilsSEM$Response == 'log(prim.consumption)' & results.soilsSEM$Predictor == 'log(S)'],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  
  SR.indirect.pred <- results.soilsSEM$`Std. Estimate`[results.soilsSEM$Response == 'log(MaxTL)' & results.soilsSEM$Predictor == 'log(S)'] *
    results.soilsSEM$`Std. Estimate`[results.soilsSEM$Response == 'log(second.consumption)' & results.soilsSEM$Predictor == 'log(MaxTL)'],
  MTL.indirect.pred <- 0,
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- 0,
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0
)
eff.table_soils <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot
soils.effects <- ggplot(eff.table_soils, 
                        aes(x = flux, y = std.effect, fill = flux, pattern = effect.type)) + 
  ylab('Effect size') + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5) +
  geom_bar_pattern(stat = "identity", position = "stack", pattern_color = 'White', pattern_fill = "white", width = .98,
                   pattern_alpha = 0.7, pattern_angle = 55, pattern_density = 0.75, pattern_spacing = 0.1, pattern_size=0) +
  scale_fill_manual(values = c("#C257579E", "#3A67AE9E")) + 
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  coord_flip(ylim=c(-0.5, 0.5)) + scale_y_continuous(labels = c(-0.5, '',0.0, '',0.5)) +
  facet_grid(fct_relevel(FW_prop,'Taxon richness', 'Max TL', 'Trophic dissim.')~., scales = "free_y", labeller = label_wrap_gen(width=10)) +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), 
        strip.clip = "off", strip.background = element_blank(), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

ggsave("Soils effects.svg", soils.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### STREAMS ####
ggplot(meta.Streams, aes(x = NPP.scale, y = log(S))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Streams, aes(x = NPP.scale, y = log(MaxTL))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Streams, aes(x = NPP.scale, y = logit(sim.prim.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Streams, aes(x = NPP.scale, y = logit(sim.sec.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Streams, aes(x = NPP.scale, y = log(prim.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Streams, aes(x = NPP.scale, y = log(second.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2))

mod0 = lme(log(S) ~ NPP.scale, random=~1|study_ID, data=meta.Streams, method="ML")
  plot(mod0)
  qqnorm(mod0)
  summary(mod0)
mod0a = lme(log(S) ~ NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams, method='ML') # Best model
  plot(mod0a)
  qqnorm(mod0a)
  summary(mod0a)
mod0b = lme(log(S) ~ poly(NPP.scale, 2), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams, method='ML')
  plot(mod0b)
  qqnorm(mod0b)
  summary(mod0b)
anova(mod0, mod0a, mod0b)


mod1 = lme(log(MaxTL) ~ log(S), random=~1|study_ID, data=meta.Streams, method='ML')
  plot(mod1)
  qqnorm(mod1)
  summary(mod1)
mod1a = lme(log(MaxTL) ~ log(S), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams, method='ML') # Best model
  plot(mod1a)
  qqnorm(mod1a)
  summary(mod1a)
anova(mod1, mod1a) 
  

mod2 = lme(logit(sim.prim.cons) ~ log(MaxTL), random=~1|study_ID, data=meta.Streams, method='ML')
  plot(mod2)
  qqnorm(mod2)
  summary(mod2)
mod2a = lme(logit(sim.prim.cons) ~ log(MaxTL), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams, method='ML') # Best model
  plot(mod2a)
  qqnorm(mod2a)
  summary(mod2a)
anova(mod2, mod2a)   
  

mod3 = lme(logit(sim.sec.cons) ~ log(S), random=~1|study_ID, data=meta.Streams, method='ML')
  plot(mod3)
  qqnorm(mod3)
  summary(mod3)
mod3a = lme(logit(sim.sec.cons) ~ log(S), weights=varExp(form = ~S), random=~1|study_ID, data=meta.Streams, method='ML') # Best model
  plot(mod3a)
  qqnorm(mod3a)
  summary(mod3a)
anova(mod3, mod3a)  
  

mod4 = lme(log(prim.consumption) ~ log(S), random=~1|study_ID, data=meta.Streams)
  plot(mod4, which=1)
  qqnorm(mod4)
  summary(mod4)
  
  
mod5 = lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), random=~1|study_ID, data=meta.Streams)
  plot(mod5, which=1)
  qqnorm(mod5)
  summary(mod5)
  

### Maximal model
SEM.Streams2 <- psem(
  
  lme(log(S) ~ NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|study_ID, weights=varExp(), data=meta.Streams),
  lme(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), random=~1|study_ID, data=meta.Streams),
  lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), random=~1|study_ID, data=meta.Streams),

  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Streams2)


### Min adequate model
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

summary(SEM.Streams3)

results.streamsSEM <- summary(SEM.Streams3)$coefficients[c(1:7, 11),c(1:5, 8, 7)]
names(results.streamsSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
streamsSEM_table <- nice_table(results.streamsSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")


#### Create results dataframe for summary boxplots 
std.effect <- c(
  SR.direct.pred <- 0,
  MTL.direct.pred <- results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(second.consumption)' & results.streamsSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.pred <- results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(second.consumption)' & results.streamsSEM$Predictor == 'logit(sim.sec.cons)'],
  SR.direct.prim <- results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(prim.consumption)' & results.streamsSEM$Predictor == 'log(S)'],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,

  SR.indirect.pred <- results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'logit(sim.sec.cons)' & results.streamsSEM$Predictor == 'log(S)'] *
    results.streamsSEM$`Std. Estimate`[results.streamsSEM$Response == 'log(second.consumption)' & results.streamsSEM$Predictor == 'logit(sim.sec.cons)'],
  MTL.indirect.pred <- 0,
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- 0,
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0
)
eff.table_streams <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot
streams.effects <- ggplot(eff.table_streams, 
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

ggsave("Streams effects.svg", streams.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### Lakes ####
ggplot(meta.Lakes, aes(x = NPP.scale, y = log(S))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Lakes, aes(x = NPP.scale, y = log(MaxTL))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Lakes, aes(x = NPP.scale, y = logit(sim.prim.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Lakes, aes(x = NPP.scale, y = logit(sim.sec.cons))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Lakes, aes(x = NPP.scale, y = log(prim.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
ggplot(meta.Lakes, aes(x = NPP.scale, y = log(second.consumption))) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", formula = y ~ x + I(x^2))

mod0 = gls(log(S) ~ NPP.scale, data=meta.Lakes, method='ML') # Best model
  plot(mod0)
  qqnorm(mod0)
  summary(mod0)
mod0a = gls(log(S) ~ NPP.scale, weights=varPower(form=~NPP.scale), data=meta.Lakes, method='ML')
  plot(mod0a)
  qqnorm(mod0a)
  summary(mod0a)
mod0b = gls(log(S) ~ poly(NPP.scale, 2), weights=varPower(form=~NPP.scale), data=meta.Lakes, method='ML')
  plot(mod0b)
  qqnorm(mod0b)
  summary(mod0b)
anova(mod0, mod0a, mod0b)

mod1 = gls(log(MaxTL) ~ log(S) + NPP.scale, data=meta.Lakes, method='ML')
  plot(mod1)
  qqnorm(mod1)
  summary(mod1)
mod1a = gls(log(MaxTL) ~ log(S) + NPP.scale, weights=varPower(form=~S), data=meta.Lakes, method='ML')
  plot(mod1a)
  qqnorm(mod1a)
  summary(mod1a)
mod1b = gls(log(MaxTL) ~ log(S) + poly(NPP.scale, 2), weights=varPower(form=~S), data=meta.Lakes, method='ML') # Best model
  plot(mod1b)
  qqnorm(mod1b)
  summary(mod1b)
anova(mod1, mod1a, mod1b)

mod2 = gls(logit(sim.prim.cons) ~ log(S) + NPP.scale, data=meta.Lakes, na.action = na.omit, method='ML')
  plot(mod2)
  qqnorm(mod2)
  summary(mod2)
mod2a = gls(logit(sim.prim.cons) ~ log(S) + NPP.scale, weights=varExp(form=~S), data=meta.Lakes, na.action = na.omit, method='ML') # Best model
  plot(mod2a)
  qqnorm(mod2a)
  summary(mod2a)
mod2b = gls(logit(sim.prim.cons) ~ log(S) + poly(NPP.scale, 2), weights=varExp(form=~S), data=meta.Lakes, na.action = na.omit, method='ML')
  plot(mod2b)
  qqnorm(mod2b)
  summary(mod2b)
anova(mod2, mod2a, mod2b)

mod3 = gls(logit(sim.sec.cons) ~ log(MaxTL), data=meta.Lakes, method='ML')
  plot(mod3, which=1)
  qqnorm(mod3)
  summary(mod3)
mod3a = gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varExp(form=~MaxTL), data=meta.Lakes, method='ML') # Best model
  plot(mod3a, which=1)
  qqnorm(mod3a)
  summary(mod3a)
mod3b = gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varComb(varExp(form=~S), varExp(form=~MaxTL)), data=meta.Lakes, method='ML')
  plot(mod3b, which=1)
  qqnorm(mod3b)
  summary(mod3b)
anova(mod3, mod3a, mod3b)

mod4 = gls(log(prim.consumption) ~ log(MaxTL) + log(S), data=meta.Lakes, method='ML')
  plot(mod4, which=1)
  qqnorm(mod4)
  summary(mod4)
mod4a = gls(log(prim.consumption) ~ log(MaxTL) + log(S), weights=varPower(), data=meta.Lakes, method='ML') # Best model
  plot(mod4a, which=1)
  qqnorm(mod4a)
  summary(mod4a)
anova(mod4, mod4a)

mod5 = gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S),  data=meta.Lakes, method='ML') # Best model
  plot(mod5, which=1)
  qqnorm(mod5)
  summary(mod5)
mod5a = gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S), weights=varExp(form=~sim.sec.cons),  data=meta.Lakes, method='ML')
  plot(mod5a, which=1)
  qqnorm(mod5a)
  summary(mod5a)
anova(mod5, mod5a)

### Maximal model
SEM.Lakes2 <- psem(

  gls(log(S) ~ NPP.scale, data=meta.Lakes),
  gls(log(MaxTL) ~ log(S) + NPP.scale, weights=varPower(form=~S), data=meta.Lakes),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, weights=varExp(form=~S), data=meta.Lakes),
  gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, weights=varExp(form=~log(MaxTL)), data=meta.Lakes),
  gls(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), weights=varPower(), data=meta.Lakes),
  gls(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), data=meta.Lakes),

  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Lakes2)


### Min adequate model
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


summary(SEM.Lakes3)

results.lakesSEM <- summary(SEM.Lakes3)$coefficients[c(1:10, 14),c(1:5, 8, 7)]
names(results.lakesSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
lakesSEM_table <- nice_table(results.lakesSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")


####calculate std. effect size for quadratic variables sensu Henseler et al. 2012 (https://doi.org/10.1057/ejis.2011.36)
#Refit SEM without NPP.scale quadratic variable
SEM.Lakes3_no.NPP_R2 <- rsquared(psem(
  
  gls(log(S) ~ NPP.scale, data=meta.Lakes),
  gls(log(MaxTL) ~ log(S), weights=varPower(form=~S), data=meta.Lakes),
  gls(logit(sim.prim.cons) ~ log(S), weights=varExp(form=~S), data=meta.Lakes),
  gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varExp(form=~log(MaxTL)), data=meta.Lakes),
  gls(log(prim.consumption) ~ log(MaxTL) + log(S), weights=varPower(), data=meta.Lakes),
  gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S), data=meta.Lakes),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
))
SEM.Lakes3_R2 <- rsquared(SEM.Lakes3)

#Calculate Cohen's f² = (R²full - R²reduced) / (1 - R²full)
SEM.Lakes3_Q.Std.Est <- data.frame(SEM.Lakes3_R2[,1],
                                   (SEM.Lakes3_R2[,5] - SEM.Lakes3_no.NPP_R2[,5]) / (1 - SEM.Lakes3_R2[,5]),
                                   row.names=NULL)
colnames(SEM.Lakes3_Q.Std.Est) <- c("response", "F^2")


## Create results daframe for summary boxplots
std.effect <- c(
  SR.direct.pred <- 0,
  MTL.direct.pred <- 0,
  SIM.direct.pred <- results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(second.consumption)' & results.lakesSEM$Predictor == 'logit(sim.sec.cons)'],
  SR.direct.prim <- results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(prim.consumption)' & results.lakesSEM$Predictor == 'log(S)'],
  MTL.direct.prim <- results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(prim.consumption)' & results.lakesSEM$Predictor == 'log(MaxTL)'],
  SIM.direct.prim <- 0,

  SR.indirect.pred <- 
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(MaxTL)' & results.lakesSEM$Predictor == 'log(S)'] *
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'logit(sim.sec.cons)' & results.lakesSEM$Predictor == 'log(MaxTL)'] *
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(second.consumption)' & results.lakesSEM$Predictor == 'logit(sim.sec.cons)'],
  MTL.indirect.pred <- 
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'logit(sim.sec.cons)' & results.lakesSEM$Predictor == 'log(MaxTL)'] *
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(second.consumption)' & results.lakesSEM$Predictor == 'logit(sim.sec.cons)'],
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- 
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(MaxTL)' & results.lakesSEM$Predictor == 'log(S)'] *
    results.lakesSEM$`Std. Estimate`[results.lakesSEM$Response == 'log(prim.consumption)' & results.lakesSEM$Predictor == 'log(MaxTL)'],
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0
)
eff.table_lakes <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot
lakes.effects <- ggplot(eff.table_lakes, 
       aes(x = flux, y = std.effect, fill = flux, pattern = effect.type)) + 
  ylab('Effect size') + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5) +
  geom_bar_pattern(stat = "identity", position = "stack", pattern_color = 'White', pattern_fill = "white", width = .98,
                   pattern_alpha = 0.7, pattern_angle = 55, pattern_density = 0.75, pattern_spacing = 0.1, pattern_size=0) +
  scale_fill_manual(values = c("#C257579E", "#3A67AE9E")) + 
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  coord_flip(ylim=c(-0.5, 0.5)) + scale_y_continuous(labels = c(-0.5, '',0.0, '',0.5)) +
  facet_grid(fct_relevel(FW_prop,'Taxon richness', 'Max TL', 'Trophic dissim.')~., scales = "free_y", labeller = label_wrap_gen(width=10)) +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), 
        strip.clip = "off", strip.background = element_blank(), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

ggsave("Lakes effects.svg", lakes.effects, width = 6, height = 9, units = "cm", bg='transparent')




#################################################
#### Map of food web sampling sites - Fig. 2 ####
#################################################

world <- map_data("world")

all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                  levels = c("Marine", "Soils", "Streams", "Lakes"))

ggplot() +
  geom_map(data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.1) + 
  xlim(-130, 75) + ylim(-52, 80) +
  geom_point(data = all_data, size=2, stroke = 1,
    aes(lon, lat, color = ecosystem.type),
    alpha = 0.4, shape=21) +
  geom_jitter(width = 25) +
# scale_size_continuous(range=c(4,5)) +
  theme_void()+
  labs(colour = ("Ecosystem type")) +
  theme(panel.background = element_rect(fill = "white", colour = NA),legend.key=element_blank(),
    legend.position = c(.2, .4), legend.justification = c("right", "top"), 
    legend.box.just = "left", legend.margin = margin(6, 6, 6, 6)) + 
  guides(colour = guide_legend(override.aes = list(size=2.5, stroke = 1.5, shape=21, alpha=0.9)), shape=19)
ggsave("map_food webs.png", width = 17, height = 11, units = "cm")

