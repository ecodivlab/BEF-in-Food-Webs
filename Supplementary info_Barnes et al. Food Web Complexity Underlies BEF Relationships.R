#=========================================================================================================================
#     
# Code for Barnes et al. 'Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning'
# Main text analyses and figure
#
#=========================================================================================================================

# This code generates bivariate and SEM models testing for the effect of taxa richness, maximum trophic level, and
# trophic dissimilarity on total primary consumption and predation rates (derived from estimated energy fluxes) in
# 319 food webs, spanning marine, lake, stream, and soil food webs. 
# 
# All models presented in Barnes et al. are run according to the described methods in the manuscript.
# The data accompanying this script are available on figshare DOI: 10.6084/m9.figshare.28646129

## Load packages ##
library(tidyverse); library(sjPlot); library(ggeffects); library(gridExtra); library(piecewiseSEM);
library(patchwork); library(nlme); library(grid); library(car); library(rempsyc); library(ggpattern);
library(ggh4x); library(scales); library(ggtext); library(ggrain)

## Clear environment and read ecosystem-specific food web data sets ##
rm(list=ls())

setwd("C:\\Users\\barnesa\\OneDrive - The University of Waikato\\FuSED\\BEF-in-Food-Webs")

NPP.proxy <- read.csv("NDVI and Chlorophyll-a/data/proxy-npp.csv") ## import NDVI & Chl-a data
meta.Marine <- read.csv('meta.Marine.csv')  
  meta.Marine <- meta.Marine %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name")
  colnames(meta.Marine)[19] <- "NPP.proxy"
  meta.Marine <- meta.Marine %>% mutate(clean_var = ifelse(is.na(NPP.proxy), 0.5, NPP.proxy),  #Temporary - replace BS NAs with 0.5
                                        NPP.scale = logit((clean_var - min(clean_var, na.rm = TRUE)) /
    (max(clean_var, na.rm = TRUE) - min(clean_var, na.rm = TRUE))))
  # Calculate inverse of the leading eigenvalue to represent stability (so positive values = greater resilience) 
  meta.Marine$stability <- -scale(log(meta.Marine$stability))[,1] # scale and log-transform for model fitting
  
meta.Soils <- read.csv('meta.Soils.csv')
  meta.Soils <- meta.Soils %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name")
  colnames(meta.Soils)[c(7,19)] <- c("temperature_C", "NPP.proxy")
  meta.Soils <- meta.Soils %>% mutate(NPP.scale = logit(NPP.proxy))
  # Calculate inverse of the leading eigenvalue to represent stability (so positive values = greater resilience)
  meta.Soils$stability <- -scale(log(meta.Soils$stability))[,1] # scale and log-transform for model fitting
                                      
meta.Streams <- read.csv('meta.Streams.csv')
  meta.Streams <- meta.Streams %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name")
  colnames(meta.Streams)[c(5,20)] <- c("temperature_C", "NPP.proxy")
  meta.Streams <- meta.Streams %>% mutate(clean_var = ifelse(NPP.proxy < 0, 0, NPP.proxy), # replace negative value at Cananeia SP6 with 0
                  NPP.scale = logit(clean_var))
  # Calculate inverse of the leading eigenvalue to represent stability (so positive values = greater resilience)
  meta.Streams$stability <- -scale(log(meta.Streams$stability))[,1] # scale and log-transform for model fitting
  
meta.Lakes <- read.csv('meta.Lakes.csv')
  meta.Lakes <- meta.Lakes %>% left_join(NPP.proxy %>% select(FW_name, metric, avg), by = "FW_name")
  colnames(meta.Lakes)[19] <- "NPP.proxy"
  meta.Lakes <- meta.Lakes %>% mutate(NPP.scale = logit(NPP.proxy))
  # Calculate inverse of the leading eigenvalue to represent stability (so positive values = greater resilience)
  meta.Lakes$stability <- -scale(log(meta.Lakes$stability))[,1] # scale and log-transform for model fitting

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

## Total flux
S.flux_Global <- lme(log10(flux) ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data, na.action=na.omit)
plot(S.flux_Global, which=1)
qqnorm(S.flux_Global)
summary(S.flux_Global)

S.flux_Globala=update(S.flux_Global, random = ~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID)) #Best model
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

all_data$ecosystem.type <- factor(all_data$ecosystem.type, 
                                  levels = c("Marine", "Soils", "Streams", "Lakes"))


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


## Stability
S.stability_Global <- lme(stability ~ log10(S), random = ~1|ecosystem.type/study_ID, data=all_data, na.action=na.omit)
plot(S.stability_Global, which=1)
qqnorm(S.stability_Global)
summary(S.stability_Global)

S.stability_Globala=update(S.stability_Global, random = ~1|ecosystem.type/study_ID, weights=varIdent(form=~1|study_ID)) #Best model
anova(S.stability_Global,S.stability_Globala)
plot(S.stability_Globala)
qqnorm(S.stability_Globala)
summary(S.stability_Globala)

S.stability_Globalb=update(S.stability_Global, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)))
anova(S.stability_Global,S.stability_Globala,S.stability_Globalb)
plot(S.stability_Globalb)
qqnorm(S.stability_Globalb)
summary(S.stability_Globalb)

S.stability_Global = S.stability_Globalb
## Graph BEF for total flux 
Stability_Global=ggpredict(S.stability_Global, terms = "S")
S.total.sjp_global <- ggplot(Stability_Global, aes(x, predicted)) + 
  ylab('Stability') + xlab("Multitrophic taxon richness") + 
  geom_line(linewidth=1.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), linetype=2, alpha = 0.15) +
  geom_point(data = all_data, aes(x = S, y = stability, colour=ecosystem.type), size=2, alpha=0.4) +
  labs(colour = ("Ecosystem type")) +
  #scale_y_continuous(breaks=breaks_log(n = 6), labels = label_log(digits = 2)) +
  scale_x_continuous(breaks=breaks_log(n=7), expand = expansion(mult = c(0.05, .1)))+
  coord_trans(x = "log10")+
  theme(legend.position="top", plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"), text=element_text(size=14)) 


#### Supplementary Fig 2. Raincloud plot of food web properties ####

max.TL <- ggplot(all_data, aes(x = ecosystem.type, y = MaxTL, fill = 	ecosystem.type)) +
  geom_rain(alpha = .5, cov='ecosystem.type',
            boxplot.args.pos = list(width = .1, position = position_nudge(x = -.2))) + 
  ylab('Maximum trophic level') + xlab('Ecosystem type') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt")) 

prim.similarity <- ggplot(all_data, aes(x = ecosystem.type, y = sim.prim.cons, fill = 	ecosystem.type)) +
  geom_rain(alpha = .5, cov='ecosystem.type',
            boxplot.args.pos = list(width = .1, position = position_nudge(x = -.2))) + 
  ylab('1\u00B0 consumer dissimilarity') + xlab('Ecosystem type') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"))

pred.similarity <- ggplot(all_data, aes(x = ecosystem.type, y = sim.sec.cons, fill = 	ecosystem.type)) +
  geom_rain(alpha = .5, cov='ecosystem.type',
            boxplot.args.pos = list(width = .1, position = position_nudge(x = -.2))) + 
  ylab('Predator dissimilarity') + xlab('Ecosystem type') +
  theme(legend.position="none", 
        plot.margin = unit(c(5.5, 5.5, 6, 5.5), "pt"))

FW.Properties <- Richness_main <- grid.arrange(patchworkGrob(max.TL | prim.similarity | pred.similarity)) 

#ggsave(filename = "FW.Properties.svg", FW.Properties, width = 200, height = 75, unit = "mm", dpi = 300)


##################################################################
#### Structural Equation Models including food web stability #####
##################################################################


## Create vectors for effect size plotting ##
flux <- c(rep(c(rep("stability", 3), rep("predation", 3), rep("primary consumption",3)),2))
FW_prop <- rep(c("Taxon richness","Max TL","Trophic dissim."), 6)
effect.type <- factor(c(rep("direct", 9),rep("indirect", 9)), levels=c('indirect', 'direct'))


#### Cross-ecosystem SEM ####

modStability = lme(stability ~ logit(sim.sec.cons) + logit(sim.prim.cons) + logit(sim.prim.cons) + log(MaxTL), random=~1|study_ID, data=all_data) 
  plot(modStability, which=1)
  qqnorm(modStability)
modStabilitya = lme(stability ~ logit(sim.sec.cons) + logit(sim.prim.cons) + logit(sim.prim.cons) + log(MaxTL), weights=varPower(form=~S), random=~1|ecosystem.type/study_ID, data=all_data) 
  plot(modStabilitya, which=1)
  qqnorm(modStabilitya)
modStabilityb = lme(stability ~ log(prim.consumption) + log(second.consumption), weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), random=~1|ecosystem.type/study_ID, data=all_data) # Best model
  plot(modStabilityb, which=1)
  qqnorm(modStabilityb)
anova(modStability, modStabilitya, modStabilityb)
summary(modStabilityb)


### Maximal model
SEM.all <- psem(
  lme(log(S) ~ NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~NPP.scale)), data=all_data),
  lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~NPP.scale)), data=all_data),
  lme(logit(sim.prim.cons) ~ log(S), random=~1|ecosystem.type/study_ID, weights=varComb(varPower(form=~NPP.scale), varExp(form=~S)), data=all_data),
  lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varPower(form=~NPP.scale), data=all_data),
  lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), data=all_data),
  lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=all_data),
  lme(stability ~ log(prim.consumption) + log(second.consumption), weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), random=~1|ecosystem.type/study_ID, data=all_data),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.all)


### Min adequate model
SEM.all2 <- psem(
  lme(log(S) ~ NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~NPP.scale)), data=all_data),
  lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~NPP.scale)), data=all_data),
  lme(logit(sim.prim.cons) ~ log(S), random=~1|ecosystem.type/study_ID, weights=varComb(varPower(form=~NPP.scale), varExp(form=~S)), data=all_data),
  lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varPower(form=~NPP.scale), data=all_data),
  lme(log(prim.consumption) ~ log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), data=all_data),
  lme(log(second.consumption) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|ecosystem.type/study_ID, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=all_data),
  lme(stability ~ log(prim.consumption) + log(second.consumption) + log(S) + logit(sim.prim.cons) + log(MaxTL), weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), random=~1|ecosystem.type/study_ID, data=all_data),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.all2)
results.allSEM <- summary(SEM.all2)$coefficients[c(1:17),c(1:5, 8, 7)]
names(results.allSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
fun <- function(x) {
  formatC(x, format = "f", digits = 3)
}
allSEM_table <- nice_table(results.allSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(allSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/allSEM_table.docx")


## Create results dataframe for summary boxplots
std.effect <- c(
  
  SR.direct.prim <- results.allSEM[8,6],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  SR.direct.stab <- results.allSEM[15,6],
  MTL.direct.stab <- results.allSEM[17,6],
  SIM.direct.stab <- results.allSEM[16,6],
  SR.direct.pred <- results.allSEM[11,6],
  MTL.direct.pred <- results.allSEM[10,6],
  SIM.direct.pred <- 0,

  SR.indirect.prim <- 0,
  MTL.indirect.prim <- 0,
  SIM.indirect.prim <- 0,
  SR.indirect.stab <- (results.allSEM[8,6]*results.allSEM[13,6] + 
                         results.allSEM[4,6]*results.allSEM[16,6] +
                         results.allSEM[2,6]*results.allSEM[17,6] +
                         results.allSEM[2,6]*results.allSEM[10,6]*results.allSEM[17,6]+
                         results.allSEM[11,6]*results.allSEM[14,6]),
  MTL.indirect.stab <- results.allSEM[10,6]*results.allSEM[14,6],
  SIM.indirect.stab <- 0,
  SR.indirect.pred <- (results.allSEM[2,6]*results.allSEM[10,6]),
  MTL.indirect.pred <- 0,
  SIM.indirect.pred <- 0
)
eff.table_all <- data.frame(FW_prop, flux, effect.type, std.effect)

## Produce SEM effects plot for Fig. 3 ##
global.effects <- ggplot(eff.table_all, 
                         aes(x = flux, y = std.effect, fill = flux, pattern = effect.type)) + 
  ylab('Effect size') + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5) +
  geom_bar_pattern(stat = "identity", position = "stack", pattern_color = 'White', pattern_fill = "white", width = .98,
                   pattern_alpha = 0.7, pattern_angle = 55, pattern_density = 0.75, pattern_spacing = 0.1, pattern_size=0) +
  scale_fill_manual(values = c("#F5C2699E", "#C257579E", "#3A67AE9E")) + 
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  coord_flip(ylim=c(-0.5, 0.5)) + scale_y_continuous(labels = c(-0.5, '',0.0, '',0.5)) +
  facet_grid(fct_relevel(FW_prop,'Taxon richness', 'Max TL', 'Trophic dissim.')~., scales = "free_y", labeller = label_wrap_gen(width=10)) +
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), 
        strip.clip = "off", strip.background = element_part_rect(side = "l", linewidth = .5), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))


ggsave("Global effects stability.svg", global.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### Construct ecosystem-specific SEMs as in Fig. 4 ####


#### MARINE ####

mod1 = gls(log(MaxTL) ~ log(S), data=meta.Marine)
plot(mod1, which=1)
qqnorm(mod1)
summary(mod1)

mod2 = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Marine)
plot(mod2, which=1)
qqnorm(mod2)
summary(mod2)

mod2a = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Marine) #Best model
plot(mod2a, which=1)
qqnorm(mod2a)
summary(mod2a)

#anova(mod2, mod2a)

mod3 = gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S), data=meta.Marine) #Best model
plot(mod3, which=1)
qqnorm(mod3)
summary(mod3)

# mod3a = gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Marine)
# plot(mod3a, which=1)
# qqnorm(mod3a)
# summary(mod3a)
# 
# anova(mod3, mod3a)

mod4 = gls(log(prim.consumption) ~ logit(sim.prim.cons), data=meta.Marine) #Best model
plot(mod4, which=1)
qqnorm(mod4)
summary(mod4)

# mod4a = gls(log(prim.consumption) ~ logit(sim.prim.cons), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Marine)
# plot(mod4a, which=1)
# qqnorm(mod4a)
# summary(mod4a)
# 
# anova(mod4, mod4a)

mod5 = gls(log(second.consumption) ~ logit(sim.sec.cons), data=meta.Marine) # Best model
plot(mod5, which=1)
qqnorm(mod5)
summary(mod5)
# 
# mod5a = gls(log(second.consumption) ~ logit(sim.sec.cons), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Marine)
# plot(mod5a, which=1)
# qqnorm(mod5a)
# summary(mod5a)
# 
# anova(mod5, mod5a)

### Maximal model
SEM.Marine2 <- psem(
  
  gls(log(S) ~ NPP.scale, data=meta.Marine),
  gls(log(MaxTL) ~ log(S) + NPP.scale, data=meta.Marine),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Marine),
  gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Marine),
  gls(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), data=meta.Marine),
  gls(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), data=meta.Marine),

  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Marine2)

### Min adequate model
SEM.Marine3 <- psem(
  
  gls(log(S) ~ NPP.scale, data=meta.Marine),
  gls(log(MaxTL) ~ log(S) + NPP.scale, data=meta.Marine),
  gls(logit(sim.prim.cons) ~ log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Marine),
  gls(logit(sim.sec.cons) ~ log(S) + log(MaxTL) + NPP.scale, data=meta.Marine),
  gls(log(prim.consumption) ~ log(S) + NPP.scale + log(MaxTL), data=meta.Marine),
  gls(log(second.consumption) ~ log(MaxTL) + NPP.scale, data=meta.Marine),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Marine3)
results.marineSEM <- summary(SEM.Marine3)$coefficients[c(1:8,12),c(1:5, 8, 7)]
names(results.marineSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
fun <- function(x) {
  formatC(x, format = "f", digits = 3)
}
marineSEM_table <- nice_table(results.marineSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(marineSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/marineSEM_table.docx")


## Create results dataframe for summary boxplots
std.effect <- c(
  SR.direct.pred <- 0,
  MTL.direct.pred <- 0,
  SIM.direct.pred <- results.marineSEM[8,6],
  SR.direct.prim <- results.marineSEM[5,6],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  
  SR.indirect.pred <- ((results.marineSEM[1,6]*results.marineSEM[4,6]) + (results.marineSEM[3,6]*results.marineSEM[8,6])),
  MTL.indirect.pred <- (results.marineSEM[4,6]*results.marineSEM[8,6]),
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
        strip.clip = "off", strip.background = element_part_rect(side = "r", linewidth = .5), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

#ggsave("Marine effects.svg", marine.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### SOIL ####
mod1 = gls(log(MaxTL) ~ log(S), data=meta.Soils) #Best model
plot(mod1, which=1)
qqnorm(mod1)
summary(mod1)

# mod1a = gls(log(MaxTL) ~ log(S), weights=varIdent(form=~1|study_ID), data=meta.Soils)
# plot(mod1, which=1)
# qqnorm(mod1)
# summary(mod1)
# 
# anova(mod1, mod1a)

mod2 = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Soils)
plot(mod2, which=1)
qqnorm(mod2)
summary(mod2)


mod3 = gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S), data=meta.Soils)
plot(mod3, which=1)
qqnorm(mod3)
summary(mod3)

# mod3a = gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Soils)
# plot(mod3a, which=1)
# qqnorm(mod3a)
# summary(mod3a)
# 
# anova(mod3, mod3a)

mod4 = gls(log(prim.consumption) ~ logit(sim.prim.cons) + log(MaxTL), data=meta.Soils)
plot(mod4, which=1)
qqnorm(mod4)
summary(mod4)

mod4a = gls(log(prim.consumption) ~ logit(sim.prim.cons) + log(MaxTL), 
            weights=varComb(varIdent(form=~1|study_ID)), data=meta.Soils) #Best model
plot(mod4a, which=1)
qqnorm(mod4a)
summary(mod4a)

anova(mod4, mod4a)

mod5 = gls(log(second.consumption) ~ logit(sim.sec.cons) + log(MaxTL), data=meta.Soils)
plot(mod5, which=1)
qqnorm(mod5)
summary(mod5)

mod5a = gls(log(second.consumption) ~ logit(sim.sec.cons) + log(MaxTL), weights=varComb(varIdent(form=~1|study_ID)), data=meta.Soils) 
plot(mod5a, which=1)
qqnorm(mod5a)
summary(mod5a)

mod5b = gls(log(second.consumption) ~ logit(sim.sec.cons) + log(MaxTL), 
            weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Soils) #Best model
plot(mod5a, which=1)
qqnorm(mod5a)
summary(mod5a)

anova(mod5, mod5a, mod5b)


### Maximal model
SEM.Soils2 <- psem(

  gls(log(S) ~ NPP.scale, data=meta.Soils),
  gls(log(MaxTL) ~ log(S) + NPP.scale, data=meta.Soils),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Soils),
  gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, data=meta.Soils),
  gls(log(prim.consumption) ~ logit(sim.prim.cons) + log(MaxTL), weights=varComb(varIdent(form=~1|study_ID)), data=meta.Soils),
  gls(log(second.consumption) ~ logit(sim.sec.cons) + log(MaxTL), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Soils),

  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Soils2)


### Min adequate model
SEM.Soils3 <- psem(
  
  gls(log(S) ~ NPP.scale, data=meta.Soils),
  gls(log(MaxTL) ~ log(S), data=meta.Soils),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Soils),
  gls(logit(sim.sec.cons) ~ NPP.scale, data=meta.Soils),
  gls(log(prim.consumption) ~ log(S) + NPP.scale, weights=varComb(varIdent(form=~1|study_ID)), data=meta.Soils),
  gls(log(second.consumption) ~ log(S) + log(MaxTL), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Soils),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)

summary(SEM.Soils3)
results.soilsSEM <- summary(SEM.Soils3)$coefficients[c(1:7, 11),c(1:5, 8, 7)]
names(results.soilsSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
soilsSEM_table <- nice_table(results.soilsSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(soilsSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/soilsSEM_table.docx")


## Create results dataframe for summary boxplots
std.effect <- c(
  SR.direct.pred <- results.soilsSEM[6,6],
  MTL.direct.pred <- results.soilsSEM[7,6],
  SIM.direct.pred <- 0,
  SR.direct.prim <- results.soilsSEM[5,6],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,
  
  SR.indirect.pred <- (results.soilsSEM[1,6]*results.soilsSEM[7,6]),
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
        strip.clip = "off", strip.background = element_part_rect(side = "l", linewidth = .5), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

#ggsave("Soils effects.svg", soils.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### STREAMs ####

mod1 = lme(log(MaxTL) ~ log(S), random=~1|study_ID, data=meta.Streams)
plot(mod1, which=1)
qqnorm(mod1)
summary(mod1)

mod2 = lme(logit(sim.prim.cons) ~ log(MaxTL) + log(S), weights=varIdent(form=~1|study_ID), random=~1|study_ID, data=meta.Streams)
plot(mod2, which=1)
qqnorm(mod2)
summary(mod2)

mod3 = lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), random=~1|study_ID, data=meta.Streams)
plot(mod3, which=1)
qqnorm(mod3)
summary(mod3)

mod4 = lme(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), random=~1|study_ID, data=meta.Streams)
plot(mod4, which=1)
qqnorm(mod4)
summary(mod4)

mod5 = lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), random=~1|study_ID, data=meta.Streams)
plot(mod5, which=1)
qqnorm(mod5)
summary(mod5)


### Maximal model
SEM.Streams2 <- psem(
  
  lme(log(S) ~ NPP.scale, random=~1|study_ID, data=meta.Streams),
  lme(log(MaxTL) ~ log(S) + NPP.scale, random=~1|study_ID, data=meta.Streams),
  lme(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, random=~1|study_ID, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Streams),
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
  
  lme(log(S) ~ NPP.scale, random=~1|study_ID, data=meta.Streams),
  lme(log(MaxTL) ~ log(S), random=~1|study_ID, data=meta.Streams),
  lme(logit(sim.prim.cons) ~ log(MaxTL), random=~1|study_ID, weights=varIdent(form=~1|study_ID), data=meta.Streams),
  lme(logit(sim.sec.cons) ~ log(S) + log(MaxTL) + NPP.scale, random=~1|study_ID, weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Streams),
  lme(log(prim.consumption) ~ log(S) , random=~1|study_ID, data=meta.Streams),
  lme(log(second.consumption) ~ log(MaxTL) + logit(sim.sec.cons), random=~1|study_ID, data=meta.Streams),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)


summary(SEM.Streams3)
results.streamsSEM <- summary(SEM.Streams3)$coefficients[c(1:8, 12),c(1:5, 8, 7)]
names(results.streamsSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
streamsSEM_table <- nice_table(results.streamsSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(streamsSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/streamsSEM_table.docx")

## Create results dataframe for summary boxplots
std.effect <- c(
  SR.direct.pred <- 0,
  MTL.direct.pred <- results.streamsSEM[7,6],
  SIM.direct.pred <- results.streamsSEM[8,6],
  SR.direct.prim <- results.streamsSEM[5,6],
  MTL.direct.prim <- 0,
  SIM.direct.prim <- 0,

  SR.indirect.pred <- (results.streamsSEM[3,6]*results.streamsSEM[8,6]),
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
        strip.clip = "off", strip.background = element_part_rect(side = "r", linewidth = .5), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

#ggsave("Streams effects.svg", streams.effects, width = 6, height = 9, units = "cm", bg='transparent')




#### Lakes ####
mod0 = gls(log(S) ~ NPP.scale, weights=varExp(form=~S), data=meta.Lakes)
plot(mod0, which=1)
qqnorm(mod0)
summary(mod0)

mod1 = gls(log(MaxTL) ~ log(S), data=meta.Lakes)
plot(mod1, which=1)
qqnorm(mod1)
summary(mod1)

mod1a = gls(log(MaxTL) ~ log(S) + NPP.scale, weights=varPower(form=~S), data=meta.Lakes)
plot(mod1a, which=1)
qqnorm(mod1a)
summary(mod1a)

# mod2 = gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S), data=meta.Lakes, na.action = na.omit)
# plot(mod2, which=1)
# qqnorm(mod2)
# summary(mod2)

mod2a = gls(logit(sim.prim.cons) ~ log(S) + NPP.scale, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), data=meta.Lakes, na.action = na.omit)
plot(mod2a, which=1)
qqnorm(mod2a)
summary(mod2a)

mod3 = gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varComb(varExp(form=~S), varExp(form=~MaxTL)), data=meta.Lakes)
plot(mod3, which=1)
qqnorm(mod3)
summary(mod3)

mod4 = gls(log(prim.consumption) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Lakes)
plot(mod4, which=1)
qqnorm(mod4)
summary(mod4)

mod5 = gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S),  data=meta.Lakes)
plot(mod5, which=1)
qqnorm(mod5)
summary(mod5)


### Maximal model
SEM.Lakes2 <- psem(

  gls(log(S) ~ NPP.scale, data=meta.Lakes),
  gls(log(MaxTL) ~ log(S) + NPP.scale, weights=varPower(form=~S), data=meta.Lakes),
  gls(logit(sim.prim.cons) ~ log(MaxTL) + log(S) + NPP.scale, weights=varComb(varPower(form=~S)), data=meta.Lakes),
  gls(logit(sim.sec.cons) ~ log(MaxTL) + log(S) + NPP.scale, weights=varPower(form=~log(MaxTL)), data=meta.Lakes),
  gls(log(prim.consumption) ~ log(MaxTL) + logit(sim.prim.cons), data=meta.Lakes),
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
  gls(log(MaxTL) ~ log(S) + NPP.scale, weights=varPower(form=~S), data=meta.Lakes),
  gls(logit(sim.prim.cons) ~ log(S) + NPP.scale, weights=varComb(varIdent(form=~1|study_ID), varExp(form=~S)), data=meta.Lakes),
  gls(logit(sim.sec.cons) ~ log(MaxTL), weights=varPower(form=~log(MaxTL)), data=meta.Lakes),
  gls(log(prim.consumption) ~ log(MaxTL) + log(S), weights=varComb(varIdent(form=~1|study_ID), varPower(form=~S)), data=meta.Lakes),
  gls(log(second.consumption) ~ logit(sim.sec.cons) + log(S), data=meta.Lakes),
  
  logit(sim.prim.cons) %~~% logit(sim.sec.cons),
  log(prim.consumption) %~~% logit(sim.sec.cons),
  log(second.consumption) %~~% logit(sim.prim.cons),
  log(second.consumption) %~~% log(prim.consumption)
)


summary(SEM.Lakes3)
results.lakesSEM <- summary(SEM.Lakes3)$coefficients[c(1:8, 12),c(1:5, 8, 7)]
names(results.lakesSEM) <- c("Response", "Predictor", "Estimate", "SE", "df", "Std. Estimate", "p")
lakesSEM_table <- nice_table(results.lakesSEM, col.format.custom = c(3:4,6:7), format.custom = "fun")

#flextable::save_as_docx(lakesSEM_table, path = "C:/Users/barnesa/OneDrive - The University of Waikato/FuSED/Data/lakesSEM_table.docx")

## Create results daframe for summary boxplots
std.effect <- c(
  SR.direct.pred <- results.lakesSEM[7,6],
  MTL.direct.pred <- 0,
  SIM.direct.pred <- results.lakesSEM[8,6],
  SR.direct.prim <- results.lakesSEM[5,6],
  MTL.direct.prim <- results.lakesSEM[6,6],
  SIM.direct.prim <- 0,

  SR.indirect.pred <- (results.lakesSEM[1,6]*results.lakesSEM[4,6]*results.lakesSEM[8,6]+results.lakesSEM[3,6]*results.lakesSEM[8,6]),
  MTL.indirect.pred <- (results.lakesSEM[4,6]*results.lakesSEM[8,6]),
  SIM.indirect.pred <- 0,
  SR.indirect.prim <- results.lakesSEM[1,6]*results.lakesSEM[6,6],
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
        strip.clip = "off", strip.background = element_part_rect(side = "l", linewidth = .5), strip.text.y = element_text(size = 13), 
        legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), 
        axis.title.y=element_blank(), axis.text=element_text(size=13), axis.title=element_text(size=13))

#ggsave("Lakes effects.svg", lakes.effects, width = 6, height = 9, units = "cm", bg='transparent')


