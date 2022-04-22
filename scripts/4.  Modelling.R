### 4. Modelling ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Date created: 29th June 2020 ###

#### LOAD PACKAGES ####

library(ggregplot)
library(tidyverse)
library(esquisse)
library(rbenchmark)
library(DataCombine) 
library(ggeffects)
library(corrplot)
library(brms)
library(rstanarm)
library(sjPlot)
library(sjstats)
library(sjmisc)
library(modelr)
library(tidybayes)
library(effects)
library(ncar)
library(ggpubr)
library(clickR)
library(viridis)    

#### SET UP WORKSPACE ####
#Import Data
tea_full <- read.csv(file="data/fullpropertieslong.csv") # read in soil data

#remove 'Inf' values from Green colums
tea_full = filter(tea_full, Green != "-Inf")

#### CORRELATION MATRIX ####
colnames(tea_full) # view column numbers
tea_cor <- tea_full[,c(1,2,3,5,10,11,23,24,25,26)] # select numeric variables
res <- cor(tea_cor, use = "complete.obs")
round(res, 2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
cor(tea_full$thermsum,tea_full$Elevation,  method = c("pearson"), use = "complete.obs")
p.mat <- cor.mtest(res)$p
round(p.mat, 2)
head(p.mat[, 1:5])

#### Run Bayesian models ####
## Scale variables 
tea_full$Soilmoist_scaled <- scale(tea_full$Soilmoist, center = 0)
tea_full$ALT_scaled <- scale(tea_full$ALT, center = 0)
tea_full$meantemp_scaled <- scale(tea_full$meantemp, center = 0)
tea_full$Green_scaled <- scale(tea_full$Green, center = 0)
tea_full$Rooibos_scaled <- scale(tea_full$Rooibos, center = 0)
tea_full$elevation_scaled <- scale(tea_full$elevation, center = 0)
tea_full$slope_scaled <- scale(tea_full$slope, center = 0)
tea_full$stab_scaled <- scale(tea_full$stab, center = 0)
tea_full$rate_scaled <- scale(tea_full$rate, center = 0)
tea_full$aspect_scaled <- scale(tea_full$aspect, center = 0)
tea_full$therm_scaled <- scale(tea_full$thermsum, center = 0)


scaling_attributes <- data.frame(predictor = c("Soilmoist",
                                               "ALT",
                                               "meantemp",
                                               "Green",
                                               "Rooibos",
                                               "elevation",
                                               "slope",
                                               "stab",
                                               "rate",
                                               "aspect",
                                               "therm"),
                                 scale = c(
                                   attributes(tea_full$Soilmoist_scaled)$"scaled:scale",
                                   attributes(tea_full$ALT_scaled)$"scaled:scale",
                                   attributes(tea_full$meantemp_scaled)$"scaled:scale",
                                   attributes(tea_full$Green_scaled)$"scaled:scale",
                                   attributes(tea_full$Rooibos_scaled)$"scaled:scale",
                                   attributes(tea_full$elevation_scaled)$"scaled:scale",
                                   attributes(tea_full$slope_scaled)$"scaled:scale",
                                   attributes(tea_full$stab_scaled)$"scaled:scale",
                                   attributes(tea_full$rate_scaled)$"scaled:scale",
                                   attributes(tea_full$aspect_scaled)$"scaled:scale",
                                   attributes(tea_full$therm_scaled)$"scaled:scale"
                                 ))

# aspect as factor
tea_full$aspect <- as.integer(tea_full$aspect)

# group aspect by cardinal directions
tea_full <- tea_full %>% 
  mutate(aspect_class = case_when(
    between(tea_full$aspect, 45, 134) ~ "East",
    between(tea_full$aspect, 135, 224) ~ "South",
    between(tea_full$aspect, 225, 314) ~ "West",
    between(tea_full$aspect, 315, 360) ~ "North",
    between(tea_full$aspect, 0, 44) ~ "North"
  ))


# 1 - Soil Property Models with default priors
# green mass loss with ALT interactive term
alttherm_green_m <- brm(bf(Green_scaled ~  ALT_scaled*therm_scaled + (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 2, cores = 2,
                       warmup = 2000)
plot_model(alttherm_green_m)
summary(alttherm_green_m)

alttherm_green_alt <- brm(bf(Green_scaled ~  ALT_scaled+therm_scaled + (1|Plot),
                           family = brmsfamily('Gaussian')), data = tea_full,
                        iter = 8000,
                        chains = 2, cores = 2,
                        warmup = 2000)
plot_model(alttherm_green_alt)
summary(alttherm_green_alt)


# Make null model
alttherm_green_null <- brm(bf(Green_scaled ~  ALT_scaled + therm_scaled,
                             family = brmsfamily('Gaussian')), data = tea_full,
                          iter = 8000,
                          chains = 2, cores = 2,
                          warmup = 2000)


# save model output
save(alttherm_green_m, file = "models/scaled/alttherm_green_m.RData")
save(alttherm_green_null, file = "models/scaled/alttherm_green_null.RData")
save(alttherm_green_alt, file = "models/scaled/alttherm_green_alt.RData")

# green mass loss with soil interactive term
soiltherm_green_m <- brm(bf(Green_scaled ~  Soilmoist_scaled*therm_scaled + (1|Plot),
                           family = brmsfamily('Gaussian')), data = tea_full,
                        iter = 8000,
                        chains = 2, cores = 2,
                        warmup = 2000)

plot(soiltherm_green_m)
summary(soiltherm_green_m)


soiltherm_green_alt <- brm(bf(Green_scaled ~  Soilmoist_scaled + therm_scaled + (1|Plot),
                            family = brmsfamily('Gaussian')), data = tea_full,
                         iter = 8000,
                         chains = 2, cores = 2,
                         warmup = 2000)

plot_model(soiltherm_green_m_alt)
summary(soiltherm_green_m_alt)


soiltherm_green_null <- brm(bf(Green_scaled ~  Soilmoist_scaled + therm_scaled,
                               family = brmsfamily('Gaussian')), data = tea_full,
                            iter = 8000,
                            chains = 2, cores = 2,
                            warmup = 2000)

# save model output
save(soiltherm_green_m, file = "models/scaled/soiltherm_green_mm.RData")
save(soiltherm_green_null, file = "models/scaled/soiltherm_green_null.RData")
save(soiltherm_green_alt, file = "models/scaled/soiltherm_green_alt.RData")



# therm version
alttherm_red_m <- brm(bf(Rooibos_scaled ~  ALT_scaled*therm_scaled+ (1|Plot),
                         family = brmsfamily('Gaussian')), data = tea_full,
                      iter = 8000,
                      chains = 2, cores = 2,
                      warmup = 2000)

plot_model(alttherm_red_m)

alttherm_red_alt <- brm(bf(Rooibos_scaled ~  ALT_scaled+therm_scaled+ (1|Plot),
                         family = brmsfamily('Gaussian')), data = tea_full,
                      iter = 8000,
                      chains = 2, cores = 2,
                      warmup = 2000)

plot_model(alttherm_red_alt)

alttherm_red_null <- brm(bf(Rooibos_scaled ~  ALT_scaled+therm_scaled,
                            family = brmsfamily('Gaussian')), data = tea_full,
                         iter = 8000,
                         chains = 2, cores = 2,
                         warmup = 2000)

# save model output
save(alttherm_red_m, file = "models/scaled/alttherm_red_mm.RData")
save(alttherm_red_null, file = "models/scaled/alttherm_red_null.RData")
save(alttherm_red_alt, file = "models/scaled/alttherm_red_alt.RData")

# red mass loss with soil interactive term

# therm version
soiltherm_red_m <- brm(bf(Rooibos_scaled ~  Soilmoist_scaled*therm_scaled+ (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 2, cores = 2,
                       warmup = 2000)

soiltherm_red_alt <- brm(bf(Rooibos_scaled ~  Soilmoist_scaled+therm_scaled+ (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 2, cores = 2,
                       warmup = 2000)

soiltherm_red_null <- brm(bf(Rooibos_scaled ~  Soilmoist_scaled+therm_scaled,
                             family = brmsfamily('Gaussian')), data = tea_full,
                          iter = 8000,
                          chains = 2, cores = 2,
                          warmup = 2000)

# save model output
save(soiltherm_red_m, file = "models/scaled/soiltherm_red_mm.RData")
save(soiltherm_red_null, file = "models/scaled/soiltherm_red_null.RData")
save(soiltherm_red_alt, file = "models/scaled/soiltherm_red_alt.RData")

# decomp rate with ALT interactive term

# therm sum
alttherm_rate_m <- brm(bf(rate_scaled ~  ALT_scaled*therm_scaled + (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 4, cores = 2,
                       warmup = 2000)

alttherm_rate_alt <- brm(bf(rate_scaled ~  ALT_scaled + therm_scaled + (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 4, cores = 2,
                       warmup = 2000)

alttherm_rate_null<- brm(bf(rate_scaled ~  ALT_scaled+therm_scaled,
                            family = brmsfamily('Gaussian')), data = tea_full,
                         iter = 8000,
                         chains = 2, cores = 2,
                         warmup = 2000)

# save model output
save(alttherm_rate_m, file = "models/scaled/alttherm_rate_mm.RData")
save(alttherm_rate_null, file = "models/scaled/alttherm_rate_null.RData")
save(alttherm_rate_alt, file = "models/scaled/alttherm_rate_alt.RData")


# decomp rate with soil interactive term
# therm sum
soiltherm_rate_m <- brm(bf(rate_scaled ~  Soilmoist_scaled*therm_scaled+ (1|Plot),
                           family = brmsfamily('Gaussian')), data = tea_full,
                        iter = 8000,
                        chains = 2, cores = 2,
                        warmup = 2000)

soiltherm_rate_alt <- brm(bf(rate_scaled ~  Soilmoist_scaled+therm_scaled+ (1|Plot),
                           family = brmsfamily('Gaussian')), data = tea_full,
                        iter = 8000,
                        chains = 2, cores = 2,
                        warmup = 2000)

soiltherm_rate_null <- brm(bf(rate_scaled ~  Soilmoist_scaled+therm_scaled,
                              family = brmsfamily('Gaussian')), data = tea_full,
                           iter = 8000,
                           chains = 2, cores = 2,
                           warmup = 2000)


# save model output
save(soiltherm_rate_m, file = "models/scaled/soiltherm_rate_mm.RData")
save(soiltherm_rate_null, file = "models/scaled/soiltherm_rate_null.RData")
save(soiltherm_rate_alt, file = "models/scaled/soiltherm_rate_alt.RData")

# stabilisation rate with ALT interactive term

# therm version
alttherm_stab_m <- brm(bf(stab_scaled ~  ALT_scaled*therm_scaled+ (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 2, cores = 2,
                       warmup = 2000)

alttherm_stab_alt <- brm(bf(stab_scaled ~  ALT_scaled+therm_scaled+ (1|Plot),
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 2, cores = 2,
                       warmup = 2000)

summary(alttherm_stab_alt)

alttherm_stab_null <- brm(bf(stab_scaled ~  ALT_scaled+therm_scaled,
                             family = brmsfamily('Gaussian')), data = tea_full,
                          iter = 8000,
                          chains = 2, cores = 2,
                          warmup = 2000)

# save model output

save(alttherm_stab_m, file = "models/scaled/alttherm_stab_m.RData")
save(alttherm_stab_null, file = "models/scaled/alttherm_stab_null.RData")
save(alttherm_stab_alt, file = "models/scaled/alttherm_stab_alt.RData")


# stabilisation rate with soil interactive term
# therm version
soiltherm_stab_m <- brm(bf(stab_scaled ~  Soilmoist_scaled*therm_scaled + (1|Plot),
                           family = brmsfamily('Gaussian')), data = tea_full,
                        iter = 8000,
                        chains = 2, cores = 2,
                        warmup = 2000)

soiltherm_stab_alt <- brm(bf(stab_scaled ~  Soilmoist_scaled+therm_scaled + (1|Plot),
                           family = brmsfamily('Gaussian')), data = tea_full,
                        iter = 8000,
                        chains = 2, cores = 2,
                        warmup = 2000)

soiltherm_stab_null <- brm(bf(stab_scaled ~  Soilmoist_scaled+therm_scaled,
                              family = brmsfamily('Gaussian')), data = tea_full,
                           iter = 8000,
                           chains = 2, cores = 2,
                           warmup = 2000)

# save model output
save(soiltherm_stab_m, file = "models/scaled/soiltherm_stab_mm.RData")
save(soiltherm_stab_null, file = "models/scaled/soiltherm_stab_null.RData")
save(soiltherm_stab_alt, file = "models/scaled/soiltherm_stab_alt.RData")


# 2 - Topography Models - scaled because otherwise we see non-convergence!
topo_green_m <- brm(bf(Green_scaled ~ elevation_scaled + slope_scaled + aspect_class + (1|Plot),
                       family = brmsfamily('Gaussian')), data = tea_full,
                    iter = 8000,
                    chains = 2, cores = 2,
                    warmup = 2000)

topo_green_null <- brm(bf(Green_scaled ~ elevation_scaled + slope_scaled + aspect_class,
                          family = brmsfamily('Gaussian')), data = tea_full,
                       iter = 8000,
                       chains = 2, cores = 2,
                       warmup = 2000)

topo_red_m <- brm(bf(Rooibos_scaled ~ elevation_scaled + slope_scaled + aspect_class + (1|Plot),
                     family = brmsfamily('Gaussian')), data = tea_full,
                  iter = 8000,
                  chains = 2, cores = 2,
                  warmup = 2000)

topo_red_null <- brm(bf(Rooibos_scaled ~ elevation_scaled + slope_scaled + aspect_class,
                        family = brmsfamily('Gaussian')), data = tea_full,
                     iter = 8000,
                     chains = 2, cores = 2,
                     warmup = 2000)

topo_rate_m <- brm(bf(rate_scaled ~ elevation_scaled + slope_scaled + aspect_class + (1|Plot),
                      family = brmsfamily('Gaussian')), data = tea_full,
                   iter = 8000,
                   chains = 2, cores = 2,
                   warmup = 2000)

topo_rate_null <- brm(bf(rate_scaled ~ elevation_scaled + slope_scaled + aspect_class,
                         family = brmsfamily('Gaussian')), data = tea_full,
                      iter = 8000,
                      chains = 2, cores = 2,
                      warmup = 2000)

topo_stab_m <- brm(bf(stab_scaled ~ elevation_scaled + slope_scaled + aspect_class + (1|Plot),
                      family = brmsfamily('Gaussian')), data = tea_full,
                   iter = 8000,
                   chains = 2, cores = 2,
                   warmup = 2000)

topo_stab_null <- brm(bf(stab_scaled ~ elevation_scaled + slope_scaled + aspect_class,
                         family = brmsfamily('Gaussian')), data = tea_full,
                      iter = 8000,
                      chains = 2, cores = 2,
                      warmup = 2000)

# save model output
save(topo_green_m, file = "models/scaled/topo_green_m.RData")
save(topo_red_m, file = "models/scaled/topo_red_m.RData")
save(topo_stab_m, file = "models/scaled/topo_stab_m.RData")
save(topo_rate_m, file = "models/scaled/topo_rate_m.RData")
save(topo_green_null, file = "models/scaled/topo_green_null.RData")
save(topo_red_null, file = "models/scaled/topo_red_null.RData")
save(topo_stab_null, file = "models/scaled/topo_stab_null.RData")
save(topo_rate_null, file = "models/scaled/topo_rate_null.RData")

# examine topo model results
summary(topo_green_m) # negligible / slight negative all
summary(topo_red_m) # negligible / slight negative all
summary(topo_rate_m) # negligible 
summary(topo_stab_m) # negligible 

#### Load Bayesian Models ####
load("models/topo_green_m 2.Rdata")
load("models/topo_red_m 2.Rdata")
load("models/topo_rate_m 2.Rdata")
load("models/topo_stab_m 2.Rdata")
load("models/scaled/soiltherm_green_mm.Rdata")
load("models/scaled/soiltherm_red_mm.Rdata")
load("models/scaled/soiltherm_rate_mm.Rdata")
load("models/scaled/soiltherm_stab_mm.Rdata")
load("models/scaled/alttherm_green_m.Rdata")
load("models/scaled/alttherm_red_mm.Rdata")
load("models/scaled/alttherm_rate_mm.Rdata")
load("models/scaled/alttherm_stab_m.Rdata")


##### PLOT: ELEVATION #####

ele_predictions <- ggpredict(topo_green_m, terms = c("elevation_scaled"))

(ele_plot <- ggplot() +
    geom_line(data = ele_predictions, aes(y = predicted*37.38227087, x = x*45.60150226),
              size = 2, colour = "#E00D5E") +
    geom_ribbon(data = ele_predictions, aes(ymin = conf.low*37.38227087, ymax = conf.high*37.38227087, 
                                            x = x*45.60150226), fill = "#E845BA", alpha = 0.1, colour = "#E845BA") +
    geom_point(data = tea_full, aes(x = elevation_scaled*45.60150226,  y = Green_scaled*37.38227087, shape = Plot),
               alpha = 1) +
    scale_shape_manual(values=seq(0,8)) +
    theme_classic() +
    ggtitle("a)") +
    theme_classic(base_size = 16) +
    labs(y = 'Green Tea Mass Loss %', x = "Elevation (m)"))
(ele_plot <- ele_plot  + theme(legend.position = "none"))


ele_stab_predictions <- ggpredict(topo_stab_m, terms = c("elevation_scaled"))


(elestab_plot <- ggplot() +
    geom_line(data = ele_stab_predictions, aes(y = predicted*0.56594405, x = x*45.60150226),
              size = 2, colour = "#E00D5E") +
    geom_ribbon(data = ele_stab_predictions, aes(ymin = conf.low*0.565944059, ymax = conf.high*0.56594405, 
                                                 x = x*45.60150226), fill = "#E845BA", alpha = 0.1, colour = "#E845BA") +
    geom_point(data = tea_full, aes(x = elevation_scaled*45.60150226,  y = stab_scaled*0.56594405, colour = Plot),
               alpha = 1) +
    scale_color_viridis_d(option = "plasma") +
    theme_classic(base_size = 16) +
    theme_classic() +
    labs(y = 'Stabilisation Factor (S)', x = "Elevation (m)"))



ele_rate_predictions <- ggpredict(topo_rate_m, terms = c("elevation_scaled"))


(elerate_plot <- ggplot() +
    geom_line(data = ele_rate_predictions, aes(y = predicted*0.03314069, x = x*45.60150226),
              size = 2, colour = "#E00D5E") +
    geom_ribbon(data = ele_rate_predictions, aes(ymin = conf.low*0.03314069, ymax = conf.high*0.03314069, 
                                                 x = x*45.60150226), fill = "#E845BA", alpha = 0.1, colour = "#E845BA") +
    geom_point(data = tea_full, aes(x = elevation_scaled*45.60150226,  y = rate_scaled*0.03314069, shape = Plot),
               alpha = 1) +
    scale_shape_manual(values=seq(0,8)) +
    theme_classic() +
    ggtitle("c)") +
    theme_classic(base_size = 16) +
    labs(y = 'Decomposition Rate (k)', x = "Elevation (m)"))
(elerate_plot <- elerate_plot  + theme(legend.position = "none"))


ele_stab_predictions <- ggpredict(topo_stab_m, terms = c("Elevation_scaled"))


(elev_modelplots <- ggarrange(ele_plot,  
                              elerate_plot,  nrow=2))



#### PLOT ASPECT ####
asp_green_predictions <- ggpredict(topo_green_m, terms = c("aspect_class"))
asp_stab_predictions <- ggpredict(topo_stab_m, terms = c("aspect_class"))
asp_rate_predictions <- ggpredict(topo_rate_m, terms = c("aspect_class"))

green_asp <- plot(asp_green_predictions)
(green_asp <- green_asp + theme_classic(base_size = 16) + 
    labs(y = 'Green Tea Mass Loss (scaled)', x = "Cardinal Aspect",title="b)"))

stab_asp <- plot(asp_stab_predictions)
stab_asp <- stab_asp + theme_classic(base_size = 16) + 
  labs(y = 'Stabilisation Factor (scaled)', x = "Cardinal Aspect",title=NULL)

rate_asp <- plot(asp_rate_predictions)
rate_asp <-rate_asp + theme_classic(base_size = 16) +  ggtitle("d)") +
  labs(y = 'Decomposition Rate (scaled)', x = "Cardinal Aspect",title="d)")

(asp_modelplots <- ggarrange(green_asp,
                             rate_asp,
                             nrow=2))

(topo_modelplots <- ggarrange(elev_modelplots, 
                              asp_modelplots, 
                              nrow=1, ncol = 2))

ggsave(topo_modelplots, filename = "figures/manuscript/elev_asp_panel.png", height = 7, width = 7)



#### SCALING PLOTS ####

tea_full$Plot <- as.factor(tea_full$Plot)


tea_mean <- tea_full %>% group_by(Plot) %>% summarise(Green_mean = mean(Green,na.rm = T), Soilmoist_mean = mean(Soilmoist), 
                                                      ALT_mean = mean(ALT), meantemp_mean = mean(meantemp, na.rm = T), 
                                                      thermsum_mean = mean(thermsum, na.rm = T),
                                                      Green_SE = sd(Green, na.rm = T)/sqrt(length(na.omit(Green))), 
                                                      Soilmoist_SE = sd(Soilmoist)/sqrt(length(Soilmoist)), 
                                                      ALT_SE = sd(ALT)/sqrt(length(ALT)), 
                                                      meantemp_SE = sd(meantemp, na.rm = T)/sqrt(length(na.omit(meantemp))),
                                                      thermsum_SE = sd(thermsum,na.rm = T)/sqrt(length((thermsum))))

tea_mean$Plot <- as.factor(tea_mean$Plot)
#c("FO", "FP", "DG", "ECS", "ECW", "ECN", "CHKO", "CHHE", "WET", "CONTROL")
tea_full$Plot <- as.factor(tea_full$Plot)
#c("FO", "FP", "DG", "ECS", "ECW", "ECN", "CHKO", "CHHE", "WET", "CONTROL")

# reorder factor levels for plot
tea_full$Plot <- factor(tea_full$Plot, levels =  c("FO", "DG", "FP",
                                                   "ECS", "ECW", "ECN",
                                                   "CHKO","CHHE"))

tea_mean$Plot <- factor(tea_mean$Plot, levels =  c("FO", "DG", "FP",
                                                   "ECS", "ECW", "ECN",
                                                   "CHKO","CHHE"))

# figure of soil moisture relationship 
(Soilmoist_fig <- ggplot(data = tea_full) +
    geom_point(aes(y = Green, x = Soilmoist, colour = Plot)) +
    geom_smooth(aes(y = Green, x = Soilmoist, colour = Plot, group = Plot), method = "lm", se = F, size = 0.5) +
    geom_smooth(aes(y = Green, x = Soilmoist), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    xlim(15, 50) +
    ylim(30, 43) +
    labs(y = 'Green Tea Mass Loss %', x = "Soil Moisture %") +
    theme_classic())
(Soilmoist_fig <- Soilmoist_fig +  theme(legend.position = "none"))


ggsave(Soilmoist_fig, file = "figures/Soilmoist.png", height = 3, width = 4)

(ALT_fig <- ggplot(data = tea_full) +
    geom_point(aes(y = Green, x = ALT, colour = Plot)) +
    geom_smooth(aes(y = Green, x = ALT, colour = Plot, group = Plot), method = "lm", se = F, size = 0.5) +
    geom_smooth(aes(y = Green, x = ALT), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    xlim(15, 50) +
    ylim(30, 43) +
    labs(y = 'Green Tea Mass Loss %', x = "Active Layer Thickness (cm)") +
    theme_classic())

(ALT_fig <- ALT_fig +  theme(legend.position = "none"))

ggsave(ALT_fig, file = "figures/ALT.png", height = 3, width = 4)


(thermsum_fig <- ggplot(data = tea_full) +
    geom_point(aes(y = Green, x = thermsum, colour = Plot)) +
    geom_smooth(aes(y = Green, x = thermsum, colour = Plot, group = Plot), method = "lm", se = F, size = 0.5) +
    geom_smooth(aes(y = Green, x = thermsum), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    labs(y = 'Green Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)") +
    theme_classic())
(thermsum_fig <- thermsum_fig +  theme(legend.position = "none"))


ggsave(thermsum_fig, file = "figures/thermsum.png", height = 3, width = 4)

tea_mean$Plot <- factor(tea_mean$Plot, levels =  c("FO", "DG", "FP",
                                                   "ECS", "ECW", "ECN",
                                                   "CHKO","CHHE"))
lower_y <- tea_mean$Green_mean - tea_mean$Green_SE
upper_y <- tea_mean$Green_mean + tea_mean$Green_SE
lower_x <- tea_mean$Soilmoist_mean - tea_mean$Soilmoist_SE
upper_x <- tea_mean$Soilmoist_mean + tea_mean$Soilmoist_SE

(Soilmoist_mean_fig <- ggplot(data = tea_mean) +
    geom_point(aes(y = Green_mean, x = Soilmoist_mean, colour = Plot)) +
    geom_errorbar(aes(y = Green_mean, x = Soilmoist_mean, ymin = lower_y, ymax = upper_y, width=0, colour = Plot)) +
    geom_errorbarh(aes(y = Green_mean, xmin = lower_x, xmax = upper_x, colour = Plot)) +
    geom_smooth(aes(y = Green_mean, x = Soilmoist_mean), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    labs(y = 'Green Tea Mass Loss %', x = "Soil Moisture %") +
    xlim(15, 50) +
    ylim(30, 43) +
    theme_classic())

ggsave(Soilmoist_mean_fig, file = "figures/Soilmoist_mean.png", height = 3, width = 4)


lower_x <- tea_mean$ALT_mean - tea_mean$ALT_SE
upper_x <- tea_mean$ALT_mean + tea_mean$ALT_SE

(ALT_mean_fig <- ggplot(data = tea_mean) +
    geom_point(aes(y = Green_mean, x = ALT_mean, colour = Plot)) +
    geom_errorbar(aes(y = Green_mean, x = ALT_mean, ymin = lower_y, ymax = upper_y, width=0, colour = Plot)) +
    geom_errorbarh(aes(y = Green_mean, xmin = lower_x, xmax = upper_x, colour = Plot)) +
    geom_smooth(aes(y = Green_mean, x = ALT_mean), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    labs(y = 'Green Tea Mass Loss %', x = "Active Layer Thickness (cm)") +
    xlim(0, 65) +
    ylim(30, 43) +
    theme_classic())

ggsave(ALT_mean_fig, file = "figures/ALT_mean.png", height = 3, width = 4)

lower_y <- tea_mean$meantemp_mean - tea_mean$meantemp_SE
upper_y <- tea_mean$meantemp_mean + tea_mean$meantemp_SE
lower_x <- tea_mean$meantemp_mean - tea_mean$meantemp_SE
upper_x <- tea_mean$meantemp_mean + tea_mean$meantemp_SE

(meantemp_mean_fig <- ggplot(data = tea_mean) +
    geom_point(aes(y = Green_mean, x = meantemp_mean, colour = Plot)) +
    geom_errorbar(aes(y = Green_mean, x = meantemp_mean, ymin = lower_y, ymax = upper_y, width=0, colour = Plot)) +
    geom_errorbarh(aes(y = Green_mean, xmin = lower_x, xmax = upper_x, colour = Plot)) +
    geom_smooth(aes(y = Green_mean, x = meantemp_mean), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    labs(y = 'Green Tea Mass Loss %', x = "Mean Surface Temperature (\u00B0C)") +
    xlim(6.5, 7.25) +
    ylim(30, 43) +
    theme_classic())

ggsave(meantemp_mean_fig, file = "figures/meantemp_mean.png", height = 3, width = 4)

lower_y <- tea_mean$thermsum_mean - tea_mean$thermsum_SE
upper_y <- tea_mean$thermsum_mean + tea_mean$thermsum_SE
lower_x <- tea_mean$thermsum_mean - tea_mean$thermsum_SE
upper_x <- tea_mean$thermsum_mean + tea_mean$thermsum_SE

(thermsum_mean_fig <- ggplot(data = tea_mean) +
    geom_point(aes(y = Green_mean, x = thermsum_mean, colour = Plot)) +
    geom_errorbar(aes(y = Green_mean, x = thermsum_mean, ymin = lower_y, ymax = upper_y, width=0, colour = Plot)) +
    geom_errorbarh(aes(y = Green_mean, xmin = lower_x, xmax = upper_x, colour = Plot)) +
    geom_smooth(aes(y = Green_mean, x = thermsum_mean), method = "lm", colour = "black") +
    scale_color_viridis_d(option = "plasma") +
    labs(y = 'Green Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)") +
    ylim(30, 43) +
    theme_classic())

ggsave(thermsum_mean_fig, file = "figures/thermtemp_mean.png", height = 3, width = 4)



#### Bayes lots for panel ####

predictions_alt <- ggpredict(alttherm_green_alt, terms = c("ALT_scaled"))

(poster_alt <- ggplot() +
    geom_line(data = predictions_alt, aes(y = predicted*37.38227087, x = x*42.66834713), col = "#FF7F00",
              size = 1) +
    geom_ribbon(data = predictions_alt, aes(ymin = conf.low*37.38227087, ymax = conf.high*37.38227087, 
                                            x = x*42.66834713), fill = "#FF7F00",  alpha = 0.4) +
    geom_point(data = tea_full, aes(x = ALT,  y = Green),col = "#FF7F00", 
               alpha = 0.8) +
    xlim(0, 65) +
    theme_classic() +
    labs(y = 'Green Tea Mass Loss %', x = "Active Layer Thickness (cm)"))

(poster_alt <- poster_alt +  theme(legend.position = "none"))



predictions_temp <- ggpredict(alttherm_green_alt, terms = c("therm_scaled"))

(poster_temp <- ggplot() +
    geom_line(data = predictions_temp, aes(y = predicted*37.38227087, x = x*132.57754955), col = "#FF3E96",
              size = 1) +
    geom_ribbon(data = predictions_temp, aes(ymin = conf.low*37.38227087, ymax = conf.high*37.38227087, 
                                             x = x*132.57754955), fill = "#FF3E96",  alpha = 0.4) +
    geom_point(data = tea_full, aes(x = thermsum,  y = Green),col = "#FF3E96", 
               alpha = 0.8) +
    theme_classic() +
    labs(y = 'Green Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)"))
(poster_temp <- poster_temp +  theme(legend.position = "none"))

predictions_moi <- ggpredict(soiltherm_green_alt, terms = c("Soilmoist_scaled"))

(poster_moi <- ggplot() +
    geom_line(data = predictions_moi, aes(y = predicted*37.38227087, x = x*34.68700926), col = "#66CDAA",
              size = 1) +
    geom_ribbon(data = predictions_moi, aes(ymin = conf.low*37.38227087, ymax = conf.high*37.38227087, 
                                            x = x*34.68700926), fill = "#66CDAA",  alpha = 0.4) +
    geom_point(data = tea_full, aes(x = Soilmoist,  y = Green),col = "#66CDAA", 
               alpha = 0.8) +
    theme_classic() +   
    labs(y = 'Green Tea Mass Loss %', x = "Soil Moisture (%)"))
(poster_moi <- poster_moi +  theme(legend.position = "none"))



(poster_bayesplot <- ggarrange(poster_temp, poster_moi,poster_alt,nrow =1))

ggsave(poster_bayesplot, filename = "figures/manuscript/env_bayes_scatter.png", height = 3, width = 9)


(local_scaleplot <- ggarrange(thermsum_fig, Soilmoist_fig,ALT_fig,nrow =1))

ggsave(local_scaleplot, filename = "figures/manuscript/local_bayes_scatter.png", height = 3, width = 9)

#### Mega monster table - Full ####
# using code adapted from Mariana: https://github.com/ShrubHub/ShrubHub/blob/master/scripts/users/mgarciacriado/traits_vs_ranges/scripts/8_model_table.R

library(broom)
library(stargazer)

# Let's break this down since otherwise R crashes
# Load a few models, run the table, and then close R and clean the environment
# Rinse and repeat with the rest of models

# Jonathan Chang's function
p_summarize <- function(model) {
  brms::posterior_summary(model) %>% 
    as_tibble(rownames = "parameter")
}


# add model objects to list
models.list <- list(soiltherm_green_alt,soiltherm_red_alt,
                    soiltherm_rate_alt,soiltherm_stab_alt,
                    alttherm_green_alt,alttherm_red_alt,
                    alttherm_rate_alt,alttherm_stab_alt,
                    topo_green_m,topo_red_m,
                    topo_rate_m,topo_stab_m)

# compile model titles
model_names <- c(
  "Soil moisture & thermal sum vs Green Tea Mass Loss",
  "Soil moisture & thermal sum vs Rooibos Tea Mass Loss",
  "Soil moisture & thermal sum vs Decomposition Rate",
  "Soil moisture & thermal sum vs Stabilisation Factor",
  "Active Layer & thermal sum vs Green Tea Mass Loss",
  "Active Layer & thermal sum vs Rooibos Tea Mass Loss",
  "Active Layer & thermal sum vs Decomposition Rate",
  "Active Layer & thermal sum vs Stabilisation Factor",
  "Topography vs Green Tea Mass Loss",
  "Topography vs Rooibos Tea Mass Loss",
  "Topography vs Decomposition Rate",
  "Topography vs Stabilisation Factor")

# number the models 1 through 12
model_number <- 1:12 
# bind the model tables together
mod.df <- data.frame(model_number, model_names)

# Extract parameters
mod.table <- lapply(models.list, p_summarize) %>% 
  bind_rows(.id = "model_number") 

# Add model name to table
mod.table$model_number <- as.integer(mod.table$model_number)
mod.table.final <- left_join(mod.table, mod.df, by = "model_number")

# Clean model parameters
mod.table.final <- mod.table.final %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final$model_names[duplicated(mod.table.final$model_names)] <- "  "
mod.table.final$model_names <- as.character(mod.table.final$model_names)
mod.table.final$model_number[duplicated(mod.table.final$model_number)] <- "  "

colnames(mod.table.final) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final <- mod.table.final[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final <- mod.table.final %>% mutate_if(is.numeric, round, digits = 3)


# Save in csv
write.csv(mod.table.final, "models/scaled/non_int_table.csv")

# remove underscores which are fudging html table production
mod.table.final$Term <- sapply(mod.table.final$Term , function(x) gsub("_", "",  x))

# Convert to table
stargazer(mod.table.final, type = "html",  summary = FALSE, out = "models/scaled/non_int_table.html")

#### Mega monster table - Null ####

# Let's break this down since otherwise R crashes
# Load a few models, run the table, and then close R and clean the environment
# Rinse and repeat with the rest of models

# Jonathan Chang's function
p_summarize <- function(model) {
  brms::posterior_summary(model) %>% 
    as_tibble(rownames = "parameter")
}


# add model objects to list
models.list.null <- list(soiltherm_green_null,soiltherm_red_null,
                         soiltherm_rate_null,soiltherm_stab_null,
                         alttherm_green_null,alttherm_red_null,
                         alttherm_rate_null,alttherm_stab_null,
                         topo_green_null,topo_red_null,
                         topo_rate_null,topo_stab_null)

# compile model titles
model_names.null <- c(
  "Soil moisture & surface temperature vs Green Tea Mass Loss",
  "Soil moisture & surface temperature vs Rooibos Tea Mass Loss",
  "Soil moisture & surface temperature vs Decomposition Rate",
  "Soil moisture & surface temperature vs Stabilisation Factor",
  "Active Layer & surface temperature vs Green Tea Mass Loss",
  "Active Layer & surface temperature vs Rooibos Tea Mass Loss",
  "Active Layer & surface temperature vs Decomposition Rate",
  "Active Layer & surface temperature vs Stabilisation Factor",
  "Topography vs Green Tea Mass Loss",
  "Topography vs Rooibos Tea Mass Loss",
  "Topography vs Decomposition Rate",
  "Topography vs Stabilisation Factor")

# number the models 1 through 12
model_number_null <- 1:12 
# bind the model tables together
mod.df.null <- data.frame(model_number_null, model_names.null)

# Extract parameters
mod.table.null <- lapply(models.list.null, p_summarize) %>% 
  bind_rows(.id = "model_number_null") 

# Add model name to table
mod.table.null$model_number_null <- as.integer(mod.table.null$model_number_null)
mod.table.final.null <- left_join(mod.table.null, mod.df.null, by = "model_number_null")

# Clean model parameters
mod.table.final.null <- mod.table.final.null %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final.null$model_names[duplicated(mod.table.final.null$model_names_null)] <- "  "
mod.table.final.null$model_names <- as.character(mod.table.final.null$model_names_null)
mod.table.final.null$model_number_null[duplicated(mod.table.final.null$model_number_null)] <- "  "

colnames(mod.table.final.null) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final.null <- mod.table.final.null[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final.null <- mod.table.final.null %>% mutate_if(is.numeric, round, digits = 3)


# Save in csv
write.csv(mod.table.final.null, "models/sclaed/bayes_table_null.csv")

# remove underscores which are fudging html table production
mod.table.final.null$Term <- sapply(mod.table.final.null$Term , function(x) gsub("_", "",  x))

# Convert to table
stargazer(mod.table.final.null, type = "html",  summary = FALSE, out = "models/scaled/bayes_table_null.html")

#### Mega monster table - Thermal Sum ####
# using code adapted from Mariana: https://github.com/ShrubHub/ShrubHub/blob/master/scripts/users/mgarciacriado/traits_vs_ranges/scripts/8_model_table.R

library(broom)
library(stargazer)

# Let's break this down since otherwise R crashes
# Load a few models, run the table, and then close R and clean the environment
# Rinse and repeat with the rest of models

# Jonathan Chang's function
p_summarize <- function(model) {
  brms::posterior_summary(model) %>% 
    as_tibble(rownames = "parameter")
}


# add model objects to list
models.list <- list(soiltherm_green_m,soiltherm_red_m,
                    soiltherm_rate_m,soiltherm_stab_m,
                    alttherm_green_m,alttherm_red_m,
                    alttherm_rate_m,alttherm_stab_m,
)

# compile model titles
model_names <- c(
  "Soil moisture & thermal sum vs Green Tea Mass Loss",
  "Soil moisture & thermal sum vs Rooibos Tea Mass Loss",
  "Soil moisture & thermal sum vs Decomposition Rate",
  "Soil moisture & thermal sum vs Stabilisation Factor",
  "Active Layer & thermal sum vs Green Tea Mass Loss",
  "Active Layer & thermal sum vs Rooibos Tea Mass Loss",
  "Active Layer & thermal sum vs Decomposition Rate",
  "Active Layer & thermal sum vs Stabilisation Factor")

# number the models 1 through 8
model_number <- 1:8 
# bind the model tables together
mod.df <- data.frame(model_number, model_names)

# Extract parameters
mod.table <- lapply(models.list, p_summarize) %>% 
  bind_rows(.id = "model_number") 

# Add model name to table
mod.table$model_number <- as.integer(mod.table$model_number)
mod.table.final <- left_join(mod.table, mod.df, by = "model_number")

# Clean model parameters
mod.table.final <- mod.table.final %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final$model_names[duplicated(mod.table.final$model_names)] <- "  "
mod.table.final$model_names <- as.character(mod.table.final$model_names)
mod.table.final$model_number[duplicated(mod.table.final$model_number)] <- "  "

colnames(mod.table.final) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final <- mod.table.final[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final <- mod.table.final %>% mutate_if(is.numeric, round, digits = 3)


# Save in csv
write.csv(mod.table.final, "models/scaled/thermalsum_bayes_table.csv")

# remove underscores which are fudging html table production
mod.table.final$Term <- sapply(mod.table.final$Term , function(x) gsub("_", "",  x))

# Convert to table
stargazer(mod.table.final, type = "html",  summary = FALSE, out = "models/scaled/thermal sum_bayes_table.html")






#### Soil Plot Panel - Thermal Sum Version - KEEP ####
#Put data in data frame 
predictions_altgreen <- ggpredict(alttherm_green_m, terms = c("therm_scaled", "ALT_scaled"))
predictions_altgreen$depth <- factor(predictions_altgreen$group, levels = c("0.56", "0.93", "1.29"),
                                     labels = c("Shallow", "Medium", "Deep"))

# plot alt and temp interaction predictions
(alttemp_green_plot <- ggplot() +
    geom_line(data = predictions_altgreen, aes(y = predicted*37.38227087, x = x*132.57754955
, colour = depth),
              size = 2) +
    geom_ribbon(data = predictions_altgreen, aes(ymin = conf.low*37.38227087, ymax = conf.high*37.38227087, 
                                                 x = x*132.57754955, fill = depth), alpha = 0.2) +
    geom_point(data = tea_full, aes(x = therm_scaled*132.57754955,  y = Green_scaled*37.38227087, shape = Plot), 
               alpha = 0.7) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#F2F200","#E38D1C","#D41613")) +
    scale_colour_manual(values = c("#F2F200","#E38D1C","#D41613")) +
    theme_classic() +
    ggtitle("a)") +
    labs(y = 'Green Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)"))
(alttemp_green_plot <- alttemp_green_plot +  theme(legend.position = "none"))

#Put data in data frame - soil edition
predictions_soilgreen <- ggpredict(soiltherm_green_m, terms = c("therm_scaled", "Soilmoist_scaled"))
predictions_soilgreen$moist <- factor(predictions_soilgreen$group, levels = c("0.53", "0.91", "1.3"),
                                      labels = c("Dry", "Medium", "Moist"))


# plot alt and temp interaction predictions
(soiltemp_green_plot <- ggplot() +
    geom_line(data = predictions_soilgreen, aes(y = predicted*37.38227087, x = x*132.57754955, colour = moist),
              size = 2) +
    geom_ribbon(data = predictions_soilgreen, aes(ymin = conf.low*37.38227087, ymax = conf.high*37.38227087, 
                                                  x = x*132.57754955, fill = moist),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = thermsum,  y = Green, shape = Plot), 
               alpha = 0.7) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    scale_colour_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    theme_classic() +
    ggtitle("b)") +
    labs(y = 'Green Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)"))
(soiltemp_green_plot <- soiltemp_green_plot  + theme(legend.position = "none"))

#### rate ####
#Put data in data frame 
predictions_altrate <- ggpredict(alttherm_rate_m, terms = c("therm_scaled", "ALT_scaled"))
predictions_altrate$depth <- factor(predictions_altrate$group, levels = c("0.55", "0.92", "1.29"),
                                    labels = c("Shallow", "Medium", "Deep"))

# plot alt and temp interaction predictions
(alttemp_rate_plot <- ggplot() +
    geom_line(data = predictions_altrate, aes(y = predicted*0.03314069, x = x*132.57754955, colour = depth),
              size = 2) +
    geom_ribbon(data = predictions_altrate, aes(ymin = conf.low*0.03314069, ymax = conf.high*0.03314069, 
                                                x = x*132.57754955, fill = depth),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = therm_scaled*132.57754955,  y = rate_scaled*0.03314069, shape = Plot), 
               alpha = 0.7) +
    # ylim(0,0.085) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#EAF041","#E38D1C","#D41613")) +
    scale_colour_manual(values = c("#EAF041","#E38D1C","#D41613")) +
    theme_classic() +
    ggtitle("c)") +
    labs(y = 'Decomposition Rate (k)', x = "28 Day Thermal Sum (\u00B0C)",color = "Active Layer Depth"))

(alttemp_rate_plot <- alttemp_rate_plot + theme(legend.position = "none"))

#Put data in data frame - soil edition
predictions_soilrate <- ggpredict(soiltherm_rate_m, terms = c("therm_scaled", "Soilmoist_scaled"))
predictions_soilrate$moist <- factor(predictions_soilrate$group, levels = c("0.54", "0.92", "1.31"),
                                     labels = c("Dry", "Medium", "Moist"))


# plot alt and temp interaction predictions
(soiltemp_rate_plot <- ggplot() +
    geom_line(data = predictions_soilrate, aes(y = predicted*0.03314069, x = x*132.57754955, colour = moist),
              size = 2) +
    geom_ribbon(data = predictions_soilrate, aes(ymin = conf.low*0.03314069, ymax = conf.high*0.03314069, 
                                                 x = x*132.57754955, fill = moist),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = therm_scaled*132.57754955,  y = rate_scaled*0.03314069, shape = Plot), 
               alpha = 0.7) +
    #ylim(0,0.085) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    scale_colour_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    theme_classic() +
    ggtitle("d)") +
    labs(y = 'Decomposition Rate (k)', x = "28 Day Thermal Sum (\u00B0C)",color = "Moisture Index"))
(soiltemp_rate_plot <- soiltemp_rate_plot + theme(legend.position = "none"))



#### Rooibos ####
#Put data in data frame 
predictions_altred <- ggpredict(alttherm_red_m, terms = c("therm_scaled", "ALT_scaled"))
predictions_altred$depth <- factor(predictions_altred$group, levels = c("0.55", "0.92", "1.29"),
                                   labels = c("Shallow", "Medium", "Deep"))

# plot alt and temp interaction predictions
(alttemp_red_plot <- ggplot() +
    geom_line(data = predictions_altred, aes(y = predicted*14.20226993, x = x*132.57754955, colour = depth),
              size = 2) +
    geom_ribbon(data = predictions_altred, aes(ymin = conf.low*14.20226993, ymax = conf.high*14.20226993, 
                                               x = x*132.57754955, fill = depth),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = thermsum,  y = Rooibos, shape = Plot), 
               alpha = 0.7) +
    ylim(0,22) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#EAF041","#E38D1C","#D41613")) +
    scale_colour_manual(values = c("#EAF041","#E38D1C","#D41613")) +
    theme_classic() +
    ggtitle("e)") +
    labs(y = 'Rooibos Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)",color = "Active Layer Depth"))

(alttemp_red_plot <- alttemp_red_plot + theme(legend.position = "none"))


#Put data in data frame - soil edition
predictions_soilred <- ggpredict(soiltherm_red_m, terms = c("therm_scaled", "Soilmoist_scaled"))
predictions_soilred$moist <- factor(predictions_soilred$group, levels = c("0.53", "0.92", "1.3"),
                                    labels = c("Dry", "Medium", "Moist"))


# plot alt and temp interaction predictions
(soiltemp_red_plot <- ggplot() +
    geom_line(data = predictions_soilred, aes(y = predicted*14.20226993, x = x*132.57754955, colour = moist),
              size = 2) +
    geom_ribbon(data = predictions_soilred, aes(ymin = conf.low*14.20226993, ymax = conf.high*14.20226993, 
                                                x = x*132.57754955, fill = moist),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = thermsum,  y = Rooibos, shape = Plot), 
               alpha = 0.7) +
    ylim(0,22) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    scale_colour_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    theme_classic() +
    ggtitle("f)") +
    labs(y = 'Rooibos Tea Mass Loss %', x = "28 Day Thermal Sum (\u00B0C)",color = "Moisture Index"))
(soiltemp_red_plot <- soiltemp_red_plot + theme(legend.position = "none"))


#### stab ####
#Put data in data frame 
predictions_altstab <- ggpredict(alttherm_stab_m, terms = c("therm_scaled", "ALT_scaled"))
predictions_altstab$depth <- factor(predictions_altstab$group, levels = c("0.56", "0.93", "1.29"),
                                    labels = c("Shallow", "Medium", "Deep"))

# plot alt and temp interaction predictions
(alttemp_stab_plot <- ggplot() +
    geom_line(data = predictions_altstab, aes(y = predicted*0.56594405, x = x*132.57754955, colour = depth),
              size = 2) +
    geom_ribbon(data = predictions_altstab, aes(ymin = conf.low*0.56594405, ymax = conf.high*0.56594405, 
                                                x = x*132.57754955, fill = depth),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = thermsum,  y = stab, shape = Plot), 
               alpha = 0.7) +
    ylim(0,0.8) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#EAF041","#E38D1C","#D41613")) +
    scale_colour_manual(values = c("#EAF041","#E38D1C","#D41613")) +
    theme_classic() +
    ggtitle("g)") +
    labs(y = 'Stabilisation Factor (S)', x = "28 Day Thermal Sum (\u00B0C)",color = "Active Layer Depth"))

(alttemp_stab_plot <- alttemp_stab_plot + guides(fill=guide_legend(title="Active Layer Thickness"),
                                                 color=guide_legend(title="Active Layer Thickness"))
  + theme(legend.position = c(.8,.2)) + guides(shape=FALSE))


#Put data in data frame - soil edition
predictions_soilstab <- ggpredict(soiltherm_stab_m, terms = c("therm_scaled", "Soilmoist_scaled"))
predictions_soilstab$moist <- factor(predictions_soilstab$group, levels = c("0.53", "0.91", "1.3"),
                                     labels = c("Dry", "Medium", "Moist"))


# plot alt and temp interaction predictions
(soiltemp_stab_plot <- ggplot() +
    geom_line(data = predictions_soilstab, aes(y = predicted*0.56594405, x = x*132.57754955, colour = moist),
              size = 2) +
    geom_ribbon(data = predictions_soilstab, aes(ymin = conf.low*0.56594405, ymax = conf.high*0.56594405, 
                                                 x = x*132.57754955, fill = moist),  alpha = 0.2) +
    geom_point(data = tea_full, aes(x = thermsum,  y = stab, shape = Plot), 
               alpha = 0.7) +
    ylim(0,0.8) +
    scale_shape_manual(values=seq(0,8)) +
    scale_fill_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    scale_colour_manual(values = c("#FAD60DCC","#72CC8A","#5E51AD")) +
    theme_classic() +
    ggtitle("h)") +
    labs(y = 'Stabilisation Factor (S)', x = "28 Day Thermal Sum (\u00B0C)",color = "Moisture Index"))

(soiltemp_stab_plot <- soiltemp_stab_plot + guides(fill=guide_legend(title="Moisture Index"),
                                                   color=guide_legend(title="Moisture Index")) + theme(legend.position = c(.8,.2)) 
  +guides(shape=FALSE))





#### merge plot together ####

(soil_modelplots <- ggarrange(soiltemp_green_plot, 
                              soiltemp_rate_plot, 
                              soiltemp_red_plot,
                              soiltemp_stab_plot,nrow=4))

(alt_modelplots <- ggarrange(alttemp_green_plot,    
                             alttemp_rate_plot,
                             alttemp_red_plot,
                             alttemp_stab_plot,nrow=4))



(belowground_bayesplot<- ggarrange(alt_modelplots,soil_modelplots,  nrow=1,common.legend = FALSE))


ggsave(belowground_bayesplot, filename = "figures/manuscript/thermsum_bayes_scatter.png", height = 16, width = 10 )



