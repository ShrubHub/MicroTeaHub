### 3. Semivariograms ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Date created: 3rd Feb 2021 ###

#### LOAD PACKAGES ####

library(tidyverse)
library(esquisse)
library(gstat)
library(readxl)
library(sp)
library(gridExtra)
library(grid)
library(openair)
library(cowplot)
library(ggpubr)

#### SET UP WORKSPACE ####
teapoint_full <- read.csv(file="data/fullpropertieslong.csv") # read in soil data
colnames(teapoint_full)

#### THERMAL SUM SEMIVAR ####
plot(teapoint_full$X, teapoint_full$Y)

semivar_therm <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$thermsum) 

# Remove NA rows
semivar_therm <- na.omit(object = semivar_therm)

semivar_therm <- variogram(rate~1, data = semivar_therm, locations=~x+y)

plot(semivar_therm)


semivar_fit_therm <- fit.variogram(semivar_therm,
                                 model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_therm, semivar_fit_therm)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_therm, maxdist = max(semivar_therm$dist))

(thermsemi <- ggplot() +
    geom_point(data = semivar_therm, aes(x = dist, y = gamma),size = 1.6, col = "#dd1c77") +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma),size = 1.6, col = "#dd1c77") +
    theme_classic(base_size = 15) +
    labs(x = "Distance (m)", y = "Semivariance (Y)") +
    ggtitle("Thermal Sum"))

#### ACTIVE LAYER SEMIVAR ####
plot(teapoint_full$X, teapoint_full$Y)

semivar_ALT <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                         rate = teapoint_full$ALT) 

# Remove NA rows
semivar_ALT <- na.omit(object = semivar_ALT)

semivar_ALT <- variogram(rate~1, data = semivar_ALT, locations=~x+y)

plot(semivar_ALT)


semivar_fit_ALT <- fit.variogram(semivar_ALT,
                                 model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_ALT, semivar_fit_ALT)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_ALT, maxdist = max(semivar_ALT$dist))

(ALTsemi <- ggplot() +
  geom_point(data = semivar_ALT, aes(x = dist, y = gamma),size = 1.6, col = "#d95f0e") +
  geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma),size = 1.6, col = "#d95f0e") +
  theme_classic(base_size = 15) +
  labs(x = "Distance (m)", y = "Semivariance (Y)") +
  ggtitle("Active Layer Thickness"))

#### RATE  SEMIVAR ####
semivar_RATE <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$rate) 

# Remove NA rows
semivar_RATE <- na.omit(object = semivar_RATE)

semivar_RATE <- variogram(rate~1, data = semivar_RATE, locations=~x+y)

plot(semivar_RATE)

semivar_fit_RATE <- fit.variogram(semivar_RATE,
                             model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_RATE, semivar_fit_RATE)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_RATE, maxdist = max(semivar_RATE$dist))

(RATEsemi <- ggplot() +
    geom_point(data = semivar_RATE, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic(base_size = 15) +
    labs(x = "Distance (m)", y = "Semivariance (k)") +
    ggtitle("Decomposition Rate"))


#### MOISTURE  SEMIVAR ####
semivar_MOI <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$Soilmoist) 

# Remove NA rows
semivar_MOI <- na.omit(object = semivar_MOI)

semivar_MOI <- variogram(rate~1, data = semivar_MOI, locations=~x+y)

plot(semivar_MOI)

semivar_fit_MOI <- fit.variogram(semivar_MOI,
                                  model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_MOI, semivar_fit_MOI)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_MOI, maxdist = max(semivar_MOI$dist))

(MOIsemi <- ggplot() +
    geom_point(data = semivar_MOI, aes(x = dist, y = gamma), size = 1.6, col = "#1c9099") +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma), size = 1.6, col = "#1c9099") +
    theme_classic(base_size = 15) +
    labs(x = "Distance (m)", y = "Semivariance (Y)") +
    ggtitle("Soil Moisture"))


#### GREEN TEA MASS LOSS  SEMIVAR ####
semivar_GREEN <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$Green) 

# Remove NA rows
semivar_GREEN <- na.omit(object = semivar_GREEN)
# Remove -Inf rows
semivar_GREEN <- semivar_GREEN[!is.infinite(rowSums(semivar_GREEN)),]

semivar_GREEN <- variogram(rate~1, data = semivar_GREEN, locations=~x+y)

plot(semivar_GREEN)

semivar_fit_GREEN <- fit.variogram(semivar_GREEN,
                                 model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_GREEN, semivar_fit_GREEN)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_GREEN, maxdist = max(semivar_GREEN$dist))

(GREENsemi <- ggplot() +
    geom_point(data = semivar_GREEN, aes(x = dist, y = gamma),size = 1.6, colour = "#4de1a1") +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma),size = 1.6, colour = "#4de1a1") +
    theme_classic(base_size = 15) +
    labs(x = "Distance (m)", y = "Semivariance (Y)") +
    ggtitle("Green Tea % Mass Loss"))

#### RED TEA MASS LOSS  SEMIVAR ####
semivar_RED <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                            rate = teapoint_full$Rooibos) 

# Remove NA rows
semivar_RED <- na.omit(object = semivar_RED)

semivar_RED <- variogram(rate~1, data = semivar_RED, locations=~x+y)

plot(semivar_RED)

semivar_fit_RED <- fit.variogram(semivar_RED,
                                   model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_RED, semivar_fit_RED)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_RED, maxdist = max(semivar_RED$dist))

(RED_semi <- ggplot() +
    geom_point(data = semivar_RED, aes(x = dist, y = gamma),size = 1.6, colour = "#F2275D") +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma),size = 1.6, colour = "#F2275D") +
    theme_classic(base_size = 15) +
    labs(x = "Distance (m)", y = "Semivariance (Y)") +
    ggtitle("Rooibos Tea % Mass Loss"))

#### ELEVATION  SEMIVAR ####
semivar_ELEV <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$elevation) 

# Remove NA rows
semivar_ELEV <- na.omit(object = semivar_ELEV)

semivar_ELEV <- variogram(rate~1, data = semivar_ELEV, locations=~x+y)

plot(semivar_ELEV)

semivar_fit_ELEV <- fit.variogram(semivar_ELEV,
                                 model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_ELEV, semivar_fit_ELEV)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_ELEV, maxdist = max(semivar_ELEV$dist))

(ELEVsemi <- ggplot() +
    geom_point(data = semivar_ELEV, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic(base_size = 15) +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("DSM Elevation"))

#### Stablisation  SEMIVAR ####
semivar_stab <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                           rate = teapoint_full$stab) 

# Remove NA rows
semivar_stab <- na.omit(object = semivar_stab)

semivar_stab <- variogram(rate~1, data = semivar_stab, locations=~x+y)

plot(semivar_stab)

semivar_fit_stab <- fit.variogram(semivar_stab,
                                  model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_stab, semivar_fit_stab)

semivar_fit_fort <- variogramLine(semivar_fit_stab, maxdist = max(semivar_stab$dist))

(STABsemi <- ggplot() +
    geom_point(data = semivar_stab, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic(base_size = 18) +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("Stabilisation Rate"))

#### Slope  SEMIVAR ####
semivar_SLO <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                           rate = teapoint_full$slope) 

# Remove NA rows
semivar_SLO <- na.omit(object = semivar_SLO)

semivar_SLO <- variogram(rate~1, data = semivar_SLO, locations=~x+y)

plot(semivar_SLO)

semivar_fit_SLO <- fit.variogram(semivar_SLO,
                                  model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_SLO, semivar_fit_SLO)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_SLO, maxdist = max(semivar_SLO$dist))

(SLOsemi <- ggplot() +
    geom_point(data = semivar_SLO, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic(base_size = 18) +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("Slope"))

#### Aspect  SEMIVAR ####
semivar_ASP <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$aspect) 

# Remove NA rows
semivar_ASP <- na.omit(object = semivar_ASP)

semivar_ASP <- variogram(rate~1, data = semivar_ASP, locations=~x+y)

plot(semivar_ASP)

semivar_fit_ASP <- fit.variogram(semivar_ASP,
                                 model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_ASP, semivar_fit_ASP)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_ASP, maxdist = max(semivar_ASP$dist))

(ASPsemi <- ggplot() +
    geom_point(data = semivar_ASP, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic(base_size = 18) +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("Aspect"))



#### Hillshade SEMIVAR ####
semivar_HILL <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                           rate = teapoint_full$Hillshade) 

# Remove NA rows
semivar_HILL <- na.omit(object = semivar_HILL)

semivar_HILL <- variogram(rate~1, data = semivar_HILL, locations=~x+y)

plot(semivar_HILL)

semivar_fit_HILL <- fit.variogram(semivar_HILL,
                                  model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_HILL, semivar_fit_HILL)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_HILL, maxdist = max(semivar_HILL$dist))

(HILLsemi <- ggplot() +
    geom_point(data = semivar_HILL, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic() +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("Hillshade"))

#### Meantemp SEMIVAR ####
semivar_MEAN <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                           rate = teapoint_full$meantemp) 

# Remove NA rows
semivar_MEAN <- na.omit(object = semivar_MEAN)

semivar_MEAN <- variogram(rate~1, data = semivar_MEAN, locations=~x+y)

plot(semivar_MEAN)

semivar_fit_MEAN <- fit.variogram(semivar_MEAN,
                                  model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_MEAN, semivar_fit_MEAN)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_MEAN, maxdist = max(semivar_MEAN$dist))

(MEANsemi <- ggplot() +
    geom_point(data = semivar_MEAN, aes(x = dist, y = gamma),size = 1.6, colour = "#dd1c77") +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma),size = 1.6, colour = "#dd1c77") +
    theme_classic() +
    labs(x = "Distance (m)", y = "Semivariance (Y)") +
    ggtitle("Mean Temperature"))

#### Mintemp SEMIVAR ####
semivar_MIN <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                           rate = teapoint_full$mintemp) 

# Remove NA rows
semivar_MIN  <- na.omit(object = semivar_MIN )

semivar_MIN  <- variogram(rate~1, data = semivar_MIN , locations=~x+y)

plot(semivar_MIN )

semivar_fit_MIN <- fit.variogram(semivar_MIN ,
                                  model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_MIN , semivar_fit_MIN)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_MIN, maxdist = max(semivar_MIN$dist))

(MINsemi <- ggplot() +
    geom_point(data = semivar_MIN, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic() +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("Min Temperature"))

#### Maxtemp SEMIVAR ####
semivar_MAX <- data.frame(x = teapoint_full$X, y = teapoint_full$Y, 
                          rate = teapoint_full$maxtemp) 

# Remove NA rows
semivar_MAX  <- na.omit(object = semivar_MAX )

semivar_MAX   <- variogram(rate~1, data = semivar_MAX, locations=~x+y)

plot(semivar_MAX )

semivar_fit_MAX <- fit.variogram(semivar_MAX ,
                                 model=vgm(c("Exp","Sph", "Gau", "Mat")))

plot(semivar_MAX, semivar_fit_MAX)

#grab values for semivar and semivarfit for ggplot
semivar_fit_fort <- variogramLine(semivar_fit_MAX, maxdist = max(semivar_MAX$dist))

(MAXsemi <- ggplot() +
    geom_point(data = semivar_MAX, aes(x = dist, y = gamma)) +
    geom_line(data = semivar_fit_fort, aes(x = dist, y = gamma)) +
    theme_classic() +
    labs(x = "Distance (m)", y = "Semivariance (\u03B3)") +
    ggtitle("Max Temperature"))

#### GROUP & PANEL SEMIVARIOGRAMS ####
topography_panel <- grid.arrange(ELEVsemi,SLOsemi,ASPsemi,HILLsemi,ncol=1)
ggsave(topography_panel, filename = "figures/semivariograms_topography.png", height = 8, width = 5) 

soilproperties_panel <- grid.arrange(ALTsemi,MOIsemi,ncol=1)
ggsave(soilproperties_panel, filename = "figures/semivariograms_soilproperties.png", height = 8, width = 5)

decomposition_panel <- grid.arrange(RATEsemi,GREENsemi,REDsemi,STABsemi,ncol=1)
ggsave(decomposition_panel, filename = "figures/semivariograms_decomposition.png", height = 8, width = 5)

climate_panel <- grid.arrange(MEANsemi,MINsemi,MAXsemi,ncol=1)
ggsave(climate_panel, filename = "figures/semivariograms_climate.png", height = 8, width = 5)



#### MAP PANEL ####
##### plots with maps? ####
(green_map <- (ggplot(teapoint_full) +
                 aes(x = X, y = Y, colour = Green) +
                 geom_point(size = 2.62) +
                 scale_color_viridis_c(option = "viridis") +
                 labs(color = "Green Tea Mass Loss %") +
                 labs(x = "Latitude", y = "Longitude") +
                 ggtitle("d)") +
                 theme_classic(base_size = 13)))
(green_map <- green_map+ theme(legend.position = c(.8,.2)))

(green_plot <- green_map + annotation_custom(ggplotGrob(GREENsemi), xmin = 581900, xmax = 582500, 
                       ymin = 7719800, ymax = 7720100))


(red_map <- (ggplot(teapoint_full) +
                 aes(x = X, y = Y, colour = Rooibos) +
                 geom_point(size = 2.62) +
               scale_color_distiller(palette = "YlOrRd") +                 
               labs(color = "Rooibos Tea Mass Loss %") +
                 labs(x = "Latitude", y = "Longitude") +
                 ggtitle("a)") +
                 theme_classic(base_size = 13)))
(red_map <- red_map+ theme(legend.position = c(.8,.2)))

(red_plot <- red_map + annotation_custom(ggplotGrob(RED_semi), xmin = 581900, xmax = 582500, 
                                             ymin = 7719800, ymax = 7720100))



(soilmoist_map <- (ggplot(teapoint_full) +
                     aes(x = X, y = Y, colour = Soilmoist) +
                     geom_point(size = 2.62) +
                     scale_color_viridis_c(option = "cividis") +
                     labs(color = "Soil Moisture %") +
                     labs(x = "Latitude", y = "Longitude") +
                     ggtitle("c)") +
                     theme_classic(base_size = 13)))
(soilmoist_map <- soilmoist_map + theme(legend.position = c(.8,.2)))
(soil_plot <- soilmoist_map + annotation_custom(ggplotGrob(MOIsemi), xmin = 581900, xmax = 582500, 
                                             ymin = 7719800, ymax = 7720100))


(thermsum_map <- (ggplot(teapoint_full) +
                    aes(x = X, y = Y, colour = thermsum) +
                    geom_point(size = 2.62) +
                    scale_color_viridis_c(option = "magma") +
                    labs(color = "Thermal Sum (\u00B0C)") +
                    labs(x = "Latitude", y = "Longitude") +
                    ggtitle("e)") +
                    theme_classic(base_size = 13)))
(thermsum_map <- thermsum_map + theme(legend.position = c(.8,.2)))
(temp_plot <- thermsum_map + annotation_custom(ggplotGrob(thermsemi), xmin = 581900, xmax = 582500, 
                                             ymin = 7719800, ymax = 7720100))


(ALT_map <- (ggplot(teapoint_full) +
               aes(x = X, y = Y, colour = ALT) +
               geom_point(size = 2.62) +
               scale_color_viridis_c(option = "inferno") +
               labs(color = "Active Layer Thickness (cm)") +
               labs(x = "Latitude", y = "Longitude") +
               ggtitle("b)") +
               theme_classic(base_size = 13)))
(ALT_map <- ALT_map + theme(legend.position = c(.8,.2)))

(ALT_plot <- ALT_map + annotation_custom(ggplotGrob(ALTsemi), xmin = 581900, xmax = 582500, 
                                             ymin = 7719800, ymax = 7720100))
 
semi_map_panel <- ggarrange(red_plot,  ALT_plot, soil_plot,green_plot, temp_plot,
          ncol = 2, nrow = 3)
semi_map_panel

ggsave(semi_map_panel, filename = "figures/manuscript/semivariograms_maps.png", height = 18, width = 16)


