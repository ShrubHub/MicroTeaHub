### 6. Model interpolation maps ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Date: 2nd March 2021 ###
# Interpolate between points

#### LOAD PACKAGES ####
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(tidyverse)
library(esquisse)
library(ggridges)
library(gridExtra)
library(viridis)
library(spatstat)  
library(maptools)  
library(dynatopmodel)
library(ggeffects)


#### Load full tea data ####
tea_full <- read.csv(file="data/fullpropertieslong.csv") # read in soil data

# turn data into shapefile
# remove NA rows and prepare spatial points layer
tea_fullXY <- tea_full[!is.na(tea_full$X), ]
tea_fullXY_sp <- SpatialPoints(cbind(tea_fullXY$X, tea_fullXY$Y), 
                               proj4string=CRS("+init=epsg:32607"))

# get landscape data
dem <- raster('data/spatial/microclima/landscape_DSM_crop.tif')

# crop to size of microclimate model
e <- extent(581951.1, 583164.2, 7719428, 7720069)
dem_crop <- crop(dem, e, snap = "out")
plot(dem_crop) 

#### Make contour map ####
contour <- rasterToContour(dem_crop)
class(contour)
plot(contour)
plot(dem_crop, add=TRUE)

#### Topographic Wetness Index ####
# resample raster to ensure x and y res has same res
r <- raster(ext=extent(dem_crop), res=1) # 1 metre res
dem <- raster::resample(dem_crop, r, method = "bilinear")

# calculate wetness index using dynatopmodel package
twi <- upslope.area(dem_crop, log = TRUE, atb = TRUE, deg = 0.1, fill.sinks = FALSE)
sp::plot(twi$atb, main=c("Topographic Wetness Index log (m^2/m)"))
hist(twi$atb)

#save mean twi as rasters
writeRaster(twi$atb, 'data/spatial/twi.tif', format = 'GTiff')

# code to load data
twi <- raster('data/spatial/twi.tif')

# check to see how TWI layer aligns with soil moisture data points
## extract the raster values from these locations   
cord.UTM_mean2 <- spTransform((tea_fullXY_sp), crs(twi))

wetness_pred_point <- raster::extract(twi, cord.UTM_mean2, method='simple',df=TRUE)

## need to write CSV file to join to the master datasheet with all varibles 
combinePointValue=cbind(tea_fullXY, wetness_pred_point)

# rename new column
colnames(combinePointValue)[26] <- "wetness_pred_point"

# compare predicted values with real green tea mass loss values
ggplot(combinePointValue) +
  aes(x = Soilmoist, y = wetness_pred_point, colour = Plot) +
  geom_point(size = 3L) +
  scale_color_hue() +
  labs(x = "Real Soil Moisture Readings  % ", y = "Topographic Wetness Index log (m^2/m) ") +
  theme_classic()

# looks like lots of variaion in soil moist within site, but again each site is defined by soil moist characteristics
# confirm with twi boxplot
combinePointValue %>%
  filter(!(Plot %in% "WET CONTROL")) %>%
  ggplot() +
  aes(x = Plot, y = wetness_pred_point, fill = Plot, colour = Plot) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "inferno") +
  scale_color_viridis_d(option = "inferno") +
  theme_minimal()

#### Floodplain classifications #####
plot(twi, main="Topographic Wetness Index log (m^2/m)")
points(tea_full$X,tea_full$Y, pch=10, cex=.1, col = 1)

# calculate floodplain (local low) vs non-floodplain (local high) as rough soil moist binary
# high wetness index = > 8
# low wetness index = < 8
hist(twi)

# divide TWI into floodplain vs non-floodplain zones
# create classification matrix
reclass_df <- c(0, 2, NA,
                2, 8, 1,  # 0-9 = class 1 (dry) 
                8, Inf, 3) # 9+ = class 2 (wet) # 3, as floodplain sites have are roughtly 3 units higher than non-floodplain
reclass_df

# reshape the object into a matrix with columns and rows
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)

reclass_m

# reclassify the raster using the reclass object - reclass_m
twi_classified <- reclassify(twi,
                             reclass_m)

# view reclassified data
barplot(twi_classified,
        main = "Number of pixels in each class")
# assign all pixels that equal 0 to NA or no data value
twi_classified[twi_classified == 0] <- NA

# plot reclassified data
plot(twi_classified)
plot(twi)


#### Green tea mass loss prediction map ####
# get thermsum tif data
thermsum <- raster('data/spatial/microclima/thermsum.tif')

# crop twi to meantemp extent
e <- extent(581951.1, 583164.2, 7719428, 7720069)
twi_crop <- crop(twi_classified, e, snap = "out")

# load model results
load("models/soiltherm_green_mm.Rdata")
summary(soiltherm_green_m)
soiltemp_green_m
# mean temp effect size on green tea mass loss is -0.10
# intercept iss 45.17 

# resample meantemp model to be the same resolution as TWI
temp_crop <- resample(thermsum, twi_crop)
plot(temp_crop)


# raster calculaton (* by the effect size) for thermsum only
r <- 45.17 + (temp_crop*-0.10) + twi_crop
plot(r)

#plot mean min and max temp outputs with spiral points
green_temp_predict <- plot(r,main = "Predicted Green Tea Mass Loss (%)") 
green_temp_predict <- green_temp_predict + points(tea_full$X,tea_full$Y, pch=10, cex=.1, col = 1)

# reorder plot names
tea_full$Plot <- factor(tea_full$Plot, levels =  c("FO", "DG", "FP",
                                                     "ECS", "ECW", "ECN",
                                                     "CHKO","CHHE"))
# green tea plot
(ridge_green <- ggplot(tea_full) +                                     
  aes(x = Green, fill = Plot, color = Plot, alpha = 0.9) +
  geom_density(adjust = 1L) +
    scale_fill_viridis_d(option = "inferno") +
    scale_color_viridis_d(option = "inferno") +
  labs(x = "Observed Green Tea Mass Loss (%)", y = "Frequency", fill = "Plot") +
  theme_classic())

(ridge_green <- (ridge_green + guides(alpha = FALSE, color = FALSE)))

# soil moist plot
(ridge_soil <- ggplot(tea_full) +                                     
    aes(x = Soilmoist, fill = Plot, color = Plot, alpha = 0.9) +
    geom_density(adjust = 1L) +
    scale_fill_viridis_d(option = "inferno") +
    scale_color_viridis_d(option = "inferno") +
    labs(x = "Observed Soil Moisture (%)", y = "Frequency", fill = "Plot") +
    theme_classic())

(ridge_soil <- (ridge_soil + guides(alpha = FALSE, color = FALSE)))

#save surf temp plots
ggsave(green_temp_predict, filename = "figures/manuscript/thermsum_predict.png", height = 7, width = 9)


