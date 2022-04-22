### 1. Cleaning teabag weight data ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Date Created: 21st January 2020 ###

#### LOAD PACKAGES ####

library(tidyverse)
library(esquisse)
library(ggridges)
library(gridExtra)
library(forcats)
library(readxl)
library(lme4)
library(ggeffects)
library(MCMCglmm)
library(viridis)
library(gstat)
library(viridis)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)

#### SET UP WORKSPACE ####

tea <- read.csv(file="data/drone_weights_clean.csv") # read in clean weight data

#### RENAME PLOT LEVELS ####
levels(tea$Plot) # check level order
levels(tea$Plot) <- c("","?","Collison Head HE",
                      "Collison Head KO","Dry Grass",
                      "Dry Control","East Creek North",
                      "East Creek South", "East Creek Wet",
                      "Flooded Orca", "Flower Plots", "Wet Control")

#### RECALCULATE MASS LOSS PERCENTAGE ####

tea$Tea_final <- as.numeric(as.character(tea$Tea_final)) # final weights converted to numeric
tea$Tea_init <- as.numeric(as.character(tea$Tea_init)) # initial weights converted to numeric
tea$Tea_loss <- 1-(tea$Tea_final/tea$Tea_init) # mass loss calculcated from final & initial weights
tea$Tea_loss <- tea[!is.na(tea$Tea_loss),1] # remove the NAs forced by coersion
tea$Plot <- na.omit(tea$Plot) # remove the NAs forced by coersion
tea$Tea_loss <- na.omit(tea$Tea_loss)

# Further cleaning
tea <- tea %>%
  filter(!(Plot %in% c("", "?"))) %>%
  filter(!(Tea_loss %in% ""))  %>%
  filter(!Tea_loss %in% "-Inf") %>%
  filter(!Plot %in% "Dry Control")


# Simple Red & Green Mass loss density plot
(p1 <- ggplot(tea) +                                     
 aes(x = Tea_loss, fill = Tea_Type, color = Tea_Type, alpha = 0.9) +
 geom_density(adjust = 1L) +
 scale_fill_brewer(palette = "Dark2") +
 scale_color_brewer(palette = "Dark2") +
 labs(x = "Mass Loss (g)", y = "# Teabags", fill = "Teabag Type") +
 theme_linedraw())


# Geom ridge density plot by site
(p2 <-ggplot(tea, mapping = aes(y = Plot, x = Tea_loss*100, fill = Tea_Type, ordered=TRUE)) +
  geom_density_ridges(alpha = 0.8,
                      jittered_points = TRUE, point_alpha=0.2) +
  scale_fill_manual(values= c("#00AFBB","#FC4E07")) +
  scale_color_manual(values= c("#00AFBB","#FC4E07")) +
  labs(x = "Mass Loss %", y = "QHI Sites", fill = "Teabag Type") +
  theme_linedraw())


#### CALCULATE STABILISATION RATE & DECOMPOSITION RATE ####

# define fixed hr and hg values
hr <- 0.552
hg <- 0.842

# total teabag burial days - overall length of experiment
days <- 28 

# need to clean tea_id so i can pair green and rooibos
tea$Tea_ID<-as.character(tea$Tea_ID) # convert Tea ID to character
tea$Tea_ID <- substr(tea$Tea_ID,1,nchar(tea$Tea_ID)-1) # remove last value (G or R) from character

# simplify and elongate dataframe  to pair green and rooibos by ID
tea_long <- tea %>%
  group_by(Plot,Tea_Type) %>% # remove empty variables
  dplyr::select(Tea_ID,Tea_loss)  %>%
  spread(Tea_Type,Tea_loss) # individual Green and Red mass loss columns

# calculate stabilisation rate
tea_long$stab <- 1-(tea_long$Green/hg)
tea_long$ar<-hr*(1-tea_long$stab)
tea_long$rate<-log(tea_long$ar/((1-tea_long$Rooibos)-(1-tea_long$ar)))/days #Calculate k

# Geom ridge density plot by site - RATE

(p3 <- ggplot(tea_long, aes(x = rate, y = Plot, fill = Plot)) +
  geom_density_ridges(alpha = 0.8,jittered_points = TRUE, point_alpha=0.2) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Decomposition rate (k)", y = "QHI Sites") +
  theme_linedraw())
p3 <- (p3 + theme(legend.position = "none"))

# Geom ridge density plot by site - STABILISATION

(p4 <- ggplot(tea_long, aes(x = stab, y = Plot, fill = Plot)) +
  geom_density_ridges(alpha = 0.8,jittered_points = TRUE, point_alpha=0.2) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Stabilisation Factor (S)", y = "QHI Sites") +
  theme_linedraw())

p4 <- (p4 + theme(legend.position = "none"))
       
#### MERGE ACTIVE LAYER, SOIL MOISTURE AND SOIL PROPERTIES TO LONG DATA ####

#load full soil property compilation for soil moisture and active layer depth readings
soilprop <- read_excel("data/Drone teabag experiment 2017 9th Aug.xlsx",
                       skip=1)

#change ID column name
names(soilprop)[names(soilprop)=="Teabag N"]<- "Tea_ID"

#reorder data frame so ID comes first & remove redundant columns
soilprop <- soilprop[,c(8,2,3,4,5,9,11,12,13)]

#prepare to join tables by converting Tea_ID to numeric
tea_long <- tea_long %>% 
  mutate(Tea_ID = as.numeric(Tea_ID))

#join tables together
tea_full <- left_join(tea_long, soilprop)

#rename soil variables
names(tea_full)[names(tea_full)=="Active layer depth (cm)"]<- "ALT"
names(tea_full)[names(tea_full)=="Soil Moisture"]<- "Soilmoist"

#convert mass loss to percentage values
tea_full$Green <- tea_full$Green*100
tea_full$Rooibos <- tea_full$Rooibos*100
tea_full$Green <- as.numeric(tea_full$Green)

# Plot soil moisture properties by site
soilplot <- ggplot(tea_full) +
  aes(x = Plot, y = Soilmoist, fill = Plot, colour = Plot) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_viridis_d(option = "plasma") +
  labs(x = "Plot", y = "Soil moisture %") +
  theme_classic() 
soilplot <- soilplot + theme(legend.position = "none") 

# Plot active layer properties by site
altplot <- ggplot(tea_full) +
   aes(x = Plot, y = ALT, fill = Plot, colour = Plot) +
   geom_boxplot() +
   scale_fill_viridis_d(option = "plasma") +
   scale_color_viridis_d(option = "plasma") +
   labs(x = "Plot", y = "Active Layer Depth (cm)") +
   theme_classic() 
altplot <- altplot+ theme(legend.position = "none") 

#Save soil prop as its own CSV
write.csv(tea_full, file = "data/soilpropertieslong.csv", row.names = FALSE)


#### MERGE SPATIAL POINT DATA ####
#load in file containing GPS values
topoprop <- read.csv(file="data/spatial/tea_topo_extract.csv")

#prepare to join tables by converting Tea_ID to factor
tea_full <- tea_full %>% 
  mutate(Tea_ID = as.factor(Tea_ID))

#reorder data frame so ID comes first & remove redundant columns
topoprop <- topoprop[,c(5,7,8,9,10,1,2)]
tea_full <- tea_full[,c(2,3,4,5,6,7,1,8,9,10,11,12,13,14,15)]

#join tables together
tea_full <- left_join(tea_full,topoprop,copy=TRUE)

#reorder data frame so ID comes first & remove redundant columns
tea_full <- tea_full[,c(1,2,3,4,5,6,7,8,
                        9,10,11,12,13,14,15,16,17,20,21)]


#### MERGE MICROCLIMA AND TOPO DATA ####
#load raster files from microclima output, 1m resolution for study area
dir <- "data/spatial/microclima/"
files <- list.files(path = dir, pattern = ".tif")
rasters <- lapply(paste0(dir, files), raster)
meantemp <- raster('data/spatial/microclima/meantemp.tif')
mintemp <- raster('data/spatial/microclima/mintemp.tif')
maxtemp <- raster('data/spatial/microclima/maxtemp.tif')
thermsum <- raster('data/spatial/microclima/thermsum.tif')
elev <- raster('data/spatial/elevation.tif')
aspect <- raster('data/spatial/aspect.tif')
slope <- raster('data/spatial/slope.tif')

#plot mean min and max temp outputs with spiral points
meanplot <- plot(meantemp, main = "Mean Surface Temperature (\u00B0C)", col = inferno(100))
thermplot <- plot(thermsum, main = "28-Day Thermal Sum (\u00B0C)", col = inferno(100))

 #save surf temp plots
ggsave(thermplot, filename = "figures/manuscript/thermplot.png", height = 7, width = 9)

#extract climate values for the teapoints
str(tea_full)
plot(thermsum,col=inferno(100), main="Mean temperature")
points(tea_full$X,tea_full$Y, pch=10, cex=.1, col = 3)

# Remove NA rows and prepare spatial points layer
tea_fullXY <- tea_full[!is.na(tea_full$X), ]
tea_fullXY_sp <- SpatialPoints(cbind(tea_fullXY$X, tea_fullXY$Y), 
                            proj4string=CRS("+init=epsg:32607"))
cord.UTM_mean <- spTransform((tea_fullXY_sp), crs(meantemp))
cord.UTM_max <- spTransform((tea_fullXY_sp), crs(maxtemp))
cord.UTM_min <- spTransform((tea_fullXY_sp), crs(mintemp))
cord.UTM_therm <- spTransform((tea_fullXY_sp), crs(thermsum))
cord.UTM_elev <- spTransform((tea_fullXY_sp), crs(elev))
cord.UTM_slope <- spTransform((tea_fullXY_sp), crs(slope))
cord.UTM_aspect <- spTransform((tea_fullXY_sp), crs(aspect))


## Now we just need to extract the raster values from these locations                
mean_point <- raster::extract(meantemp, cord.UTM_mean, method='simple',df=TRUE)
max_point <- raster::extract(maxtemp, cord.UTM_max, method='simple',df=TRUE)
min_point <- raster::extract(mintemp, cord.UTM_min, method='simple',df=TRUE)
therm_point <- raster::extract(thermsum, cord.UTM_min, method='simple',df=TRUE)
elev_point <- raster::extract(elev, cord.UTM_min, method='simple',df=TRUE)
slope_point <- raster::extract(slope, cord.UTM_min, method='simple',df=TRUE)
asp_point <- raster::extract(aspect, cord.UTM_min, method='simple',df=TRUE)

## need to write CSV file to join to the master datasheet with all varibles 
combinePointValue=cbind(tea_fullXY,mean_point,min_point,
                        max_point,therm_point, elev_point,
                        slope_point,asp_point)
str(combinePointValue)
#remove waste columns 
combinePointValue$ID...20 <- NULL
combinePointValue$ID...22 <- NULL
combinePointValue$ID...24 <- NULL
combinePointValue$ID...26 <- NULL
combinePointValue$ID...28 <- NULL
combinePointValue$ID...29 <- NULL
combinePointValue$ID...30 <- NULL
combinePointValue$ID...31 <- NULL
combinePointValue$ID...32 <- NULL
combinePointValue$ID...33 <- NULL
combinePointValue$Aspect <- NULL



#### Save full prop dataset as its own CSV ####
write.csv(combinePointValue, file = "data/fullpropertieslong.csv", row.names = FALSE)


