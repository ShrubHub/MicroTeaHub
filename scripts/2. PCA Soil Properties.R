### 2. PCA using soil properties data ###
### Elise Gallois, elise.gallois94@gmail.com ###
### Date created: 9th April 2020 ###

#### LOAD PACKAGES ####

library(tidyverse)
library(esquisse)
library(ggpubr)
library(forcats)
library(FactoMineR)
library(factoextra)
library(corrplot)

#### SET UP WORKSPACE ####
soil_long <- read.csv(file="data/fullpropertieslong.csv") # read in soil data

# choose active -aka numeric- columns
#reorder data frame so ID comes first & remove redundant columns
colnames(soil_long)
soil_long <- soil_long[,c(7,6,1,2,3,5,10,11,23,24,25,26)]
# reorder plot names
soil_long$Plot <- factor(soil_long$Plot, levels =  c("FO", "DG", "FP",
                                                     "ECS", "ECW", "ECN",
                                                     "CHKO","CHHE"))

# create dataframes with and without relevant variables - and with no NAs
soil.active <- soil_long[3:12]
soil.active <- na.omit(soil.active, na.action = "omit")
soil.inactive <- na.omit(soil_long, na.action = "omit")

# rename soil active variables
soil.active <- soil.active %>% 
  dplyr::rename(
    "Green tea mass loss" = Green, 
    "Rooibos tea mass loss" = Rooibos,
    "Soil moisture" = Soilmoist,
    "Stabilisation Factor" = stab,
    "Decomposition rate" = rate,
    "Thermal Sum" = thermsum,
    "Elevation" = elevation,
    "Slope" = slope,
    "Aspect" = aspect) 

# get variables factor map
soil.pca <- PCA(soil.active, graph = TRUE)

# generate scree plot - what is relative importance of each dimension?
fviz_eig(soil.pca, addlabels = TRUE, ylim = c(0, 50)) + 
  ggtitle("c)") +
  theme_classic(base_size = 16) 
  
var <- get_pca_var(soil.pca)

# Coordinates
head(var$coord)

# Cos2: quality on the factor map
head(var$cos2)

# Contributions to the principal components
head(var$contrib)

# Coordinates of variables
head(var$coord, 4)

# Examine dimensions 1 and 2
fviz_pca_var(soil.pca, col.var = "black") 

#correlation plot - weighting of variable by dimension
(correlation_subplot <- corrplot(var$cos2,tl.col = 'black',is.corr=TRUE)) + ggtitle("d)")  


# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(soil.pca, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(soil.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Change the transparency by cos2 values
fviz_pca_var(soil.pca, alpha.var = "cos2")
head(var$contrib, 4)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(134)
soil.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(soil.km$cluster)

# Color variables by groups
fviz_pca_var(soil.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster",
             repel = TRUE)

#look at individual contributions
ind <- get_pca_ind(soil.pca)
ind

# Coordinates of individuals
head(ind$coord)

# Quality of individuals
head(ind$cos2)

# Contributions of individuals
head(ind$contrib)
fviz_pca_ind(soil.pca)
fviz_pca_ind(soil.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(soil.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)


#contribution of individuals
fviz_cos2(soil.pca, choice = "ind")

# Total contribution on PC1 and PC2
fviz_contrib(soil.pca, choice = "ind", axes = 1:2)

# Create a random continuous variable of length 136,
# Same length as the number of active individuals in the PCA
set.seed(136)
my.cont.var <- rnorm(136)

# Color individuals by the continuous variable
fviz_pca_ind(soil.pca, col.ind = my.cont.var,
             gradient.cols = NULL,
             legend.title = "Cont.Var")

#colour by groups
head(soil.active, 8)

# PCA - cluster by Plot?
fviz_pca_ind(soil.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = soil.inactive$Plot , # color by groups
             addEllipses = TRUE, ellipse.type = "confidence",# Concentration ellipses,
             legend.title = "Groups")


# PCA - cluster by soil?
fviz_pca_ind(soil.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = soil.inactive$Soil , # color by groups
             palettepalette = "pal",
             addEllipses = TRUE, ellipse.type = "confidence",# Concentration ellipses,
             legend.title = "Groups")

# PCA - cluster by vegetation?
fviz_pca_ind(soil.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = soil.inactive$Vegetation , # color by groups
             palettepalette = "pal",
             addEllipses = TRUE, ellipse.type = "confidence",# Concentration ellipses,
             legend.title = "Groups")


(biplotdim12 <- fviz_pca_biplot(soil.pca, axes = c(1,2),
                # Individuals
                geom.ind = "point",
                fill.ind = soil.inactive$Plot, col.var="black",
                col.ind = soil.inactive$Plot, 
                pointshape = 21, pointsize = 2,alpha = 0.9,
                addEllipses = TRUE, ellipse.type = "confidence",
                ellipse.alpha = 0.6,
                palettepalette = "viridis",
                  repel = TRUE) +
    theme_classic(base_size = 16) +
    scale_fill_viridis(discrete=TRUE, guide=FALSE,option = "plasma",soil.inactive$Plot, labels = c("Orca Floodplain",
                                                                                                  "Dry Grass",
                                                                                                  "Flower Plots",
                                                                                                  "East Creek South",
                                                                                                  "East Creek Wet",
                                                                                                  'East Creek North',
                                                                                                  "Collison Head Komakuk",
                                                                                                  "Collison Head Herschel")) +
    scale_color_viridis(discrete=TRUE, guide=FALSE,option = "plasma") +
    guides(fill=guide_legend(title="Plot (West to East)"))+
    ggtitle("a)") +
    guides(color = FALSE))

(biplotdim23 <- fviz_pca_biplot(soil.pca, axes = c(2,3),
                                # Individuals
                                geom.ind = "point",
                                fill.ind = soil.inactive$Plot, col.var="black",
                                col.ind = soil.inactive$Plot, 
                                pointshape = 21, pointsize = 2,alpha = 0.9,
                                addEllipses = TRUE, ellipse.type = "confidence",
                                ellipse.alpha = 0.6,
                                palettepalette = "viridis",
                                repel = TRUE) +
    theme_classic(base_size = 16) +
    scale_fill_viridis(discrete=TRUE, guide=FALSE,option = "plasma",soil.inactive$Plot, labels = c("Flooded Orca",
                                                                                                   "Dry Grass",
                                                                                                   "Flower Plots",
                                                                                                   "East Creek South",
                                                                                                   "East Creek Wet",
                                                                                                   'East Creek North',
                                                                                                   "Collison Head Komakuk",
                                                                                                   "Collison Head Herschel")) +
    scale_color_viridis(discrete=TRUE, guide=FALSE,option = "plasma") +
    guides(fill=guide_legend(title="Plot (West to East)"))+
    ggtitle("b)") +
    guides(color = FALSE))


ggsave(biplotdim12, filename = "figures/PCA_biplot1v2.png", height = 7, width = 9)
ggsave(biplotdim23, filename = "figures/PCA_biplot2v3.png", height = 7, width = 9)


simple_biplot <- fviz_pca_biplot(soil.pca, 
                          # Individuals
                          geom.ind = "point",
                          fill.ind = soil.inactive$Plot, col.ind = "black",
                          pointshape = 21, pointsize = 2,
                          addEllipses = TRUE, ellipse.type = "confidence",
                          ellipse.alpha = 0.5,
                          palettepalette = "viridis",
                          # Variables  
                          repel = TRUE)

ind.p <- fviz_pca_biplot(soil.pca, geom = "point", col.ind = soil.inactive$Plot)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              fill.ind = soil.inactive$Plot, col.ind = "black",
              subtitle = "Decomposition Characteristics",
              addEllipses = TRUE, ellipse.type = "confidence",
              xlab = "PC1", ylab = "PC2", 
              legend.title = "Location", legend.position = "top",
              ggtheme = theme_linedraw())






