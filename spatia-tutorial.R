# Author: Mako Shibata (s2471259@ed.ac.uk)
# Date created: 07/10/2025
# Aim: Tutorial for Spatial Analysis in R
# Satellite data available from 
# https://scihub.copernicus.eu/

# set working directory ----
setwd("/Users/Owner/Library/CloudStorage/OneDrive-UniversityofEdinburgh/R/Tutorial10_SpatialAnalysis/SpatialAnalysis_tutorial/CC-spatial")

# load functions ----
theme.LPI <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=12, angle=0, hjust=0.5, vjust=1, colour="black"), 
          axis.text.y=element_text(size=12, colour="black"), 
          axis.title.x=element_text(size=14, face="plain"), 
          axis.title.y=element_text(size=14, face="plain"), 
          panel.grid.major.x=element_blank(), 
          panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), 
          panel.grid.major.y=element_blank(), 
          plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), units= , "cm"), 
          plot.title = element_text(size=20, vjust=1, hjust=0.5), 
          legend.text = element_text(size=12, face="italic"), 
          legend.title = element_blank(), 
          legend.position=c(0.1,0.9), 
          legend.background = element_rect(fill='transparent'))
}
# load packages ----
install.packages("sp")
install.packages("sf")
install.packages("terra")
install.packages("raster")
install.packages("viridis")
install.packages("rasterVis")
library(sp)
library(sf)
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(rasterVis)

# import tif ----
tay <- raster("taycrop.tif")
tay

# create individual raster layer 
b1 <- raster("taycrop.tif", band=1)
b2 <- raster("taycrop.tif", band=2)
b3 <- raster("taycrop.tif", band=3)
b4 <- raster("taycrop.tif", band=4)
b5 <- raster("taycrop.tif", band=5)
b6 <- raster("taycrop.tif", band=6)
b7 <- raster("taycrop.tif", band=7)
b8 <- raster("taycrop.tif", band=8)
b9 <- raster("taycrop.tif", band=9)
b10 <- raster("taycrop.tif", band=10)
b11 <- raster("taycrop.tif", band=11)
b12 <- raster("taycrop.tif", band=12)

# compare two bands to see if the following matches. 
## number of rows/columns
## resolution
## origin 
compareRaster(b2, b3)
# checking the followings are important: 
## coordinate reference system 
## extents 

plot(b8) # only plots 100,000 pixels
image(b8) # stretches the view
zoom(b8) # run this and click on the plot viewer 
#          to set the extent 
e <- drawExtent()
cropped_tay <- crop(b7, e)
plot(cropped_tay)
image(cropped_tay)

# Visualise spectral bands ----
png('tayplot.png', width = 4, height = 4, units = "in", res = 300) # to save plot
image(b8, col= viridis_pal(option="D")(10), main="Sentinel 2 image of Loch Tay")
dev.off()   

image(b8, col= viridis_pal(option="D")(10), main="Sentinel 2 image of Loch Tay")
# option = "D" ...inferno
# (10) ...10 colours from the palette, equally spaced and in gradient

# RGB - multilayered objects 
png("RGB.png", width = 5, height = 5, units = "in", res = 300)
tayRGB <- stack(list(b4, b3, b2)) 

plotRGB(tayRGB, axes = TRUE, stretch = "lin", 
        main = "Sentinel RGB colour composite")
dev.off()

plotRGB(tayRGB, axes = TRUE, stretch = "lin", # <- "lin" stretches the pixel values to 0-225 range 
        # computer's displayable extent
        main = "Sentinel RGB colour composite")

# (help(plotRGB))? to get more arguments on composites 

# Create a FCC of the Loch Tay area using a raster stack.
tayFCC <- stack(list(b8, b3, b2))
plotRGB(tayFCC, axes = TRUE, stretch = "lin", 
        main = "Sentinel false colour composite")

# Visualisation using rasterVis 
# levelplot() ...allows level and contour plots to be made
#                of raster objects with ELEVATION data (e.g., LiDAR)
# plot3D() ...3D mapping (e.g., LiDAR)
# gplot() ...plot of a uni or a multivariate raster object

gplot(b8) + 
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_colour_viridis_c() + 
  coord_quickmap() + 
  ggtitle("West of Loch tay, raster plot") + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("ggtay.png", scale = 1.5, dpi = 300)

# visualise all bands in one figure. 
t <- stack(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11)

b1
gplot(t) + 
  geom_raster(aes(x = x, y = y, fill = value)) +  # values equals pixel values
  scale_fill_viridis_c() + # scale_fill!! not scale_colour
  facet_wrap(~variable) + # variable equals raster layer 
  coord_quickmap() +
  ggtitle("West of Loch Tay, Sentinel-2 All bands") + 
  labs(x = "longitude", y = "latitude") + 
  theme_classic() + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, size = 10, vjust = 0.5), 
        axis.text.y = element_text(size = 10))
  
ggsave("allbands.png", scale = 1.5, dpi = 300)
dev.off()  

# for a quick visualisation, use "brick"
s_tay <- brick('taycrop.tif')
plot(s_tay)
dev.off()
# Manipulate rasters: NDVI and KMN  ----
## Make a function ---- 
VI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi)/(bk + bi)
}

# img <- A raster object with multiple layers inside. 
# img[1] <- A raster object with one layer inside.Still within a multi-layer container.
# img[[1]] <- A raster layer 

# Sentinel-2 NDVI band 8 - band 4 
s_tay <- brick("taycrop.tif")
ndvi_s2 <- VI(s_tay, 8, 4) 
# brick() creates a multi-layer raster object. 

# plot NDVI on maps
ndvi_s2
png("ndviplot.png", width = 4, height = 4, units = "in", res = 300)
plot(ndvi_s2, col = rev(terrain.colors(10)), main = "Sentinel 2, Loch Tay-NDVI")
dev.off()

# create a histogram of NDVI to see distribution
png("ndvihist.png", width = 4, height = 4, units = "in", res = 300)
hist(ndvi_s2, 
     main = "NDVI Distribution Loch Tay", 
     xlab = "NDVI", 
     ylab = "Frequency", 
     col = "lightgray", 
     xlim = c(-0.5, 1), 
     breaks = 30, 
     xaxt = "n") # <- surpresses the default x axis 
axis(side = 1, at = seq(-0.5, 1, 0.05), labels = seq(-0.5, 1, 0.05))
dev.off()

# mask pixels where NDVI < 0.4 to be NA 

png("ndvimasked.png", width = 4, height = 4, units = "in", res = 300)
veg <- reclassify(ndvi_s2, cbind(-Inf, 0.4, NA))
plot(veg, col = rev(terrain.colors(10)), main = "Sentinel 2, Loch Tay-NDVI")
dev.off()
# We are reclassifying our object and making all values between
# negative infinity and 0.4 be NAs

# mask pixels where NDVI < 0.4 to be NA, 0.2 < x < 0.5 to be 1, and < 0.5 to be 2
m <- cbind(
  c(-Inf, 0.2, NA),
  c(0.2,  0.5, 0.5),
  c(0.5,    1,  1)
)

veg2 <- reclassify(ndvi_s2, m)
png("ndvimasked2.png", width = 4, height = 4, units = "in", res = 300)
plot(veg2, col = rev(terrain.colors(10)), main = "Sentinel 2, Loch Tay-NDVI")
dev.off()
print(veg2)

# Save rasters ----
writeRaster(x = ndvi_s2, 
            filename = "/Users/Owner/Library/CloudStorage/OneDrive-UniversityofEdinburgh/R/Tutorial10_SpatialAnalysis/SpatialAnalysis_tutorial/CC-spatial/ndvi.tif", 
            format = "GTiff", # GeoTIFF file 
            datatype = "INT2S")  # save as integer rather than a float

# Kmeans algorithm to run Unsupervised Classification ----
