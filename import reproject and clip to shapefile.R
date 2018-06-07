#Ancilleno Davis
#Geo 560 rcode

#This r-code has been developed by Ancilleno Davis in partial completion of the 
#requirements for PhD in Ecology Evolution and Environmental Science at Miami University in Oxford, Ohio.

#The purpose of this Rcode is to import habitat classification data 
#generated in Google Earth Engine from LAndsat 8 OLI reflectance data,
#and trim that classification data to the contiguous terrestrial area of Grand Bahama Island
#This will allow the use of more efficient summary and analysis of the remaining habitat data

#### Initial workspace parameters####
#set seed
set.seed(1981)

#Set my working directory and load necessary libraries
setwd("C:/Users/davisao2/Desktop/open source GIS/eBird r code")
#note this working directory should include:
#1: the geotiff of the raster data created using the Google Earth Engine habitat classification

#2: the shapefile of the outline of the terrestrial area of Grand Bahama and all associated files


#Load libraries
require(caret)#this has confusionMatrix for determining raster accuracy
require(cooccur)
require(data.table)#to allow rearrangment of columns by their headings
require(EcoSimR)
require("FactoMineR")#used for the Multiple Correspondence Analysis of the raster datasets
require("factoextra")#used for the Multiple Correspondence Analysis of the raster datasets
require(fmsb)
require(geosphere)
require(ggplot2)
require(ggpubr)
require(gstat)
require(lme4)
require(lmerTest)
require(plyr)
require(psych)
require(raster)#work with raster items
require(readr)
require(readxl)
require(reshape2)#contains dcast for reorganizing data into presence absence etc.
require(rgdal)#read in shapefiles
require(rgeos)
require(RStoolbox)
require(sp) #work with spatial items

dev.off()#this resets the graphics window so that margins and other adjustments are returned to default

#Assign the colors I would like to represent my habitat classes####


col7class=c(#include all 7 colors including water. This will be used to plot the raster image.
  "blue", #ocean 
  "dark green", #Pine Forest
  "brown", #Wetlands 
  "seashell2",#sand
  "gray", #Urban
  "Light green", #grass
  "yellow")#high water table plant communities
##Assign label names for the habitat classes####
legend=c("Water", 
         "Pine", 
         "Wetland", 
         "Sand", 
         "Urban", 
         "Grass", 
         "HWTC")
#import the Outline of Grand Bahama Island and reproject to WGS84####
GBOutline<-readOGR(".","Grand_Bahama_Outline")
print(proj4string(GBOutline)) #details of the projection for GBOutline
plot(GBOutline, 
     axes=T, #remove the lat long coordinates with F and add them with T
     border="black")

#Reproject GBOutline to WGS84
GBOutlineWGS84 <- 
  spTransform(GBOutline, #this command reprojects a spatial item  
              CRS("+proj=longlat +datum=WGS84"))#the new projection is WGS 84

plot(GBOutlineWGS84, 
     axes=TRUE, 
     border="black",
     las=1,
     #ylab="Latitude",
     # xlab="Longitude",
     main="Grand Bahama Contiguous Area")
###Import Rasters of habitat classification ####

#Random Forests classification map created using
#6 Landsat 8 OLI bands 2,3,4,5,6,7
#collected during 2017-01-01 through 2017-12-31

RF7Classes2017<-raster("RFclass6bands30m7classes2017.tif")
RF7Classes2017@crs #finds the coordinate reference system which is already WGS84
RF7Classes2017@extent #gives me the boundaries of the raster

plot(RF7Classes2017, 
     las=1, #rotates the y-axis labels horizontal
     ylab="Latitude",
     xlab="Longitude",
     main = "Grand Bahama Habitat Map")
#overlay the reprojected shapefile on the raster####
plot(GBOutlineWGS84, add=T)
#contiguous area habitat

#Clip area to the contiguous area of Grand Bahama Island####
GBcontiguoushab<-#this vector will hold our clipped raster
  mask(RF7Classes2017#This is the raster image to be clipped
                      , GBOutlineWGS84)#This is the shape you are clipping to

plot(GBcontiguoushab); plot(GBOutlineWGS84, add=T) #Plot the trimmed raster and the outline
summary(GBcontiguoushab) #note all the clipped area is now NA
##Save the newly clipped raster to a file you can use in other software####
writeRaster(GBcontiguoushab,#our new rater image
            "GBcontiguoushab.tif"#the file name we want to save it to
            , format="GTiff") #the raster format we are using.
##display the clipped raster####
GBcontiguoushab.range<-cellStats(GBcontiguoushab, range)
dev.off()
plot(GBcontiguoushab,#plot this raster file
     col=col7class, #a list of colors make sure they are in the order of the classes you are plotting
     axis.args=list(
      at= #how do you want the axis labels spread out?
       seq(#sequential separations from the 
        GBcontiguoushab.range[1],#lowest end of the range
        GBcontiguoushab.range[2], #to the highest end of the range
      1),#with one unit of separation between each tick mark
      labels=legend, #at each tick mark, what do you want the label to say?
     cex.axis=1),
      axes=T ,
     main = "Grand Bahama Island, Bahamas habitat map", #title of the plot 
     ylab="Latitude", #left side 
     xlab="Longitude", #x axis label
     las=1 #rotate the y axis labels to horizontal
)
plot(GBOutlineWGS84, add=T)
legend('bottomright', legend , fill=col7class, border="black",
        col=col7class, bty='n', cex=1.5)
#summarize the number of pixels and the geographic area in each class####
GBhabitatdf<-as.data.frame(GBcontiguoushab)
habitatpixels<-count(GBhabitatdf, "RFclass6bands30m7classes2017")
habitatareakm2<- #a vector of the total area in each habitat type
  (habitatpixels$freq[1:7] #1:7 selects the columns with the classes 0-6 and ignores the NA values
                 )*0.0009 #0.0009 is the number of km2 in a 30 x 30 m pixel (Landsat 8 imagery)
habitatareakm2
totalpixels<- #vector for the total pixels in the study area
  sum(habitatpixels$freq[1:7])
totalarea<-#vector of total area in study area
  sum(habitatareakm2)
percentarea<-habitatareakm2/totalarea #%of total area in each habitat type

habitat<-data.frame(legend,habitatpixels$freq[1:7],habitatareakm2,percentarea)

##Display barplots of the pixel area in each habitat type.####
BpGBHabitatPixels<-barplot(GBcontiguoushab, 
                     col=col7class,
                     axes=F,
                     main="Pixels per terrestrial habitat class on Grand Bahama 2017", 
                     xlab="Habitat type",
                     ylab="Number of Pixels",
                     horiz=FALSE,
                     ylim=c(0,400000),
                     names.arg=c("Water", 
                                 "Pine", 
                                 "Wetland", 
                                 "Sand", 
                                 "Urban", 
                                 "Grass", 
                                 "HWTC"),
                     las=1)
text(BpGBHabitatPixels, 
     habitatpixels$freq[1:7], 
     label=habitatpixels$freq[1:7], 
     pos = 3, 
     xpd = NA)
##Display barplots of the km2 area in each habitat type.####

BpGBHabitatkm2<-barplot(habitatareakm2, 
                     col=col7class,
                     axes=F,
                     main="Square Km per terrestrial habitat class on Grand Bahama 2017", 
                     xlab="Habitat type",
                     ylab="Area (sq.Km)",
                     horiz=FALSE,
                     ylim=c(0,400),
                     names.arg=c("Water", 
                                 "Pine", 
                                 "Wetland", 
                                 "Sand", 
                                 "Urban", 
                                 "Grass", 
                                 "HWTC"))
text(BpGBHabitatkm2, 
     habitatareakm2, 
     label=habitatareakm2, 
     pos = 3, 
     xpd = NA)

