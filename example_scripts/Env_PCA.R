
setwd("C:/Users/User/Desktop/2020_Corona_Fun/Coding/Rm1_ToadRDA/Code")

library(factoextra)
library(FD)
library(tibble)
library("ggpubr")

######CANE TOAD MASTER SHEET#############

ToadData_D3 <- read.csv("Cane_Toad_Master_Sheet_SC_NoNSW.CSV",stringsAsFactors=TRUE,sep=",")
str(ToadData_D3)

ToadVars_D3 <- ToadData_D3[c(10:19)]
ToadVars_D3 <- as.data.frame(lapply(ToadVars_D3, as.numeric))
str(ToadVars_D3)

PopCoords_D3 <- ToadData_D3[c(5,6)]
PopCoords_D3 <- as.data.frame(lapply(PopCoords_D3, as.numeric))
str(PopCoords_D3)

State <- ToadData_D3[c(3)]
State <- as.data.frame(lapply(State, as.factor))
str(State)

#####

#library(raster)
#climdata <- getData('worldclim',download=TRUE,var='bio',res=5)

points_D3 <- SpatialPoints(PopCoords_D3, proj4string=climdata@crs)
values_D3 <- extract(climdata,points_D3)
envdata_D3_NA <- cbind.data.frame(PopCoords_D3,values_D3)
str(envdata_D3_NA)
envdata_D3 <- envdata_D3_NA[rowSums(is.na(envdata_D3_NA)) == 0,]
str(envdata_D3)


colnames(envdata_D3) <- c("long","lat",
                          "AnnualMeanTemp","MeanDiurnalRange",
                          "Isothermality","TempSeasonality",
                          "MaxTempofWarmestMonth","MinTempofColdestMonth",
                          "TempAnnualRange","MeanTempofWettestQuarter",
                          "MeanTempofDriestQuarter","MeanTempofWarmestQuarter",
                          "MeanTempofColdestQuarter","AnnualPrecipitation",
                          "PrecipitationofWettestMonth","PrecipitationofDriestMonth",
                          "PrecipitationSeasonality","PrecipitationofWettestQuarter",
                          "PrecipitationofDriestQuarter","PrecipitationofWarmestQuarter",
                          "PrecipitationofColdestQuarter")
str(envdata_D3)

rowSums(is.na(envdata_D3_NA))

PopCoords_D3 <- PopCoords_D3[rowSums(is.na(envdata_D3_NA)) == 0,] ##REMOVE THE NA's that resulted from RASTER package


####

#install.packages("rgdal")

#library(rgdal)

#extracting aridity index
raster_ai <- raster("global_ai_et0/ai_et0/ai_et0.tif")
points_ai_D3 <- SpatialPoints(PopCoords_D3, proj4string=raster_ai@crs)
identicalCRS(raster_ai,points_ai_D3)
aridity_index_D3 <- extract(raster_ai,points_ai_D3)
str(aridity_index_D3)
aridity_index_D3

#extracting PET
raster_pet <- raster("global-et0_annual.tif/et0_yr/et0_yr.tif")
points_pet_D3 <- SpatialPoints(PopCoords_D3, proj4string=raster_pet@crs)
identicalCRS(raster_pet,points_pet_D3)
potential_evap_D3 <- extract(raster_pet,points_pet_D3) #Potential Evapotranspiration
str(potential_evap_D3)

#calculating Daylength
#install.packages("geosphere")
library("geosphere")
day_length_D3 <- NULL
counter <- 1
for (val in PopCoords_D3$Latitude) {
  day_length_D3[[counter]] <-   print(
    mean(
      daylength(val, 1:365)
    )
  )
  names(day_length_D3[[counter]]) = names(val)
  counter <- counter + 1
}

str(day_length_D3)


#Combining all extra env variables
State_NoNA <- as.data.frame(State[rowSums(is.na(envdata_D3_NA)) == 0,])
str(State_NoNA)

envdata_D3 <- cbind.data.frame(State_NoNA,envdata_D3,aridity_index_D3,potential_evap_D3,day_length_D3)
str(envdata_D3)

envdata_D3_NoNA <- envdata_D3[rowSums(is.na(envdata_D3)) == 0,] #Remove NA's that resulted from 3 extra env. measures
str(envdata_D3_NoNA)

names(envdata_D3_NoNA)[names(envdata_D3_NoNA) == "State[rowSums(is.na(envdata_D3_NA)) == 0, ]"] <- "Sample_State"
str(envdata_D3_NoNA)

##PCA all env variables

wdbc.pr <- prcomp(envdata_D3_NoNA[,c(4:25)], center = TRUE, scale = TRUE)
wdbc.pr

summary(wdbc.pr)

#library("factoextra")
fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = envdata_D3_NoNA$Sample_State, 
             col.ind = "black", 
             palette = c("#AA3377","#CCBB44", "#228833","#EE6677","#66CCEE"), 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "State") +
  theme(plot.title = element_text(color="white", hjust = 0.5))+
  xlab("Axis 1 (58.3%)") +
  ylab("Axis 1 (18.0%)") +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0.80))  +
  theme(legend.title = element_text(colour="black", size=12, 
                                    face="bold")) +
  theme(legend.text = element_text(colour="black", size=12, 
                                   face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
