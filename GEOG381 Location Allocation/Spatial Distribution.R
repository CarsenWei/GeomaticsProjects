#set my current working directory 
#setwd()

# ensure the necessary packages
library(rgdal)  #including readOGR() to import shape files, writeOGR() to make shape files
library(sp)     #including functions for creating and visualizing a patial point
library(rgeos)  #including gBuffer() 
library(spatstat) #including quadratcount() 
library(maptools) #including as()

#import data (filename, layer)
OHD <- readOGR("OntarioHomeDepot.shp", "OntarioHomeDepot") 
ONP <- readOGR(dsn="OntarioProvince.shp", layer="OntarioProvince") 

#mean coordinate points and standard distance 
OHD_X <- OHD$X
OHD_X_mean <- mean(OHD_X) 
OHD_Y <- OHD$Y
OHD_Y_mean <- mean(OHD_Y) 
  print (paste("x coordinate mean = ", format(OHD_X_mean, nsmall=2))) 
  print (paste("Y coordinate mean = ", format(OHD_Y_mean, nsmall=2))) 
xResSqrd <- (OHD_X - OHD_X_mean) ^2 
yResSqrd <- (OHD_Y - OHD_Y_mean) ^2
xSumSqrdRes <- sum(xResSqrd) 
ySumSqrdRes <- sum(yResSqrd) 
  print (paste("sum of the squared X residuals = ", format(xSumSqrdRes, nsmall=2))) 
  print (paste("sum of the squared Y residuals = ", format(ySumSqrdRes, nsmall=2))) 
stdDistance <-sqrt((xSumSqrdRes+ySumSqrdRes)/nrow(OHD))
  print (paste("standard distance= ", format(stdDistance, nsmall=2))) 

#creating and visualizing a patial point
meanCentre <- data.frame(OHD_X_mean, OHD_Y_mean) 
meanPointLoc <- cbind(meanCentre$OHD_X_mean, meanCentre$OHD_Y_mean) 
coordinates(meanCentre) <- meanPointLoc 
  print(meanCentre) 
OHD_proj <- proj4string(OHD) 
  print(OHD_proj) 
proj4string(meanCentre) <- OHD_proj 
  writeOGR(meanCentre, "meanCentre.shp", "meanCentre", driver="ESRI Shapefile", overwrite_layer = TRUE) 

#visualizing
plot(ONP, col="light gray") 
plot(OHD, col="orange", add=TRUE) 
plot(meanCentre, col="red", add=TRUE) 

#number of Home Depot stores within each standard distance 
stdist1 <- stdDistance
  print (format(stdist1, nsmall=2))
stdist2 <- stdist1 *2
  print (format(stdist2, nsmall=2))
stdist3 <- stdist1 *3
  print (format(stdist3, nsmall=2))
stdist1_buf <- gBuffer(meanCentre, width=stdist1)
  plot (stdist1_buf, add=TRUE) 
stdist2_buf <- gBuffer(meanCentre, width=stdist2)
  plot (stdist2_buf, add=TRUE) 
stdist3_buf <- gBuffer(meanCentre, width=stdist3)
  plot (stdist3_buf, add=TRUE) 
stdist1_HD_v <- over(OHD, stdist1_buf) 
  print(stdist1_HD_v) 
stdist2_HD_v <- over(OHD, stdist2_buf) 
  print(stdist2_HD_v) 
stdist3_HD_v <- over(OHD, stdist3_buf) 
  print(stdist3_HD_v) 
stdist1_HD_c <- sum(stdist1_HD_v, na.rm = TRUE) 
  print (stdist1_HD_c) 
stdist2_HD_c <- sum(stdist2_HD_v, na.rm = TRUE) 
  print (stdist2_HD_c) 
stdist3_HD_c <- sum(stdist3_HD_v, na.rm = TRUE) 
  print (stdist3_HD_c) 

#convert OHD to ppp
  ONP_minX <- ONP@bbox[1] 
  ONP_minY <- ONP@bbox[2] 
  ONP_maxX <- ONP@bbox[3] 
  ONP_maxY <- ONP@bbox[4] 
OHD_ppp <- ppp(OHD$X, OHD$Y, c(ONP_minX,ONP_maxX), c(ONP_minY, ONP_maxY)) 
  print (OHD_ppp) 
#specify the window
  ONP_win <- owin(xrange=c(ONP@bbox[1], ONP@bbox[3]), yrange=c(ONP@bbox[2], ONP@bbox[4]), poly=NULL,units="m") 
#vectors identifying where the quadrats should be crated 
  xQuadratBreaks <- c() 
  x_break <- ONP_minX 
  while (x_break < ONP_maxX + 50000) 
    {   xQuadratBreaks <- c(xQuadratBreaks, x_break)   
        x_break <- x_break + 50000 
    } 
  yQuadratBreaks <- c() 
  y_break <- ONP_minY
  while (y_break < ONP_maxY + 50000) 
  {   yQuadratBreaks <- c(yQuadratBreaks, y_break)   
      y_break <- y_break + 50000 
  } 
# create and plot quadrats
OHD_quadrats <- quadrats(ONP_win, xbreaks = xQuadratBreaks, ybreaks = yQuadratBreaks) 
#plot(OHD_quadrats, add=TRUE) #skip this plotting for later
#convert quadrats to spatial objects; assign georeference information
OHD_quadrats_sp <- as(OHD_quadrats, "SpatialPolygons") 
proj4string(OHD_quadrats_sp) <- proj4string(ONP) 
#subset the quadrats intersecting Ontario
OHD_quadrat_subset <- OHD_quadrats_sp[ONP] 
plot(OHD_quadrat_subset, add=TRUE) 

#quadrat counts (the number of Home Depot stores within each Ontario quadrat)
OHD_quadrat_subset_count <- length(OHD_quadrat_subset) 
#create a new spatial data frame where Home Depot stores are attatched with quadrat ID
HDstore_quadratID <- over(OHD, OHD_quadrat_subset) 
#table() counts the # of events in quadrats with events
HDStore_quadrat_tabulation <- table(HDstore_quadratID) 
  print(HDStore_quadrat_tabulation)
#the number of quadrats with and without events 
OHD_quadrat_subset_w_store_count <- nrow(HDStore_quadrat_tabulation)
OHD_quadrat_subset_wo_store_count = OHD_quadrat_subset_count - OHD_quadrat_subset_w_store_count
  print(paste("number of total quadrats = ", OHD_quadrat_subset_count, " ; number of quadrats with stores = ", OHD_quadrat_subset_w_store_count)) 

#to acquire the distribution of quadrat values  
HDStore_quadrat_distribution <- table(HDStore_quadrat_tabulation) 
  print(HDStore_quadrat_distribution)
HDStore_no_of_events_k <- as.vector(as.numeric(names(table(table(HDstore_quadratID))))) 
#frequency of quadrats with different event 
HDStore_no_of_quadrats_x <- as.vector(HDStore_quadrat_distribution)   
#c(A,B): COMBINE A TO THE FRONT OF B
HDStore_no_of_quadrats_x <- c(OHD_quadrat_subset_wo_store_count,HDStore_no_of_quadrats_x) 
HDStore_no_of_events_k <- c(0, HDStore_no_of_events_k) 
#the number of quadrats 
  print (sum(HDStore_no_of_quadrats_x)) 
#number of events 
  print (sum(HDStore_no_of_quadrats_x * HDStore_no_of_events_k)) 
mu <- nrow(OHD)/sum(HDStore_no_of_quadrats_x) 
  print(mu) 
