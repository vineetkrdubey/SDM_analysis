c <- raster("D:/Dahal/Ochotona.sikimaria_BIOMOD/Projection_Current/bio15.asc")
f <- raster("D:/Dahal/Ochotona.sikimaria_BIOMOD/2070_45/bio15.asc")
pts = rasterToPoints(f)

res(f)
res(c)
alt.crop
windows()
plot(c)
list()
P_url<- "https://geodata.ucdavis.edu/cmip6/30s/MIROC6/ssp585/wc2.1_30s_bioc_MIROC6_ssp585_2061-2080.tif"
download.file(P_url, destfile="ssp585_new1.tif")
unzip("wc2.1_2.5m_bio.zip")

library(rgeos)
library(rgdal)
library(sp)
library(raster)
tmax.crop <- stack("D:/Dahal/ssp585_new1.tif")
tmax.crop <- stack(tmax.crop)
tmax.crop <- mean(tmax.crop)
writeRaster(tmax.crop, "tmax.mean.tif")


