#SDM with BIOMOD2 Edited by Vineet K Dubey

# load the required packages
remove.packages("biomod2")
devtools::install_github('biomodhub/biomod2')
library(biomod2)
install.packages("corrplot")
library(corrplot)
install.packages("raster")
library(raster)
install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)
install.packages("sp")
library(sp)
install.packages("rgdal")
library(rgdal)
install.packages("dismo")
library(dismo)
install.packages("dplyr")
library(dplyr)
install.packages("sdm")
library(sdm)
install.packages("basemaps")
library("basemaps")
install.packages("magick")
library(magick)
install.packages("tidyr") 
library(tidyr) 
install.packages("mapview")
library(mapview)
install.packages("terra") 
library(terra)

# setwd
#Ailurus fulgens
gbif("Ailurus","fulgens",download= F)
sp <- gbif("Ailurus","fulgens",download= T)
class(sp)
dim(sp)
write.csv(sp, "pandagbif.csv")
table(sp$basisOfRecord)
sp <- sp %>% 
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION","OBSERVATION"))
nrow(sp)

spg <- sp %>% select(species,lon,lat,)
head(spg)
#spg$species <- 1
spg <- spg %>% drop_na()
nrow(spg)
class(spg)
coordinates(spg) <- c('lon','lat')
write.csv(spg, "ailurus.csv")
#------------


# Load species occurrences 

panda <- read.csv('ailurus.csv')
head(panda)
########################## Select the name of the studied species
myRespName <- 'Ailurus.fulgens'

############################ Get corresponding presence/absence data
myResp <- as.numeric(panda[, myRespName])

############################# Get corresponding XY coordinates
myRespXY <- panda[, c('long', 'lat')]

# Load environmental variables

##################
# download the bioclim data:
bio <- list.files("E:/woldclim/wc10", full.names=T)
bio <- raster::getData('worldclim', var='bio', res=10)
bio
names(bio)
windows()
plot(bio[[1]])  
points(spg)
e <- drawExtent()
spg <- crop(spg, e)
bioc <- crop(bio, e)
plot(bioc[[1]])
points(spg,col='red')

#------------
ex <- raster::extract(bioc,spg)
head(ex)
#######################
corr_enviroscaled = cor(ex)
corrplot(corr_enviroscaled, method = 'number')
predictors <- stack(predictors)
predictors <-  dropLayer(predictors,c("bio7","bio10","bio11","bio16","bio17","bio14", "bio15", "bio18", "bio19" ))

#####################################################################################
myExpl <-stack(predictors)
#Prepare data & parameters
#Format data (observations & explanatory variables)
# # Transform true absences into potential pseudo-absences
myResp.PA <- ifelse(myResp == 1, 1, NA)
# 
#Format Data with pseudo-absences : random method
myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.nb.rep = 4,
                                       PA.nb.absences = 1000,
                                       PA.strategy = 'random')

myBiomodData.r
plot(myBiomodData.r)
#Parameterize modeling options
# Print default modeling options
bm_DefaultModelingOptions()

# Create default modeling options
myBiomodOptions <- BIOMOD_ModelingOptions()
myBiomodOptions

#pCluster(cl)

myBiomodCV <- BIOMOD_CrossValidation(bm.format = myBiomodData.r)
head(myBiomodCV)



# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData.r,
                                    bm.options = myBiomodOptions,
                                    modeling.id = 'AllModels',
                                    nb.rep = 2,
                                    data.split.perc = 80,
                                    data.split.table = myBiomodCV,
                                    var.import = 3,
                                    metric.eval = c('TSS','ROC'),
                                    do.full.models = FALSE)
# seed.val = 123)
# nb.cpu = 8)
myBiomodModelOut


# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
get_variables_importance(myBiomodModelOut)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodModelOut)

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'median')


# Model ensemble models
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      var.import = 3,
                                      metric.eval = c('TSS', 'ROC'),
                                      prob.mean = TRUE,
                                      prob.median = TRUE,
                                      prob.cv = TRUE,
                                      prob.ci = TRUE,
                                      prob.ci.alpha = 0.05,
                                      committee.averaging = TRUE,
                                      prob.mean.weight = TRUE,
                                      prob.mean.weight.decay = 'proportional')
myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)


# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'median')


# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj)
windows()
current <- get_predictions(myBiomodProj)
names(current)
RFstack <- raster::subset(current, grep('RF', names(current), value = TRUE))



# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'CurrentEM',
                                             new.env = myExpl,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
myBiomodEMProj
plot(myBiomodEMProj)
windows()

current1 <- raster("E:/aconyx/Ailurus.fulgens/proj_Current/proj_Current_Ailurus.fulgens.tif")
current <- raster("E:/aconyx/Ailurus.fulgens/proj_Current/proj_Current_Ailurus.fulgens_ROCbin.tif")
windows()
plot(current)






#Compare range sizes
# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data("bioclim_future")
setwd("E:/aconyx")
options(timeout = 600)
biof <- raster::getData('CMIP5', var='bio', res=10, rcp=85, model='BC', year=70)
getData('worldclim', var='tmin', res=0.5, lon=5, lat=45)
names(biof)
biocf <- crop(biof, e)
predictors <- stack(biof)
myExplFuture <-  dropLayer(predictors,c("bio7","bio10","bio11","bio16","bio17","bio14", "bio15", "bio18", "bio19" ))

#myExplFuture = rast(predictors)
#data("bioclim_future")
myExplFuture = rast(bioclim_future)
# Project onto future conditions
myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                              proj.name = 'Future',
                                              new.env = myExplFuture,
                                              models.chosen = 'all',
                                              metric.binary = 'TSS',
                                              build.clamping.mask = TRUE)

# Load current and future binary projections
CurrentProj <- rast("GuloGulo/proj_Current/proj_Current_GuloGulo_TSSbin.grd")
FutureProj <- rast("GuloGulo/proj_Future/proj_Future_GuloGulo_TSSbin.grd")

# Compute differences
myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj)

myBiomodRangeSize$Compt.By.Models
plot(myBiomodRangeSize$Diff.By.Pixel)

# Represent main results 
gg = bm_PlotRangeSize(bm.range = myBiomodRangeSize, 
                      do.count = TRUE,
                      do.perc = TRUE,
                      do.maps = TRUE,
                      do.mean = TRUE,
                      do.plot = TRUE,
                      row.names = c("Species", "Dataset", "Run", "Algo"))



