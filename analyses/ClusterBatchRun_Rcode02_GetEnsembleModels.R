
#################################################
########## Read env data file ################
#################################################
rm(list = ls())
setwd('/bettik/primam/03_GetResistanceMaps')
#install.packages('biomod2', lib = '../local_lib', repos = "http://cran.us.r-project.org")
#install.packages('glmnet', lib = '../local_lib', repos = 'http://cran.us.r-project.org')
library(glmnet, lib.loc = "../local_lib")
library(biomod2, lib.loc = '../local_lib')

env <- terra::rast(here::here('data/derived-data/inputSDM/EnvironmentalVariables_France_Res1000m.tif'))
names(env) <- gsub('rr.', '', names(env), fixed = T)

# Params for pseudo-absence generation
PA.strat <- 'random'
disk.min <- 1000
disk.max <- 50000
nrep.PA <- 3 

# Params for modelling
lst.mod <- c("XGBOOST", "RF", "ANN",  "MAXNET")
nCrossVal <-5
cv.perc <- 0.8
nperm.var.imp <- 3
ens.calc <- c('EMca', 'EMwmean')

args <- commandArgs(trailingOnly=TRUE) ## get group parameters
g <- as.character(args[1]) 
c <- as.numeric(args[2]) 

occur <- terra::vect(here::here(paste0('data/derived-data/inputSDM/OccurrenceData_France_', g, '_GroupID_K=', c, '_Res1000_2010-2020.gpkg')))
n.abs <- ifelse(nrow(occur) < 5000, 5000, nrow(occur))
n.abs <- ifelse(nrow(occur) < 500, 500, n.abs)
nrep.PA <- ifelse(nrow(occur) < 500, 50, nrep.PA)

# Format Data with pseudo-absences 
myBiomodData <- biomod2::BIOMOD_FormatingData(resp.name = paste0(g, '_GroupID_', c),
                                              resp.var = occur,
                                              expl.var = env,
                                              PA.nb.rep = nrep.PA,
                                              PA.strategy = PA.strat,
                                              PA.dist.min = disk.min,
                                              PA.dist.max = disk.max,
                                              PA.nb.absences = n.abs,
                                              dir.name = 'data/derived-data/outputSDM',
                                              filter.raster = T)

myBiomodOptions <- biomod2::BIOMOD_ModelingOptions()
myBiomodOptions@XGBOOST$max.depth <- 3
myBiomodOptions@RF$nodesize <- round(nrow(occur)/10)

myBiomodModelOut <- try(biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
                                                 modeling.id = 'AllModels',
                                                 models = lst.mod,
                                                 bm.options = myBiomodOptions,
                                                 CV.strategy = 'random',
                                                 CV.nb.rep = nCrossVal,
                                                 CV.perc = cv.perc,
                                                 metric.eval = c('TSS','ROC'), 
                                                 var.import = nperm.var.imp))
saveRDS(myBiomodData, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodData')))
saveRDS(myBiomodModelOut, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodModelOut')))


# Ensemble model
myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                               models.chosen = 'all',
                                               em.by = 'all',
                                               em.algo = ens.calc ,
                                               metric.select = c('TSS'),
                                               metric.select.thresh = c(0.4),
                                               metric.eval = c('TSS', 'ROC'),
                                               var.import = nperm.var.imp)
saveRDS(myBiomodEM, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodEM')))

# Project ensemble models (building single projections)
myBiomodEMProj <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                      proj.name = 'CurrentEM',
                                                      new.env = env,
                                                      models.chosen = 'all',
                                                      metric.binary = 'all',
                                                      metric.filter = 'all')
saveRDS(myBiomodEMProj, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodEMProj')))
