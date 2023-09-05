#################################################
########## Read env data file ################
#################################################
rm(list = ls())

env <- terra::rast(here::here('data/derived-data/inputSDM/EnvironmentalVariables_France_Res1000m.tif'))
# env <- terra::subset(env, names(env)[1:2])
names(env) <- gsub('rr.', '', names(env), fixed = T)

combi.doable <- openxlsx::read.xlsx(here::here('data/raw-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))

# Params for pseudo-absence generation
PA.strat <- 'random'
disk.min <- 1000
disk.max <- 50000
nrep.PA <- 3 

# Params for modelling
lst.mod <- c("XGBOOST", "RF", "MAXNET", "ANN")
# lst.mod <- c("MAXNET", "MAXENT")
nCrossVal <- 5
cv.perc <- 0.8
nperm.var.imp <- 3
cpu = 8
ens.calc <- c('EMca', 'EMwmean')

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]

  for (c in c(1:k)) {
    
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
                                       dir.name = './data/derived-data/outputSDM',
                                       filter.raster = T)
  saveRDS(myBiomodData, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodData')))
  
  # biomod2::plot(myBiomodData)
  plot.new()
  corrplot::corrplot(cor(myBiomodData@data.env.var))
  
  myBiomodOptions <- biomod2::BIOMOD_ModelingOptions()
  
  myBiomodOptions@XGBOOST$nthread <- cpu
  myBiomodOptions@XGBOOST$max.depth <- 3
  myBiomodOptions@RF$nodesize <- round(nrow(occur)/10)
  # myBiomodOptions@RF$mtry <- 0.2
  
  tic = Sys.time()
  myBiomodModelOut <- try(biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
                                               modeling.id = 'AllModels',
                                               models = lst.mod,
                                               bm.options = myBiomodOptions,
                                               CV.strategy = 'random',
                                               CV.nb.rep = nCrossVal,
                                               CV.perc = cv.perc,
                                               metric.eval = c('TSS','ROC'), 
                                               var.import = nperm.var.imp, 
                                               nb.cpu = cpu))
                          
  Sys.time() - tic
  eval <- myBiomodModelOut@models.evaluation@val
  eval <- eval[eval$metric.eval == 'TSS',]
  saveRDS(myBiomodModelOut, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodModelOut')))

  # Get evaluation scores & variables importance
  # biomod2::get_evaluations(myBiomodModelOut)
  # biomod2::get_variables_importance(myBiomodModelOut)
  # head(biomod2::get_predictions(myBiomodModelOut))
  # head(biomod2::get_predictions(myBiomodModelOut))
  
  # Represent evaluation scores
  biomod2::bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
  biomod2::bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
  biomod2::bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'calibration', group.by = c('algo', 'algo'))
  biomod2::bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = 'validation', group.by = c('algo', 'algo'))

  # # Represent variables importance
  biomod2::bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))

  # # Represent response curves
  mods <- biomod2::get_built_models(myBiomodModelOut)
  biomod2::bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                        models.chosen = mods,
                        fixed.var = 'median') #fixed.var = valeur des autres variables fixées 

  # # Project single models
  # myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
  #                                   proj.name = 'Current',
  #                                   new.env = env,
  #                                   models.chosen = 'all',
  #                                   build.clamping.mask = TRUE,
  #                                   nb.cpu = cpu)
  # biomod2::plot(myBiomodProj)
  
  
  # Ensemble model
  tic = Sys.time()
  
  myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                 models.chosen = 'all',
                                                 em.by = 'all',
                                                 em.algo = ens.calc ,
                                                 metric.select = c('TSS'),
                                                 metric.select.thresh = c(0.4),
                                                 metric.eval = c('TSS', 'ROC'),
                                                 var.import = nperm.var.imp, 
                                                 nb.cpu = cpu)
  Sys.time() - tic
  
  saveRDS(myBiomodEM, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodEM')))
  
  # # Project ensemble models (from single projections)
  # myBiomodEMProj <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
  #                                              bm.proj = myBiomodProj,
  #                                              models.chosen = 'all',
  #                                              metric.binary = 'all',
  #                                              metric.filter = 'all', 
  #                                              nb.cpu = cpu)
  
  # Project ensemble models (building single projections)
  myBiomodEMProj <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                               proj.name = 'CurrentEM',
                                               new.env = env,
                                               models.chosen = 'all',
                                               metric.binary = 'all',
                                               metric.filter = 'all', 
                                               nb.cpu = cpu)
  myBiomodEMProj
  biomod2::plot(myBiomodEMProj)
  saveRDS(myBiomodEMProj, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodEMProj')))
  
  ################################
  # Get SDM summary
  ################################
  
  rmarkdown::render(here::here('outputs/HabSuitModTest.Rmd'), output_format = 'pdf_document', 
                      output_file = sprintf("HabSuitMod-for_%s_K=%s.pdf", g, c), 
                      params = list(group = g, Nclus = k, mode = combi.doable$mode[i], clus.id = c))
  
  }
}



################################
# Get SDM summary
################################
combi.doable <- openxlsx::read.xlsx(here::here('data/raw-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  
  for (c in c(1:k)) {
  rmarkdown::render(here::here('outputs/HabSuitModTest.Rmd'), output_format = 'pdf_document', 
                    output_file = sprintf("HabSuitMod-for_%s_K=%s.pdf", g, c), 
                    params = list(group = g, Nclus = k, mode = combi.doable$mode[i], clus.id = c))
  }
}
