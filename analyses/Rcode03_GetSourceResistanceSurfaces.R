#################################################
########## Read env data file ################
#################################################
rm(list = ls())

combi.doable <- openxlsx::read.xlsx(here::here('data/raw-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))
coef <- c(0.25, 1, 4, 8, 16, 32)

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  
  for (c in c(1:k)) {
    
    hab.suit <- terra::rast(here::here(paste0('data/derived-data/outputSDM/', g,'.GroupID.', c, '/proj_CurrentEM/proj_CurrentEM_', g,'.GroupID.', c,'_ensemble.tif')))
    ly.tokeep <- names(hab.suit)[grep('EMwmean', names(hab.suit))]
    hab.suit <- terra::subset(x = hab.suit, ly.tokeep)
    hab.suit <- hab.suit/1000
    res <- 100 - 99*(1-exp(-coef*hab.suit))/(1-exp(-coef))
    names(res) <- paste0('tranfo.coef.', coef)
    terra::writeRaster(res, here::here(paste0('outputs/ResistanceSurfaces/ResistanceSurface_', g, '_GroupID_', c, '_TransfoCoef_', coef, '.tif')))
    
    # thre <- quantile(na.omit(terra::values(hab.suit)), probs = 0.9)
    thre <- 0.8
    srce <- hab.suit
    srce[srce >= thre] <- 1
    srce[srce < thre] <- 0 
    terra::writeRaster(res, here::here(paste0('outputs/SourceLayers/SourceLayer_', g, '_GroupID_', c, '.tif')))
    
    }
}