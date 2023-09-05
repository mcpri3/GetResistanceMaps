#################################################
########## Prep. occurrence file ################
#################################################
# install.packages('biomod2')

combi.doable <- openxlsx::read.xlsx(here::here('data/raw-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))
grid.fr <- sf::st_read(here::here('data/raw-data/Grids/ReferenceGrid_France_bin_1000m.gpkg'))
occur <- readRDS(here::here('data/raw-data/SNAP-Vertebrate-Species-GBIF-INPN-IUCNOccurrenceData_France_Res1000m_2010-2020'))

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  m <- combi.doable$mode[i]
  
# Read list of SNAP species 
lst.sp <- openxlsx::read.xlsx(here::here(paste0('data/raw-data/FunctionalGroups/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_', 
g, '_GroupID_K=', k, '_Weight-', m,'.xlsx')))

for (c in unique(lst.sp$cluster.id)) {
  
  sub.lst <- lst.sp[lst.sp$cluster.id %in% c, ]
  sp <- gsub(' ', '_', sub.lst$LB_NOM_VALIDE_SPE_LEVEL_SYN)
  occur.c <- occur[, sp] 
  if (class(occur.c) == 'numeric') {
    grid.fr$Nocc <- occur.c
  } else {
    occur.c <- apply(occur.c, 1, sum)
    occur.c[occur.c>1] <- 1
    grid.fr$Nocc <- occur.c
    }
  
  grid.fr2 <- grid.fr[grid.fr$Nocc > 0, ]
  grid.fr2 <- grid.fr2[, 'Nocc']
  grid.fr2 <- sf::st_centroid(grid.fr2)
  print(nrow(grid.fr2))
  sf::st_write(grid.fr2, here::here(paste0('data/derived-data/inputSDM/OccurrenceData_France_', g, '_GroupID_K=', c, '_Res1000_2010-2020.gpkg')), driver = 'GPKG', delete_layer = T)
  
}
}

#################################################
########## Prep. env data file ################
#################################################

grid.fr <- sf::st_read(here::here('data/raw-data/Grids/ReferenceGrid_France_bin_1000m.gpkg'))
env <- readRDS(here::here('data/raw-data/EnvironmentalVariables_France_Res1000m'))
env <- env[, colnames(env) != 'climatic.prsd']
env <- env[, colnames(env) != 'climatic.cdd']
env <- env[, colnames(env) != 'climatic.fd']
env <- env[, colnames(env) != 'topo.SDalti']
env <- env[, colnames(env) != 'climatic.gdd5']
env <- env[, colnames(env) != 'climatic.swe']
env <- env[, colnames(env) != 'climatic.scd']
env <- env[, colnames(env) != 'land.compo.LandSystem']
# tokeep <- env[, c("land.compo.LandSystem", "land.compo.WAW_prop")]
# env <- env[, - grep('land.compo', colnames(env))]
# env <- cbind(env, tokeep)

grid.fr <- cbind(grid.fr, env)
sf::st_write(grid.fr, here::here('data/derived-data/inputSDM/EnvironmentalVariables_France_Res1000m.gpkg'), driver = 'GPKG', delete_layer = T)

# create raster stack
rr <- terra::rast(here::here('data/raw-data/Grids/ReferenceGrid_Europe_bin_1000m.tif'))
france <- sf::st_read(here::here('data/raw-data/Grids/CasestudyOutlines.gpkg'), layer = 'France')
france <- sf::st_transform(france, sf::st_crs(rr))
rr <- terra::crop(rr, france)

todo <- colnames(grid.fr)
todo <- todo[grep('.', todo, fixed = T)]
for (col in todo) {
  print(col)
rr <- terra::rasterize(grid.fr, rr, field = col)
assign(paste0('rr.', col), rr)
}

lys <- ls()[grep('rr.', ls(), fixed = T)]
assign('rr.stack', mget(x = lys))
rr.stack <- terra::rast(rr.stack)
terra::writeRaster(rr.stack, 'data/derived-data/inputSDM/EnvironmentalVariables_France_Res1000m.tif', overwrite = T)
