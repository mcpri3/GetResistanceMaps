
combi.doable <- openxlsx::read.xlsx(here::here('data/raw-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))
lst.param <- data.frame()

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]

  lst.param <- rbind(lst.param , data.frame(gp = g, clus = c(1:k)))
}
lst.param <- as.matrix(lst.param)
write.table(lst.param, here::here('data/derived-data/inputSDM/list-of-params-for-batchrun.txt'), row.names = F, col.names = F)
