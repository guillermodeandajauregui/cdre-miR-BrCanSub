#libraries
library(tidyverse)
library(igraph)
#read data

my_analyzed <- list.files(path = "results", pattern = "analyzed.graphml", full.names = TRUE)[c(3,4,1,2)]
names(my_analyzed) <- c("luma", "lumb", "basal", "her2")

lista_analyzed <- lapply(X = my_analyzed, FUN = function(i){
  read_graph(file = i, format = "graphml")
})

lista_de_analisis <- lapply(X = lista_analyzed, FUN = function(i){get.data.frame(i, "vertices")})
allSubtypeAnalysis <- bind_rows(lista_de_analisis, .id = "subtype")

#save as RDS

saveRDS(object = lista_analyzed, "List.Analyzed.Networks.RDS")
saveRDS(object = allSubtypeAnalysis, "allSubtype.DF.RDS")
