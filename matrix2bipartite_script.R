library(data.table)
library(tidyverse)
library(igraph)

my_files <- list.files("results", full.names = TRUE, pattern = "mi.tsv")[c(3,4,1,2)]
my_names <- c("luma", "lumb", "basal", "her2")


mi_matrices <- 
lapply(X = my_files, 
       FUN = function(i){
            fread(i)
            }
       )

names(mi_matrices) <- my_names

#filter by 0.9999 
#heuristic; 0.9990 to 0.9995 too lax 

my_quantile <- 0.9999

my_matrices_filtered <-
lapply(X = mi_matrices, FUN = function(i){
  r1 <- as.matrix(i[,-1])
  rownames(r1) <- i[[1]]
  my_quant <- quantile(r1, my_quantile)
  print(my_quant)
  r2 <- ifelse(r1 < my_quant, 0, r1)
})

#convert to networks 

my_g <-
  lapply(X = my_matrices_filtered, FUN = function(i){
    i %>% graph_from_incidence_matrix()
  })


lapply(X = my_g, FUN = function(i){
  components(i)$csize %>% table
})

#export for bipartite analyses in networkx 

sapply(1:4, FUN = function(i){
  write.graph(graph = my_g[[i]], 
              file = paste0("results/", my_names[i], "_9999.graphml"),
              format = "graphml")
})

