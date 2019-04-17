#libraries
library(tidyverse)
library(igraph)

#read data

networkList <- readRDS(file = "List.Analyzed.Networks.RDS")
networkDF   <- readRDS(file = "allSubtype.DF.RDS")

GO_data     <- readRDS(file = "GO_symbol.RDS")
GO_dict     <- readRDS(file = "goterms_dict.RDS")

#set my universe 
my_universe <- allSubtypeAnalysis %>% filter(type == FALSE) %>% pull(name)

# enrichment function 

enrich.neighbors.cdre <- function(nw){
  #get data frame
  analysis_df <- igraph::get.data.frame(nw, "vertices")
  
  #get cdre nodes
  my_cdre <- analysis_df %>% 
    dplyr::filter(type==TRUE,
           grado >= 100, 
           Redundancy <= 0.5) %>% 
    dplyr::pull(name)
  
  names(my_cdre) <- my_cdre
  
  #get neighbors 
  cdre_neighbors <- lapply(X = my_cdre, FUN = function(i){
    neighbors(nw, v = i) %>% names
  })
  
  #enrich
  
   my_results <-
   lapply(cdre_neighbors, FUN = function(j){
     r1 <- HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = GO_ALL_SYMBOL,
                                    universe = my_universe,
                                    hits = j,
                                    minGeneSetSize = 10,
                                    pAdjustMethod = "BH"
     )
     r2 <- r1 %>% 
       as.data.frame %>% rownames_to_column
     
     r2 <- left_join(x = r2, y = GO_dict, by = c("rowname" = "GO_ID"))
     
     return(r2)
       
   })
   
   return(my_results)
}

#Filter Function

FilterMyList <- function(my_list){
lapply(my_list, FUN = function(i){
  my_idx <- seq_along(i)
  names(my_idx) <- names(i)
  
  lapply(X = my_idx, FUN = function(j){
    i[[j]] %>% 
      filter(Adjusted.Pvalue < 1e-3, 
             category == "BP"
      ) 
  })
  
    })
  }


##########################################
#enrich by subtype 

EnrichmentList <- lapply(X = networkList, FUN = function(my_nw){
  enrich.neighbors.cdre(my_nw)
})

saveRDS(EnrichmentList, file = "EnrichmentOutput.RDS")

FilteredList <- FilterMyList(EnrichmentList)

FilteredList$basal$`hsa-mir-136` %>% pull(rowname)
