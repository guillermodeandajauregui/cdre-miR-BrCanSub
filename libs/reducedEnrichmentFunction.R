#Enrichment Reduction function
#using the GOSemSim package

library(GOSemSim)



reducedEnrichment <- function(rich.df, 
                              GO.ID.row   = "rowname",
                              sim.measr   = "Wang", 
                              SemData     =  godata('org.Hs.eg.db', ont = "BP"), #hsGO <- godata('org.Hs.eg.db', ont = "BP")
                              dist.mthd   = "euclidean",
                              clustr.mthd = "ward.D2",
                              no.groups   = 10,
                              cluster.nom = "my_groups"
){
  
  #get the GO IDs of the interesting GO terms 
  my_terms <- rich.df %>% dplyr::pull(GO.ID.row)
  
  #make a similarity matrix using GO sem sim
  simatrix <- GOSemSim::mgoSim(GO1 = my_terms, 
                               GO2 = my_terms, 
                               semData = SemData, 
                               measure = sim.measr, 
                               combine = NULL)
  
  #make a distance structure
  dist.struc = dist(x = simatrix, method = dist.mthd)
  
  #cluster 
  my_clust <- hclust(d = dist.struc, method = clustr.mthd)
  
  #get N groups
  my_clustered <- cutree(my_clust, k = no.groups) %>% as.data.frame 
  
  #reshape the my_clustered data frame
  
  my_clustered <- my_clustered %>% 
    rownames_to_column() %>% 
    rename(. = "slash") %>% #that dot thing is unfortunate; but this way I don't lose my "Slash" easter egg
    rename(slash = cluster.nom) %>% 
    arrange(!! rlang::sym(cluster.nom))
  
  #join with original data 
  
  my_join <- left_join(x = rich.df, 
                       y = my_clustered
  )
  
  #group by my_clustered 
  
  my_join <- 
    my_join %>% 
    group_by(!! rlang::sym(cluster.nom)) 
  
  #keep the row with the lowest adjusted p.value 
  my_join <- 
    my_join %>%
    filter(Adjusted.Pvalue == min(Adjusted.Pvalue)) %>% 
    filter(stringr::str_length(goterms) == min(stringr::str_length(goterms)) ) #to avoid repeated values
    
  
  #return the most informative columns
  
  my_join <- 
    my_join %>% 
    dplyr::select(rowname, goterms, Adjusted.Pvalue, !! rlang::sym(cluster.nom))
  
  #return
  return(my_join)
}


#####
#Get groups for all enriched values
#####

library(GOSemSim)



reducedEnrichment_groups <- function(rich.df, 
                              GO.ID.row   = "rowname",
                              sim.measr   = "Wang", 
                              SemData     =  godata('org.Hs.eg.db', ont = "BP"), #hsGO <- godata('org.Hs.eg.db', ont = "BP")
                              dist.mthd   = "euclidean",
                              clustr.mthd = "ward.D2",
                              no.groups   = 10,
                              cluster.nom = "my_groups"
){
  
  #get the GO IDs of the interesting GO terms 
  my_terms <- rich.df %>% dplyr::pull(GO.ID.row)
  
  #make a similarity matrix using GO sem sim
  simatrix <- GOSemSim::mgoSim(GO1 = my_terms, 
                               GO2 = my_terms, 
                               semData = SemData, 
                               measure = sim.measr, 
                               combine = NULL)
  
  #make a distance structure
  dist.struc = dist(x = simatrix, method = dist.mthd)
  
  #cluster 
  my_clust <- hclust(d = dist.struc, method = clustr.mthd)
  
  #get N groups
  my_clustered <- cutree(my_clust, k = no.groups) %>% as.data.frame 
  
  #reshape the my_clustered data frame
  
  my_clustered <- my_clustered %>% 
    rownames_to_column() %>% 
    rename(. = "slash") %>% #that dot thing is unfortunate; but this way I don't lose my "Slash" easter egg
    rename(slash = cluster.nom) %>% 
    arrange(!! rlang::sym(cluster.nom))
  
  #join with original data 
  
  my_join <- left_join(x = rich.df, 
                       y = my_clustered
  )
  return(my_join)
}
