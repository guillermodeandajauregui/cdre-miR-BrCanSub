FilteredList$basal$`hsa-mir-136` %>% pull(rowname)

FilteredList$basal$`hsa-mir-136` %>% head

FilteredList$basal$`hsa-mir-136`$category %>% table
library(GOSemSim)

hsGO <- godata('org.Hs.eg.db', ont = "BP")
hsGO

my_try <- FilteredList$basal$`hsa-mir-136` %>%  pull(rowname)

my_test2 <- 
GOSemSim::mgoSim(GO1 = my_try, 
                 GO2 = my_try, 
                 semData = hsGO, 
                 measure = "Wang", 
                 combine = NULL)


apply(X = my_test2, MARGIN = 1, FUN = function(i){
  r1 <- i[which(i != 1)]
  r2 <- max(r1)
})

#make a dist
dist(my_test2)

#make a clust

my_clust <- hclust(d = dist(my_test2), "ward.D2")
my_clustered <- cutree(my_clust, k = 10) %>% as.data.frame 
my_clustered <- my_clustered %>% rownames_to_column() %>% rename(. = "slash") %>% arrange(slash)

my_ressy <- FilteredList$basal$`hsa-mir-136`

my_joiny <- 
left_join(x = FilteredList$basal$`hsa-mir-136`, 
          y = my_clustered)

?summarize
my_joiny %>% group_by(slash) %>% filter(Adjusted.Pvalue == min(Adjusted.Pvalue)) %>% dplyr::select(rowname, goterms, Adjusted.Pvalue, slash)


my_try_other <- FilteredList$basal$`hsa-mir-139` %>%  pull(rowname)

my_test_other <- 
  GOSemSim::mgoSim(GO1 = my_try_other, 
                   GO2 = my_try_other, 
                   semData = hsGO, 
                   measure = "Wang", 
                   combine = NULL)


my_clust_other <- hclust(d = dist(my_test_other), "ward.D2")
my_clustered_other <- cutree(my_clust_other, k = 10) %>% as.data.frame 
my_clustered_other <- my_clustered_other %>% rownames_to_column() %>% rename(. = "slash") %>% arrange(slash)

my_ressy_other <- FilteredList$basal$`hsa-mir-139`

my_joiny_other <- 
  left_join(x = FilteredList$basal$`hsa-mir-139`, 
            y = my_clustered_other)

?summarize
my_joiny_other %>% group_by(slash) %>% filter(Adjusted.Pvalue == min(Adjusted.Pvalue)) %>% dplyr::select(rowname, goterms, Adjusted.Pvalue, slash)


my_first_reduced <- my_joiny %>% group_by(slash) %>% filter(Adjusted.Pvalue == min(Adjusted.Pvalue)) %>% dplyr::select(rowname, goterms, Adjusted.Pvalue, slash)
my_second_reduced <- my_joiny_other %>% group_by(slash) %>% filter(Adjusted.Pvalue == min(Adjusted.Pvalue)) %>% dplyr::select(rowname, goterms, Adjusted.Pvalue, slash)

