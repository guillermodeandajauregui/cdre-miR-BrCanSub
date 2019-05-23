
nodos_i <- 
get.data.frame(lista_analyzed$luma, "vertices") %>% 
  filter(type==TRUE,
         grado >= 100, 
         Redundancy <= 0.5) %>% 
  pull(name)

neighbors(graph = lista_analyzed$luma, v = nodos_i[1])  %>% names %>% datapasta::dmdclip()
neighbors(graph = lista_analyzed$luma, v = nodos_i[2])  %>% names


probando <- 
allSubtypeAnalysis %>% 
  filter(type==TRUE,
         grado >= 100, 
         Redundancy <= 0.5, 
         subtype == "luma") %>% 
  pull(name) 

probando1 <-
neighbors(lista_analyzed$luma, v = probando[2]) %>% names

probando1
my_universe <- allSubtypeAnalysis %>% filter(type == FALSE) %>% pull(name)


probando3 <- 
HTSanalyzeR::multiHyperGeoTest(collectionOfGeneSets = GO_ALL_SYMBOL, 
                               universe = my_universe, 
                               hits = probando1, 
                               minGeneSetSize = 10, 
                               pAdjustMethod = "BH"
                               )


probando3 %>% as.data.frame %>% rownames_to_column %>%  dplyr::arrange(Adjusted.Pvalue) %>% head %>%  
  left_join(y = goterms, by = c("rowname" = "GO_ID"))
