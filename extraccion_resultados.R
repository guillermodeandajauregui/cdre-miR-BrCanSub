library(igraph)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(graphlayouts)

lista_analyzed <- readRDS("List.Analyzed.Networks.RDS")
List_reducedEnrichment <- readRDS("List_reducedEnrichment.RDS")
List_reducedEnrichmentGroups <- readRDS("List_reducedEnrichmentGroups.RDS")
EnrichmentList <- readRDS("EnrichmentOutput.RDS")
lista_de_analisis <- lapply(X = lista_analyzed, FUN = function(i){get.data.frame(i, "vertices")})
allSubtypeAnalysis <- bind_rows(lista_de_analisis, .id = "subtype")

lista_analyzed %>% lapply(FUN = function(i){
  #components(i)$csize %>% max
  #E(i) %>% length
  #which(components(i)$csize!=1) %>% length
  components(i)$csize
}
)


nueva_lista_analisis <- lapply(X = lista_analyzed, FUN = function(i){
  V(i)$componente <- components(i)$membership
  return(i)
})

nueva_lista_analisis <- lapply(X = nueva_lista_analisis, FUN = function(i){
  get.data.frame(i, "vertices")
})

nueva_lista_analisis %>% lapply(FUN = function(i){
  #table(i[["type"]]) %>% print()
  #i %>% filter(grado != 0) %>% pull("type")  %>% table 
  #i %>% E() %>% length
  #i %>% filter()
  le_top <- i %>% pull(componente) %>% table %>% sort(decreasing = TRUE) %>% .[1] %>% names()
  i %>% filter(componente == le_top) %>% pull(type) %>% table
})


EnrichmentList %>% lapply(FUN = function(i){
  i %>% lapply(FUN = function(j){
    j %>% nrow
  })
})


List_reducedEnrichment


prueba <- as_tbl_graph(lista_analyzed$luma)

p.dopest <- ggraph(prueba) +
  geom_node_point() + 
  geom_edge_link()+
  #scale_edge_color_manual(values = c("lavenderblush4", "goldenrod1")) +
  #scale_edge_alpha_manual(values = c(0.5, 0.25)) +
  #geom_node_text(aes(label=name), repel=T)+
  #ggtitle(i) +
  theme_graph() 

p.dopest

########
allSubtypeAnalysis %>% filter(grado == 1, Redundancy == 0) %>% pull(subtype) %>% table
allSubtypeAnalysis %>% filter(grado == 0) %>% pull(subtype) %>% table

#####
#for network visualization of enrichment
mir_function_nw <-
lapply(List_reducedEnrichmentGroups, function(i){
  r1 <- i %>% bind_rows(.id = "miR")
  }) %>% bind_rows(.id = "subtype") %>% 
  dplyr::select(subtype, miR, rowname, my_groups) %>% 
  mutate(subtype_miR = glue::glue("{subtype}_{miR}"))

#write_delim(x = mir_function_nw, path = "for_presentation/mir_function_nw.tsv", delim = "\t", col_names = TRUE)  

#
# mir_function_nw %>% filter(subtype == "luma") %>% write_delim(path = "for_presentation/luma_dict.tsv", 
#                                                               delim = "\t", 
#                                                               col_names = TRUE
#                                                               )  
# 
# mir_function_nw %>% filter(subtype == "lumb") %>% write_delim(path = "for_presentation/lumb_dict.tsv", 
#                                                               delim = "\t", 
#                                                               col_names = TRUE
# )  
# 
# mir_function_nw %>% filter(subtype == "basal") %>% write_delim(path = "for_presentation/basl_dict.tsv", 
#                                                               delim = "\t", 
#                                                               col_names = TRUE
# )  

mir_function_nw %>% filter(subtype == "basal") %>% pull(var = my_groups) %>% table
mir_function_nw %>% filter(subtype == "luma") %>% pull(var = my_groups) %>% table
mir_function_nw %>% filter(subtype == "lumb") %>% pull(var = my_groups) %>% table

#####

#get another dictionary of functions across all conditions 

my_prefunction_list <- 
lapply(List_reducedEnrichmentGroups, FUN = function(subtype){
  lapply(subtype, FUN = function(mir){
    as_tibble(mir %>% dplyr::select(rowname, goterms, my_groups))
  })
})

my_prefunction_bind <-
lapply(my_prefunction_list, FUN = function(subtype){
  bind_rows(subtype, .id = "miR")
}) %>% bind_rows(, .id = "subtype")

my_prefunction_bind<-
my_prefunction_bind %>% 
  mutate(sub_mir_code = glue::glue("{subtype}_{miR}")) %>% 
  mutate(my_groups = str_pad(my_groups, width = 2, side = "left", pad = "0")) %>% 
  select(-c(subtype, miR)) %>% 
  spread(key = sub_mir_code, value = my_groups, fill = "00") %>% 
  unite(col = my_final, 3:8, sep = "", remove = FALSE)

my_prefunction_bind %>% write_delim(path = "for_presentation/dict_all.tsv", delim = "\t", quote_escape = "double")

my_prefunction_bind_preunbind <-
  lapply(my_prefunction_list, FUN = function(subtype){
    bind_rows(subtype, .id = "miR")
  }) %>% bind_rows(, .id = "subtype")

#get neighbors of each commodore 
allSubtypeAnalysis %>% head
my_cdres <- allSubtypeAnalysis %>% filter(grado >= 100, Redundancy <= 0.5) %>% pull(name) %>% unique
names(my_cdres) <- my_cdres
neighbors(lista_analyzed$luma, v = my_cdres[1]) %>% names

cdre_neighborhood_list<-
lapply(lista_analyzed, FUN = function(subtype){
  lapply(my_cdres, FUN = function(cdre){
    neighbors(subtype, v = cdre) %>% names
  })
})

cdre_neighborhood_list_unlisted<-unlist(cdre_neighborhood_list, recursive = FALSE)

neighborhood_jaccard <- 
sapply(cdre_neighborhood_list_unlisted, FUN = function(i){
  sapply(cdre_neighborhood_list_unlisted, FUN = function(j){
    the_intersection <- intersect(i,j)
    the_union <- union(i,j) 
    the_jaccard <- length(the_intersection)/length(the_union)
  })
})

max(neighborhood_jaccard[which(neighborhood_jaccard!=1)])

write_delim(as.data.frame(neighborhood_jaccard), path = "for_presentation/jaccard_neighborhoods.txt", delim = "\t")

allSubtypeAnalysis %>% filter(name%in%my_cdres) %>% arrange(name)
allSubtypeAnalysis %>% filter(name=="hsa-mir-190b")

my_prefunction_bind_preunbind %>% head
my_prefunction_bind_preunbind %>% group_by(subtype, miR) %>% summarise(count = n()) %>% xtable::xtable()

#systematic joining of group characteristic names and group membership, to get the count of GO terms by group  

#example
my_prefunction_list$luma$`hsa-mir-139`
left_join(x = my_prefunction_list$luma$`hsa-mir-139`, 
          y = List_reducedEnrichment$luma$`hsa-mir-139`, 
          by = "my_groups", 
          suffix = c("", ".exemplar")) %>% 
  group_by(my_groups, rowname.exemplar, goterms.exemplar) %>% 
  summarise(count = n())
  
my_subtypes_no <- seq_along(my_prefunction_list)
names(my_subtypes_no) <- names(my_prefunction_list)

my_functions_with_characteristic_term <- 
lapply(X = my_subtypes_no, FUN = function(subtype){
  
  mirs_by_subtype <- seq_along(my_prefunction_list[[subtype]])
  names(mirs_by_subtype) <- names(my_prefunction_list[[subtype]])
  lapply(mirs_by_subtype, FUN = function(mir){
    left_join(x = my_prefunction_list[[subtype]][[mir]],
              y = List_reducedEnrichment[[subtype]][[mir]],
              by = "my_groups",
              suffix = c("", ".exemplar")
              ) %>% 
      select(-"Adjusted.Pvalue")
    }) 
  })

df_functions_with_char_term_count <- 
lapply(my_functions_with_characteristic_term, FUN = function(subtype){
  lapply(subtype, FUN = function(mir){
    mir %>% group_by(my_groups, rowname.exemplar, goterms.exemplar) %>% 
      summarise(count = n())
  }) %>% bind_rows(.id = "mir")
}) %>% bind_rows(.id = "subtype")

my_subtypes <- unique(df_functions_with_char_term_count$subtype)
names(my_subtypes) <- my_subtypes
lapply(X = my_subtypes, FUN = function(subtipo){
  #print(subtipo)
  df_functions_with_char_term_count%>% 
    dplyr::filter(subtype == subtipo)
})

df_functions_with_char_term_count %>% head

df_functions_with_char_term_count %>% 
  rename("group number" = my_groups, 
         "characteristic GO-BP" = rowname.exemplar, 
         "characteristic GO-BP name" = goterms.exemplar 
         ) %>% 
  dplyr::ungroup() %>% group_by(subtype, mir) %>% 
#group_walk(~ write.csv(.x, file = paste0("testing", .y$subtype, .y$mir, ".csv")))
  group_walk(~ print(xtable::xtable(x = .x), 
                          file = paste0("for_presentation/latex_tables/", 
                                       .y$subtype, "_", .y$mir, 
                                       "_xlatex.txt"),
                     include.rownames = FALSE
                          )
             )

#write out as one big file
df_functions_with_char_term_count %>% 
  rename("group number" = my_groups, 
         "characteristic GO-BP" = rowname.exemplar, 
         "characteristic GO-BP name" = goterms.exemplar 
  ) %>% 
  dplyr::ungroup() %>% group_by(subtype, mir) %>% write_delim(path = "results/aggregated_functions_all.txt", delim = "\t")

allSubtypeAnalysis %>% filter(name%in%my_cdres) %>% arrange(name) %>% 
  select(-c(ClusteringRobbinsAlexander, ClusteringDot, type, id))

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

allSubtypeAnalysis %>% mutate(subtype = as_factor(subtype)) %>% pull(subtype)
  
allSubtypeAnalysis %>% filter(name%in%my_cdres) %>% 
  mutate(subtype = as_factor(subtype)) %>% 
  arrange(name, subtype) %>% 
  rename(degree = grado) %>% 
  select(-c(ClusteringRobbinsAlexander, ClusteringDot, type, id)) %>% 
  clipr::write_clip()

#las viejas buenas distribuciones de grado

my_network <- lista_analyzed$luma
my_network %>% get.data.frame("vertices") %>% head
my_network %>% 
  get.data.frame("vertices") %>% 
  dplyr::rename(degree = grado) %>% 
  dplyr::filter(type==TRUE) %>% 
  group_by(degree) %>% 
  tally() %>% 
  mutate(resultado = 1 - (cumsum(n)/sum(n)),
         mocos = sum(n))
  
  
  pull(degree) %>% 
  table %>% 
  as.data.frame %>% 
  head

my_dgree_plot <- function(my_network){
  
  #mir
  #dg.top <- my_network %>% get.data.frame("vertices") %>% dplyr::filter(type == TRUE)  %>% pull(grado) %>% table %>% as.data.frame
  #gene
  #dg.bot <- my_network %>% get.data.frame("vertices") %>% dplyr::filter(type == FALSE) %>% pull(grado) %>% table %>% as.data.frame 
}

my_dgree_plot <- function(my_network){
  
  #mir
  #dg.top <- my_network %>% get.data.frame("vertices") %>% dplyr::filter(type == TRUE)  %>% pull(grado) %>% table %>% as.data.frame
  #gene
  #dg.bot <- my_network %>% get.data.frame("vertices") %>% dplyr::filter(type == FALSE) %>% pull(grado) %>% table %>% as.data.frame 
  #colnames(dg.top)[1] <- "grado"
  #colnames(dg.bot)[1] <- "grado"
  
  dg.top <- get.data.frame(my_network, "vertices")%>%dplyr::filter(type==TRUE)%>%
    group_by(grado)%>%tally()%>%
    mutate(resultado = 1 - (cumsum(n)/sum(n)))%>%
    mutate(masaca = c(1, resultado[1:n()-1]))
  
  dg.bot <- get.data.frame(my_network, "vertices")%>%dplyr::filter(type==FALSE)%>%
    group_by(grado)%>%tally()%>%
    mutate(resultado = 1 - (cumsum(n)/sum(n)))%>%
    mutate(masaca = c(1, resultado[1:n()-1]))
  
  
  dg.list = list(mir  = dg.top, 
                 gene = dg.bot)
  
  bind_rows(dg.list, .id = "type") %>% 
    mutate(type = as_factor(type)) %>% 
    group_by(type)
  
  
  
  #dg.all <- full_join(x = dg.top, 
  #                    y = dg.bot, by = "grado", 
  #                    suffix = c(".top", ".bottom")) %>% 
  #  arrange(grado)
  
  #return(dg.all)
  # p <-
  # dg.all %>%
  #   ggplot(mapping = aes(grado, masaca.top, color = "red")) +
  #   geom_line() +
  #   geom_point() +
  #   theme(text = element_text(size=20)) +
  #   ylab(label = "P(K)") +
  #   xlab(label = "K") +
  #   geom_line(aes(grado, masaca.bottom, color = "blue")) +
  #   geom_point(aes(grado, masaca.bottom, color = "blue")) +
  #   theme_minimal() +
  #   theme(legend.position="none")
  # return(p)
}

my_dgree_list <- lapply(X = lista_analyzed, FUN = my_dgree_plot)

my_grouped_dg <-
my_dgree_list %>% bind_rows(.id = "subtype") %>% 
  mutate(subtype = as_factor(subtype)) %>% 
  group_by(subtype, type)

my_grouped_dg %>% 
  dplyr::filter(type == "mir") %>% 
  ggplot(aes(x = grado, y = masaca, color = type)) +
  xlab(label = "K") +
  ylab(label = "1 - P(K)") + 
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  scale_color_manual(values = c("blueviolet", "darkgoldenrod1")) +
  facet_wrap(. ~ subtype) + 
  theme_minimal()



my_grouped_dg %>% 
  dplyr::filter(type == "gene") %>% 
  ggplot(aes(x = grado, y = masaca, color = type)) +
  xlab(label = "K") +
  ylab(label = "1 - P(K)") + 
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  scale_color_manual(values = c("darkgoldenrod1")) +
  facet_wrap(. ~ subtype) + 
  theme_minimal()

#### write out enrichment 

lapply(X = EnrichmentList, FUN = function(subtype){
  bind_rows(subtype, .id = "mir")
  }) %>% bind_rows(.id = "subtype") %>% write_delim(path = "results/enrichment_results_all.txt", delim = "\t")


#### write out big groups representative