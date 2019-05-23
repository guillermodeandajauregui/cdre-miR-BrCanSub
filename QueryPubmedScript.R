library(rentrez)
#count PMIDs

testing_reduced_list <- List_reducedEnrichment$luma$`hsa-mir-139`

my_mirna <- "hsa-mir-139"
my_gos   <- testing_reduced_list %>% pull(goterms) %>% as.character()

my_query_list <- lapply(X = my_gos, FUN = function(gogo){
  my_query = glue::glue("{my_mirna} \"{gogo}\" ")
  my_requests = rentrez::entrez_search(db = "pubmed", term = my_query)
  Sys.sleep(1)
  return(my_requests)
})



testing_reduced_list <- List_reducedEnrichment$luma


mirna_gobp_list <- 
lapply(X = List_reducedEnrichment, FUN = function(subtype){
  safely_entrecize = safely(.f = entrez_search)
  #list mirnas 
  my_mirnames <- names(subtype)
  names(my_mirnames) <- my_mirnames
  
  #iterate over mirnames 
  lapply(X = my_mirnames, FUN = function(mir){
    my_gos <- subtype[[mir]] %>% pull(goterms) %>% as.character()
    names(my_gos) <- subtype[[mir]] %>% pull(rowname) %>% as.character()
    #iterate over go terms 
    my_query_list <- lapply(X = my_gos, FUN = function(gogo){
      my_query = glue::glue("{mir} \"{gogo}\" ")
      my_requests = safely_entrecize(db = "pubmed", term = my_query)
      print("will nap now for 0.1 seconds")
      Sys.sleep(0.1)
      print("on it!")
      return(my_requests)
    })
    
  })
})


sapply(X = mirna_gobp_list$luma$`hsa-mir-139`, FUN = function(i){
  i[["result"]]$count
})

lapply(X = mirna_gobp_list$luma, FUN = function(j){
  sapply(X = j, FUN = function(i){
    i[["result"]]$count
  })
})

lapply(X = mirna_gobp_list, FUN = function(subtipo){
  lapply(X = subtipo, FUN = function(j){
    sapply(X = j, FUN = function(i){
      i[["result"]]$ids
    })
  })
  })
  
mirna_pubmed_list <- 
  lapply(X = List_reducedEnrichment, FUN = function(subtype){
    safely_entrecize = safely(.f = entrez_search)
    #list mirnas 
    my_mirnames <- names(subtype)
    names(my_mirnames) <- my_mirnames
    
    #iterate over mirnames 
    lapply(X = my_mirnames, FUN = function(mir){
      my_query = glue::glue("{mir}")
      my_requests = safely_entrecize(db = "pubmed", term = my_query)
      return(my_requests)
      
    })
  })

lapply(mirna_pubmed_list, FUN = function(i){
  sapply(i, FUN = function(j){
    j[["result"]]$count
  })
})






  ####################################

mirna_gobp_list$luma$`hsa-mir-139`$`GO:0001525`$result$count

mirna_gobp_list$luma$`hsa-mir-139` %>% lapply(function(i){
  r1 <- i[["result"]][["count"]]
  r2 <- as.data.frame(r1)
  }) %>% 
  bind_rows(.id = "goterm")
  
  
mirna_gobp_list$luma %>% 
lapply(FUN = function(i){
  lapply(i, function(j){
    r1 <- j[["result"]][["count"]]
    r2 <- as.data.frame(r1)
  })%>% bind_rows(.id = "goterm")
}) %>% bind_rows(.id = "mir")
   

mir_gobp_pubmed_df <- 
mirna_gobp_list%>% 
  lapply(FUN = function(k){
    lapply(X = k, FUN = function(i){
      lapply(i, function(j){
        r1 <- j[["result"]][["count"]]
        r2 <- as.data.frame(r1)
      })%>% bind_rows(.id = "goterm")
    }) %>% bind_rows(.id = "mir")
  }) %>% bind_rows(.id = "subtype")

mir_pubmed_df <- 
mirna_pubmed_list%>% 
  lapply(FUN = function(i){
    
      lapply(i, function(j){
        r1 <- j[["result"]][["count"]]
        r2 <- as.data.frame(r1)
      })%>% 
      bind_rows(.id = "mir")
    }) %>% 
  bind_rows(.id = "subtype")

mir_gobp_pubmed_df %>% head
mir_gobp_pubmed_df %>% filter(r1!=0)

mir_gobp_pubmed_df %>% filter(r1!=0) %>% 
  #group_by(subtype) %>% 
  filter(subtype == "basal") %>% 
  ggplot(aes(x = goterm, y = mir, fill = r1)) +
  
  #mutate(subtype_mir = glue::glue("{subtype}_{mir}")) %>% 
  #ggplot(aes(x = goterm, y = subtype_mir, fill = r1)) +
  
  geom_raster(aes(fill = r1)) + 
  geom_text(aes(label = r1)) + 
  theme_minimal() +
  scale_fill_distiller(palette = "Purples", direction = 1)

?scale_fill_distiller


mir_gobp_pubmed_df %>% filter(r1!=0)

rbind(as.matrix(cbind(mir_gobp_pubmed_df[,c(1,2)],1)),
      as.matrix(mir_gobp_pubmed_df[,c(2:4)])
      )


mir_gobp_pubmed_df %>% 
  filter(r1!=0) %>% 
  left_join(y = GO_dict, by = c("goterm" = "GO_ID")) %>% 
  write_tsv(path = "results/resultados_validados.tsv")
  
