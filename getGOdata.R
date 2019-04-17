#GO list from HTSanalyzer
library(HTSanalyzeR)
library(org.Hs.eg.db)
library(GO.db)

#get GO terms
GO_MF <- GOGeneSets(species="Hs", ontologies=c("MF"))
GO_BP <- GOGeneSets(species="Hs", ontologies=c("BP"))
GO_CC <- GOGeneSets(species="Hs", ontologies=c("CC"))

GO_ALL <- c(GO_MF, GO_BP, GO_CC)

#make map 
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
xx[GO_MF$`GO:0000010`] %>% unlist()

#convert all GO to gene symbol
length(GO_ALL)
length(GO_BP)
length(GO_MF)
length(GO_CC)

GO_ALL_SYMBOL <- lapply(X = GO_ALL, FUN = function(i){
  xx[i] %>% unlist()
})



#make terms dictionary 
goterms <- Term(GOTERM)
goterms <- as.data.frame(goterms) 
goterms <- goterms %>% rownames_to_column()
colnames(goterms)[1] <- "GO_ID"
goterms <- goterms %>% 
  mutate("category" = ifelse(test = GO_ID%in%names(GO_BP), 
                             yes  = "BP", 
                             no   = ifelse(test = GO_ID%in%names(GO_CC), 
                                           yes = "CC", 
                                           no = ifelse(test = GO_ID%in%names(GO_MF), 
                                                       yes = "MF", 
                                                       no = "X")
                                           ) 
                               )
         )

goterms <- goterms %>% filter(category != "X")

#saves 
saveRDS(object = GO_ALL_SYMBOL, file = "GO_symbol.RDS")
saveRDS(object = goterms, file = "goterms_dict.RDS")
