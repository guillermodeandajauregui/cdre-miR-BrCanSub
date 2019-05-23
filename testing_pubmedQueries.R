library(RISmed)
library(rentrez)
library(fulltext)
library(tidytext)

#testing rismed
testingshit <- EUtilsQuery("myeloma[ti] smith[au]")

testingshit
summary(testingshit)

#didnt like. Also why do we need kovalchik's email in the query? There must be a reason

#testing rentrez

#get the unique ids
katipo_search <- rentrez::entrez_search(db="popset", term="Latrodectus katipo[Organism]")
katipo_search$ids
##how many? 
katipo_search$count

#get the summaries
katipo_summs <- entrez_summary(db="popset", id=katipo_search$ids)
katipo_summs$`167843256`$title

entrez_dbs()

my_search <- rentrez::entrez_search(db="pubmed", term = "Guillermo de Anda-Jauregui[Au")
my_summs  <- rentrez::entrez_summary(db = "pubmed", id = my_search$ids)

sapply(my_summs, FUN = function(i){i[["title"]]})


my_test_search <- rentrez::entrez_search(db = "pubmed", term = "cancer", use_history = TRUE)
my_summs  <- rentrez::entrez_summary(db = "pubmed", web_history = my_test_search$web_history, retmax = 25)
my_summs$`31018252`$authors

#search for a micro RNA 

my_mirna_search <- rentrez::entrez_search(db = "pubmed", term = " \"hsa-mir-136\"[All Fields] AND cancer[All Fields]", use_history = TRUE)
my_mirna_search
my_mirna_search$count
my_mirna_search$QueryTranslation
my_mirna_search$web_history

entrez_db_searchable(db = "pubmed")

my_summs <- rentrez::entrez_summary(db = "pubmed", web_history = my_mirna_search$web_history)
lapply(my_summs, function(i){i[["title"]]})


my_mirna_1 <- "hsa-mir-136"
my_first_reduced$goterms %>% as.character()

my_mirna_2 <- "hsa-mir-139"
my_second_reduced$goterms %>% as.character()

my_list_results <- 
lapply(X = my_second_reduced$goterms %>% as.character(), FUN = function(i){
  #my query 
  #my_query = paste0(my_mirna_1, " ", "(",i,")")
  my_query = glue::glue("{my_mirna_1} \"{i}\" ")
  
  rentrez::entrez_search(db = "pubmed", term = my_query)
  
  
})

my_list_results[[1]]$ids
my_list_results[[2]]$ids

my_summaries_1 <- entrez_summary(db="pubmed", id=my_list_results[[1]]$ids)
my_summaries_2 <- entrez_summary(db="pubmed", id=my_list_results[[2]]$ids)


rentrez::extract_from_esummary(esummaries = my_summaries_1, elements = "title")

