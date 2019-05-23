library(tidyverse)

normexpresion <- function(expmatrix){
  r1 <- rowMeans(expmatrix[,-1]) / apply(expmatrix[,-1], MARGIN = 1, FUN = max)
  names(r1)<- expmatrix %>% pull(1)
  r1 <- rownames_to_column(as.data.frame(r1), var = "name")
  colnames(r1)[2]  <- "avgExp"
  r1
}

luma_fpkm  <- vroom::vroom(file = "luma/luma_FPKM.tsv")
luma_rpmmm <- vroom::vroom(file = "luma/luma_mirna_rpmmm.tsv")

lumb_fpkm  <- vroom::vroom(file = "lumb/lumb_FPKM.tsv")
lumb_rpmmm <- vroom::vroom(file = "lumb/lumb_mirna_rpmmm.tsv")

basal_fpkm  <- vroom::vroom(file = "basal/basal_FPKM.tsv")
basal_rpmmm <- vroom::vroom(file = "basal/basal_mirna_rpmmm.tsv")

her2_fpkm  <- vroom::vroom(file = "her2/her2_FPKM.tsv")
her2_rpmmm <- vroom::vroom(file = "her2/her2_mirna_rpmmm.tsv")

bind_rows(normexpresion(luma_fpkm), normexpresion(luma_rpmmm)) %>% write_delim(path = "results/avgExp_luma.txt", delim = "\t")
bind_rows(normexpresion(lumb_fpkm), normexpresion(lumb_rpmmm)) %>% write_delim(path = "results/avgExp_lumb.txt", delim = "\t")
bind_rows(normexpresion(basal_fpkm), normexpresion(basal_rpmmm)) %>% write_delim(path = "results/avgExp_basal.txt", delim = "\t")
bind_rows(normexpresion(her2_fpkm), normexpresion(her2_rpmmm)) %>% write_delim(path = "results/avgExp_her2.txt", delim = "\t")