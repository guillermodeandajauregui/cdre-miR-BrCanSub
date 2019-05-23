########
#testing
########

#libraries 
source("libs/functions_mi.R")

#paths
path_mir <- "basal/basal_FPKM.tsv"
path_rna <- "basal/basal_mirna_rpmmm.tsv"

#read data 

mir <- as.data.frame(readr::read_tsv(path_mir))
rna <- as.data.frame(readr::read_tsv(path_rna))

#discretizing 

tempus <- proc.time()
d.mir <- par_discretizer(mir, korez = 10)
tempus <- proc.time() - tempus
print(tempus)

tempus <- proc.time()
d.rna <- par_discretizer(rna, korez = 10)
tempus <- proc.time() - tempus
print(tempus)

#mi calculating 

tempus <- proc.time()
mirXrna <- par_mi_calc(sources = d.mir, 
                       targets = d.rna, 
                       korez = 10)
tempus <- proc.time() - tempus
print(tempus)

mi_matrix <- bind_rows(!!!mirXrna, #explicit splicing
                       .id = "mirna/gen")

write_tsv(mi_matrix, "basal_mi.tsv")
