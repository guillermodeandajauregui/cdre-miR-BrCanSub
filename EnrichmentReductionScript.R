#get reduced enrichment for cdre-miR on all breast cancer subtypes

source("libs/enrichmentFunctions.R")
source("libs/reducedEnrichmentFunction.R")

#read data
EnrichmentList <- readRDS("EnrichmentOutput.RDS")
#filter according to the ad-hoc criteria
EnrichmentList <- FilterMyList(EnrichmentList)

#prepare the gene ontology data
hsGO <- godata('org.Hs.eg.db', ont = "BP") #make sure you filtered out all Non BP categories

#for each cdre-mir in each subtype
List_reducedEnrichment <- lapply(X = EnrichmentList, FUN = function(subtype){
  lapply(X = subtype, FUN = function(cdre.rich){
                            reducedEnrichment(rich.df = cdre.rich, 
                                              GO.ID.row   = "rowname",
                                              sim.measr   = "Wang", 
                                              SemData     =  hsGO, 
                                              dist.mthd   = "euclidean",
                                              clustr.mthd = "ward.D2",
                                              no.groups   = 10,
                                              cluster.nom = "my_groups")
  })
})

saveRDS(object = List_reducedEnrichment, file = "List_reducedEnrichment.RDS")
