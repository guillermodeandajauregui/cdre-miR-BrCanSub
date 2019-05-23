library(igraph)
library(tidyverse)

prueboide <- read.graph(file = "results/basal_9999.analyzed.graphml", "graphml")
warnings()

get.data.frame(prueboide, "vertices") %>% head

get.data.frame(prueboide, "vertices") %>% pull(Redundancy)

get.data.frame(prueboide, "vertices") %>% 
  filter(grado!=0) %>% 
  ggplot(aes(x = Redundancy, y = grado)) + 
  geom_vline(mapping = aes(xintercept = 0.5)) +
  geom_hline(mapping = aes(yintercept = 100)) +
  geom_point()

get.data.frame(prueboide, "vertices") %>% pull(grado) %>% table
get.data.frame(prueboide, "vertices") %>% pull(grado) %>% mean
get.data.frame(prueboide, "vertices") %>% pull(grado) %>% median

get.data.frame(prueboide, "vertices") %>% filter(grado != 0) %>% pull(grado) %>% median
get.data.frame(prueboide, "vertices") %>% filter(grado != 0) %>% pull(grado) %>% mean

get.data.frame(prueboide, "vertices") %>% pull(Redundancy) %>% table()


pruebota <- read.graph(file = "results/basal_9999.analyzed.graphml", "graphml")
pruebaz <- read.graph(file = "probando.graphml", "graphml")

get.data.frame(pruebaz, "vertices") %>% 
  filter(grado!=0) %>% 
  ggplot(aes(x = Redundancy, y = grado)) + 
  geom_vline(mapping = aes(xintercept = 0.5)) +
  geom_hline(mapping = aes(yintercept = 100)) +
  geom_point()

my_analyzed <- list.files(path = "results", pattern = "analyzed.graphml", full.names = TRUE)[c(3,4,1,2)]
names(my_analyzed) <- c("luma", "lumb", "basal", "her2")
lista_analyzed <- lapply(X = my_analyzed, FUN = function(i){
  read_graph(file = i, format = "graphml")
})

get.data.frame(lista_analyzed$her2, "vertices") %>% 
  filter(type == TRUE) %>% 
  ggplot(aes(x = Redundancy, y = grado)) + 
  geom_vline(mapping = aes(xintercept = 0.5)) +
  geom_hline(mapping = aes(yintercept = 100)) +
  geom_point()

lista_cdrePlots <- 
lapply(lista_analyzed, FUN = function(i){
  get.data.frame(i, "vertices") %>% 
    filter(type == TRUE) %>% 
    ggplot(aes(x = Redundancy, y = grado)) + 
    geom_vline(mapping = aes(xintercept = 0.5)) +
    geom_hline(mapping = aes(yintercept = 100)) +
    geom_point()
})



get.data.frame(lista_analyzed$her2, "vertices") %>% 
  filter(type == TRUE) %>% 
  filter(grado != 0) %>% 
  ggplot(aes(x = Redundancy, y = grado)) + 
  geom_vline(mapping = aes(xintercept = median(Redundancy))) +
  geom_hline(mapping = aes(yintercept = median(grado))) +
  geom_point()

get.data.frame(lista_analyzed$her2, "vertices") %>% 
  filter(type == TRUE) %>%
  filter(grado != 0) %>% 
  ggplot(aes(x = Redundancy, y = grado)) + 
  geom_vline(mapping = aes(xintercept = mean(Redundancy))) +
  geom_hline(mapping = aes(yintercept = mean(grado))) +
  geom_point()


get.data.frame(lista_analyzed$luma, "vertices") %>% 
  filter(type==TRUE,
         grado >= 100, 
         Redundancy <= 0.5)

get.data.frame(lista_analyzed$lumb, "vertices") %>% 
  filter(type==TRUE,
         grado >= 100, 
         Redundancy <= 0.5)

get.data.frame(lista_analyzed$basal, "vertices") %>% 
  filter(type==TRUE,
         grado >= 100, 
         Redundancy <= 0.5)

get.data.frame(lista_analyzed$her2, "vertices") %>% 
  filter(type==TRUE,
         grado >= 100, 
         Redundancy <= 0.5)

####

lista_de_analisis <- lapply(X = lista_analyzed, FUN = function(i){get.data.frame(i, "vertices")})
allSubtypeAnalysis <- bind_rows(lista_de_analisis, .id = "subtype")
allSubtypeAnalysis %>% head

allSubtypeAnalysis %>%
  filter(type == TRUE) %>% 
    ggplot(aes(x = Redundancy, y = grado, color = as.factor(subtype))) + 
    geom_vline(mapping = aes(xintercept = 0.5)) +
    geom_hline(mapping = aes(yintercept = 100)) +
    geom_point() +
    scale_y_log10() +
    ggrepel::geom_text_repel(aes(label = ifelse((grado >= 100 & Redundancy <= 0.5), 
                                 yes = as.character(name), 
                                 no = "")
                  )
              )

allSubtypeAnalysis %>% filter(name %in% c("hsa-mir-139", "hsa-mir-136"))

#################
allSubtypeAnalysis %>%
  filter(type == TRUE) %>%
  filter(grado != 0) %>% 
  ggplot(aes(x = Redundancy, y = grado, color = as_factor(subtype))) + 
  #scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "yellow")) +
  scale_color_manual(values = c("firebrick3", "#00BA38", "#619CFF", "goldenrod3")) +
  geom_vline(mapping = aes(xintercept = 0.5)) +
  geom_hline(mapping = aes(yintercept = 100)) +
  geom_point(alpha = 0.5) +
  scale_y_log10() +
  ggrepel::geom_text_repel(aes(label = ifelse((grado >= 100 & Redundancy <= 0.5), 
                                              yes = as.character(name), 
                                              no = "" 
                                              )
  ), show.legend = FALSE
  )+
  xlab("Redundancy") + 
  ylab("Degree")+
  theme_minimal() + 
  labs(colour = "Subtype")
  
  
allSubtypeAnalysis %>%
  filter(type == TRUE) %>%
  filter(grado != 0) %>% 
  mutate(subtype = as.factor(subtype)) %>% 
  pull(subtype)
