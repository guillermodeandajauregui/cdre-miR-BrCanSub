####################
#discretization and 
#mutual information 
#calculation
#functions
####################

#packrat unbundle
packrat::on()
#packrat::unbundle(bundle = "packrat/bundles/cdre-miR-BrCanSub-2019-04-10.tar.gz", where = ".")
packrat::restore()

#libraries
library("infotheo")
library("dplyr")
library("readr")

#functions 

#discretization

discretizer <- function(expmatrix, diskz = "equalfreq"){
  my_index <- 1:nrow(expmatrix)
  names(my_index) <- dplyr::pull(expmatrix, 1)
  my_discrete <- lapply(X = my_index, FUN = function(i){
    my_row <- expmatrix[i,-1]
    my_row <- infotheo::discretize(X = my_row, disc = diskz)
    
    
  })
  return(my_discrete)
}

par_discretizer <- function(expmatrix, korez = 2, diskz = "equalfreq"){
  my_index <- 1:nrow(expmatrix)
  names(my_index) <- dplyr::pull(expmatrix, 1)
  my_discrete <- parallel::mclapply(X = my_index, mc.cores = korez, FUN = function(i){
    my_row <- expmatrix[i,-1]
    my_row <- infotheo::discretize(X = my_row, disc = diskz)
    
  })
  return(my_discrete)
}

#MI calculation 

mi_calc <- function(sources, targets){
  lapply(X = sources, FUN = function(i){
    sapply(X = targets, FUN = function(j){
      infotheo::mutinformation(X = i, 
                               Y = j)
    })
  })
}

par_mi_calc <- function(sources, targets, korez = 2){
  parallel::mclapply(X = sources, mc.cores = korez, FUN = function(i){
    sapply(X = targets, FUN = function(j){
      infotheo::mutinformation(X = i, 
                               Y = j)
    })
  })
}