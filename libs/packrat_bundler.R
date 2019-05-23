#packages used in this project
#install.packages("tidyverse")
#install.packages("data.table")
packrat::on()

install.packages("infotheo")
install.packages("dplyr")
install.packages("readr")

packrat::snapshot()
packrat::bundle(project = NULL, file = NULL, include.lib = TRUE, overwrite = TRUE)
