## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4, tibble.print_max = 4)

## ----load library, message=FALSE, warning=FALSE-------------------------------
library(devtools)  
install_github("fabio-morea/CCD")
library(CCD)
library(igraph)
library(tidyverse)

## ----echo=TRUE----------------------------------------------------------------
comm <- consensus_community_detection(roc45bc,
                              t=100, 
                              method='LV', 
                              gamma_lim = 0.6) 

## -----------------------------------------------------------------------------
plot(comm,roc45bc )


## -----------------------------------------------------------------------------
plot(roc45bc, vertex.color = if_else(comm$gamma > .01,"red", "green"))

## ----echo=TRUE----------------------------------------------------------------
comm <- consensus_community_detection(roc84bc,
                              t=100, 
                              method='LV', 
                              gamma_lim = 0.6)
comm$gamma
plot(comm, roc84bc)
plot(roc84bc, vertex.color = if_else(comm$gamma > .01,"red", "green"))


