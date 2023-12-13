## ----setup, echo=FALSE, message=FALSE-----------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4, tibble.print_max = 4)

## ----load library, message=FALSE, warning=FALSE-------------------------------
library(devtools)  
devtools::install_github("fabio-morea/CCD")
library(CCD)
library(igraph)
library(tidyverse)

## -----------------------------------------------------------------------------
plot(roc45bc)

## ----c1, echo=TRUE------------------------------------------------------------
comm <- cluster_louvain(roc45bc) 
print(comm$membership)
plot(comm, roc45bc)
 


## -----------------------------------------------------------------------------
comm <- consensus_community_detection(roc45bc,
                              t=1000, 
                              method='LV', 
                              p = 0.6, 
                              shuffle = T)
print(comm$gamma)
plot(comm, roc45bc)
plot(roc45bc, 
     vertex.size = 30,
     vertex.label = comm$gamma,
     vertex.color = if_else(comm$gamma > .2,"yellow", "green"))

## ----c2, echo=TRUE------------------------------------------------------------
comm <- consensus_community_detection(roc84bc,
                              t=1000, 
                              method='LV', 
                              p = 0.55, 
                              shuffle = T)
print(comm$gamma)



## -----------------------------------------------------------------------------
plot(comm, roc84bc, vertex.label = NA)

## -----------------------------------------------------------------------------
plot(roc84bc, 
     vertex.size = 30,
     vertex.label = comm$gamma,
     vertex.color = if_else(comm$gamma > .2,"yellow", "green"))


## -----------------------------------------------------------------------------
#setting the parameters
V(roc84bc)$community <- comm$membership
E(roc84bc)$w <- E(roc84bc)$weight
# build the community network
gc <- make_community_network(roc84bc)
#plot
plot(gc, vertex.size = V(gc)$size*10, 
     vertex.label = V(gc)$size)

