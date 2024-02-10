#' Basic function to find communities
#' @export
find_communities <-
    function(g,
             method,
             r = NA,
             s = c(10),
             IMtrials = 1,
             verbose = FALSE) {
        require(igraph)
                gu <- as.undirected(g, mode = 'each')
        method = substr(method, 1, 2)
        if (method == "LV") {
            if (all(is.na(r))){r<- c(1.0)}
            comms <- cluster_louvain(gu, resolution = sample(r, 1))
        } else if (method == "ML") {
            if (all(is.na(r))){r<- c(1.0)}
            comms <- multilevel.community(gu, resolution = sample(r, 1))
        } else if (method == "LD") {
            if (all(is.na(r))){r<- c(quantile(strength(g))[2] / (gorder(g) - 1))}
            comms <- cluster_leiden(gu, resolution_parameter = sample(r, 1))
        } else if (method == "FG") {
            comms <- fastgreedy.community(gu)
        } else if (method == "IM") {
            comms <- infomap.community(gu, nb.trials = IMtrials) 
        } else if (method == "LP") {
            comms <- label.propagation.community(gu)
        } else if (method == "WT") {
            comms <- walktrap.community(gu, steps = sample(s, 1))
        } else if (method == "LE") {
            comms <- leading.eigenvector.community(gu)
        } else {
            print("No valid method")
            stop
        }
        comms$algorithm = method
        
        if (verbose == TRUE) {
            print(paste("Community detection with ", method, "completed."))
        }
        return(comms)
    }