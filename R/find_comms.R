#' Basic function to find communities
#' @export
find_communities <-
    function(g,
             method,
             r = c(1.0),
             s = c(10),
             verbose = FALSE) {
        require(igraph)
                gu <- as.undirected(g, mode = 'each')
        method = substr(method, 1, 2)
        if (method == "LV") {
            comms <- cluster_louvain(gu, resolution = sample(r, 1))
        } else if (method == "ML") {
            comms <- multilevel.community(gu, resolution = sample(r, 1))
        } else if (method == "LD") {
            comms <-
                cluster_leiden(gu, resolution_parameter = quantile(strength(g))[2] / (gorder(g) - 1))
        } else if (method == "FG") {
            comms <- fastgreedy.community(gu)
        } else if (method == "IM") {
            comms <- infomap.community(gu)
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