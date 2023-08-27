

#' Consensus community detection
#' 
#' Performs a consensus community detection algorithm on a given network g, using one of a set of algorithms from the R iGraph library. The result is a more stable than a single-trial, and includes an estimate of uncertainty associated with each community label. 
#' @param g: the network to be analysed. It must be an iGraph object with a node attribute V(g)$id as an integer value, and an edge attribute E(g)$weight as a numeric value. If g is unweighted, E(g)$weight must be set to 1.0. The graph is treated as undirected.
#' @param t: the number of independent trials 
#' @param method: the method chosen for community detection. Possible values are:  "ML" multilevel.community(), "LD" cluster_leiden(), "FG" fastgreedy.community(), "IM" infomap.community(), "LP" label.propagation.community(), "WT" walktrap.community() and "LE" leading.eigenvector.community(g). 
#' @param resolution: the resolution parameter of LV algorithm
#' @param gamma_lim the threshold for gamma parameter used to calculate concensus, controlling the formation of single-node communities. Nodes that are assigned the same community label more than t*gamma_lim times are clusterd together; otherwise they will form a single-node community. Possible values between 0.0 and 1.0. Typical values are gamma_lim  = 0.5 (large consensus communities) or gamma_lim = 0.9 (smaller, sharper consensus communities, with a larger number of single-node communities)
#' @param shuffle: a boolean parameter. If TRUE the network vertices are randomly permuted before each trial of community detection. It allows to obtain results that are no dependent on the order of nodes and edges within the network. (default value = TRUE). 
#' 
#' @returns returns a community object, that stores the community labels as $membership and the uncertainty coefficients as $uncertainty
#' @export
#'  
consensus_community_detection <- function(g, t, method='LV', gamma_lim, resolution=c(1.0), shuffle=TRUE) {
    
    require(igraph)
    require(tidyverse)
    M <- find_communities_repeated(g,
                                   n_trials=t,
                                   method = method,
                                   shuffle = shuffle,
                                   resolution = resolution,#for Louvain
                                   verbose = FALSE)
    
    nco <- normalized_co_occurrence(M)
    
    CC <- consensus_communities(nco,gamma_lim=gamma_lim)
    
    cons_communities <- make_clusters(g, array(as.numeric(CC$cons_comm_label)))
    cons_communities$gamma<-CC$gamma
    return(cons_communities)
}

########################
########################
########################
########################
find_communities <- function(g,method, r = c(1.0),  verbose = FALSE) {
    #' community detection,
    #' retunrs a community object
    #' including the algorithm
    #' applies selected method
    #' applies resolution for LV and LD
    #' undirected(g) for LV and LD
    
    method = substr(method, 1, 2)
    if (method == "LV") {
        comms <- cluster_louvain(as.undirected(g, mode = 'each'), resolution = sample(r, 1))
    } else if (method == "ML") {
        comms <- multilevel.community(as.undirected(g, mode = 'each'), resolution = sample(r, 1))
    } else if (method == "LD") {
        comms <-
            cluster_leiden(as.undirected(g, mode = 'each'), resolution_parameter = quantile(strength(g))[2] / (gorder(g) - 1))
    } else if (method == "FG") {
        comms <- fastgreedy.community(g)
    } else if (method == "IM") {
        comms <- infomap.community(g)
    } else if (method == "LP") {
        comms <- label.propagation.community(g)
    } else if (method == "WT") {
        comms <- walktrap.community(g)
    } else if (method == "LE") {
        comms <- leading.eigenvector.community(g)
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

######################################################################

find_communities_repeated <- function(g,
                                      n_trials,
                                      method = method,
                                      shuffle = TRUE,
                                      resolution = c(1.0),
                                      verbose = FALSE) {
    
    #' community detection, repeated t times
    #' 
    membership_table <- data.frame(name = V(g)$name)

    for (i in 1:n_trials) {
        if (verbose){print(i)}
        if (shuffle == TRUE) {
            gs <- igraph::permute(g, sample(1:vcount(g),size = vcount(g),replace = FALSE))
        } else {
            gs <- g
        }
        comms <-find_communities(gs, method = method, r = resolution)
        comm_labeled <-data.frame(name = V(gs)$name, memb = comms$membership)
        membership_table <-inner_join(membership_table ,  comm_labeled, by = 'name')
        colnames(membership_table) <- c('name', seq(1:i))
    }
    membership_table <- membership_table %>% arrange(name)
    return(membership_table)
}



######################################################################

find_communities_N_times <- function(g,
                                     n_trials,
                                     methods,
                                     shuffle = FALSE,
                                     resolution = c(1.0),
                                     filename_summary = '',
                                     filename_membership = '') {
    
    #' community detection, repeated (THIS FUNCTION IS DEPRECATED)
    results_n_trials <- data.frame()
    membership_matrix <- c()


    for (i in 1:n_trials) {
        for (met in methods) {
            tmp_comms <- find_communities(g, resolution = resolution, method = met, shuffle = shuffle)
            V(g)$community <- tmp_comms$membership
            results_single <- analyse_communities(g, tmp_comms)
            results_n_trials = rbind(results_n_trials, results_single)
            membership_matrix <- rbind(membership_matrix, t(tmp_comms$membership))


        }
    }

    membership_n_trials <- data.frame(membership_matrix)
    if (filename_membership != '') {
        membership_n_trials %>% write_csv(filename_membership)
        results_n_trials %>% write_csv(filename_summary)
    }
    return(membership_n_trials)


}


# ################# analysis of communities
# returns modularity, number of communities and NMI against true labels
# communities <- cons_communities_2




normalized_co_occurrence <- function(M) {
    #' calculates normalized co occurence matrix

    names <- M$name
    M<-as.matrix(M %>% select(-name))
    n_trials <- ncol(M)
    n_nodes <- nrow(M)
    CO <-matrix(0, nrow = n_nodes,ncol =n_nodes)
    colnames(CO) <- names
    rownames(CO) <- names
    for (t in (1:n_trials)) {
        nclusters <- max(M[, t])
        for (k in 1:nclusters) {
            samecluster <- (which(M[, t] == k))
            nc <- length(samecluster)
            for (i in 1:nc) {
                for (j in 1:nc) {
                    CO[samecluster[j], samecluster[i]] <-
                        CO[samecluster[j], samecluster[i]] + 1
                }
            }
        }
    }
    #' X is a matrix of x_ij counts of co-occurrence
    #' X_normalized is a matrix of gamma_ij normalized coefficients
    #' expressing the probability that i and j are in the same community
    #' X_normalized_ij == 1.0 means that i and j have ALWAYS been classified in the same community
    #' X_normalized_ij == 0.0 means that i and j have NEVER been classified in the same community
    X_normalized <- CO / n_trials
    return (X_normalized)
}




consensus_communities <- function(nco, gamma_lim){
    #' calculates consensus based on co-occurrence matrix
    #' 
    results <- data.frame(name = as.numeric(colnames(nco)))
    results$done <- FALSE
    results$cons_comm_label <- 0

    # GAMMA
    coeffs <- nco
    diag(coeffs)<-0.0
    results$gamma <-  round( 1 -  apply(coeffs, 1, max), 4)

    community_label <- 0
    nodes_to_process = nrow(results)
    while ( nodes_to_process >= 1 )  {
        community_label <- community_label + 1
        row_test <- nco[ which.max(results$done == FALSE), ]
        nodes_above_threshold <-  (row_test > gamma_lim )
        results$cons_comm_label[ nodes_above_threshold ] <- community_label
        results$done[ nodes_above_threshold] <-  TRUE
        nodes_to_process <- sum( results$done == FALSE)
    }

    return(results %>% arrange(name))
}




########################



