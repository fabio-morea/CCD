


#' Consensus community detection
#'
#' Performs a consensus community detection algorithm on a given network g, using one of a set of algorithms from the R iGraph library. The result is a more stable than a single-trial, and includes an estimate of uncertainty associated with each community label.
#' @param g: the network to be analysed. It must be an iGraph object with a node attribute V(g)$id as an integer value, and an edge attribute E(g)$weight as a numeric value. If g is unweighted, E(g)$weight must be set to 1.0. The graph is treated as undirected.
#' @param t: the number of independent trials
#' @param method: the method chosen for community detection. Possible values are:  "ML" multilevel.community(), "LD" cluster_leiden(), "FG" fastgreedy.community(), "IM" infomap.community(), "LP" label.propagation.community(), "WT" walktrap.community() and "LE" leading.eigenvector.community(g).
#' @param resolution: the resolution parameter of LV algorithm
#' @param p the threshold for gamma parameter used to calculate concensus, controlling the formation of single-node communities. Nodes that are assigned the same community label more than t*p times are clusterd together; otherwise they will form a single-node community. Possible values between 0.0 and 1.0. Typical values are p  = 0.5 (large consensus communities) or p = 0.9 (smaller, sharper consensus communities, with a larger number of single-node communities)
#' @param shuffle: a boolean parameter. If TRUE the network vertices are randomly permuted before each trial of community detection. It allows to obtain results that are no dependent on the order of nodes and edges within the network. (default value = TRUE).
#'
#' @returns returns a community object, that stores the community labels as $membership and the uncertainty coefficients as $uncertainty
#'




#' @export
consensus_community_detection <- function(g,
             t=100,
             method = 'LV',
             p=0.6,
             group_outliers = FALSE,
             resolution = c(1.0),
             steps = c(5),
             IMtrials = 1,
             shuffle = TRUE) {
        require(igraph)
        require(tidyverse)
        M <- CCD::find_communities_repeated(
            g,
            n_trials = t,
            method = method,
            shuffle = shuffle,
            resolution = resolution,
            #for Louvain
            steps = steps,
            #for walwtrap
            verbose = FALSE
        )
        
        D <- CCD::normalized_co_occurrence(M)
        
        CC <- CCD::consensus_communities(D, p = p,group_outliers = group_outliers )
        
        cons_communities <-make_clusters(g, array(as.numeric(CC$cons_comm_label)))
        cons_communities$gamma <- CC$gamma
        cons_communities$name <- CC$name
        cons_communities$comm_size <- CC$comm_size
        return(cons_communities)
    }

#' @export

find_communities_repeated <- function(g,
                                      n_trials,
                                      method = method,
                                      shuffle = TRUE,
                                      resolution = c(1.0),
                                      steps = c(10),
                                      IMtrials = 1,
                                      verbose = FALSE) {
    membership_table <- data.frame(name = V(g)$name)
    
    for (i in 1:n_trials) {
        if (verbose) {
            print(i)
        }
        if (shuffle == TRUE) {
            gs <-
                igraph::permute(g, sample(
                    1:vcount(g),
                    size = vcount(g),
                    replace = FALSE
                ))
        } else {
            gs <- g
        }
        comms <-
            find_communities(gs,
                             method = method,
                             r = resolution,
                             s = steps)
        comm_labeled <-
            data.frame(name = V(gs)$name,
                       memb = comms$membership)
        membership_table <-
            inner_join(membership_table ,  comm_labeled, by = 'name')
        colnames(membership_table) <- c('name', seq(1:i))
    }
    return(membership_table)
}




#' @export
#'
normalized_co_occurrence <- function(M) {
    # calculates normalized co occurence matrix
    
    names <- M$name
    M <- as.matrix(M %>% select(-name))
    n_trials <- ncol(M)
    n_nodes <- nrow(M)
    CO <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    colnames(CO) <- names
    rownames(CO) <- names
    
    for (t in (1:n_trials)) {
        nclusters <- max(M[, t])
        for (k in 1:nclusters) {
            samecluster <- (which(M[, t] == k))
            nc <- length(samecluster)
            for (i in 1:nc) {
                for (j in (i+1):nc) {
                    CO[samecluster[j], samecluster[i]] <- CO[samecluster[j], samecluster[i]] + 1
                    CO[samecluster[i], samecluster[j]] <- CO[samecluster[j], samecluster[i]]
                }
            }
        }
    }
    diag(CO)<-1
    X_normalized <- CO / n_trials
    return (X_normalized)
}



#' @export

consensus_communities <- function(D, p, group_outliers = FALSE, verbose = FALSE, save_results=FALSE) {
    
    # definition of community: block within D in which dij > p 
    # this definition includes single node communities (outliers)
    
    # definition of uncertainty coefficient gamma: 
    #     (1-MEAN of di) over all nodes that are at least once in the same community
    
    
    
    results <- data.frame(name = colnames(D))
    results$done <- FALSE
    results$tmp_comm_label <- NA
    results$gamma <- NA
    results$comm_size <- NA
    results$single <- FALSE
    community_label <- 0
    nodes_to_process = nrow(results)
     
        # definition of community: block within D in which dij > p
        # this definition includes single node communities (outliers)
        
        # definition of uncertainty coefficient gamma:
        # (1-MEAN of di) over all nodes that are at least once in the same community
        
        results <- data.frame(name = colnames(D))
        results$done <- FALSE
        results$tmp_comm_label <- NA
        results$gamma <- NA
        results$comm_size <- NA
        results$single
        community_label <- 0
        nodes_to_process = nrow(results)
        
        while (nodes_to_process > 0)  {
            community_label <- community_label + 1
            
            #select a block with respect to threshold p, first row not done
            nodes_internal <- (D[which.max(results$done == FALSE),] > p)
            
            # calculate gamma for eachnode in the block
            gammas <- D[nodes_internal,]
            # ignore nodes that are never in the same community
            gammas[gammas == 0] <- NA
            
            if (sum(nodes_internal) > 1) {
                # a proper block
                results$gamma[nodes_internal] <- 1 - apply(gammas, 1, mean, na.rm = T)
            } else {
                # a single node
                results$gamma[nodes_internal] <- 1 - mean(gammas,  na.rm = T)
            }
            
            results$tmp_comm_label[nodes_internal] <- community_label
            results$comm_size[nodes_internal] <- sum(nodes_internal)
            results$done[nodes_internal] <-  TRUE
            nodes_to_process <- sum(results$done == FALSE)
        }
        
        results$gamma[is.na(results$gamma)] <- 0.0
        results$single[results$comm_size == 1] <- TRUE
        
        if (group_outliers) { results$tmp_comm_label[results$single] <- 0 }
        
        x <- results %>%
            group_by(tmp_comm_label) %>%
            summarize(n = n()) %>% arrange(-n) %>%
            mutate(cons_comm_label = row_number())
        
        results <- results %>%
            inner_join (x, by = 'tmp_comm_label') %>%
            select(name, cons_comm_label, gamma, comm_size, single)
        
        
        if (save_results) {
            results %>% write_csv('results.csv')
        }
        
        return(results)
    }

     
     



########################
