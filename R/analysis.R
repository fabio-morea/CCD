#' This function builds a network (Gc) whose nodes represent communities within the original network (G). 
#' Gc can be utilized for in-depth analysis and visualization of the relationships among these communities.
#' 
#' @param g: The network to be analyzed. It should be an iGraph object with specific attributes:
#'   - V(g)$community: An integer value representing the community assignment for each node.
#'   - E(g)$w: A numeric vector containing edge weights. If the network is unweighted, set E(g)$w to 1.0.
#' 
#' @returns A network object (Gc) with the following attributes:
#'   - $membership: Stores the community labels.
#'   - E(Gc)$w: The sum of corresponding edge weights in the original network (G).
#'   - V(Gc)$size: the number of nodes of the original network G that are represented in a node of Gc.
#'   
#' @export

make_community_network <- function (g) {
    edges_list <- g %>%
        as_long_data_frame() %>%
        select(from_community, to_community, w) %>%
        group_by(from_community, to_community) %>%
        summarize(weight = sum(w))
    
    gc <- graph_from_data_frame(edges_list, directed = FALSE)
    V(gc)$id <- (1:vcount(gc))
    
    comms <- data.frame(label = V(gc)$name)
    comms$size <- 0
    for (i in 1:length(comms$label)) {
        comms$size[i] <-
            length(V(g)$community[V(g)$community == comms$label[i]])
    }
    V(gc)$size <- comms$size
    
    return(gc)
}




analyse_communities <- function(g, communities, verbose = FALSE) {
    #' Analyse communities function
    #' Prints summary information about the communities
    #'
   
    method <- communities$algorithm
    c_membership <- communities$membership
    # modularity: need to use c_membership + 1 to handle community label 0
    mod <- round(modularity (g,  c_membership + 1), digits = 4)

    #number of communities
    nc <- length(table(c_membership))
    nc_builtin <- length(table(V(g)$comm_built_in))
    nc_norm <- nc / nc_builtin

    mu_built_in <- round(empirical_mu(g), 4)

    # NMI against "Built In Communities"
    nmi = round(aricode::NMI(as.factor(V(g)$comm_built_in), as.factor(c_membership)), 3)

    #empirical value of mu on communities found
    V(g)$community <- c_membership
    mu_emp <- round(empirical_mu(g), 4)

    if (verbose == TRUE) {
        print(paste("Muxing parameter Mu (empirical value): ", mu_emp))
        print(paste("Communities found: ", nc))
        print(paste("Modularity: ", mod))
        print(paste("Normalized Mutual Information between C.D. method and built-in communities:",nmi)
        )
    }

    return(data.frame(method , mu_built_in, mu_emp, mod, nc, nc_norm, nmi))
}


# ############# Calculate empirical value of mixing parameter MU

empirical_mu <- function(g) {
    #' Empirical Mu function
    #' Calculates the value of mixing parameter "mu"
    #'
    gdf<-as_long_data_frame(g)
    gdf$inter_comm <- (gdf$from_community != gdf$to_community)
    inter_community_links <- sum(gdf$weight[ gdf$inter_comm == TRUE])
    mu = sum(inter_community_links) / sum(gdf$weight)
    return(mu)
}

 



# find neighbours

find_neighbours <- function(g, node_name, order) {
    selected_node_id <- which(V(g)$name == node_name)
    list_neis <- ego(g, order = order, selected_node_id)
    nei_ids <- c()
    for (nn in list_neis) {
        for (x in nn) {
            nei_ids <- append(nei_ids, x)
        }
    }
    nei_ids <- c(nei_ids, selected_node_id)
    nei_ids <- unique(nei_ids)
    return(nei_ids)
}
