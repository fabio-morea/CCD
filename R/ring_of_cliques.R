
#' Functions to make a "ring of cliques"
#' @export
make_clique <- function(clique_size, comm_label) {
    G <- graph.empty(n = clique_size)
    edges <- t(combn(1:clique_size, 2))
    for (e in 1:nrow(edges)) {
        G <- add_edges(G, edges[e, ])
    }
    V(G)$community <- comm_label
    
    return(as.undirected(G))
}
 
#' @export
make_ring_of_cliques <- function(num_cliques,
                                 clique_size,
                                 add_center = TRUE,
                                 add_bridges = TRUE,
                                 filename = 'ring_of_cliques') {
    G <- as.undirected(graph.empty())
    
    for (i in 1:num_cliques) {
        next_clique <- make_clique(clique_size, comm_label = paste0("C", i))
        G <- G + next_clique
    }
    
    
    b <- vcount(G)
    print(b)
    if (add_bridges) {
        G <- add_vertices(G, num_cliques)
    }
    
    
    for (j in 1:(num_cliques)) {
        b <- b + 1
        b_start <- (j-1) * clique_size +1
        b_end <- b_start + clique_size +1
        if (b_end > (clique_size * num_cliques)) {b_end <- 2}
        if (add_bridges) {
            print(vcount(G))
            print(paste(b, b_start, b_end))
            G <- add_edges(G, c(b_start, b))
            G <- add_edges(G, c(b, b_end))
            V(G)$community[b] <- paste0("B", j)
        } else {
            G <- add_edges(G, c(b_start, b_end))
        }
        
    }
    
    if (add_center) {
        G <- add_vertices(G, 1)
        id_center <- vcount(G)
        V(G)$community[id_center] <- "A"
        for (j in 1:(num_cliques)) {
            c_start <- (j-1) * clique_size +3
            G <- add_edges(G, c(c_start , id_center))
        }
        
    }
    
    E(G)$weight <- 1.0
    
    write_graph(G, paste0(filename, '.gml'), format = 'gml')
    
    return(G)
    
}

 