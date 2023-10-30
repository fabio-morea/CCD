
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
                                 add_bridges = TRUE ) {
    f <- as.undirected(graph.empty())
    
    for (i in 1:num_cliques) {
        next_clique <- make_clique(clique_size, comm_label = paste0("C", i))
        f <- f + next_clique
    }
    
    b <- vcount(f)
    
    
    if (add_bridges) {
        f <- add_vertices(f, num_cliques)
    }
    
    
    for (j in 1:(num_cliques)) {
        b <- b + 1
        b_start <- (j-1) * clique_size +1
        b_end <- b_start + clique_size +1
        if (b_end > (clique_size * num_cliques)) {b_end <- 2}
        if (add_bridges) {
            f <- add_edges(f, c(b_start, b))
            f <- add_edges(f, c(b, b_end))
            V(f)$community[b] <- paste0("B", j)
        } else {
            f <- add_edges(f, c(b_start, b_end))
        }
        
    }
    
    if (add_center) {
        f <- add_vertices(f, 1)
        id_center <- vcount(f)
        V(f)$community[id_center] <- "A"
        for (j in 1:(num_cliques)) {
            c_start <- (j-1) * clique_size +3
            f <- add_edges(f, c(c_start , id_center))
        }
        
    }
    
    E(f)$weight <- 1.0
    
    V(g)$name <- 1:vcount(g)
    V(g)$id <-   1:vcount(g)
    
    
    return(g)
    
}

 