library("dplyr")
library("ontologyPlot")
library("ontologyIndex")
library("Rgraphviz")

#' MAIN FUNCTION
#' Get a DAG and node ordering of nodes of interest. Attempt 2, by starting on
#' root and then traversing the tree in a depth-first manner and mantaining
#' left-to-right order of siblings.
#' 
#' @param CL ontology object loaded with ontologyIndex::get_ontology.
#' @param cell_types vector of strings indicating nodes of interest.

get_cell_type_order <- function(CL, cell_types, root_term="CL:0000000") {
  
  cell_type_ancestors <- get_ancestors_ontology(cell_types, CL, root_term)
  
  # Get tree and graph coodinates
  tree <- get_tree(cell_type_ancestors, cell_types, CL)
  graph <- agopen_ontology_plot(tree)
  graph_df <- get_node_coordinates(graph)
  
  targets <- set_order(term = root_term, ontology_object = CL, 
                       graph_df = graph_df, all_terms = cell_type_ancestors, 
                       terms_of_interest = cell_types)
  
  ordered_nodes <- targets$ordered_terms
  rownames(graph_df) <- graph_df$name
  return(list(onto_plot = tree, ordered_nodes=graph_df[ordered_nodes,]))
  
}

# This function is in "ontologyPlot" but it's not exported 
# This is a wrapper to agopen from "Rgraphiz"

agopen_ontology_plot <- function(x) {
  ont_graph <- new(
    "graphAM", 
    adjMat=x[["adjacency_matrix"]], 
    edgemode="directed"
  )
  
  result <- agopen(graph=ont_graph, nodeAttrs=x[["node_attributes"]], 
                   name="ontological_plot") 
  if (length(result@AgEdge) > 0)
    for (i in 1:length(result@AgEdge)) {
      for (aai in 1:length(x[["edge_attributes"]])) {
        slot(result@AgEdge[[i]], names(x[["edge_attributes"]])[aai]) <- x[["edge_attributes"]][[aai]]
      }
    }
  return(result)
}


#' Get node coordinates from an agopen object obtained from agopen_ontology_plot
#' 
#' @param x agopen object from with agopen_ontology_plot
#'
#' @return A data.frame with x, and y coordinates. Columns are 'name', 'label', 
#' 'x', and 'y'

get_node_coordinates <- function(x) {
  
  nodes <- x@AgNode
  results <- list()
  for (i in 1:length(nodes)) {
    node <- nodes[[i]]
    results[[i]] <- data.frame(name = node@name, 
                               label = node@txtLabel@labelText, 
                               x=node@center@x, y=node@center@y, 
                               stringsAsFactors=F)
  }
  
  do.call(bind_rows, results)
  
}


#' Get all ancestors of a set of terms using an ontology
#' 
#' @param terms vector, terms of interest. 
#' @param ontology ontology object loaded with ontologyIndex.
#' @param root_term character, indicate up to what term the ancestors should be
#' obtained.
#' 
#' @return vector of strings with all ancestor terms.

get_ancestors_ontology <- function(terms, ontolgoy_object, root_term) {
  
  tree_terms <- get_descendants(CL, root_term)
  
  # All ancestors of our targets 
  ancestors <- c(get_ancestors(ontolgoy_object, terms))
  ancestors <- ancestors[ancestors %in% tree_terms]
  ancestors <- unique(ancestors)
  
  return(ancestors)
  
}


#' Get tree object for a set of terms from an ontology.
#' 
#' @param terms vector, all terms to draw.
#' @param target_terms vector, terms of interest will be assigned a different.
#' color
#' @param ontology_object the ontology object.
#' 
#' @return tree object from onto_plot 
 
get_tree <- function(terms, target_terms, ontology_object, 
                     light_color = "#C9EBF7", strong_color = "#FF9D1F") {
  
  node_colors <- rep(light_color, length(terms))
  node_colors[terms %in% target_terms] <- strong_color
  tree <- onto_plot(ontology_object, terms=terms, fillcolor = node_colors, 
                    edge_attributes = list(lwd = 0.5))
  
  return(tree)
}


#' From a node distance data.frame get a horizontal distance matrix 
#' for nodes of interest.
#' 
#' @param graph_dfd data.frame with node coordinates from get_node_coordinates.
#'
#' @return A distance matrix, diagonal set to Inf .

get_graph_horizontal_distances <- function(graph_df) {
  
  x_dist <- as.matrix(dist(graph_df[,c("x")], diag = FALSE, upper=TRUE))
  x_dist <- x_dist / max(x_dist)
  x_dist[x_dist == 0] <- Inf
  
  return(x_dist)
}


#' Order terms from left to right based on node coordinates.
#' 
#' @param terms vector. Terms matching to the column 'name' of graph_df.
#' @param graph_df data.frame. With node coordinates from get_node_coordinates.
#'
#' @return same as terms but ordered.

lef_to_right_order <- function(terms, graph_df){  
  
  graph_df <- filter(graph_df, name %in% terms)
  node_dist <- get_graph_horizontal_distances(graph_df)
  
  # getting ordering
  all_nodes <- 1:nrow(graph_df)
  current_node <- which.min(graph_df$x)
  ordered_nodes<-numeric()
  while(length(all_nodes) > 0) {
    
    ordered_nodes <- c(ordered_nodes, current_node)
    all_nodes <- all_nodes[!all_nodes == current_node]
    
    current_order <- order(node_dist[current_node,,drop=T])
    current_node <- current_order[!current_order %in% ordered_nodes][1]
  } 
  
  return(graph_df[ordered_nodes, "name", drop=T])
}


#' Get all children for a given term.
#' 
#' @param ontology_object ontology_object from ontologyIndex::get_ontology.
#' @param term character. Query term.
#' @param all_terms vector. Subset results to only include terms in this vector.
#'
#' @return vector will children terms.
 
get_children <- function(ontology_object, term, all_terms = NULL) {
  
  children <- get_term_property(ontology = ontology_object, 
                                property_name = "children", 
                                term = term)
  
  if (length(children > 0) & !is.null(all_terms)) {
    children <- children[children %in% all_terms]
  }
  
  return(children)
}


#' Seed function for recursion. Get terms in order using a tree traversion
#' algorithm. Starts with a root and traverses the tree until finding the bottom
#' leftmost node, then traverses back and prints siblings along the way in a
#' left-to-right order.  
#'  
#' @param term character. Root term to start with
#' @param ontology_object ontology_object from ontologyIndex::get_ontology.
#' @param graph_df data.frame. With node coordinates from get_node_coordinates.
#' @param all_terms vector. Only traverse nodes in this vector
#' @param terms_of_interest vector. Only order nodes in this vector
#'
#' @return environment. An environment with two object, $terms and 
#' $ordered_terms. $terms contain the original terms_of_interest, $ordered_terms
#' is the same but ordered.

set_order <- function(term, ontology_object, graph_df, all_terms, 
                      terms_of_interest) {
  
  targets <- new.env()
  targets$terms <- terms_of_interest
  targets$ordered_terms <- vector(mode = "character")
  
  set_order_recursion(term = term, ontology_object = ontology_object, 
                      graph_df = graph_df, all_terms = all_terms, targets = targets)

  return(targets)
  
}

#' Recursive function of set_order

set_order_recursion <- function(term, ontology_object, graph_df, all_terms, targets) {
  
  # term if is target
  if (term %in% targets$terms){
    targets$terms <- targets$terms[ targets$terms != term]
    targets$ordered_terms <- c(targets$ordered_terms, term)
  }
  
  # Go to children if there are any
  children <- get_children(ontology_object, term, all_terms)
  if (length(children) > 0) {
    
    # Order children starting with left most
    children <- lef_to_right_order(children, graph_df)  
    
    for(child in children) {
      set_order_recursion(term = child, ontology_object = ontology_object, 
                graph_df = graph_df, all_terms = all_terms, targets = targets)
    }
  }
  
  return (NULL)
  
}


#' Get a DAG and node ordering of nodes of interest
#' 
#' @param CL ontology object loaded with ontologyIndex
#' @param cell_types vector of strings indicating nodes of interest
#'
#' @return A list with to elements, the first "onto_plot" is a drawable object 
#' as obtained by `ontologyPlot::onto_plot`, the second "ordered_nodes" a 
#' data.frame with the nodes of interest ordered.

get_cell_type_order_old <- function(CL, cell_types) {
  
  # All terms that are under "animal cell"
  CL_tree_terms <- get_descendants(CL, "CL:0000000")
  
  # All ancestors of our cell types of interest
  cell_type_ancestors <- c(get_ancestors(CL, cell_types))
  
  cell_type_ancestors <- cell_type_ancestors[cell_type_ancestors %in% 
                                               CL_tree_terms]
  cell_type_ancestors <- unique(cell_type_ancestors)
  
  # Get tree
  light_color <- "#C9EBF7"
  strong_color <- "#FF9D1F"
  node_colors <- rep(light_color, length(cell_type_ancestors))
  node_colors[cell_type_ancestors %in% cell_types] <- strong_color
  tree <- onto_plot(CL, terms=cell_type_ancestors, 
                    fillcolor = node_colors, edge_attributes = list(lwd = 0.5))
  
  # Obtain the graph object 
  graph <- agopen_ontology_plot(tree)
  
  #browser()
  graph_df <- get_node_coordinates(graph)
  graph_df <- filter(graph_df, name %in% cell_types)
  
  # Get y distances between all of the
  node_dist <- as.matrix(dist(graph_df[,c("x","y")], diag = FALSE, upper=TRUE))
  y_dist <- as.matrix(dist(graph_df[,c("x")], diag = FALSE, upper=TRUE))
  
  y_dist <- y_dist / max(y_dist)
  node_dist <- (node_dist / max(node_dist)) * 1.3
  full_dist <- node_dist + y_dist
  full_dist[full_dist == 0] <- Inf
  
  # getting ordering
  all_nodes <- 1:nrow(graph_df)
  current_node <- which.min(graph_df$x)
  ordered_nodes<-numeric()
  while(length(all_nodes) > 0) {
    
    ordered_nodes <- c(ordered_nodes, current_node)
    all_nodes <- all_nodes[!all_nodes == current_node]
    
    current_order <- order(node_dist[current_node,,drop=T])
    current_node <- current_order[!current_order %in% ordered_nodes][1]
  } 
  
  return(list(onto_plot = tree, ordered_nodes=graph_df[ordered_nodes,]))
  
}
