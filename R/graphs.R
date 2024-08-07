#' Creating graphs and related operations
#'
#' @description
#' The function enable to create graph objects using as input regulon objects returned by
#' `pruneRegulon` or `addWeights`. Both weighted and unweighted graphs can be
#' created that can further be visualized using dedicated functions.
#'
#' @details
#' \code{buildGraph} function creates a directed graph based on the output of
#' the \code{getRegulon} function. Four modes are available: (1) `tg` in which
#' connections are made directly between transcription factor and target genes. Even if
#' the same tf-tg pair is connected in the original regulon object through many
#' regulatory elements then only one edge is created. In such a case, when weighted
#' graph is created, weights are summarized by the aggregating function (by default
#' the maximum absolute value with the sign of the original value). Similarly, aggregation
#' is made in the `re` mode leaving only unique transcription factor-regulatory element pairs.
#' In `tripartite` mode edges connect transcription factors with regulatory elements and
#' regulatory elements with target genes. The same weights are used for both edges
#' that correspond to the single row in the regulon data frame (tf-re and re-tg). Note
#' that the original regulon structure is not fully preserved because each row is now
#' represented by two edges which are independent from each other. Thus they can be
#' coupled with different edges connected to the same regulatory element building the
#' path from transcription factor to the target gene of another transcription factor
#' through the shared regulatory element.
#'
#' \code{buildDiffGraph} a graph difference by subtracting the edges of \code{graph_obj_2}
#' from those of the \code{graph_obj_1}. If \code{weighted} is set to \code{TRUE} then for each
#' ordered pair of vertices (nodes) the difference in number of edges between \code{graph_obj_1}
#' and \code{graph_obj_1} is calculated. The result is used to set the number of
#' corresponding edges in output graph. Note that unless \code{abs_diff} is set to
#' \code{TRUE} any non-positive difference will translate into lack of the edges
#' for a corresponding ordered pair of vertices in the output graph (equivalent
#' to 0 value in the respective position in adjacency matrix). In case of
#' weighted graphs, the weight of the output graph is calculated as a difference
#' of the corresponding weights between input graphs.
#'
#' \code{addCentrality} calculates degree centrality for each vertex using
#' \code{igraph::strength}.
#'
#' With \code{normalizeCentrality} function the normalized values of centrality
#' are calculated from the original ones divided by
#' \code{FUN}(total number of non-zero edges associated with each node).
#'
#' \code{rankTfs} assign ranks to transcription factors according to degree
#' centrality of their vertices
#'
#' @param regulon an object returned by the getRegulon or addWeights function.
#' @param mode a character specifying which type of graph will be built. In \code{'tg'} mode
#' a bipartite graph is built by connecting transcription factors directly to the target genes
#' and ignoring information about mediating regulatory elements; in \code{'pairs'} mode
#' transcription factors are connected to unique target gene-regulatory element pairs;
#' in \code{'tripartite'} mode the network is made up of three types of vertices (nodes):
#' transcription factors, regulatory elements and target genes; here the path from
#' target gene to regulatory element always contains a regulatory element; in
#' \code{'re'} mode data in the target genes is dropped and only connections are
#' between transcription factors and regulatory elements.
#' @param graph,graph_obj_1,graph_obj_2  an igraph object.
#' @param weights a character specifying which variable should be used to assign
#' weights to edges.
#' @param cluster a character specifying the name of the cluster column which should be used
#' to retrieve weight values from regulon object. Using this argument makes sense
#' only with combination with \code{weights} parameter when it points to the regulon column
#' that is a matrix.
#' @param aggregation_function a function used to aggregate weights of duplicated edges,
#'  which might appear due to the many transcription factor converging at the same regulatory
#'  element; starting from this point each transcription factor is supposed to have a separate
#'  connection to the target gene, perhaps the same one across several connections. In
#'  \code{tripartite} mode this might result in many edges in the same node pair, however weights might
#'  differ since they are inherited from different tf-re-tg triplets (rows) in the
#'  \code{regulon} object. Similarly, duplicated edges are generated by one
#'  transcription factor using a regulatory element multiple times to reach different
#'  target genes. In \code{tg} mode the edges became duplicated if one transcription
#'  factor reaches the same target genes through many regulatory elements.
#' @param weighted a logical indicating whether weighted graphs are used; in \code{tripartite}
#' mode tf-re-tg triplet is decomposed into two edges corresponding to tf-re and re-tg pairs,
#' and both edges inherit the same weight, which was originally assigned to the parent triplet.
#' @param abs_diff a logical indicating whether absolute difference in the number
#' edges or their weights will be calculated.
#' @param FUN a function used for normalization. The input to this
#' function is be the number of edges connected with each node (incident edges).
#' @param type_attr a character corresponding to the name of the vertex attribute
#' which indicate the type of vertex.
#' @param na_replace a logical indicating whether NA values for weights should be
#' replaced with zeros.
#' @param keep_original_names A logical indicating whether gene names should be used as
#' node names in the output graph. Note that this might lead to the duplicated node
#' names if the same gene is present in two layers (transcription factors and target genes).
#' @param filter_edges A numeric defining the cutoff weight used for filtering out edges
#' which have weights equal or greater than cutoff. The isolated vertices are removed then
#' from the graph. Defaults to NULL in which case no filtering is applied.
#' @return an igraph object.
#' \code{rankTfs} returns a data.frame with transcription factors sorted according
#' to the value of the \code{centrality} attribute.
#' @examples
#' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[1:10], LETTERS[1:10]),1,  paste, collapse = '')
#' regulon <- DataFrame(tf = sample(tf_set, 5e3, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[1:10], LETTERS[1:10], LETTERS[1:10])
#' gene_set <- apply(gene_set,1,function(x) paste0(x,collapse=''))
#' regulon$target <- sample(gene_set, 5e3, replace = TRUE)
#' regulon$idxATAC <- 1:5e3
#' regulon$corr <- runif(5e3)*0.5+0.5
#' regulon$weights <- matrix(runif(15000), nrow=5000, ncol=3)
#' colnames(regulon$weights) <- c('all','cluster1', 'cluster2')
#' graph_tripartite <- buildGraph(regulon, cluster='all', mode = 'tripartite')
#'
#' # build bipartite graph using regulatory element-target gene pairs
#' graph_pairs_1 <- buildGraph(regulon, cluster = 'cluster1', mode = 'pairs')
#' graph_pairs_2 <- buildGraph(regulon, cluster = 'cluster2', mode = 'pairs')
#' graph_diff <- buildDiffGraph(graph_pairs_1, graph_pairs_2)
#' graph_diff <- addCentrality(graph_diff)
#' graph_diff <- normalizeCentrality(graph_diff)
#' tf_ranking <- rankTfs(graph_diff)
#'
#' @importFrom igraph V E graph_from_data_frame vcount delete_edges delete_vertices degree graph_from_adjacency_matrix as_adjacency_matrix V<- incident_edges strength list.edge.attributes list.vertex.attributes
#' @export
buildGraph <- function(regulon,
                       mode = c("tg", "tripartite", "re", "pairs"),
                       weights = "weights",
                       cluster = "all",
                       aggregation_function = function(x) x[which.max(abs(x))],
                       na_replace = TRUE,
                       keep_original_names = TRUE,
                       filter_edges = NULL) {
  
  if (!is.null(weights) && !weights %in%
      colnames(regulon))
    stop(sprintf("%s column should be present in the regulon",
                 weights))
  mode <- match.arg(mode)
  # give names to the peaks and target genes which will be easy to extract
  regulon$idxATAC <- paste0(as.character(regulon$idxATAC),
                            "_peak")
  regulon$target <- paste0(regulon$target,
                           "_target_gene")
  vertex_columns <- switch(mode,
                           re = c("tf", "idxATAC"),
                           pairs = c("tf", "idxATAC","target"),
                           tripartite = c("idxATAC","target"),
                           tg = c("tf","target"))
  
  graph_data <- regulon[, vertex_columns]
  
  if (!is.null(weights)) {
    if (is.matrix(regulon[,weights])) {
      stopifnot(cluster %in% colnames(regulon[,weights]))
      weights_df <- data.frame(regulon[, weights][, cluster])
    } else weights_df <- data.frame(regulon[,weights])
    if (any(is.na(weights_df)) & na_replace) {
      message("Replacement of na values for weights with 0")
      weights_df[[1]][is.na(weights_df[[1]])] <- 0
    }
    colnames(weights_df) <- weights
    graph_data <- cbind(graph_data,
                        weights_df)
  }
  message(sprintf("Building graph using %s as edge weights",
                  weights))
  if (mode == "tripartite") {
    # add tf-re data
    colnames(graph_data) <- c("from", "to", weights)
    graph_data_tf_re <- data.frame(from = regulon$tf,
                                   to = regulon$idxATAC)
    if (!is.null(weights)) {
      graph_data_tf_re[, weights] <- weights_df[[1]]
    }
    graph_data <- rbind(graph_data,
                        graph_data_tf_re)
    rm(graph_data_tf_re)
    vertex_columns <- c("from", "to")
  }
  
  if (mode == "pairs") {
    # create node names corresponding to re-tg pairs
    graph_data$target <- paste(graph_data$idxATAC,
                               graph_data$target, sep = "@")
    graph_data <- graph_data[,c("tf", "target", weights)]
    vertex_columns <- c("tf","target")
  }
  
  if (is.null(weights)) {
    # avoid duplicated edges in the case of unweighted graph
    graph_data <- unique(graph_data)
  } else {
    colnames(graph_data)[colnames(graph_data) ==
                           weights] <- "weight"
    graph_data <- aggregate_edges(graph_data,
                                  aggregation_function)
  }
  
  epiregulon_graph <- graph_from_data_frame(graph_data)
  if (mode == "tripartite") {
    layer_numb <- rep(1, vcount(epiregulon_graph))
    layer_numb[grepl("_target_gene$",
                     V(epiregulon_graph)$name)] <- 2
    layer_numb[grepl("_peak$",
                     V(epiregulon_graph)$name)] <- 3
    V(epiregulon_graph)$layer <- layer_numb
  }
  # set 'type' attribute for vertices required by bipartite graphs
  vertex_type <- rep("transcription factor",
                     vcount(epiregulon_graph))
  vertex_type[grepl("_peak$",
                    V(epiregulon_graph)$name)] <- "regulatory element"
  vertex_type[grepl("_target_gene$",
                    V(epiregulon_graph)$name)] <- "target gene"
  V(epiregulon_graph)$type <- vertex_type
  # transform character constants to numeric values for later use by graphics
  # functions
  V(epiregulon_graph)$type.num <- match(V(epiregulon_graph)$type,
                                        c("transcription factor",
                                          "peak", "target gene"))
  
  # restore original names
  if (keep_original_names) {
    V(epiregulon_graph)$name <- gsub("_target_gene$|_peak$",
                                     "", V(epiregulon_graph)$name)
  }
  if (!is.null(filter_edges) &&
      !is.null(weights)) {
    epiregulon_graph <- delete_edges(epiregulon_graph,
                                     E(epiregulon_graph)[E(epiregulon_graph)$weight <=
                                                           filter_edges])
    epiregulon_graph <- delete_vertices(epiregulon_graph,
                                        V(epiregulon_graph)[degree(epiregulon_graph) ==
                                                              0])
  }
  epiregulon_graph
}

#' @rdname buildGraph
#' @export
buildDiffGraph <- function(graph_obj_1, graph_obj_2, weighted = TRUE, abs_diff = TRUE) {
  checkmate::assertClass(graph_obj_1, "igraph")
  checkmate::assertClass(graph_obj_2, "igraph")
  if (!identical(V(graph_obj_1)$name, V(graph_obj_2)$name)) {
    stop("The nodes should be the same in both graphs")
  }
  transformation_function <- ifelse(abs_diff, abs, identity)
  if (weighted) {
    res <- graph_from_adjacency_matrix(transformation_function(as_adjacency_matrix(graph_obj_1, attr = "weight") -
                                                                 as_adjacency_matrix(graph_obj_2, attr = "weight")), weighted = TRUE)
    # remove zero-weight edges
    res <- delete_edges(res, E(res)[E(res)$weight == 0])
  } else {
    res <- graph_from_adjacency_matrix(abs(as_adjacency_matrix(graph_obj_1) - as_adjacency_matrix(graph_obj_2)),
                                       weighted = FALSE)
  }
  
  if (!identical(V(graph_obj_1)$type, V(graph_obj_2)$type)) {
    warning("Types of nodes differ between graphs. Only those from the first graph are used.")
  }
  V(res)$type <- V(graph_obj_1)$type
  V(res)$type.num <- V(graph_obj_1)$type.num
  
  # remove nodes with no edges
  edge_numbers <- vapply(incident_edges(res, V(res), mode = "all"), length, FUN.VALUE = numeric(1))
  res <- delete_vertices(res, V(res)[edge_numbers == 0])
  res
}

#' @rdname buildGraph
#' @export
addCentrality <- function(graph) {
  checkmate::assertClass(graph, "igraph")
  V(graph)$centrality <- strength(graph)
  graph
}

#' @rdname buildGraph
#' @export
normalizeCentrality <- function(graph, FUN = sqrt, weighted = TRUE) {
  checkmate::assertClass(graph, "igraph")
  if (!"centrality" %in% list.vertex.attributes(graph))
    stop("Vertices do not have 'centrality' attribute")
  if (!"weight" %in% list.edge.attributes(graph) & weighted)
    stop("Set 'weight' attribute to edges or use with 'weighted = FALSE'")
  
  # calculate number of edges for each node
  edge_numbers <- vapply(incident_edges(graph, V(graph), mode = "all"),
                         length, FUN.VALUE = numeric(1))
  
  V(graph)$centrality <- V(graph)$centrality/FUN(edge_numbers)
  graph
}

#' @rdname buildGraph
#' @export
rankTfs <- function(graph, type_attr = "type") {
  checkmate::assertClass(graph, "igraph")
  rank_df <- data.frame(tf = V(graph)$name[order(V(graph)$centrality[vertex_attr(graph,type_attr) == "transcription factor"],
                                                 decreasing = TRUE)],
                        centrality = sort(V(graph)$centrality[vertex_attr(graph, type_attr) == "transcription factor"],
                                          decreasing = TRUE))
  rank_df$rank <- base::rank(-rank_df$centrality)
  rank_df
}

#' Plot a graph build based on \code{getRegulon} output
#'
#' This function takes an input an igraph object created by any of the following:
#' \code{buildGraph}, \code{addCentrality}, \code{igraph::strength}, \code{normalizeCentrality}.
#' It makes a force-directed layout plot to visualize it at a high level.
#'
#' @param graph an igraph object
#' @param layout a layout specification. Any values that are valid for
#' \link[ggraph]{ggraph} or \link[ggraph]{create_layout} will work. Defaults to
#' 'stress'. Consider also trying 'mds', 'nicely', and 'fr' while you experiment.
#' @param label_size an integer indicating how large the labels of highlighted
#' transcription factors should be
#' @param tfs_to_highlight a character vector specifying which TFs in the plot
#' should be highlighted. Defaults to NULL (no labels).
#' @param edge_alpha a numeric value between 0 and 1 indicating the level of
#' transparency to use for the edge links in the force-directed layout. Defaults
#' to 0.02.
#' @param point_size a numeric value indicating the size of nodes in the force-directed layout
#' @param point_border_size a numeric value indicating the size of point
#' borders for nodes in the force-directed layout
#' @param label_alpha a numeric value between 0 and 1 indicating the level of
#' transparency to use for the labels of highlighted nodes
#' @param label_nudge_x a numeric value indicating the shift of the labels
#' along the x axis that should be used in the force-directed layout
#' @param label_nudge_y A numeric value indicating the shift of the labels
#' along the y axis that should be used in the force-directed layout.
#' @param ... optional additional arguments to pass to \link[ggraph]{create_layout}
#' @return a ggraph object
#' @importFrom ggraph create_layout ggraph geom_edge_link geom_node_label
#' @importFrom ggplot2 aes theme_void labs
#' @author Timothy Keyes, Tomasz Wlodarczyk
#' @examples
#' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[seq_len(5)], LETTERS[seq_len(5)]),1,  paste, collapse = '')
#' regulon <- data.frame(tf = sample(tf_set, 5e2, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[seq_len(5)], LETTERS[seq_len(5)], LETTERS[seq_len(5)])
#' gene_set <- apply(gene_set,1,function(x) paste0(x,collapse=''))
#' regulon$target <- sample(gene_set, 5e2, replace = TRUE)
#' regulon$idxATAC <- seq_len(5e2)
#' regulon$corr <- runif(5e2)*0.5+0.5
#' regulon$weights <- runif(500)
#' #create igraph object
#' graph_tripartite <- buildGraph(regulon, mode = 'tripartite')
#' plotEpiregulonNetwork(graph_tripartite, tfs_to_highlight = sample(unique(tf_set),3),
#' edge_alpha = 0.2)
#' @export
plotEpiregulonNetwork <- function(graph,
                                  layout = "stress", label_size = 3,
                                  tfs_to_highlight = NULL, edge_alpha = 0.02,
                                  point_size = 1, point_border_size = 0.5,
                                  label_alpha = 0.8, label_nudge_x = 0.2,
                                  label_nudge_y = 0.2, ...) {
  checkmate::assertClass(graph,
                         "igraph")
  my_layout <- create_layout(graph,
                             layout = layout, ...)
  highlighted <- my_layout[my_layout$name %in%
                             tfs_to_highlight, ]
  my_plot <- ggraph(graph = my_layout) +
    geom_edge_link(alpha = edge_alpha) +
    geom_node_point(aes(fill = type),
                    shape = 21, size = point_size,
                    stroke = point_border_size) +
    geom_node_label(aes(label = name),
                    data = highlighted, alpha = label_alpha,
                    nudge_x = label_nudge_x,
                    nudge_y = label_nudge_y,
                    size = label_size) +
    geom_node_point(aes(fill = type),
                    shape = 21, data = highlighted,
                    size = 3, stroke = point_border_size) +
    theme_void() + labs(fill = NULL)
  
  return(my_plot)
}

#' Plot graph according to grouping factor
#'
#' Plot graph with separate weights for different levels of the grouping factor
#'
#' @param regulon an object returned by the getRegulon or addWeights function
#' @param cutoff a numeric used to select values of the variables passed in `clusters`
#' parameter. Values greater than `cutoff` are retained and used as
#' graph edge weights.
#' @param tf a character vector storing the names of transcription factors to be
#' included in the graph
#' @param weight a string indicating the name of the column in the regulon to be used as
#' the weight of the edges
#' @param clusters a character vector indicating the clusters to be plotted
#' @param layout a layout specification. Any values that are valid for
#' \link[ggraph]{ggraph} or \link[ggraph]{create_layout} will work.
#' @author Xiaosai Yao, Tomasz Wlodarczyk
#' @return a ggraph object
#' @examples
#' #' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[1:10], LETTERS[1:10]),1,  paste, collapse = '')
#' regulon <- S4Vectors::DataFrame(tf = sample(tf_set, 5e3, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[1:10], LETTERS[1:10], LETTERS[1:10])
#' gene_set <- apply(gene_set,1,function(x) paste0(x,collapse=''))
#' regulon$target <- sample(gene_set, 5e3, replace = TRUE)
#' regulon$idxATAC <- 1:5e3
#' regulon$weight <- cbind(data.frame(C1 = runif(5e3), C2 = runif(5e3),
#' C3 = runif(5e3)))
#' plotDiffNetwork(regulon, tf = unique(tf_set)[1:3],
#' clusters = c('C1', 'C2', 'C3'), cutoff = 0.2)
#' @export

plotDiffNetwork <- function(regulon, cutoff = 0.01, tf = NULL,
                            weight = "weight", clusters, layout = "stress") {
  regulon.tf <- list()
  for (cluster in clusters) {
    # apply cutoff
    idx <- which(regulon$tf %in% tf & regulon[, weight][,cluster] > cutoff)
    regulon_cluster <- regulon[idx, c("tf", "target", weight)]
    
    # rename colnames as weight to be consistent across all groups
    regulon_cluster$weight <- regulon_cluster[, weight][,cluster]
    
    # rename tf to be tf_cluster
    regulon_cluster[, "tf"] <- paste0(regulon_cluster[,"tf"], "_", cluster)
    regulon.tf[[cluster]] <- regulon_cluster
  }
  
  combined.regulon <- do.call("rbind", regulon.tf)
  
  combined.graph <- buildGraph(combined.regulon, mode = "tg",
                               weights = "weight")
  
  plotEpiregulonNetwork(combined.graph, layout = layout,
                        tfs_to_highlight = unique(combined.regulon$tf), label_nudge_x = 0.1,
                        label_nudge_y = 0.1)
}

aggregate_edges <- function(graph_data, FUN) {
  grouping_factors <- paste(setdiff(colnames(graph_data), "weight"),
                            collapse = "+")
  aggregation_formula <- eval(parse(text = paste0("weight", "~",
                                                  grouping_factors)))
  stats::aggregate(graph_data, aggregation_formula, FUN)
}

