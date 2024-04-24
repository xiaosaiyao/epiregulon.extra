#' Find interaction partners of a transcription factor of interest
#' @param graph a igraph object from `buildGraph` or `buildDiffGraph`
#' @param focal_tf character string indicating the name of the transcription factors to find
#' interaction partners of
#' @export
#' @return A list with elements corresponding to each transcription factor apart from
#' the focal one. Each list element is represented as a data frame with columns containing
#' names of all target genes shared with focal transcription factor, weights of edges
#' connecting transcription factor with target genes, equivalent weights for focal transcription
#' factor and the element wise product of both weight columns.
#' @examples
#' regulon <- data.frame(tf = sample(letters[1:4], 100, replace = TRUE), idxATAC= 1:100,
#' target = sample(letters[5:14], 100, replace = TRUE))
#' regulon$weights <- runif(100)
#' GRN_graph <- buildGraph(regulon)
#' partners <- findPartners(GRN_graph, 'a')
#' @importFrom igraph V E ego vertex_attr edge_attr delete_vertices subgraph shortest_paths similarity intersection subcomponent rewire each_edge ends
findPartners <- function(graph, focal_tf) {
  checkmate::assert_character(focal_tf, len = 1)
  checkmate::assert_class(graph, "igraph")

  if (length(unique(vertex_attr(graph, name = "type"))) != 2)
    stop("Graph should contain vertices of two types. Input graph should not be created with the use of tripartite mode")
  focal_vertex <- V(graph)[V(graph)$name == focal_tf & V(graph)$type == "transcription factor"]
  if (length(focal_vertex) == 0)
    stop(sprintf("Focal transcription factor (%s) is not present in the input graph"),
         focal_tf)
  focal_tf_targets <- ego(graph, nodes = focal_vertex, mode = "out", mindist = 1)[[1]]

  # keep other transcription factors
  redundant_vertices <- setdiff(as.numeric(V(graph)), 
                                as.numeric(V(graph)[V(graph)$type == "transcription factor"]))

  # keep focal tf's target genes
  redundant_vertices <- setdiff(redundant_vertices, as.numeric(focal_tf_targets))

  # prune network removing target genes which do not belong to the focal tf's
  # regulon
  graph <- delete_vertices(graph, redundant_vertices)
  if (any(duplicated(ends(graph, es = E(graph)))))
    stop("Duplicated edges")
  
  # find focal tf in the new graph
  focal_vertex <- V(graph)[V(graph)$name == focal_tf & V(graph)$type == "transcription factor"]
  other_tfs <- V(graph)[V(graph)$type == "transcription factor" & V(graph)$name !=
                          focal_tf]
  all_tfs <- c(focal_vertex, other_tfs)  # focal tf first on the list
  focal_tf_targets <- V(graph)[V(graph)$type == "target gene"]

  # create regulons using focal tf's target genes
  tf_targets_subgraphs <- lapply(all_tfs, function(tf, targets) 
    subgraph(graph, vids = c(tf, intersection(targets, subcomponent(graph, tf, mode = "out")))),
    focal_tf_targets)
  focal_weights <- edge_attr(tf_targets_subgraphs[[1]], index = E(tf_targets_subgraphs[[1]]), name = "weight")
  common_targets <- lapply(setdiff(seq_along(all_tfs), 1), 
                           get_common_targets, all_tfs, tf_targets_subgraphs, focal_weights, focal_tf_targets)
  do.call(c, common_targets)
}

get_common_targets <- function(idx, tfs, tf_targets_subgraphs, focal_weights, focal_tf_targets){
  res_list <- list()
  tf <- V(tf_targets_subgraphs[[idx]])[V(tf_targets_subgraphs[[idx]])$type == "transcription factor" &
                                       V(tf_targets_subgraphs[[idx]])$name == tfs[idx]$name]
  # determine targets shared between focal tf and current tf
  common_targets <- V(tf_targets_subgraphs[[idx]])[V(tf_targets_subgraphs[[idx]])$type ==
                                                   "target gene"]
  # extract edges connection shared targets with the current tf
  edges_to_tf <- shortest_paths(tf_targets_subgraphs[[idx]], from = tf, to = common_targets,
                                output = "epath")
  edges_to_tf <- do.call(c, edges_to_tf)
  # find indices of common targets in the vector of all targets of focal
  # tf
  common_targets_idx <- match(common_targets$name, focal_tf_targets$name)
  res_list[[tf$name]] <- data.frame(target = common_targets$name, focal_weight = focal_weights[common_targets_idx],
                                    other_tf_weight = edge_attr(tf_targets_subgraphs[[idx]],
                                                                        index = edges_to_tf, name = "weight"))  # extract weights from the edges
  res_list[[tf$name]]$weight_product <- res_list[[tf$name]]$focal_weight *
    res_list[[tf$name]]$other_tf_weight
  res_list
}

#' Calculate Jaccard Similarity between regulons of all transcription factors
#' @param graph a igraph object from `buildGraph` or `buildDiffGraph`
#' @importFrom igraph V similarity
#' @return A matrix with Jaccard similarity between all pairs of transcription factors.
#' @export
#' @examples
#' regulon <- data.frame(tf = sample(letters[1:4], 100, replace = TRUE), idxATAC= 1:100,
#' target = sample(letters[5:14], 100, replace = TRUE))
#' regulon$weights <- runif(100)
#' GRN_graph <- buildGraph(regulon)
#' similarity <- calculateJaccardSimilarity(GRN_graph)
calculateJaccardSimilarity <- function(graph) {
  all_tfs <- V(graph)[V(graph)$type == "transcription factor"]
  res <- similarity(graph, vids = all_tfs, method = "jaccard", mode = "out")
  rownames(res) <- colnames(res) <- V(graph)[all_tfs]$name
  res
}

#' Calculate similarity score from permuted graphs to estimate background similarity
#' @param graph an igraph object from `buildGraph` or `buildDiffGraph`
#' @param focal_tf character string indicating the name of the transcription factors to
#' calculate similarity score
#' @param n an integer indicating the number of permutations
#' @param p a scalar indicating the probability of rewiring the graphs
#' @return A matrix with Jaccard similarity between the focal transcription factor and all pairs of transcription factors
#' for n permuted graphs
#' @importFrom igraph V rewire
#' @export
#' @examples
#' regulon <- data.frame(tf = sample(letters[1:4], 100, replace = TRUE), idxATAC= 1:100,
#' target = sample(letters[5:14], 100, replace = TRUE))
#' regulon$weights <- runif(100)
#' GRN_graph <- buildGraph(regulon)
#' permuted_graph <- permuteGraph(GRN_graph, focal_tf = "a")
permuteGraph <- function(graph, focal_tf, n = 100, p = 1) {
  if (!focal_tf %in% names(V(graph)[V(graph)$type == "transcription factor"]))
    stop(focal_tf, " vertex shoud be present in the graph")
  all_tfs <- names(V(graph)[V(graph)$type == "transcription factor"])
  permute_matrix <- matrix(rep(NA, length(all_tfs) * n),
                           nrow = length(all_tfs))
  rownames(permute_matrix) <- all_tfs

  for (iteration in seq_len(n)) {
    diff_graph_permute <- rewire(graph, each_edge(prob = p))
    similarity_score <- calculateJaccardSimilarity(diff_graph_permute)
    permute_matrix[, iteration] <- similarity_score[focal_tf, all_tfs]
  }
  permute_matrix

}


