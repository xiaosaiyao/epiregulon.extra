#' Test for differential TF activity between pairs of single cell clusters/groups
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param clusters A character or integer vector of cluster or group labels for single cells
#' @param test.type String indicating the type of statistical tests to be passed to scran::findMarkers, can be "t", "wilcox". or "binom"
#' @param pval.type A string specifying how p-values are to be combined across pairwise comparisons for a given group/cluster. For more details see \link[scran]{combineMarkers}.
#' @param direction A string specifying direction of differential TF activity, can be "any", "up" or "down"
#' @param logvalues logical indicating whether activities are computed from logged gene expression or not.
#' If activity is computed from linear values of gene expression, setting logvalues to FALSE will return the difference.
#' If activity is computed from logged values of gene expression, setting logvalues to TRUE will return the log changes.
#' @param ... Further arguments to pass to scran::findMarkers
#'
#' @return A named list of dataframes containing differential TF activity test results for each cluster/group
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' score.combine <- cbind(matrix(runif(2000,0,2), 20,100), matrix(runif(2000,0,10), 20,100))
#' rownames(score.combine) <- paste0("TF",1:20)
#' colnames(score.combine) <- paste0("cell",1:200)
#' cluster <- c(rep(1,100),rep(2,100))
#' markers <- findDifferentialActivity(
#' activity_matrix = score.combine,
#' clusters = cluster,
#' pval.type = "some",
#' direction = "up",
#' test.type = "t")
#' sig.genes <- getSigGenes(markers, fdr_cutoff = 1, summary_cutoff = 0.1)
#' @author Xiaosai Yao, Shang-yang Chen

findDifferentialActivity <- function(activity_matrix,
                                     clusters,
                                     test.type = c("t", "wilcox", "binom"),
                                     pval.type = c("some", "any", "all"),
                                     direction = c("any","up","down"),
                                     logvalues = TRUE,
                                     ...){

  test.type = match.arg(test.type)
  pval.type = match.arg(pval.type)
  direction <- match.arg(direction)
  all_tfs <- rownames(activity_matrix)
  activity_matrix <- stats::na.omit(as.matrix(activity_matrix))
  filtered_tfs <- rownames(activity_matrix)
  if(!identical(all_tfs, filtered_tfs)){
    warning(sprintf("The follwing transciption factors are exclded due to the NA values in the data matrix: "),
            paste(setdiff(all_tfs, filtered_tfs), collapse=", "))
  }
  tf_markers <- scran::findMarkers(activity_matrix, clusters, test.type=test.type,
                                   pval.type=pval.type, direction=direction,
                                   sorted = FALSE, ...)


  if (!isTRUE(logvalues)){

    for (cluster in unique(clusters)) {

      # replace logFC with diff
      colnames(tf_markers[[cluster]]) <- gsub("logFC", "diff", colnames(tf_markers[[cluster]]))
    }
  }
  return(tf_markers)

}

#' Compile and summarize the output from findDifferentialActivity function
#'
#' @param da_list List of dataframes from running findDifferentialActivity
#' @param fdr_cutoff A numeric scalar to specify the cutoff for FDR value. Default is 0.05
#' @param summary_cutoff A numeric scalar to specify the cutoff for log fold change or difference
#' @param topgenes A integer scalar to indicate the number of top ordered genes to include in output
#' @param direction A string specifying direction for which differential TF activity was calculated, can be "any", "up" or "down"
#'
#' @return A compiled dataframe of TFs with differential activities across clusters/groups
#' @export
#'
#' @examples
#' set.seed(1)
#' score.combine <- cbind(matrix(runif(2000,0,2), 20,100), matrix(runif(2000,0,10), 20,100))
#' rownames(score.combine) <- paste0("TF",1:20)
#' colnames(score.combine) <- paste0("cell",1:200)
#' cluster <- c(rep(1,100),rep(2,100))
#' markers <- findDifferentialActivity(score.combine, cluster, pval.type = "some", direction = "up",
#' test.type = "t")
#' sig.genes <- getSigGenes(markers, fdr_cutoff = 1, summary_cutoff = 0.1)
#' utils::head(sig.genes)
#' @author Xiaosai Yao, Shang-yang Chen

getSigGenes <- function(da_list,
                        fdr_cutoff = 0.05,
                        summary_cutoff = NULL,
                        topgenes = NULL,
                        direction = c("any","up","down")){

  direction <- match.arg(direction)
  classes <- names(da_list)

  if (direction %in% c("up","down")) {
    direction_factor <- c(down = -1, up = 1)[direction]}

  top.list <- lapply(seq_along(da_list), function(i){
    da_genes <- as.data.frame(da_list[[i]])
    p.value.name <- colnames(da_genes)[grep("p.value", colnames(da_genes))]
    FDR.name <- colnames(da_genes)[grep("FDR", colnames(da_genes))]
    summary.name <- colnames(da_genes)[grep("summary", colnames(da_genes))]
    da_genes <- da_genes[,c(p.value.name, FDR.name,  summary.name)]

    if (is.null(summary_cutoff)){

      if (direction %in% c("up","down")) {
        summary_cutoff <- signif(stats::quantile(da_genes[,summary.name]*direction_factor, 0.95, na.rm=TRUE), digits = 2 )
      } else{
        summary_cutoff <- signif(stats::quantile(abs(da_genes[,summary.name]), 0.95, na.rm=TRUE), digits = 2)
      }
    }else {
      summary_cutoff <- summary_cutoff

    }

    message ("Using a cutoff of ", summary_cutoff, " for class ", classes[i], " for direction equal to ", direction)
    if (direction %in% c("up","down")) {
      da_genes <- da_genes[which(da_genes[,FDR.name] < fdr_cutoff &
                                   da_genes[, summary.name]*direction_factor > summary_cutoff), ]
    } else {
      da_genes <- da_genes[which(da_genes[,FDR.name] < fdr_cutoff &
                                   abs(da_genes[, summary.name]) > summary_cutoff), ]
    }


    if (nrow(da_genes) != 0){
      da_genes$class <- classes[[i]]
      da_genes$tf <- rownames(da_genes); rownames(da_genes) <- NULL
    }

    if (is.null(topgenes)){
      da_genes <- da_genes[order(da_genes[,FDR.name], -(da_genes[, 3])),]
    } else {
      da_genes <- da_genes[utils::head(order(da_genes[,FDR.name], -(da_genes[, 3])), topgenes),]
    }


    return(da_genes)
  })


  return(do.call(rbind, top.list))

}


regulonEnrich_ <- function(TF, regulon, corr, corr_cutoff, genesets) {
  message(TF)
  regulon.TF <- unique(regulon$target[which(regulon$tf ==
                                              TF & regulon[, corr] >
                                              corr_cutoff)])
  if (length(intersect(regulon.TF, genesets$genes)) < 3 | all(!genesets$genes %in% regulon.TF)) {
    results <- data.frame(p.adjust = NA,
                          Description = NA, GeneRatio = 0,
                          Odds.Ratio = NA)
  } else {
    enrichresults <- clusterProfiler::enricher(regulon.TF,
                                               TERM2GENE = genesets)

    results <- enrichresults@result
    results$GeneRatio <- (as.numeric(lapply(strsplit(results$GeneRatio, split = "/"), "[", 1)))/
      (as.numeric(lapply(strsplit(results$GeneRatio, split = "/"), "[", 2)))
    results$BgRatio <- (as.numeric(lapply(strsplit(results$BgRatio, split = "/"), "[", 1)))/
      (as.numeric(lapply(strsplit(results$BgRatio, split = "/"), "[", 2)))
    results$Odds.Ratio <- results$GeneRatio/results$BgRatio
    results <- results[order(results$p.adjust),]
    results$Description <- factor(as.character(results$Description),
                                  levels = unique(as.character(results$Description[nrow(results):1])))
  }
  return(results)
}


#' Perform geneset enrichment of user-defined regulons
#'
#' @param TF A character vector of TF names
#' @param regulon A matrix of weighted regulon consisting of tf, targets, corr and weight
#' @param weight String indicating the column name that should be used to filter target genes for geneset enrichment. Default is 'weight'.
#' @param weight_cutoff A numeric scalar to indicate the cutoff to filter on the column specified by `weight`. Default is 0.5.
#' @param genesets A dataframe with the first column being the name of the geneset and the second column being the name of the genes
#'
#' @return A dataframe showing the significantly enriched pathways
#' @export
#'
#' @examples
#' #retrieve genesets
#' H <- EnrichmentBrowser::getGenesets(org = 'hsa', db = 'msigdb',
#' cat = 'H', gene.id.type = 'SYMBOL' )
#' C6 <- EnrichmentBrowser::getGenesets(org = 'hsa', db = 'msigdb',
#' cat = 'C6', gene.id.type = 'SYMBOL' )
#'
#' #combine genesets and convert genesets to be compatible with enricher
#' gs <- c(H,C6)
#' gs.list <- do.call(rbind,lapply(names(gs), function(x) {
#' data.frame(gs=x, genes=gs[[x]])}))
#'
#' head(gs.list)
#'
#' #get regulon
#' library(dorothea)
#' data(dorothea_hs, package = 'dorothea')
#' regulon <- dorothea_hs
#' enrichment_results <- regulonEnrich(c('ESR1','AR'), regulon = regulon, weight = 'mor',
#' genesets = gs.list)
#'
#' @author Xiaosai Yao

regulonEnrich <- function(TF, regulon, weight = "weight", weight_cutoff = 0.5,
                          genesets) {
  regulonls <- lapply(TF, function(x) {
    regulonEnrich_(x, regulon, weight, weight_cutoff, genesets)
  })
  names(regulonls) <- TF
  return(regulonls)
}
