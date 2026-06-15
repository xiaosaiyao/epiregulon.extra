set.seed(102001)
n.cells <- 300
n.tf <- 50
n.groups <- 3
mu <- rnorm(n.tf*3, 1, runif(n.tf*3))
X <- unlist(lapply(mu, function(x) rnorm(100,x,runif(1))))
X.mat <- matrix(X, nrow=n.tf, ncol=n.cells, byrow=TRUE)
rownames(X.mat) <- paste0("gene_", 1:50)
da_genes <- findDifferentialActivity(X.mat, clusters=rep(1:3, each=100), direction="up")
res <- list()

for (i in 1:3){
  threshold_logFC <- round(quantile(da_genes[[i]][,"logFC"], 0.95),1)
  res[[i]] <- da_genes[[i]][da_genes[[i]][,"logFC"]> threshold_logFC & 
                              da_genes[[i]][,"log.FDR"] < log(0.05),]
  res[[i]] <- as.data.frame(res[[i]][,c("log.p.value","log.FDR","logFC")])
  res[[i]]$class <- as.character(i)
  res[[i]]$tf <- rownames(res[[i]])
  rownames(res[[i]]) <- NULL
  res[[i]] <- res[[i]][order(res[[i]]$log.FDR, -(res[[i]][, 3])),]
}

res <- do.call(rbind, res)

test_that("getSigGenes works correctly", {
  expect_identical(res, getSigGenes(da_genes, direction="up"))
})

da_genes <- findDifferentialActivity(X.mat, clusters=rep(1:3, each=100), direction="down")
res <- list()
for (i in 1:3){
  threshold_logFC <- round(quantile(da_genes[[i]][,"logFC"], 0.05),1)
  res[[i]] <- da_genes[[i]][da_genes[[i]][,"logFC"]<threshold_logFC & da_genes[[i]][,"log.FDR"]< log(0.05),]
  res[[i]] <- as.data.frame(res[[i]][,c("log.p.value","log.FDR","logFC")])
  if(nrow(res[[i]]) == 0) next
  res[[i]]$class <- as.character(i)
  res[[i]]$tf <- rownames(res[[i]])
  rownames(res[[i]]) <- NULL
  res[[i]] <- res[[i]][order(res[[i]]$log.FDR, res[[i]][, 3]),]
}

res <- do.call(rbind, res)

test_that("getSigGenes works correctly when direction is 'down'", {
  expect_identical(res, getSigGenes(da_genes, direction="down"))
})



set.seed(4563)
genesets <- data.frame(gs=c(rep("geneset_1", 50), rep("geneset_2", 70), rep("geneset_3", 30)),
                       genes=c(paste0("gene_", 1:50), paste0("gene_", sample(50, 25)),
                                 paste0("gene_", 51:95), paste0("gene_", sample(1:95,15)),
                                 paste0("gene_", 96:110)))

regulon <- data.frame(tf=c(rep("A", 130), rep("B", 200)), target=paste0("gene_", 111:440))
regulon <- rbind(regulon, data.frame(tf= rep("A", 20), target=paste0("gene_", sample(1:50, 20))))
regulon <- rbind(regulon, data.frame(tf=rep("A", 10), target=paste0("gene_", sample(96:110, 10))))
regulon <- rbind(regulon, data.frame(tf=rep("B", 40), target=paste0("gene_", sample(51:95, 40))))
regulon$weight <- runif(nrow(regulon))

enrichresults_1 <- enricher(regulon[regulon$tf=="A" & regulon$weight > 0.25,"target"], TERM2GENE=genesets)
enrichresults_2 <- enricher(regulon[regulon$tf=="B" & regulon$weight > 0.25,"target"], TERM2GENE=genesets)

res <- list ()
for (enrich_results in list(enrichresults_1, enrichresults_2)){
  results <- enrich_results@result
  results$GeneRatio <- (as.numeric(lapply(strsplit(results$GeneRatio, split="/"), "[", 1)))/
    (as.numeric(lapply(strsplit(results$GeneRatio,split="/"), "[", 2)))
  results$BgRatio <- (as.numeric(lapply(strsplit(results$BgRatio, split="/"), "[", 1)))/
    (as.numeric(lapply(strsplit(results$BgRatio,split="/"), "[", 2)))
  results$Odds.Ratio <- results$GeneRatio/results$BgRatio
  results <- results[order(results$p.adjust),]
  results$Description <- factor(as.character(results$Description),
                                levels=unique(as.character(results$Description[nrow(results):1])))
  res <- c(res, list(results))
  
}

names(res) <- c("A", "B")

test_that("regulonEnrich works correctly", {
  expect_equal(regulonEnrich(LETTERS[1:2], regulon, weight="weight", weight_cutoff=0.25,
                             genesets), res)
})


set.seed(100)
n_cells=2000
n_tfs <- 100
n_clusters <- 5
clusters <- sample(letters[1:n_clusters], 2000, replace=TRUE)
unique_clusters <- unique(clusters)
tf_names <- paste0("TF_", 1:n_tfs)
tf_names <- sample(tf_names) # random order
tf_to_cluster_matrix <- matrix(1, nrow=n_tfs, ncol=n_cells)
for (i in seq_len(n_tfs)){
  tf_to_cluster_matrix[i, which(clusters==sample(unique_clusters,1))] <- 2
}
activity_matrix <- matrix(nrow=n_tfs, ncol=n_cells)
activity_matrix[1:length(activity_matrix)] <- rnorm(length(activity_matrix), as.numeric(tf_to_cluster_matrix), 1)

colnames(activity_matrix) <- paste0("cell_", 1:n_cells)
rownames(activity_matrix) <- tf_names
tf_markers <- scran::findMarkers(activity_matrix, 
                                 clusters, 
                                 pval.type="some", 
                                 direction="up", 
                                 sorted=FALSE, 
                                 log.p=TRUE
)

tf_markers_simple <- findDifferentialActivity(activity_matrix, clusters, pval.type="some", direction="up")
test_that("findDifferentialActivity preserves the gene order", {
  expect_equal(rownames(tf_markers_simple[[1]]), tf_names)
})
test_that("findDifferentialActivity calculates the pvalues correctly", {
  expect_equal(tf_markers_simple[[1]]$log.p.value, tf_markers[[1]]$log.p.value)
})
test_that("findDifferentialActivity calculates the FDR correctly", {
  expect_equal(tf_markers_simple[[1]]$log.FDR, tf_markers[[1]]$log.FDR)
})
test_that("findDifferentialActivity calculates the logFC correctly", {
  expect_equal(tf_markers_simple[[1]]$logFC, tf_markers[[1]]$summary.logFC)
})
