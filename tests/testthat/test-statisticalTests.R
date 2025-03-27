set.seed(102001)
n.cells <- 300
n.tf <- 50
n.groups <- 3
mu <- rnorm(n.tf*3, 1, runif(n.tf*3))
X <- unlist(lapply(mu, function(x) rnorm(100,x,runif(1))))
X.mat <- matrix(X, nrow=n.tf, ncol=n.cells, byrow = TRUE)
rownames(X.mat) <- paste0("gene_", 1:50)
da_genes <- findDifferentialActivity(X.mat, clusters = rep(1:3, each=100), direction = "up")
res <- list()
for(i in 1:3){
  threshold_logFC <- round(quantile(da_genes[[i]][,"summary.logFC"], 0.95),1)
  res[[i]] <- da_genes[[i]][da_genes[[i]][,"summary.logFC"]>threshold_logFC & da_genes[[i]][,"FDR"]<0.05,]
  res[[i]] <- as.data.frame(res[[i]][,c("p.value","FDR","summary.logFC")])
  res[[i]]$class <- as.character(i)
  res[[i]]$tf <- rownames(res[[i]])
  rownames(res[[i]]) <- NULL
  res[[i]] <- res[[i]][order(res[[i]]$FDR, -(res[[i]][, 3])),]
}

res <- do.call(rbind, res)

test_that("getSigGenes works correctly", {
  expect_identical(res, getSigGenes(da_genes, direction = "up"))
})

da_genes <- findDifferentialActivity(X.mat, clusters = rep(1:3, each=100), direction = "down")
res <- list()
for(i in 1:3){
  threshold_logFC <- round(quantile(da_genes[[i]][,"summary.logFC"], 0.05),1)
  res[[i]] <- da_genes[[i]][da_genes[[i]][,"summary.logFC"]<threshold_logFC & da_genes[[i]][,"FDR"]<0.05,]
  res[[i]] <- as.data.frame(res[[i]][,c("p.value","FDR","summary.logFC")])
  if(nrow(res[[i]]) == 0) next
  res[[i]]$class <- as.character(i)
  res[[i]]$tf <- rownames(res[[i]])
  rownames(res[[i]]) <- NULL
  res[[i]] <- res[[i]][order(res[[i]]$FDR, res[[i]][, 3]),]
}

res <- do.call(rbind, res)

test_that("getSigGenes works correctly when direction is 'down'", {
  expect_identical(res, getSigGenes(da_genes, direction="down"))
})


set.seed(102001)
mu <- round(exp(runif(n.tf*3)*3)+1)
X <- unlist(lapply(mu, function(x) rpois(100,x)))
X.mat <- matrix(X, nrow=n.tf, ncol=n.cells, byrow = TRUE)
rownames(X.mat) <- paste0("gene_", 1:50)

clusters <- rep(letters[1:3], each=100)
da_genes <- findDifferentialActivity(X.mat, clusters = clusters, direction = "up")
da_genes2 <- findDifferentialActivity(X.mat, clusters = clusters, direction = "up", logvalues=FALSE)
for(i in 1:nrow(X.mat)){
  clust_means <- split(X.mat[i,], clusters) |> lapply(mean)
  for (cluster in unique(clusters)){
    diff_columns <-  da_genes[[cluster]][,grep("logFC", colnames(da_genes[[cluster]]))]
    colnames(diff_columns) <- gsub("logFC","diff", colnames(diff_columns))
    da_genes[[cluster]] <- cbind(da_genes[[cluster]], diff_columns)
    other_clusters <- setdiff(clusters, cluster)
    for(other_cluster in other_clusters){
      da_genes[[cluster]][i,paste0("logFC.", other_cluster)] <- log2(clust_means[[cluster]]/clust_means[[other_cluster]])
    da_genes[[cluster]][i,"summary.logFC"] <- log2(clust_means[[cluster]]/(mean(unlist(clust_means[other_clusters]))))
    da_genes[[cluster]] <- da_genes[[cluster]][,colnames(da_genes2[[cluster]])]
  }

  }
}

for (cluster in unique(clusters)){
  da_genes[[cluster]][,"summary.logFC"] <- setNames(da_genes[[cluster]][,"summary.logFC"], paste0("gene_", 1:50))
  other_clusters <- setdiff(clusters, cluster)
  for(other_cluster in other_clusters){
    da_genes[[cluster]][,paste0("logFC.", other_cluster)] <- setNames(da_genes[[cluster]][,paste0("logFC.", other_cluster)],paste0("gene_", 1:50))
  }
}


test_that("findDifferentialActivity calculates logFC correctly for counts on linear scale", {
  expect_equal(da_genes2, da_genes)
})

set.seed(4563)
genesets <- data.frame(gs = c(rep("geneset_1", 50), rep("geneset_2", 70), rep("geneset_3", 30)),
                       genes = c(paste0("gene_", 1:50), paste0("gene_", sample(50, 25)),
                                 paste0("gene_", 51:95), paste0("gene_", sample(1:95,15)),
                                 paste0("gene_", 96:110)))

regulon <- data.frame(tf = c(rep("A", 130), rep("B", 200)), target = paste0("gene_", 111:440))
regulon <- rbind(regulon, data.frame(tf= rep("A", 20), target = paste0("gene_", sample(1:50, 20))))
regulon <- rbind(regulon, data.frame(tf = rep("A", 10), target = paste0("gene_", sample(96:110, 10))))
regulon <- rbind(regulon, data.frame(tf = rep("B", 40), target = paste0("gene_", sample(51:95, 40))))
regulon$weight <- runif(nrow(regulon))

enrichresults_1 <- clusterProfiler::enricher(regulon[regulon$tf=="A" & regulon$weight > 0.25,"target"], TERM2GENE = genesets)
enrichresults_2 <- clusterProfiler::enricher(regulon[regulon$tf=="B" & regulon$weight > 0.25,"target"], TERM2GENE = genesets)

res <- list ()
for(enrich_results in list(enrichresults_1, enrichresults_2)){
  results <- enrich_results@result
  results$GeneRatio <- (as.numeric(lapply(strsplit(results$GeneRatio,
                                                   split = "/"), "[",
                                          1)))/(as.numeric(lapply(strsplit(results$GeneRatio,
                                                                           split = "/"), "[",
                                                                  2)))
  results$BgRatio <- (as.numeric(lapply(strsplit(results$BgRatio,
                                                 split = "/"), "[",
                                        1)))/(as.numeric(lapply(strsplit(results$BgRatio,
                                                                         split = "/"), "[",
                                                                2)))
  results$Odds.Ratio <- results$GeneRatio/results$BgRatio
  results <- results[order(results$p.adjust),
  ]
  results$Description <- factor(as.character(results$Description),
                                levels = unique(as.character(results$Description[nrow(results):1])))
  res <- c(res, list(results))

}

names(res) <- c("A", "B")

test_that("regulonEnrich works correctly", {
  expect_equal(regulonEnrich(LETTERS[1:2], regulon, weight = "weight", weight_cutoff = 0.25,
                             genesets), res)
})
