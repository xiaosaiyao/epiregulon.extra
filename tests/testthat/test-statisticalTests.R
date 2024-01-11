set.seed(102001)
n.cells <- 300
n.tf <- 50
n.groups <- 3
mu <- rnorm(n.tf*3, 1, runif(n.tf*3))
X <- unlist(lapply(mu, function(x) rnorm(100,x,runif(1))))
X.mat <- matrix(X, nrow=n.tf, ncol=n.cells, byrow = TRUE)
rownames(X.mat) <- paste0("gene_", 1:50)
da_genes <- findDifferentialActivity(X.mat, clusters = rep(1:3, each=100))
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
  expect_identical(res, getSigGenes(da_genes))
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
