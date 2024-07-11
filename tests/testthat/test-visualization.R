library(vdiffr)
set.seed(20913)
example_sce <- scuttle::mockSCE()
example_sce <- scuttle::logNormCounts(example_sce)
example_sce <- scater::runPCA(example_sce)
example_sce <- scater::runUMAP(example_sce)
set.seed(20913)
example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)


test_that("plotActivityDim works correctly",{
  expect_doppelganger(
    title = "Activity plot reduced dim",
    fig = plotActivityDim(sce = example_sce, activity = logcounts(example_sce),
                           tf = c("Gene_0001","Gene_0002"),  label = "cluster")

  )
})

activity.matrix <- matrix(c(runif(100)*1.3, runif(100)),byrow=TRUE, nrow=2)
rownames(activity.matrix) <- c("TF1", "TF2")
activity.matrix[, 1:50] <- activity.matrix[, 1:50]*(runif(100)+0.8)

test_that("plotActivityViolin works correctly",{
  expect_doppelganger(
    title = "Violin plot",
    fig = plotActivityViolin(activity.matrix, clusters = rep(c("a", "b"), each=50),
                          tf = c("TF1","TF2"))
  )
})

test_that("plotBubble works correctly",{
expect_doppelganger(
  title = "Bubble plot",
  fig = plotBubble(activity_matrix = logcounts(example_sce),
                   tf = c("Gene_0001","Gene_0002"),  clusters = example_sce$cluster)
)
})


#retrieve genesets
H <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
cat = "H", gene.id.type = "SYMBOL" )
C6 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
cat = "C6", gene.id.type = "SYMBOL" )

#combine genesets and convert genesets to be compatible with enricher
gs <- c(H,C6)
gs.list <- do.call(rbind,lapply(names(gs), function(x) {
data.frame(gs=x, genes=gs[[x]])}))

#get regulon
library(dorothea)
data(dorothea_hs, package = "dorothea")
regulon <- dorothea_hs
enrichment_results <- regulonEnrich(c("ESR1","AR"), regulon = regulon, weight = "mor",
genesets = gs.list)

test_that("enrichPlot works correctly",{
  expect_doppelganger(
    title = "Enrichment plot",
    fig = enrichPlot(results = enrichment_results)
  )
})
