library(vdiffr)
set.seed(20913)
example_sce <- scuttle::mockSCE()
example_sce <- scuttle::logNormCounts(example_sce)
example_sce <- scater::runPCA(example_sce)
example_sce <- scater::runUMAP(example_sce)
example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)


test_that("plotActivityDim works correctly",{
  expect_doppelganger(
    title = "Activity plot reduced dim",
    fig = plotActivityDim(sce = example_sce, activity = logcounts(example_sce),
                           tf = c("Gene_0001","Gene_0002"),  label = "cluster")

  )
})

test_that("plotActivityViolin works correctly",{
  expect_doppelganger(
    title = "Violin plot",
    fig = plotActivityDim(sce = example_sce, activity = logcounts(example_sce),
                          tf = c("Gene_0001","Gene_0002"),  label = "cluster")
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

regulon <- data.frame(tf=c(rep("Gene_0001",10),rep("Gene_0002",20)),
target <- sample(rownames(example_sce),30), weight = rnorm(30))

test_that("plotHeatmapRegulon works correctly",{
  expect_doppelganger(
    title = "Heatmap regulon",
    fig = plotHeatmapRegulon(example_sce, tfs=c("Gene_0001","Gene_0002"), regulon=regulon,
                             cell_attributes="cluster", col_gap = "cluster", column_title_rot = 90)
  )
})


activity_matrix <- matrix(rnorm(10*200), nrow=10, ncol=200)
rownames(activity_matrix) <- sample(rownames(example_sce),10)

test_that("plotHeatmapActivity works correctly",{
  expect_doppelganger(
    title = "Heatmap activity",
    fig = plotHeatmapActivity(activity_matrix=activity_matrix, sce=example_sce,
                              tfs=rownames(activity_matrix), cell_attributes="cluster", col_gap="cluster")
  )
})


