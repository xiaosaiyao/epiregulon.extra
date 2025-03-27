regulon <- data.frame(tf = rep(LETTERS[1:5], times = c(5,5,6,3,1)), target = NA)
regulon$target[regulon$tf=="A"] <- LETTERS[6:10]
regulon$target[regulon$tf=="B"][c(1,3,5)] <- regulon$target[regulon$tf=="A"][c(1,2,4)]
regulon$target[regulon$tf=="C"][6] <- regulon$target[regulon$tf=="A"][3]
regulon$target[regulon$tf=="E"] <- "B"
regulon$target[is.na(regulon$target)] <- sample(LETTERS[11:15], sum(is.na(regulon$target)), replace = TRUE)
regulon$weights <- runif(nrow(regulon))
regulon$weights <- regulon$weights * sample(c(1,-1), length(regulon$weights),replace = TRUE)

res <- vector(mode = "list", 4)
names(res) <- LETTERS[2:5]
res <- lapply(res, function(x) data.frame(target=character(0), focal_weight=numeric(0),
                                          other_tf_weight=numeric(0), weight_product=numeric(0)))

res$B <- rbind(res$B, data.frame(target=regulon$target[regulon$tf=="A"][c(1,2,4)],
                                 focal_weight=regulon$weights[regulon$tf=="A"][c(1,2,4)],
                                 other_tf_weight=regulon$weights[regulon$tf=="B"][c(1,3,5)],
                                 weight_product=NA))

res$C <- rbind(res$C, data.frame(target=regulon$target[regulon$tf=="A"][3],
                                 focal_weight=regulon$weights[regulon$tf=="A"][3],
                                 other_tf_weight=regulon$weights[regulon$tf=="C"][6],
                                 weight_product=NA))

res <- lapply(res, function(x) {x$weight_product <- x$focal_weight*x$other_tf_weight; x})

test_graph <- buildGraph(regulon)
partners_graph <- findPartners(test_graph, "A")

test_that("findPartners function works correctly", {
  expect_identical(partners_graph[order(names(partners_graph))], res[order(names(res))])
})

set.seed(4201)

regulon <- data.frame(tf = rep(LETTERS[1:7], times = rmultinom(1,200, prob = c(0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.28))[,1])[sample(200)],
                      target = c(rep(c("F", "G"),3), paste0("gene_", 1:30)[sample(1:30,194, replace=TRUE)])[sample(200)])

# add tfs disjoint to the previous one in terms of their regulons
regulon <- rbind(regulon, data.frame(tf = LETTERS[9:11][sample(1:3,30,replace = TRUE)],
                      target = paste0("gene_", 31:40)[sample(1:10, 30, replace = TRUE)]))

test_graph <- buildGraph(regulon, weights = NULL)
tfs <- unique(regulon$tf)
tf_n <- length(tfs)
js_matrix <- matrix(0, tf_n, tf_n, dimnames = list(unique(regulon$tf), unique(regulon$tf)))
js_matrix <- js_matrix + diag(1, tf_n, tf_n)

for(i in seq_len(tf_n-1)){
  for(j in i:tf_n)
  js_matrix[i,j] <- js_matrix[j,i] <- length(intersect(regulon[regulon$tf==tfs[i], "target"], regulon[regulon$tf==tfs[j], "target"]))/length(unique(regulon[regulon$tf %in% tfs[c(i,j)], "target"]))
}

tf_order <- V(test_graph)$name[V(test_graph)$type=="transcription factor"]
js_matrix <- js_matrix[tf_order, tf_order]

permute_matrix <- permuteGraph(test_graph, n = 90, focal_tf = V(test_graph)$name[V(test_graph)$type=="transcription factor"][1])

test_that("calculateJaccardSimilarity function works correctly", {
  expect_equal(calculateJaccardSimilarity(test_graph), js_matrix)
})

test_that("permuteGraph function works correctly", {
  expect_identical(dim(permute_matrix), c(9L,90L))
  expect_true(all(apply(permute_matrix,2,sd)[apply(permute_matrix,2,function(x) sum(x)!=0)]>0))
  expect_true(!any(is.na(permute_matrix)))
  expect_identical(rownames(permute_matrix), tfs)
})
