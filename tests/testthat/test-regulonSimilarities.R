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
