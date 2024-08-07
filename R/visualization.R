#' @importFrom ggplot2 element_text scale_color_gradient
plotActivityDim_ <- function(sce, activity_matrix, tf, dimtype,
    label, legend.label, colors, limit, ...) {

    tf.activity <- as.numeric(activity_matrix[tf, ])
    sce$activity <- tf.activity

    g <- scater::plotReducedDim(sce, dimred = dimtype, colour_by = "activity",
        text_by = label, ...)

    g <- g + scale_color_gradient(low = colors[1], high = colors[2],
        limit = limit, oob = scales::squish) + ggtitle(tf) +
        labs(color = legend.label) + theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5))

    return(g)

}

#' Plot cell-level reduced dimension results stored in a SingleCellExperiment object, colored by activities for a list of TFs
#'
#' @param sce  A SingleCellExperiment object containing dimensionality reduction coordinates
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param dimtype String indicating the name of dimensionality reduction matrix to be extracted from the SingleCellExperiment
#' @param label String corresponding to the field in the colData of sce for annotation on plot
#' @param ncol A integer to specify the number of columns in the combined plot, if combine == TRUE
#' @param nrow A integer to specify the number of rows in the combined plot, if combine == TRUE
#' @param title A string to specify the name of the combined plot
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#' @param legend.label String indicating the name of variable to be plotted on the legend
#' @param colors A vector of 2 colors for the intensity, with the first element referring to the lower value and
#' the second element referring to the higher value. Default is c('blue','yellow').
#' @param limit A vector of lower and upper bounds for the color scale. The default option is NULL and will adjust
#' to minimal and maximal values
#' @param ... Additional arguments from scater::plotReducedDim
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @importFrom SummarizedExperiment colData
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce <- scater::runPCA(example_sce)
#' example_sce <- scater::runUMAP(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' plotActivityDim(sce = example_sce, activity = logcounts(example_sce),
#' tf = c('Gene_0001','Gene_0002'),  label = 'cluster')

#' @author Xiaosai Yao, Shang-yang Chen
#'
plotActivityDim <- function(sce = NULL,
    activity_matrix, tf, dimtype = "UMAP",
    label = NULL, ncol = NULL, nrow = NULL,
    title = NULL, combine = TRUE, legend.label = "activity",
    colors = c("blue", "yellow"), limit = NULL,
    ...) {

    # give warning for genes absent in tf list
    missing <- tf[which(!tf %in% rownames(activity_matrix))]
    if (!identical(missing, character(0))) {
        writeLines(paste0(missing,
            " not found in activity matrix. Excluded from plots"))
    }

    tf <- tf[which(tf %in% rownames(activity_matrix))]


    gs <- lapply(tf, function(x) {
        suppressMessages(return(plotActivityDim_(sce,
            activity_matrix, x, dimtype,
            label, legend.label, colors,
            limit, ...)))
    })

    if (combine == TRUE) {

        gs <- patchwork::wrap_plots(gs,
            ncol = ncol, nrow = nrow) +
            patchwork::plot_annotation(title = title)

        return(gs)

    } else {

        return(gs)

    }

}


#' @importFrom ggplot2 ggplot geom_violin theme_classic ggtitle ylab theme scale_fill_manual facet_grid geom_boxplot
plotActivityViolin_ <- function(activity_matrix, tf, clusters,
    legend.label, colors, text_size, facet_grid_variable,
    boxplot) {

    tf.activity <- as.numeric(activity_matrix[tf, ])

    df <- data.frame(activity = tf.activity, clusters = clusters)

    if (!is.null(facet_grid_variable)) {
        df$facet <- facet_grid_variable
    }


    g <- ggplot(df, aes(x = clusters,
        y = activity, fill = clusters)) + geom_violin() +
        theme_classic(base_size = 12) + ggtitle(tf) + ylab(legend.label) +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 1))

    if (!is.null(colors)) {
        g <- g + scale_fill_manual(values = colors)
    }

    if (!is.null(facet_grid_variable)) {
        g <- g + facet_grid(stats::reformulate("facet",
            "."), scales = "free", space = "free")
    }

    if (boxplot) {
        g <- g + geom_boxplot(width = 0.1)
    }
    g <- g + theme(text = element_text(size = text_size))


    return(g)

}

#' Generate violin plots of inferred activities for a list of TFs grouped by cluster/group labels
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param clusters A vector of cluster or group labels for single cells
#' @param ncol A integer to indicate the number of columns in the combined plot, if `combine = TRUE`
#' @param nrow A integer to indicate the number of rows in the combined plot, if `combine = TRUE`
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#' @param legend.label String indicating the name of variable to be plotted on the legend
#' @param colors  A character vector representing the names of colors
#' @param title String indicating the title of the plot if `combine = TRUE`
#' @param text_size Scalar indicating the font size of the title
#' @param facet_grid_variable  A character vector of a secondary label to split the plots by facet_grid
#' @param boxplot logical indicating whether to add boxplot on top of violin plot
#'
#' @return A combined ggplot object or a list of ggplots if `combine = FALSE`
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' plotActivityViolin(activity_matrix = logcounts(example_sce),
#' tf = c('Gene_0001','Gene_0002'),  clusters = example_sce$cluster)
#'
#' @author Xiaosai Yao, Shang-yang Chen
plotActivityViolin <- function(activity_matrix, tf, clusters,
    ncol = NULL, nrow = NULL, combine = TRUE, legend.label = "activity",
    colors = NULL, title = NULL, text_size = 10, facet_grid_variable = NULL,
    boxplot = FALSE) {

    # give warning for genes absent in tf list
    missing <- tf[which(!tf %in% rownames(activity_matrix))]
    if (!identical(missing, character(0))) {
        message(missing, " not found in activity matrix. Excluded from plots")
    }
    tf <- tf[which(tf %in% rownames(activity_matrix))]

    gs <- lapply(tf, function(x) {
        return(plotActivityViolin_(activity_matrix, x,
            clusters, legend.label, colors, text_size,
            facet_grid_variable, boxplot))
    })

    if (combine == TRUE) {

        gs <- patchwork::wrap_plots(gs, ncol = ncol, nrow = nrow) +
            patchwork::plot_annotation(title = title,
                theme = theme(plot.title = element_text(hjust = 0.5)))

        return(gs)

    } else {

        return(gs)

    }

}

#' Generate bubble plots of relative activities across cluster/group labels for a list of TFs
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param clusters A character or integer vector of cluster or group labels for single cells
#' @param bubblesize String indicating the variable from findDifferentialActivity output to scale size of bubbles
#' by either `FDR` or `summary.logFC`. Default is `FDR`.
#' @param color.theme String indicating the color theme used for the bubble plot and corresponding to the color options
#' in `scale_color_viridis_c`
#' @param legend.label String indicating the name of legend corresponding to the color scale
#' @param x.label String indicating the x axis label
#' @param y.label String indicating the y axis label
#' @param title String indicating the title of the plot
#' @param ... Additional arguments to pass to `findDifferentialActivity`
#'
#' @return A ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c scale_size_continuous theme_classic theme labs ylab xlab ggtitle
#'
#' @examples
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' plotBubble(activity_matrix = logcounts(example_sce),
#' tf = c('Gene_0001','Gene_0002'),  clusters = example_sce$cluster)
#' @author Shang-yang Chen
plotBubble <- function(activity_matrix, tf, clusters,
    bubblesize = c("FDR", "summary.logFC"), color.theme = "viridis",
    legend.label = "relative_activity", x.label = "clusters",
    y.label = "transcription factors", title = "TF activity",
    ...) {

    bubblesize <- match.arg(bubblesize)

    # give warning for genes absent in tf list
    missing <- tf[which(!tf %in% rownames(activity_matrix))]
    if (!identical(missing, character(0))) {
        message(missing, " not found in activity matrix. Excluded from plots")
    }
    tf <- tf[which(tf %in% rownames(activity_matrix))]

    # find logFC and FDR of TFs
    markers <- findDifferentialActivity(activity_matrix,
        clusters, ...)
    markers <- suppressMessages(getSigGenes(markers,
        fdr_cutoff = 1.5, logFC_cutoff = -100))
    markers <- markers[which(markers$tf %in% tf),
        ]
    levels <- make.names(unique(tf[tf %in% markers$tf]))
    markers$tf <- make.names(markers$tf)
    # rename markers
    colnames(markers)[colnames(markers) == "class"] <- "clusters"


    # z normalize activity and compute mean by cluster
    tf.activity <- activity_matrix[tf, , drop = FALSE]
    df <- data.frame(clusters = clusters, t(as.matrix(tf.activity)))
    df.mean <- stats::aggregate(. ~ clusters, df,
        mean)
    zscores <- apply(df.mean[, -1], 2, scale)
    df.mean <- data.frame(clusters = df.mean[, 1],
        as.data.frame(zscores))
    df.plot <- suppressMessages(reshape2::melt(df.mean,
        id.variable = "clusters", variable.name = "tf",
        value.name = "relative_activity"))

    # merge logFC, FDR and mean activity
    df.plot <- merge(df.plot, markers, by = c("tf",
        "clusters"))
    df.plot$tf <- factor(as.character(df.plot$tf),
        levels = levels)

    # generate bubble plots
    if (bubblesize == "FDR") {
        logpval <- -log10(df.plot$FDR)
        max.logpval <- max(logpval[is.finite(logpval)])
        logpval <- replace(logpval, is.infinite(logpval),
            max.logpval)
        g <- ggplot(df.plot, aes(clusters,
            tf, color = relative_activity)) +
            geom_point(stat = "identity", aes(size = logpval)) +
            scale_color_viridis_c(option = color.theme) +
            scale_size_continuous("-logpval", range = c(0,
                7)) + theme_classic(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45,
                hjust = 1)) + labs(color = legend.label) +
            ylab(y.label) + xlab(x.label) + ylab(y.label) +
            xlab(x.label) + ggtitle(title)
    } else if (bubblesize == "summary.logFC") {
        g <- ggplot(df.plot, aes(clusters,
            tf, color = relative_activity)) +
            geom_point(stat = "identity", aes(size = summary.logFC)) +
            scale_color_viridis_c(option = color.theme) +
            scale_size_continuous("summary logFC", range = c(0,
                7)) + theme_classic(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45,
                hjust = 1)) + labs(color = legend.label) +
            ylab(y.label) + xlab(x.label) + ggtitle(title)
    }
    return(g)
}

#' @importFrom ggplot2 ggplot aes scale_colour_gradient geom_point 
#' @importFrom ggplot2 coord_flip theme_bw ggtitle ylab theme unit element_blank

enrichPlot_ <- function(results, title, top) {
    results$logP.adj <- -log10(results$p.adjust)
    ggplot(results[seq_len(top), ], aes(y = logP.adj,
        x = Description, color = GeneRatio)) +
        scale_colour_gradient(high = "red", low = "blue") +
        geom_point(stat = "identity", aes(size = Odds.Ratio)) +
        coord_flip() + theme_bw() + ggtitle(title) +
        ylab(expression("-log"[10] ~ "pval")) +
        theme(text = element_text(size = 10),
            axis.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45,
                hjust = 1), plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.5, "lines"))
}

#' Plot results of regulonEnrich
#'
#' @param results  Output from regulonEnrich
#' @param top An integer to indicate the number of pathways to plot ranked by significance. Default is 15.
#' @param ncol An integer to indicate the number of columns in the combined plot, if combine == TRUE. Default is 3.
#' @param title String indicating the title of the combined plot
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @export

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
#' # plot graph
#' enrichPlot(results = enrichment_results )
#'
#' @author Xiaosai Yao
enrichPlot <- function(results, top = 15, ncol = 3, title = NULL,
    combine = TRUE) {

    gs <- lapply(names(results), function(x) {
        enrichPlot_(results[[x]], x, top)
    })

    if (combine == TRUE) {

        gs <- patchwork::wrap_plots(gs, ncol = ncol) +
            patchwork::plot_annotation(title = title,
                theme = theme(plot.title = element_text(hjust = 0.5)))

        return(gs)

    } else {

        return(gs)

    }

}

#' Plot targets genes of transcription factors in regulons
#'
#' @param sce A SingleCellExperiment object containing information of cell attributes
#' @param tfs A character vector indicating the names of the transcription factors to be plotted
#' @param regulon A dataframe of regulons containing `tf`, `targets` and a column for filtering the regulons
#' @param regulon_column String indicating the column names to be used for filtering regulons
#' @param regulon_cutoff A scalar indicating the minimal value to retain the regulons for plotting
#' @param downsample Integer indicating the number of cells to sample from the matrix
#' @param scale Logical indicating whether to scale the heatmap
#' @param center Logical indicating whether to center the heatmap
#' @param color_breaks A vector indicating numeric breaks as input to `circlize::colorRamp2`
#' @param colors A vector of colors corresponding to values in `breaks` as input to `circlize::colorRamp2`
#' @param cell_attributes A character vector matching the column names of `colData(sce)` to be used for plotting
#' @param col_gap String indicating the cell attribute to split the columns of the heatmap by
#' @param exprs_values A string specifying which assay in `assays(object)` to obtain expression values from
#' @param use_raster Logical indicating whether to use rasterization to reduce image size
#' @param raster_quality Integer indicating the raster quality. The higher the value, the better the resolution
#' @param cluster_rows Logical indicating whether to cluster rows
#' @param cluster_columns Logical indicating whether to cluster columns
#' @param border Logical indicating whether to add border around heatmap
#' @param show_column_names Logical indicating whether to show column names
#' @param column_col A list specifying the colors in the columns.
#' See [here](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#simple-annotation)
#' @param row_col A list specifying the colors in the rows.
#' See [here](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#simple-annotation)
#' @param ... other arguments for `ComplexHeatmap::Heatmap`
#' @return A Heatmap-class object.
#' @export
#' @examples
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' regulon <- data.frame(tf=c(rep('Gene_0001',10),rep('Gene_0002',20)),
#' target = sample(rownames(example_sce),30), weight = rnorm(30))
#' #plot heatmap and rotate labels
#' plotHeatmapRegulon(example_sce, tfs=c('Gene_0001','Gene_0002'), regulon=regulon,
#' cell_attributes='cluster', col_gap = 'cluster', column_title_rot = 90)
#' @author Xiaosai Yao

plotHeatmapRegulon <- function(sce,
    tfs, regulon, regulon_column = "weight",
    regulon_cutoff = 0.1,
    downsample = 1000,
    scale = TRUE,
    center = TRUE,
    color_breaks = c(-2,
        0, 2),
    colors = c("blue","white", "red"),
    cell_attributes, col_gap = NULL,
    exprs_values = "logcounts",
    use_raster = TRUE,
    raster_quality = 10,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = TRUE, show_column_names = FALSE,
    column_col = NULL,
    row_col = NULL, ...) {

    downsample_seq <- seq(from = 1,
        to = ncol(sce),
        by = floor(max(1,
            ncol(sce)/downsample)))

    # keep only targets belonging to TFs and meeting cutoff
    if (is.matrix(regulon[[regulon_column]])) {

        regulon <- regulon[which(regulon$tf %in% tfs &
                             apply(regulon[[regulon_column]], 1,
                                   function(x) any(x > regulon_cutoff))), ]

    } else {
        regulon <- regulon[regulon$tf %in% tfs &
                             regulon[,regulon_column] > regulon_cutoff, ,drop=FALSE]
    }

    regulon.split <- S4Vectors::split(regulon, f <- regulon$tf)

    # remove duplicated genes from each tf
    regulon.split <- lapply(regulon.split, function(x) x[!duplicated(x$target), ])

    regulon <- do.call(rbind, as.list(regulon.split))

    # remove targets not found in sce
    regulon <- regulon[regulon$target %in% rownames(sce),]
    targets <- regulon$target

    sce <- sce[targets, downsample_seq]


    right_annotation <- data.frame(tf = regulon$tf)
    top_annotation <- data.frame(colData(sce)[cell_attributes])

    if (!is.null(col_gap)) {
        column_split <- top_annotation[col_gap]
    } else {
        column_split <- NULL
    }



    mat <- as.matrix(SummarizedExperiment::assay(sce,
        exprs_values))
    mat <- t(scale(t(mat),
        scale = scale,
        center = center))

    col_fun <- circlize::colorRamp2(color_breaks,
        colors)

    ComplexHeatmap::Heatmap(mat,
        col = col_fun,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(df = top_annotation,
            col = column_col),
        right_annotation = ComplexHeatmap::rowAnnotation(df = right_annotation,
            col = row_col),
        row_split = right_annotation,
        column_split = column_split,
        use_raster = use_raster,
        raster_quality = raster_quality,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        border = border,
        show_column_names = show_column_names,
        ...)
}

#' Plot transcription factor activity
#'
#' @param activity_matrix A matrix of values, such as TF activities inferred from calculateActivity
#' @param sce A SingleCellExperiment object containing information of cell attributes
#' @param tfs A character vector indicating the names of the transcription factors to be plotted
#' @param downsample Integer indicating the number of cells to sample from the matrix
#' @param scale Logical indicating whether to scale the heatmap
#' @param center Logical indicating whether to center the heatmap
#' @param color_breaks A vector indicating numeric breaks as input to `circlize::colorRamp2`
#' @param colors A vector of colors corresponding to values in `breaks` as input to `circlize::colorRamp2`
#' @param cell_attributes A character vector matching the column names of `colData(sce)` to be used for plotting
#' @param col_gap String indicating the cell attribute to split the columns of the heatmap by
#' @param use_raster Logical indicating whether to use rasterization to reduce image size
#' @param raster_quality Integer indicating the raster quality. The higher the value, the better the resolution
#' @param cluster_rows Logical indicating whether to cluster rows
#' @param cluster_columns Logical indicating whether to cluster columns
#' @param border Logical indicating whether to add border around heatmap
#' @param show_column_names Logical indicating whether to show column names
#' @param ... other arguments for `ComplexHeatmap::Heatmap`
#' @return A Heatmap-class object.
#' @export
#' @examples
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' activity_matrix <- matrix(rnorm(10*200), nrow=10, ncol=200)
#' rownames(activity_matrix) <- sample(rownames(example_sce),10)
#' plotHeatmapActivity(activity_matrix=activity_matrix, sce=example_sce,
#' tfs=rownames(activity_matrix), cell_attributes='cluster', col_gap='cluster')
#' @author Xiaosai Yao
plotHeatmapActivity <- function(activity_matrix, sce, tfs, downsample = 1000,
    scale = TRUE, center = TRUE, color_breaks = c(-2, 0, 2),
    colors = c("blue", "white", "red"), cell_attributes = NULL,
    col_gap = NULL, use_raster = TRUE, raster_quality = 10,
    cluster_rows = TRUE, cluster_columns = FALSE, border = TRUE,
    show_column_names = FALSE, ...) {

    tfs <- tfs[tfs %in% rownames(activity_matrix)]
    downsample_seq <- seq(from = 1, to = ncol(sce), by = floor(max(1,
        ncol(sce)/downsample)))

    sce <- sce[, downsample_seq]
    top_annotation <- data.frame(colData(sce)[cell_attributes])

    if (!is.null(col_gap)) {
        column_split <- top_annotation[col_gap]
    } else {
        column_split <- NULL
    }

    activity_matrix <- activity_matrix[tfs, downsample_seq]
    activity_matrix <- Matrix::t(scale(Matrix::t(activity_matrix),
        scale = scale, center = center))

    col_fun <- circlize::colorRamp2(color_breaks, colors)

    ComplexHeatmap::Heatmap(activity_matrix, col = col_fun,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(df = top_annotation),
        column_split = top_annotation[col_gap], use_raster = use_raster,
        raster_quality = raster_quality, cluster_rows = cluster_rows,
        cluster_columns = cluster_columns, border = border,
        show_column_names = show_column_names, ...)

}

#' @importFrom ggplot2 ggplot geom_boxplot geom_point aes ggtitle
plotDiagnostic <- function(idx, regulon, expMatrix, exp_assay,
    exp_cutoff = 1, peakMatrix, peak_assay, peak_cutoff = 0,
    clusters) {

    target <- regulon$target[idx]
    tf <- regulon$tf[idx]
    peak <- regulon$idxATAC[idx]

    tg_plot <- list()
    for (cluster in unique(clusters)) {
        target_exp <- SummarizedExperiment::assay(expMatrix,
            exp_assay)[target, clusters == cluster]
        tf_exp <- SummarizedExperiment::assay(expMatrix, exp_assay)[tf,
            clusters == cluster]
        peak_accessibility <- SummarizedExperiment::assay(peakMatrix,
            peak_assay)[peak, clusters == cluster]

        if (is.null(exp_cutoff)) {
            exp_cutoff <- mean(tf_exp)
        }

        if (is.null(peak_cutoff)) {
            peak_cutoff <- mean(peak_accessibility)
        }
        tf_exp.bi <- binarize(tf_exp, cutoff = exp_cutoff)
        peak.bi <- binarize(peak_accessibility, cutoff = peak_cutoff)
        tf_re.bi <- tf_exp.bi * peak.bi
        tf_re.bi <- factor(tf_re.bi, levels = c(0, 1))

        target_group <- data.frame(groups = as.vector(tf_re.bi),
            target = as.vector(target_exp))
        tg_plot[[cluster]] <- ggplot(target_group, aes(groups,
            target)) + geom_boxplot() + geom_point(position = "jitter") +
            ggtitle(paste0("target:", target, " tf:", tf, " peak:",
                peak, " in ", cluster, "\n", " GRN corr:",
                round(regulon$corr[idx, cluster], 2), " GRN pval:",
                round(regulon$pval[idx, cluster], 2), " weight:",
                round(regulon$weight[idx, cluster], 2)))

    }

    tg_plot


}

binarize <- function(input_vector, cutoff) {
    filtered <- as.numeric(input_vector > cutoff)
}

normalizeCols <- function(mat = NULL, scaleTo = NULL) {
    colSm <- Matrix::colSums(mat)
    if (!is.null(scaleTo)) {
        mat@x <- scaleTo * mat@x/rep.int(colSm, Matrix::diff(mat@p))
    } else {
        mat@x <- mat@x/rep.int(colSm, Matrix::diff(mat@p))
    }
    return(mat)
}

