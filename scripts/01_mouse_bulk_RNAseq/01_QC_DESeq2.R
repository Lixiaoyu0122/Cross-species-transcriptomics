#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

options(stringsAsFactors = FALSE)

run <- "/home/li/Wangxinming/FSTL1_mouse_bulkRNA/rerun_20260720_n9"

count_file <- file.path(
    run, "counts", "final", "gene_counts_raw.tsv"
)

metadata_file <- file.path(
    run, "counts", "final", "sample_metadata_deseq2.tsv"
)

annotation_file <- file.path(
    run, "counts", "final", "gene_annotation.tsv"
)

outdir <- file.path(
    run, "results", "deseq2_qc"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
    count_file,
    metadata_file,
    annotation_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop(
        "Missing input files:\n",
        paste(missing_files, collapse = "\n")
    )
}

# ============================================================
# 1. Read input data
# ============================================================

count_df <- read.delim(
    count_file,
    row.names = 1,
    check.names = FALSE
)

metadata <- read.delim(
    metadata_file,
    check.names = FALSE
)

annotation <- read.delim(
    annotation_file,
    check.names = FALSE
)

# Convert the count data frame into a numeric matrix.
# DESeq2 requires an integer matrix, not a data.frame/list.
non_numeric_columns <- colnames(count_df)[
    !vapply(count_df, is.numeric, logical(1))
]

if (length(non_numeric_columns) > 0) {
    stop(
        "Non-numeric columns detected in count matrix: ",
        paste(non_numeric_columns, collapse = ", ")
    )
}

count_mat <- as.matrix(count_df)
storage.mode(count_mat) <- "integer"

if (anyDuplicated(rownames(count_mat))) {
    stop("Duplicated gene IDs detected in count matrix.")
}

required_metadata_columns <- c(
    "sample",
    "group",
    "display_group"
)

if (!all(required_metadata_columns %in% colnames(metadata))) {
    stop(
        "Metadata must contain: ",
        paste(required_metadata_columns, collapse = ", ")
    )
}

if (anyDuplicated(metadata$sample)) {
    stop("Duplicated sample names detected in metadata.")
}

rownames(metadata) <- metadata$sample

if (!setequal(colnames(count_mat), metadata$sample)) {
    stop(
        "Count matrix and metadata sample names do not match.\n",
        "Count samples: ",
        paste(colnames(count_mat), collapse = ", "),
        "\nMetadata samples: ",
        paste(metadata$sample, collapse = ", ")
    )
}

metadata <- metadata[colnames(count_mat), , drop = FALSE]

metadata$group <- factor(
    metadata$group,
    levels = c("Sham", "MI", "MTans")
)

if (any(is.na(metadata$group))) {
    stop("Unexpected group names detected.")
}

group_counts <- table(metadata$group)

if (!all(group_counts == c(3, 3, 3))) {
    stop(
        "Incorrect experimental design: ",
        paste(
            names(group_counts),
            group_counts,
            collapse = "; "
        )
    )
}

if (any(is.na(count_mat))) {
    stop("NA values detected in count matrix.")
}

if (any(count_mat < 0)) {
    stop("Negative values detected in count matrix.")
}

# count_mat has already been converted to integer mode above

cat("Samples:", ncol(count_mat), "\n")
cat("Genes before filtering:", nrow(count_mat), "\n")
cat("Groups:\n")
print(group_counts)

# ============================================================
# 2. Construct DESeq2 object and filter genes
# ============================================================

dds_all <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = metadata,
    design = ~ group
)

# Retain genes with at least 10 counts in at least 3 samples
keep <- rowSums(DESeq2::counts(dds_all) >= 10) >= 3

dds <- dds_all[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")

if (nrow(dds) < 5000) {
    stop(
        "Too few genes retained after filtering: ",
        nrow(dds)
    )
}

filter_summary <- data.frame(
    metric = c(
        "genes_before_filter",
        "genes_after_filter",
        "genes_removed",
        "filter_rule"
    ),
    value = c(
        nrow(dds_all),
        nrow(dds),
        nrow(dds_all) - nrow(dds),
        "count >= 10 in at least 3 samples"
    )
)

write.table(
    filter_summary,
    file.path(outdir, "prefilter_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 3. DESeq2 size-factor normalization
# ============================================================

dds <- estimateSizeFactors(dds)

size_factor_table <- data.frame(
    sample = colnames(dds),
    group = as.character(colData(dds)$group),
    raw_library_size = colSums(DESeq2::counts(dds)),
    size_factor = sizeFactors(dds),
    size_factor_relative_to_median =
        sizeFactors(dds) / median(sizeFactors(dds))
)

write.table(
    size_factor_table,
    file.path(outdir, "deseq2_size_factors.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

normalized_counts <- DESeq2::counts(
    dds,
    normalized = TRUE
)

write.table(
    data.frame(
        gene_id = rownames(normalized_counts),
        normalized_counts,
        check.names = FALSE
    ),
    file.path(outdir, "normalized_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 4. VST transformation and PCA
# ============================================================

vsd <- vst(
    dds,
    blind = TRUE
)

vst_matrix <- assay(vsd)

write.table(
    data.frame(
        gene_id = rownames(vst_matrix),
        vst_matrix,
        check.names = FALSE
    ),
    file.path(outdir, "vst_expression_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pca <- prcomp(
    t(vst_matrix),
    center = TRUE,
    scale. = FALSE
)

variance_pct <- (
    pca$sdev^2 /
    sum(pca$sdev^2)
) * 100

pca_data <- data.frame(
    sample = rownames(pca$x),
    group = as.character(
        metadata[rownames(pca$x), "group"]
    ),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3]
)

write.table(
    pca_data,
    file.path(outdir, "pca_coordinates.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(
        principal_component = paste0(
            "PC",
            seq_along(variance_pct)
        ),
        variance_explained_pct = variance_pct
    ),
    file.path(outdir, "pca_variance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pca_plot <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        shape = group,
        label = sample
    )
) +
    geom_point(size = 4) +
    geom_text(
        vjust = -0.8,
        check_overlap = TRUE
    ) +
    theme_bw(base_size = 12) +
    labs(
        x = sprintf(
            "PC1 (%.1f%%)",
            variance_pct[1]
        ),
        y = sprintf(
            "PC2 (%.1f%%)",
            variance_pct[2]
        ),
        title = "PCA of VST-normalized counts",
        shape = "Group"
    )

ggsave(
    file.path(outdir, "PCA_PC1_PC2.pdf"),
    pca_plot,
    width = 7,
    height = 5
)

ggsave(
    file.path(outdir, "PCA_PC1_PC2.png"),
    pca_plot,
    width = 7,
    height = 5,
    dpi = 300
)

# ============================================================
# 5. Sample correlation and distance
# ============================================================

correlation_matrix <- cor(
    vst_matrix,
    method = "pearson"
)

distance_matrix <- as.matrix(
    dist(t(vst_matrix))
)

write.table(
    correlation_matrix,
    file.path(outdir, "sample_correlation.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

write.table(
    distance_matrix,
    file.path(outdir, "sample_distance.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

heatmap_annotation <- data.frame(
    Group = metadata$group
)

rownames(heatmap_annotation) <- metadata$sample

pheatmap(
    correlation_matrix,
    annotation_col = heatmap_annotation,
    annotation_row = heatmap_annotation,
    border_color = NA,
    filename = file.path(
        outdir,
        "sample_correlation_heatmap.pdf"
    ),
    width = 7,
    height = 6
)

pheatmap(
    distance_matrix,
    annotation_col = heatmap_annotation,
    annotation_row = heatmap_annotation,
    border_color = NA,
    filename = file.path(
        outdir,
        "sample_distance_heatmap.pdf"
    ),
    width = 7,
    height = 6
)

correlation_summary <- do.call(
    rbind,
    lapply(
        colnames(correlation_matrix),
        function(sample) {
            values <- correlation_matrix[, sample]
            values <- values[names(values) != sample]

            nearest_sample <- names(
                which.max(values)
            )

            data.frame(
                sample = sample,
                group = as.character(
                    metadata[sample, "group"]
                ),
                nearest_sample = nearest_sample,
                nearest_group = as.character(
                    metadata[nearest_sample, "group"]
                ),
                nearest_correlation = max(values),
                median_correlation = median(values),
                minimum_correlation = min(values)
            )
        }
    )
)

write.table(
    correlation_summary,
    file.path(
        outdir,
        "sample_correlation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. Fstl1 normalized expression
# ============================================================

fstl1_annotation <- annotation[
    tolower(annotation$gene_name) == "fstl1",
    ,
    drop = FALSE
]

fstl1_ids <- intersect(
    fstl1_annotation$gene_id,
    rownames(normalized_counts)
)

if (length(fstl1_ids) != 1) {
    stop(
        "Expected exactly one retained Fstl1 gene; found: ",
        paste(fstl1_ids, collapse = ", ")
    )
}

fstl1_id <- fstl1_ids[1]

fstl1_data <- data.frame(
    gene_id = fstl1_id,
    gene_name = "Fstl1",
    sample = colnames(normalized_counts),
    group = as.character(metadata$group),
    normalized_count = as.numeric(
        normalized_counts[
            fstl1_id,
            colnames(normalized_counts)
        ]
    )
)

write.table(
    fstl1_data,
    file.path(
        outdir,
        "Fstl1_normalized_counts.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

fstl1_group_summary <- do.call(
    rbind,
    lapply(
        split(
            fstl1_data$normalized_count,
            fstl1_data$group
        ),
        function(x) {
            data.frame(
                n = length(x),
                mean = mean(x),
                median = median(x),
                sd = sd(x),
                min = min(x),
                max = max(x)
            )
        }
    )
)

fstl1_group_summary$group <- rownames(
    fstl1_group_summary
)

rownames(fstl1_group_summary) <- NULL

fstl1_group_summary <- fstl1_group_summary[
    ,
    c(
        "group",
        "n",
        "mean",
        "median",
        "sd",
        "min",
        "max"
    )
]

write.table(
    fstl1_group_summary,
    file.path(
        outdir,
        "Fstl1_group_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

fstl1_plot <- ggplot(
    fstl1_data,
    aes(
        x = group,
        y = log2(normalized_count + 1)
    )
) +
    geom_boxplot(
        outlier.shape = NA,
        width = 0.55
    ) +
    geom_jitter(
        width = 0.08,
        size = 3
    ) +
    theme_bw(base_size = 12) +
    labs(
        x = NULL,
        y = "log2(normalized count + 1)",
        title = "Fstl1 normalized expression"
    )

ggsave(
    file.path(
        outdir,
        "Fstl1_normalized_expression.pdf"
    ),
    fstl1_plot,
    width = 5,
    height = 4
)

ggsave(
    file.path(
        outdir,
        "Fstl1_normalized_expression.png"
    ),
    fstl1_plot,
    width = 5,
    height = 4,
    dpi = 300
)

# ============================================================
# 7. Save objects and software information
# ============================================================

saveRDS(
    dds,
    file.path(
        outdir,
        "dds_sizefactor_estimated.rds"
    )
)

saveRDS(
    vsd,
    file.path(
        outdir,
        "vsd_blind_true.rds"
    )
)

writeLines(
    capture.output(sessionInfo()),
    file.path(
        outdir,
        "R_sessionInfo.txt"
    )
)

cat("============================================================\n")
cat("DESeq2 SAMPLE QC PASSED\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes before filtering:", nrow(dds_all), "\n")
cat("Genes after filtering:", nrow(dds), "\n")
cat("Output directory:", outdir, "\n")
cat("============================================================\n")
