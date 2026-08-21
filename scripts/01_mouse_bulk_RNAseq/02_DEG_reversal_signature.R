#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

options(stringsAsFactors = FALSE)

run <- "/home/li/Wangxinming/FSTL1_mouse_bulkRNA/rerun_20260720_n9"

dds_file <- file.path(
    run, "results", "deseq2_qc",
    "dds_sizefactor_estimated.rds"
)

annotation_file <- file.path(
    run, "counts", "final",
    "gene_annotation.tsv"
)

outdir <- file.path(
    run, "results", "deseq2_deg"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(dds_file)) {
    stop("Missing DESeq2 object: ", dds_file)
}

if (!file.exists(annotation_file)) {
    stop("Missing annotation: ", annotation_file)
}

# ============================================================
# 1. Load and validate data
# ============================================================

dds <- readRDS(dds_file)

annotation <- read.delim(
    annotation_file,
    check.names = FALSE
)

required_annotation <- c(
    "gene_id",
    "gene_name",
    "gene_biotype"
)

if (!all(required_annotation %in% colnames(annotation))) {
    stop("Gene annotation is incomplete.")
}

if (anyDuplicated(annotation$gene_id)) {
    stop("Duplicated gene IDs found in annotation.")
}

if (!"group" %in% colnames(colData(dds))) {
    stop("group column missing from DESeq2 colData.")
}

colData(dds)$group <- factor(
    as.character(colData(dds)$group),
    levels = c("Sham", "MI", "MTans")
)

if (any(is.na(colData(dds)$group))) {
    stop("Unexpected group labels.")
}

group_table <- table(colData(dds)$group)

if (!all(group_table == c(3, 3, 3))) {
    stop(
        "Expected Sham=3, MI=3, MTans=3; observed: ",
        paste(names(group_table), group_table, collapse = "; ")
    )
}

design(dds) <- ~ group

cat("Samples:", ncol(dds), "\n")
cat("Genes entering DESeq2:", nrow(dds), "\n")
cat("Group design:\n")
print(group_table)

# ============================================================
# 2. Fit DESeq2 model
# ============================================================

set.seed(20260720)

dds <- DESeq(
    dds,
    test = "Wald",
    quiet = FALSE
)

saveRDS(
    dds,
    file.path(outdir, "dds_DESeq_fitted.rds")
)

writeLines(
    resultsNames(dds),
    file.path(outdir, "DESeq2_coefficient_names.txt")
)

# ============================================================
# 3. Extract contrasts
# ============================================================

alpha_value <- 0.05

res_mi_sham <- results(
    dds,
    contrast = c("group", "MI", "Sham"),
    alpha = alpha_value
)

res_mtans_mi <- results(
    dds,
    contrast = c("group", "MTans", "MI"),
    alpha = alpha_value
)

res_mtans_sham <- results(
    dds,
    contrast = c("group", "MTans", "Sham"),
    alpha = alpha_value
)

make_result_table <- function(res, comparison) {

    df <- as.data.frame(res)

    df$gene_id <- rownames(df)

    ann <- annotation[
        match(df$gene_id, annotation$gene_id),
        required_annotation,
        drop = FALSE
    ]

    output <- data.frame(
        gene_id = df$gene_id,
        gene_name = ann$gene_name,
        gene_biotype = ann$gene_biotype,
        baseMean = df$baseMean,
        log2FoldChange = df$log2FoldChange,
        lfcSE = df$lfcSE,
        stat = df$stat,
        pvalue = df$pvalue,
        padj = df$padj,
        comparison = comparison,
        stringsAsFactors = FALSE
    )

    output$significance_FDR05 <- ifelse(
        !is.na(output$padj) &
        output$padj < 0.05 &
        output$log2FoldChange > 0,
        "Up",
        ifelse(
            !is.na(output$padj) &
            output$padj < 0.05 &
            output$log2FoldChange < 0,
            "Down",
            "NS"
        )
    )

    output$significance_FDR05_LFC05 <- ifelse(
        !is.na(output$padj) &
        output$padj < 0.05 &
        output$log2FoldChange >= 0.5,
        "Up",
        ifelse(
            !is.na(output$padj) &
            output$padj < 0.05 &
            output$log2FoldChange <= -0.5,
            "Down",
            "NS"
        )
    )

    output <- output[
        order(
            is.na(output$padj),
            output$padj,
            output$pvalue
        ),
        ,
        drop = FALSE
    ]

    rownames(output) <- NULL

    output
}

tab_mi_sham <- make_result_table(
    res_mi_sham,
    "MI_vs_Sham"
)

tab_mtans_mi <- make_result_table(
    res_mtans_mi,
    "MTans_vs_MI"
)

tab_mtans_sham <- make_result_table(
    res_mtans_sham,
    "MTans_vs_Sham"
)

write.table(
    tab_mi_sham,
    file.path(outdir, "DEG_MI_vs_Sham_all.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    tab_mtans_mi,
    file.path(outdir, "DEG_MTans_vs_MI_all.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    tab_mtans_sham,
    file.path(outdir, "DEG_MTans_vs_Sham_all.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 4. Write significant DEG subsets
# ============================================================

write_deg_subsets <- function(tab, prefix) {

    fdr05 <- tab[
        !is.na(tab$padj) &
        tab$padj < 0.05,
        ,
        drop = FALSE
    ]

    fdr05_lfc05 <- tab[
        !is.na(tab$padj) &
        tab$padj < 0.05 &
        abs(tab$log2FoldChange) >= 0.5,
        ,
        drop = FALSE
    ]

    fdr05_lfc1 <- tab[
        !is.na(tab$padj) &
        tab$padj < 0.05 &
        abs(tab$log2FoldChange) >= 1,
        ,
        drop = FALSE
    ]

    write.table(
        fdr05,
        file.path(outdir, paste0(prefix, "_FDR05.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        fdr05_lfc05,
        file.path(outdir, paste0(prefix, "_FDR05_absLFC05.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        fdr05_lfc1,
        file.path(outdir, paste0(prefix, "_FDR05_absLFC1.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

write_deg_subsets(tab_mi_sham, "DEG_MI_vs_Sham")
write_deg_subsets(tab_mtans_mi, "DEG_MTans_vs_MI")
write_deg_subsets(tab_mtans_sham, "DEG_MTans_vs_Sham")

# ============================================================
# 5. DEG count summary
# ============================================================

count_deg <- function(tab, padj_cutoff, lfc_cutoff) {

    valid <- (
        !is.na(tab$padj) &
        tab$padj < padj_cutoff &
        abs(tab$log2FoldChange) >= lfc_cutoff
    )

    up <- sum(
        valid &
        tab$log2FoldChange > 0
    )

    down <- sum(
        valid &
        tab$log2FoldChange < 0
    )

    c(
        up = up,
        down = down,
        total = up + down
    )
}

comparison_tables <- list(
    MI_vs_Sham = tab_mi_sham,
    MTans_vs_MI = tab_mtans_mi,
    MTans_vs_Sham = tab_mtans_sham
)

summary_rows <- list()

for (comparison in names(comparison_tables)) {

    tab <- comparison_tables[[comparison]]

    for (padj_cutoff in c(0.05, 0.10)) {

        for (lfc_cutoff in c(0, 0.5, 1)) {

            counts_result <- count_deg(
                tab,
                padj_cutoff,
                lfc_cutoff
            )

            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                comparison = comparison,
                padj_cutoff = padj_cutoff,
                abs_log2FC_cutoff = lfc_cutoff,
                up = counts_result["up"],
                down = counts_result["down"],
                total = counts_result["total"]
            )
        }
    }
}

deg_summary <- do.call(rbind, summary_rows)

rownames(deg_summary) <- NULL

write.table(
    deg_summary,
    file.path(outdir, "DEG_count_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. Construct M-Tans reversal table
# ============================================================

select_for_merge <- function(tab, suffix) {

    output <- tab[
        ,
        c(
            "gene_id",
            "gene_name",
            "gene_biotype",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj"
        ),
        drop = FALSE
    ]

    colnames(output)[4:7] <- paste0(
        colnames(output)[4:7],
        "_",
        suffix
    )

    output
}

reversal <- merge(
    select_for_merge(tab_mi_sham, "MI_vs_Sham"),
    select_for_merge(tab_mtans_mi, "MTans_vs_MI"),
    by = c(
        "gene_id",
        "gene_name",
        "gene_biotype"
    ),
    all = FALSE
)

reversal <- merge(
    reversal,
    select_for_merge(tab_mtans_sham, "MTans_vs_Sham"),
    by = c(
        "gene_id",
        "gene_name",
        "gene_biotype"
    ),
    all = FALSE
)

mi_lfc <- reversal$log2FoldChange_MI_vs_Sham
drug_lfc <- reversal$log2FoldChange_MTans_vs_MI

mi_padj <- reversal$padj_MI_vs_Sham
drug_padj <- reversal$padj_MTans_vs_MI

opposite_direction <- (
    (mi_lfc > 0 & drug_lfc < 0) |
    (mi_lfc < 0 & drug_lfc > 0)
)

reversal$reversal_direction <- ifelse(
    mi_lfc > 0 & drug_lfc < 0,
    "MI_up_MTans_down",
    ifelse(
        mi_lfc < 0 & drug_lfc > 0,
        "MI_down_MTans_up",
        "not_reversed"
    )
)

reversal$strict_FDR05 <- (
    !is.na(mi_padj) &
    !is.na(drug_padj) &
    mi_padj < 0.05 &
    drug_padj < 0.05 &
    opposite_direction
)

reversal$strict_FDR05_absLFC05 <- (
    reversal$strict_FDR05 &
    abs(mi_lfc) >= 0.5 &
    abs(drug_lfc) >= 0.5
)

reversal$strict_FDR05_absLFC1 <- (
    reversal$strict_FDR05 &
    abs(mi_lfc) >= 1 &
    abs(drug_lfc) >= 1
)

reversal$relaxed_FDR10_absLFC05 <- (
    !is.na(mi_padj) &
    !is.na(drug_padj) &
    mi_padj < 0.10 &
    drug_padj < 0.10 &
    opposite_direction &
    abs(mi_lfc) >= 0.5 &
    abs(drug_lfc) >= 0.5
)

reversal$reversal_score <- ifelse(
    opposite_direction,
    abs(mi_lfc) + abs(drug_lfc),
    NA_real_
)

reversal$restoration_fraction <- ifelse(
    opposite_direction &
    abs(mi_lfc) > 1e-8,
    -drug_lfc / mi_lfc,
    NA_real_
)

reversal$MTans_vs_Sham_FDR05 <- (
    !is.na(reversal$padj_MTans_vs_Sham) &
    reversal$padj_MTans_vs_Sham < 0.05
)

reversal <- reversal[
    order(
        !reversal$strict_FDR05_absLFC05,
        -reversal$reversal_score,
        reversal$padj_MI_vs_Sham,
        reversal$padj_MTans_vs_MI
    ),
    ,
    drop = FALSE
]

rownames(reversal) <- NULL

write.table(
    reversal,
    file.path(outdir, "MTans_reversal_all_genes.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

strict_reversal <- reversal[
    reversal$strict_FDR05,
    ,
    drop = FALSE
]

strict_effect_reversal <- reversal[
    reversal$strict_FDR05_absLFC05,
    ,
    drop = FALSE
]

strict_lfc1_reversal <- reversal[
    reversal$strict_FDR05_absLFC1,
    ,
    drop = FALSE
]

relaxed_reversal <- reversal[
    reversal$relaxed_FDR10_absLFC05,
    ,
    drop = FALSE
]

write.table(
    strict_reversal,
    file.path(outdir, "MTans_reversal_strict_FDR05.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    strict_effect_reversal,
    file.path(
        outdir,
        "MTans_reversal_strict_FDR05_absLFC05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    strict_lfc1_reversal,
    file.path(
        outdir,
        "MTans_reversal_strict_FDR05_absLFC1.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    relaxed_reversal,
    file.path(
        outdir,
        "MTans_reversal_relaxed_FDR10_absLFC05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

reversal_summary <- data.frame(
    signature = c(
        "strict_FDR05_direction",
        "strict_FDR05_absLFC05",
        "strict_FDR05_absLFC1",
        "relaxed_FDR10_absLFC05"
    ),
    total = c(
        nrow(strict_reversal),
        nrow(strict_effect_reversal),
        nrow(strict_lfc1_reversal),
        nrow(relaxed_reversal)
    ),
    MI_up_MTans_down = c(
        sum(
            strict_reversal$reversal_direction ==
            "MI_up_MTans_down"
        ),
        sum(
            strict_effect_reversal$reversal_direction ==
            "MI_up_MTans_down"
        ),
        sum(
            strict_lfc1_reversal$reversal_direction ==
            "MI_up_MTans_down"
        ),
        sum(
            relaxed_reversal$reversal_direction ==
            "MI_up_MTans_down"
        )
    ),
    MI_down_MTans_up = c(
        sum(
            strict_reversal$reversal_direction ==
            "MI_down_MTans_up"
        ),
        sum(
            strict_effect_reversal$reversal_direction ==
            "MI_down_MTans_up"
        ),
        sum(
            strict_lfc1_reversal$reversal_direction ==
            "MI_down_MTans_up"
        ),
        sum(
            relaxed_reversal$reversal_direction ==
            "MI_down_MTans_up"
        )
    )
)

write.table(
    reversal_summary,
    file.path(outdir, "MTans_reversal_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. Extract Fstl1 results
# ============================================================

fstl1_id <- annotation$gene_id[
    tolower(annotation$gene_name) == "fstl1"
]

fstl1_id <- intersect(
    fstl1_id,
    rownames(dds)
)

if (length(fstl1_id) != 1) {
    stop(
        "Expected one Fstl1 gene; found: ",
        paste(fstl1_id, collapse = ", ")
    )
}

fstl1_id <- fstl1_id[1]

fstl1_result <- reversal[
    reversal$gene_id == fstl1_id,
    ,
    drop = FALSE
]

write.table(
    fstl1_result,
    file.path(outdir, "Fstl1_DESeq2_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. Volcano plots
# ============================================================

make_volcano <- function(tab, title_text, filename_prefix) {

    plot_df <- tab

    plot_df$status <- plot_df$significance_FDR05_LFC05

    plot_df$neg_log10_padj <- -log10(
        pmax(plot_df$padj, 1e-300)
    )

    plot_df$neg_log10_padj[
        is.na(plot_df$neg_log10_padj)
    ] <- 0

    valid <- plot_df[
        !is.na(plot_df$padj),
        ,
        drop = FALSE
    ]

    label_df <- head(
        valid[
            order(valid$padj),
            ,
            drop = FALSE
        ],
        10
    )

    label_df$label <- ifelse(
        is.na(label_df$gene_name) |
        label_df$gene_name == "",
        label_df$gene_id,
        label_df$gene_name
    )

    p <- ggplot(
        plot_df,
        aes(
            x = log2FoldChange,
            y = neg_log10_padj,
            shape = status
        )
    ) +
        geom_point(
            alpha = 0.55,
            size = 1.2
        ) +
        geom_vline(
            xintercept = c(-0.5, 0.5),
            linetype = "dashed"
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed"
        ) +
        geom_text(
            data = label_df,
            aes(label = label),
            size = 2.7,
            check_overlap = TRUE,
            vjust = -0.5
        ) +
        theme_bw(base_size = 12) +
        labs(
            title = title_text,
            x = "log2 fold change",
            y = "-log10 adjusted P value",
            shape = "Status"
        )

    ggsave(
        file.path(
            outdir,
            paste0(filename_prefix, ".pdf")
        ),
        p,
        width = 7,
        height = 5.5
    )

    ggsave(
        file.path(
            outdir,
            paste0(filename_prefix, ".png")
        ),
        p,
        width = 7,
        height = 5.5,
        dpi = 300
    )
}

make_volcano(
    tab_mi_sham,
    "MI versus Sham",
    "Volcano_MI_vs_Sham"
)

make_volcano(
    tab_mtans_mi,
    "M-Tans versus MI",
    "Volcano_MTans_vs_MI"
)

make_volcano(
    tab_mtans_sham,
    "M-Tans versus Sham",
    "Volcano_MTans_vs_Sham"
)

# ============================================================
# 9. MA plots and dispersion plot
# ============================================================

pdf(
    file.path(outdir, "MA_plots.pdf"),
    width = 7,
    height = 6
)

plotMA(
    res_mi_sham,
    alpha = 0.05,
    ylim = c(-8, 8),
    main = "MI versus Sham"
)

plotMA(
    res_mtans_mi,
    alpha = 0.05,
    ylim = c(-8, 8),
    main = "M-Tans versus MI"
)

plotMA(
    res_mtans_sham,
    alpha = 0.05,
    ylim = c(-8, 8),
    main = "M-Tans versus Sham"
)

dev.off()

pdf(
    file.path(outdir, "DESeq2_dispersion_plot.pdf"),
    width = 7,
    height = 6
)

plotDispEsts(dds)

dev.off()

# ============================================================
# 10. Reversal-signature heatmap
# ============================================================

if (nrow(strict_effect_reversal) >= 2) {

    vsd_model <- vst(
        dds,
        blind = FALSE
    )

    top_n <- min(
        50,
        nrow(strict_effect_reversal)
    )

    heatmap_genes <- head(
        strict_effect_reversal$gene_id,
        top_n
    )

    expression <- assay(vsd_model)[
        heatmap_genes,
        ,
        drop = FALSE
    ]

    row_sd <- apply(
        expression,
        1,
        sd
    )

    expression <- expression[
        row_sd > 0,
        ,
        drop = FALSE
    ]

    expression_z <- t(
        scale(t(expression))
    )

    heatmap_labels <- strict_effect_reversal$gene_name[
        match(
            rownames(expression_z),
            strict_effect_reversal$gene_id
        )
    ]

    heatmap_labels[
        is.na(heatmap_labels) |
        heatmap_labels == ""
    ] <- rownames(expression_z)[
        is.na(heatmap_labels) |
        heatmap_labels == ""
    ]

    rownames(expression_z) <- make.unique(
        heatmap_labels
    )

    sample_annotation <- data.frame(
        Group = colData(dds)$group
    )

    rownames(sample_annotation) <- colnames(dds)

    pheatmap(
        expression_z,
        annotation_col = sample_annotation,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        border_color = NA,
        filename = file.path(
            outdir,
            "Heatmap_top50_MTans_reversal.pdf"
        ),
        width = 8,
        height = 10
    )
}

# ============================================================
# 11. Software and run records
# ============================================================

writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "R_sessionInfo.txt")
)

cat("============================================================\n")
cat("DESEQ2 DIFFERENTIAL ANALYSIS PASSED\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes tested:", nrow(dds), "\n")
cat("Strict reversal genes:", nrow(strict_reversal), "\n")
cat(
    "Strict reversal genes with |log2FC| >= 0.5:",
    nrow(strict_effect_reversal),
    "\n"
)
cat("Output directory:", outdir, "\n")
cat("============================================================\n")
