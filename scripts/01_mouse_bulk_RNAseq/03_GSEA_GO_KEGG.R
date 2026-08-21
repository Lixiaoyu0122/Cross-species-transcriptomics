#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
    library(enrichplot)
    library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260720)

run <- "/home/li/Wangxinming/FSTL1_mouse_bulkRNA/rerun_20260720_n9"

deg_dir <- file.path(run, "results", "deseq2_deg")
result_dir <- file.path(run, "results", "GO_KEGG_GSEA")
export_dir <- file.path(
    run, "export_desktop", "06_GO_KEGG_GSEA"
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

comparison_files <- c(
    MI_vs_Sham = "DEG_MI_vs_Sham_all.tsv",
    MTans_vs_MI = "DEG_MTans_vs_MI_all.tsv",
    MTans_vs_Sham = "DEG_MTans_vs_Sham_all.tsv"
)

input_files <- file.path(
    deg_dir,
    unname(comparison_files)
)

missing_files <- input_files[
    !file.exists(input_files) |
    file.info(input_files)$size == 0
]

if (length(missing_files) > 0) {
    stop(
        "Missing DEG files:\n",
        paste(missing_files, collapse = "\n")
    )
}

write_enrichment <- function(object, prefix) {

    all_file <- paste0(prefix, "_all.csv")
    sig_file <- paste0(prefix, "_FDR05.csv")

    if (is.null(object)) {
        write.csv(data.frame(), all_file, row.names = FALSE)
        write.csv(data.frame(), sig_file, row.names = FALSE)
        return(invisible(NULL))
    }

    result <- as.data.frame(object)

    write.csv(
        result,
        all_file,
        row.names = FALSE
    )

    if (
        nrow(result) > 0 &&
        "p.adjust" %in% colnames(result)
    ) {
        significant <- result[
            !is.na(result$p.adjust) &
            result$p.adjust < 0.05,
            ,
            drop = FALSE
        ]
    } else {
        significant <- result[0, , drop = FALSE]
    }

    write.csv(
        significant,
        sig_file,
        row.names = FALSE
    )
}

save_enrichment_plot <- function(
    object,
    prefix,
    title_text
) {

    if (is.null(object)) {
        return(invisible(NULL))
    }

    result <- as.data.frame(object)

    if (
        nrow(result) == 0 ||
        !"p.adjust" %in% colnames(result)
    ) {
        return(invisible(NULL))
    }

    significant_n <- sum(
        !is.na(result$p.adjust) &
        result$p.adjust < 0.05
    )

    if (significant_n == 0) {
        return(invisible(NULL))
    }

    plot_object <- enrichplot::dotplot(
        object,
        showCategory = min(15, significant_n)
    ) +
        ggtitle(title_text)

    ggsave(
        paste0(prefix, ".pdf"),
        plot_object,
        width = 9,
        height = 6
    )

    ggsave(
        paste0(prefix, ".png"),
        plot_object,
        width = 9,
        height = 6,
        dpi = 300
    )
}

map_ensembl <- function(gene_ids) {

    clean_ids <- unique(
        sub("\\..*$", "", gene_ids)
    )

    clean_ids <- clean_ids[
        !is.na(clean_ids) &
        clean_ids != ""
    ]

    mapping <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = clean_ids,
        keytype = "ENSEMBL",
        columns = c(
            "ENTREZID",
            "SYMBOL"
        )
    )

    mapping <- mapping[
        !is.na(mapping$ENTREZID),
        ,
        drop = FALSE
    ]

    unique(mapping)
}

run_go_ora <- function(
    genes,
    universe,
    ontology
) {

    genes <- unique(genes)
    universe <- unique(universe)

    if (length(genes) < 5) {
        return(NULL)
    }

    enrichGO(
        gene = genes,
        universe = universe,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
    )
}

run_kegg_ora <- function(
    genes,
    universe
) {

    genes <- unique(genes)
    universe <- unique(universe)

    if (length(genes) < 5) {
        return(NULL)
    }

    enrichKEGG(
        gene = genes,
        universe = universe,
        organism = "mmu",
        keyType = "ncbi-geneid",
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
}

build_ranked_gene_list <- function(
    result_table,
    mapping
) {

    ranked_table <- merge(
        result_table[
            ,
            c("gene_id", "stat"),
            drop = FALSE
        ],
        mapping[
            ,
            c("ENSEMBL", "ENTREZID"),
            drop = FALSE
        ],
        by.x = "gene_id",
        by.y = "ENSEMBL"
    )

    ranked_table <- ranked_table[
        !is.na(ranked_table$stat) &
        is.finite(ranked_table$stat),
        ,
        drop = FALSE
    ]

    ranked_table <- ranked_table[
        order(
            ranked_table$ENTREZID,
            -abs(ranked_table$stat)
        ),
        ,
        drop = FALSE
    ]

    ranked_table <- ranked_table[
        !duplicated(ranked_table$ENTREZID),
        ,
        drop = FALSE
    ]

    gene_list <- ranked_table$stat
    names(gene_list) <- ranked_table$ENTREZID

    sort(gene_list, decreasing = TRUE)
}

run_go_gsea <- function(
    gene_list,
    ontology
) {

    gseGO(
        geneList = gene_list,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ontology,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
    )
}

run_kegg_gsea <- function(gene_list) {

    gseKEGG(
        geneList = gene_list,
        organism = "mmu",
        keyType = "ncbi-geneid",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
    )
}

summary_list <- list()

for (comparison in names(comparison_files)) {

    cat(
        "\n============================================================\n",
        "Processing: ", comparison, "\n",
        "============================================================\n",
        sep = ""
    )

    result_table <- read.delim(
        file.path(
            deg_dir,
            comparison_files[[comparison]]
        ),
        check.names = FALSE
    )

    required_columns <- c(
        "gene_id",
        "log2FoldChange",
        "stat",
        "pvalue",
        "padj"
    )

    if (
        !all(
            required_columns %in%
            colnames(result_table)
        )
    ) {
        stop(
            "Required DEG columns missing for ",
            comparison
        )
    }

    result_table$gene_id <- sub(
        "\\..*$",
        "",
        result_table$gene_id
    )

    comparison_dir <- file.path(
        result_dir,
        comparison
    )

    ora_dir <- file.path(
        comparison_dir,
        "ORA"
    )

    gsea_dir <- file.path(
        comparison_dir,
        "GSEA"
    )

    dir.create(
        ora_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    dir.create(
        gsea_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    mapping <- map_ensembl(
        result_table$gene_id
    )

    write.csv(
        mapping,
        file.path(
            comparison_dir,
            "Ensembl_to_Entrez_mapping.csv"
        ),
        row.names = FALSE
    )

    tested_ids <- result_table$gene_id[
        !is.na(result_table$pvalue)
    ]

    universe <- unique(
        mapping$ENTREZID[
            mapping$ENSEMBL %in% tested_ids
        ]
    )

    up_ids <- result_table$gene_id[
        !is.na(result_table$padj) &
        result_table$padj < 0.05 &
        result_table$log2FoldChange >= 0.5
    ]

    down_ids <- result_table$gene_id[
        !is.na(result_table$padj) &
        result_table$padj < 0.05 &
        result_table$log2FoldChange <= -0.5
    ]

    up_entrez <- unique(
        mapping$ENTREZID[
            mapping$ENSEMBL %in% up_ids
        ]
    )

    down_entrez <- unique(
        mapping$ENTREZID[
            mapping$ENSEMBL %in% down_ids
        ]
    )

    input_summary <- data.frame(
        comparison = comparison,
        tested_Ensembl = length(unique(tested_ids)),
        background_Entrez = length(universe),
        up_Ensembl = length(unique(up_ids)),
        up_Entrez = length(up_entrez),
        down_Ensembl = length(unique(down_ids)),
        down_Entrez = length(down_entrez)
    )

    write.csv(
        input_summary,
        file.path(
            comparison_dir,
            "enrichment_input_summary.csv"
        ),
        row.names = FALSE
    )

    summary_list[[comparison]] <- input_summary

    for (direction in c("Up", "Down")) {

        genes <- if (
            direction == "Up"
        ) {
            up_entrez
        } else {
            down_entrez
        }

        for (ontology in c("BP", "CC", "MF")) {

            go_result <- run_go_ora(
                genes,
                universe,
                ontology
            )

            prefix <- file.path(
                ora_dir,
                paste0(
                    comparison,
                    "_",
                    direction,
                    "_GO_",
                    ontology
                )
            )

            write_enrichment(
                go_result,
                prefix
            )

            save_enrichment_plot(
                go_result,
                paste0(prefix, "_dotplot"),
                paste(
                    comparison,
                    direction,
                    "GO",
                    ontology
                )
            )
        }

        kegg_prefix <- file.path(
            ora_dir,
            paste0(
                comparison,
                "_",
                direction,
                "_KEGG"
            )
        )

        kegg_result <- tryCatch(
            run_kegg_ora(
                genes,
                universe
            ),
            error = function(e) {

                writeLines(
                    conditionMessage(e),
                    paste0(
                        kegg_prefix,
                        "_ERROR.txt"
                    )
                )

                NULL
            }
        )

        write_enrichment(
            kegg_result,
            kegg_prefix
        )

        save_enrichment_plot(
            kegg_result,
            paste0(
                kegg_prefix,
                "_dotplot"
            ),
            paste(
                comparison,
                direction,
                "KEGG"
            )
        )
    }

    gene_list <- build_ranked_gene_list(
        result_table,
        mapping
    )

    write.csv(
        data.frame(
            ENTREZID = names(gene_list),
            Wald_statistic = as.numeric(gene_list)
        ),
        file.path(
            gsea_dir,
            paste0(
                comparison,
                "_ranked_gene_list.csv"
            )
        ),
        row.names = FALSE
    )

    for (ontology in c("BP", "CC", "MF")) {

        go_gsea <- run_go_gsea(
            gene_list,
            ontology
        )

        prefix <- file.path(
            gsea_dir,
            paste0(
                comparison,
                "_GSEA_GO_",
                ontology
            )
        )

        write_enrichment(
            go_gsea,
            prefix
        )

        save_enrichment_plot(
            go_gsea,
            paste0(prefix, "_dotplot"),
            paste(
                comparison,
                "GSEA GO",
                ontology
            )
        )
    }

    kegg_prefix <- file.path(
        gsea_dir,
        paste0(
            comparison,
            "_GSEA_KEGG"
        )
    )

    kegg_gsea <- tryCatch(
        run_kegg_gsea(gene_list),
        error = function(e) {

            writeLines(
                conditionMessage(e),
                paste0(
                    kegg_prefix,
                    "_ERROR.txt"
                )
            )

            NULL
        }
    )

    write_enrichment(
        kegg_gsea,
        kegg_prefix
    )

    save_enrichment_plot(
        kegg_gsea,
        paste0(
            kegg_prefix,
            "_dotplot"
        ),
        paste(
            comparison,
            "GSEA KEGG"
        )
    )
}

write.csv(
    do.call(
        rbind,
        summary_list
    ),
    file.path(
        result_dir,
        "all_comparisons_input_summary.csv"
    ),
    row.names = FALSE
)

reversal_file <- file.path(
    deg_dir,
    "MTans_reversal_strict_FDR05_absLFC05.tsv"
)

if (
    file.exists(reversal_file) &&
    file.info(reversal_file)$size > 0
) {

    reversal <- read.delim(
        reversal_file,
        check.names = FALSE
    )

    if (nrow(reversal) > 0) {

        reversal$gene_id <- sub(
            "\\..*$",
            "",
            reversal$gene_id
        )

        background_table <- read.delim(
            file.path(
                deg_dir,
                "DEG_MI_vs_Sham_all.tsv"
            ),
            check.names = FALSE
        )

        background_table$gene_id <- sub(
            "\\..*$",
            "",
            background_table$gene_id
        )

        reversal_mapping <- map_ensembl(
            unique(c(
                reversal$gene_id,
                background_table$gene_id
            ))
        )

        reversal_universe <- unique(
            reversal_mapping$ENTREZID[
                reversal_mapping$ENSEMBL %in%
                background_table$gene_id
            ]
        )

        reversal_dir <- file.path(
            result_dir,
            "MTans_reversal_signature_ORA"
        )

        dir.create(
            reversal_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )

        for (direction in c(
            "MI_up_MTans_down",
            "MI_down_MTans_up"
        )) {

            selected_ids <- reversal$gene_id[
                reversal$reversal_direction ==
                direction
            ]

            selected_entrez <- unique(
                reversal_mapping$ENTREZID[
                    reversal_mapping$ENSEMBL %in%
                    selected_ids
                ]
            )

            for (ontology in c("BP", "CC", "MF")) {

                go_result <- run_go_ora(
                    selected_entrez,
                    reversal_universe,
                    ontology
                )

                prefix <- file.path(
                    reversal_dir,
                    paste0(
                        direction,
                        "_GO_",
                        ontology
                    )
                )

                write_enrichment(
                    go_result,
                    prefix
                )

                save_enrichment_plot(
                    go_result,
                    paste0(
                        prefix,
                        "_dotplot"
                    ),
                    paste(
                        direction,
                        "GO",
                        ontology
                    )
                )
            }

            kegg_prefix <- file.path(
                reversal_dir,
                paste0(
                    direction,
                    "_KEGG"
                )
            )

            kegg_result <- tryCatch(
                run_kegg_ora(
                    selected_entrez,
                    reversal_universe
                ),
                error = function(e) {

                    writeLines(
                        conditionMessage(e),
                        paste0(
                            kegg_prefix,
                            "_ERROR.txt"
                        )
                    )

                    NULL
                }
            )

            write_enrichment(
                kegg_result,
                kegg_prefix
            )

            save_enrichment_plot(
                kegg_result,
                paste0(
                    kegg_prefix,
                    "_dotplot"
                ),
                paste(
                    direction,
                    "KEGG"
                )
            )
        }
    }
}

writeLines(
    capture.output(sessionInfo()),
    file.path(
        result_dir,
        "R_sessionInfo_GO_KEGG_GSEA.txt"
    )
)

existing_export_files <- list.files(
    export_dir,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE
)

if (length(existing_export_files) > 0) {
    unlink(
        existing_export_files,
        recursive = FALSE
    )
}

result_files <- list.files(
    result_dir,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE
)

for (source in result_files) {

    relative_path <- substring(
        source,
        nchar(result_dir) + 2
    )

    destination <- file.path(
        export_dir,
        relative_path
    )

    dir.create(
        dirname(destination),
        recursive = TRUE,
        showWarnings = FALSE
    )

    if (
        !file.copy(
            source,
            destination,
            overwrite = TRUE
        )
    ) {
        stop(
            "Failed to export: ",
            source
        )
    }
}

cat("============================================================\n")
cat("GO/KEGG/GSEA PASSED\n")
cat("Comparisons:", length(comparison_files), "\n")
cat("Result directory:", result_dir, "\n")
cat("Export directory:", export_dir, "\n")
cat("============================================================\n")
