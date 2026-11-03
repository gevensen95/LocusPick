#' Read differential expression results from common tools
#'
#' Accepts the direct output of DESeq2, edgeR or limma (as R objects or saved
#' tables) and returns the standard `gene`, `lfc`, `padj` table used by
#' [defineGroups()]. The format is detected automatically and reported.
#'
#' Supported inputs:
#' * DESeq2: a `DESeqResults` object, or a table with `log2FoldChange`/`padj`.
#' * edgeR: a `TopTags` or `DGELRT`/`DGEExact` object, or a table with
#'   `logFC`/`FDR`.
#' * limma: an `MArrayLM` fit (with `coef`), or a `topTable()` table with
#'   `logFC`/`adj.P.Val`.
#' * Any `data.frame` or CSV/TSV file with recognizable columns; name them
#'   explicitly with `gene`, `lfc`, `padj` if detection fails.
#'
#' Gene identifiers can be symbols, Ensembl IDs (with or without version
#' suffix) or Entrez IDs. With `orgDb`, IDs are converted to current symbols.
#'
#' **Contrast direction.** Fold changes must be *test / reference* (e.g. KO /
#' WT, old / young). DESeq2 and edgeR objects record their contrast, which is
#' printed. Plain tables do not, so pass `direction` (e.g. `"KO/WT"`) to
#' confirm it; otherwise a warning is issued.
#'
#' @param x A results object, `data.frame`, or path to a CSV/TSV file.
#' @param gene,lfc,padj Optional column names, overriding detection. `gene`
#'   may be `"rownames"`.
#' @param orgDb Optional `OrgDb` (e.g. `org.Mm.eg.db::org.Mm.eg.db`) to convert
#'   Ensembl/Entrez IDs to symbols and harmonize outdated symbols.
#' @param direction Optional label stating the contrast, e.g. `"KO/WT"`.
#' @param coef For limma `MArrayLM` fits: the coefficient to extract.
#' @param verbose Print what was detected.
#' @return A `data.frame` with `gene`, `lfc`, `padj`, and attributes
#'   `source` (detected format) and `contrast`.
#' @examples
#' de <- read.delim(system.file("extdata", "toy_de.tsv", package = "LocusPick"))
#' tab <- data.frame(gene = de$gene, log2FoldChange = de$ko_lfc, padj = de$ko_padj)
#' readDE(tab, direction = "KO/WT")
#' @export
readDE <- function(x, gene = NULL, lfc = NULL, padj = NULL, orgDb = NULL,
                   direction = NULL, coef = NULL, verbose = TRUE) {
    contrast <- NULL
    src <- "table"
    label <- if (is.character(x) && length(x) == 1) basename(x) else
        deparse(substitute(x))[1]

    if (is(x, "DESeqResults")) {
        src <- "DESeq2 results"
        d <- as.character(S4Vectors::mcols(x)$description)
        lf <- grep("^log2 fold change", d, value = TRUE)
        if (length(lf)) contrast <- sub("^log2 fold change \\([^)]*\\):\\s*", "", lf[1])
        x <- as.data.frame(x)
    } else if (inherits(x, "TopTags")) {
        src <- "edgeR topTags"
        contrast <- x$comparison
        x <- x$table
    } else if (inherits(x, c("DGELRT", "DGEExact"))) {
        src <- paste("edgeR", class(x)[1])
        contrast <- x$comparison
        x <- .requireCall("edgeR", "topTags")(x, n = Inf)$table
    } else if (inherits(x, "MArrayLM")) {
        src <- "limma fit"
        if (is.null(coef)) {
            if (ncol(x$coefficients) > 1)
                stop("This limma fit has several coefficients; choose one with `coef`.",
                     call. = FALSE)
            coef <- 1
        }
        contrast <- if (is.character(coef)) coef else colnames(x$coefficients)[coef]
        x <- .requireCall("limma", "topTable")(x, coef = coef, number = Inf, sort.by = "none")
    } else if (is.character(x) && length(x) == 1) {
        if (!file.exists(x)) stop("DE file not found: ", x, call. = FALSE)
        sep <- if (grepl("\\.csv(\\.gz)?$", x, ignore.case = TRUE)) "," else "\t"
        x <- utils::read.table(x, header = TRUE, sep = sep, quote = "\"",
                               comment.char = "", check.names = FALSE,
                               stringsAsFactors = FALSE)
        # write.csv(row.names = TRUE) leaves an empty first header
        if (names(x)[1] %in% c("", "X", "V1")) names(x)[1] <- "rownames_col"
    } else if (is(x, "DataFrame")) {
        x <- as.data.frame(x)
    }
    if (!is.data.frame(x)) stop("Unsupported DE input of class '", class(x)[1], "'.",
                                call. = FALSE)

    cols <- names(x)
    pick <- function(user, candidates, what) {
        if (!is.null(user)) return(user)
        hit <- candidates[tolower(candidates) %in% tolower(cols)]
        if (!length(hit)) return(NULL)
        cols[match(tolower(hit[1]), tolower(cols))]
    }
    lfcCol  <- pick(lfc, c("log2FoldChange", "logFC", "log2FC", "log2_fold_change",
                           "avg_log2FC", "log2.fold_change.", "lfc"), "lfc")
    padjCol <- pick(padj, c("padj", "FDR", "adj.P.Val", "p_val_adj", "qvalue",
                            "q_value", "p.adjust", "padj_BH", "FDR.BH"), "padj")
    geneCol <- pick(gene, c("gene", "symbol", "gene_name", "gene_symbol", "SYMBOL",
                            "Gene.name", "external_gene_name", "Geneid", "gene_id",
                            "ID", "rownames_col"), "gene")
    if (is.null(lfcCol) || is.null(padjCol))
        stop("Could not find fold-change and adjusted p-value columns among: ",
             paste(cols, collapse = ", "), ". Name them with `lfc` and `padj`.",
             call. = FALSE)
    if (src == "table")
        src <- switch(lfcCol, log2FoldChange = "DESeq2-style table",
                      logFC = if (padjCol == "adj.P.Val") "limma-style table" else "edgeR-style table",
                      "table")

    ids <- if (is.null(geneCol) || identical(geneCol, "rownames")) rownames(x) else
        as.character(x[[geneCol]])
    if (is.null(ids) || identical(as.character(ids), as.character(seq_along(ids))))
        stop("No gene identifiers found; name the column with `gene`.", call. = FALSE)

    out <- data.frame(gene = ids, lfc = as.numeric(x[[lfcCol]]),
                      padj = as.numeric(x[[padjCol]]), stringsAsFactors = FALSE)
    idType <- .idType(out$gene)
    if (idType != "symbol") {
        if (is.null(orgDb))
            stop("Gene IDs look like ", idType, " IDs. Pass `orgDb` (e.g. ",
                 "org.Mm.eg.db::org.Mm.eg.db) to convert them to symbols.", call. = FALSE)
        out$gene <- .idsToSymbols(out$gene, orgDb, idType)
    } else if (!is.null(orgDb)) {
        out$gene <- harmonizeSymbols(out$gene, orgDb)
    }
    nIn <- nrow(out)
    out <- suppressMessages(deTable(out, gene = "gene", lfc = "lfc", padj = "padj"))

    if (verbose) {
        message(sprintf("%s: %s | lfc = '%s', padj = '%s', gene = %s (%s IDs) | %d genes",
                        label, src, lfcCol, padjCol,
                        if (is.null(geneCol)) "row names" else paste0("'", geneCol, "'"),
                        idType, nrow(out)))
        if (nrow(out) < nIn)
            message("  ", nIn - nrow(out), " rows dropped (missing gene/p-value or duplicates).")
        if (!is.null(contrast)) message("  contrast: ", contrast,
                                        if (!is.null(direction)) paste0("  (stated: ", direction, ")"))
        else if (!is.null(direction)) message("  contrast (stated): ", direction)
    }
    if (is.null(contrast) && is.null(direction))
        warning(label, ": the contrast direction is not recorded in a plain table. ",
                "Confirm fold changes are test/reference (e.g. KO/WT, old/young) and ",
                "pass `direction` to silence this warning.", call. = FALSE)
    attr(out, "source") <- src
    attr(out, "contrast") <- if (!is.null(contrast)) contrast else direction
    out
}

.requireCall <- function(pkg, fun) {
    .requirePkg(pkg, paste0("to read ", pkg, " objects"))
    getExportedValue(pkg, fun)
}

.idType <- function(ids) {
    ids <- ids[!is.na(ids)]
    if (!length(ids)) return("symbol")
    if (mean(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", ids)) > 0.5) return("Ensembl")
    if (mean(grepl("^[0-9]+$", ids)) > 0.5) return("Entrez")
    "symbol"
}

.idsToSymbols <- function(ids, orgDb, idType) {
    .requirePkg("AnnotationDbi", "to convert gene IDs")
    key <- if (idType == "Ensembl") sub("\\.[0-9]+$", "", ids) else ids
    kt <- if (idType == "Ensembl") "ENSEMBL" else "ENTREZID"
    sym <- suppressMessages(AnnotationDbi::mapIds(orgDb, keys = unique(key[!is.na(key)]),
                                                  column = "SYMBOL", keytype = kt,
                                                  multiVals = "first"))
    out <- unname(sym[key])
    message("  converted ", sum(!is.na(out)), "/", length(ids), " ", idType,
            " IDs to symbols")
    out
}
