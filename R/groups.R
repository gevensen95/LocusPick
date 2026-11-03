#' Standardize a differential expression table
#'
#' @param df A `data.frame` of differential expression results.
#' @param gene,lfc,padj Column names holding the gene symbol, log2 fold change
#'   and adjusted p-value.
#' @return A `data.frame` with columns `gene`, `lfc`, `padj`. Rows with a
#'   missing gene or p-value are dropped; for duplicated genes the row with the
#'   smallest adjusted p-value is kept (with a message).
#' @examples
#' deTable(data.frame(Gene = c("a", "b"), log2FoldChange = c(1, -2),
#'                    padj = c(0.01, 0.2)),
#'         gene = "Gene", lfc = "log2FoldChange", padj = "padj")
#' @export
deTable <- function(df, gene = "gene", lfc = "log2FoldChange", padj = "padj") {
    .stopIfMissingCols(df, c(gene, lfc, padj), "DE table")
    out <- data.frame(gene = as.character(df[[gene]]),
                      lfc = as.numeric(df[[lfc]]),
                      padj = as.numeric(df[[padj]]),
                      stringsAsFactors = FALSE)
    out <- out[!is.na(out$gene) & out$gene != "" & !is.na(out$padj), ]
    out <- out[order(out$padj), ]
    dup <- duplicated(out$gene)
    if (any(dup))
        message(sum(dup), " duplicated gene(s) in DE table; kept the row with ",
                "the smallest adjusted p-value.")
    out <- out[!dup, ]
    rownames(out) <- NULL
    out
}

#' Define gene groups from a perturbation and a condition contrast
#'
#' Places genes into the four combinations of *regulation by the factor*
#' (from a perturbation such as a knockout) and *direction in a condition*
#' (such as aging or treatment). For a loss-of-function perturbation, genes
#' that go down when the factor is removed are `activated` by the factor and
#' genes that go up are `repressed`.
#'
#' @param perturbation,condition DE results: tables from [deTable()] or
#'   [readDE()], or anything `readDE()` accepts (a DESeq2/edgeR/limma object
#'   or a results file), which is then read with default settings.
#' @param targets Optional character vector restricting groups to these genes
#'   (for example genes with a binding peak).
#' @param alpha Adjusted p-value cut-off, applied to both contrasts.
#' @param minLFC Minimum absolute log2 fold change, applied to both contrasts.
#'   A length-2 vector sets perturbation and condition separately.
#' @param perturbationType `"loss"` (knockout/knockdown) or `"gain"`
#'   (overexpression).
#' @param conditionName Label for the condition used in group names.
#' @return A `data.frame` with `gene`, `regulation`, `condition_dir`, `group`,
#'   and the fold changes and p-values of both contrasts.
#' @examples
#' ko  <- deTable(data.frame(gene = c("a", "b", "c"), lfc = c(-2, 1.5, -1),
#'                           padj = c(1e-5, 1e-3, 0.5)), lfc = "lfc")
#' age <- deTable(data.frame(gene = c("a", "b", "c"), lfc = c(-1, 1, 2),
#'                           padj = c(1e-3, 1e-4, 1e-6)), lfc = "lfc")
#' defineGroups(ko, age, conditionName = "age")
#' @export
defineGroups <- function(perturbation, condition, targets = NULL, alpha = 0.05,
                         minLFC = 0, perturbationType = c("loss", "gain"),
                         conditionName = "condition") {
    perturbationType <- match.arg(perturbationType)
    perturbation <- .asDE(perturbation)
    condition <- .asDE(condition)
    minLFC <- rep_len(minLFC, 2)

    sigP <- perturbation[perturbation$padj <= alpha & abs(perturbation$lfc) > minLFC[1], ]
    sigC <- condition[condition$padj <= alpha & abs(condition$lfc) > minLFC[2], ]
    if (!is.null(targets)) sigP <- sigP[sigP$gene %in% targets, ]
    m <- merge(sigP, sigC, by = "gene", suffixes = c("_perturbation", "_condition"))
    if (!nrow(m)) {
        warning("No genes are significant in both contrasts.", call. = FALSE)
    }
    actSign <- if (perturbationType == "loss") -1 else 1
    m$regulation <- ifelse(sign(m$lfc_perturbation) == actSign, "activated", "repressed")
    m$condition_dir <- ifelse(m$lfc_condition > 0, "up", "down")
    m$group <- paste(m$regulation, conditionName, m$condition_dir, sep = "_")
    lev <- paste(rep(c("activated", "repressed"), each = 2), conditionName,
                 c("down", "up"), sep = "_")
    m$group <- factor(m$group, levels = lev)
    m[order(m$group, m$gene),
      c("gene", "regulation", "condition_dir", "group", "lfc_perturbation",
        "padj_perturbation", "lfc_condition", "padj_condition")]
}

.asDE <- function(x) {
    if (is.data.frame(x) && identical(sort(intersect(names(x), c("gene", "lfc", "padj"))),
                                      c("gene", "lfc", "padj")))
        return(x)
    readDE(x)
}
