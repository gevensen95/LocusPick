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

#' Define gene groups from one or two contrasts
#'
#' With two contrasts, places genes into the four combinations of
#' *regulation by the factor* (from a perturbation such as a knockout) and
#' *direction in a condition* (such as aging or treatment). For a
#' loss-of-function perturbation, genes that go down when the factor is
#' removed are `activated` by the factor and genes that go up are `repressed`.
#'
#' Either contrast can be omitted:
#' * perturbation only: groups `activated` and `repressed`;
#' * condition only: groups `<conditionName>_up` and `<conditionName>_down`.
#'
#' To pick loci without any expression data, skip this function and call
#' [buildCandidates()] with `groups = NULL`.
#'
#' @param perturbation,condition DE results: tables from [deTable()] or
#'   [readDE()], or anything `readDE()` accepts (a DESeq2/edgeR/limma object
#'   or a results file). Either may be `NULL`, but not both.
#' @param targets Optional character vector restricting groups to these genes
#'   (for example genes with a binding peak).
#' @param alpha Adjusted p-value cut-off, applied to each contrast given.
#' @param minLFC Minimum absolute log2 fold change. A length-2 vector sets
#'   perturbation and condition separately.
#' @param perturbationType `"loss"` (knockout/knockdown) or `"gain"`
#'   (overexpression).
#' @param conditionName Label for the condition used in group names.
#' @return A `data.frame` with `gene`, `group` (a factor), and the direction,
#'   fold change and adjusted p-value columns of each contrast given
#'   (`regulation`, `lfc_perturbation`, `padj_perturbation`; `condition_dir`,
#'   `lfc_condition`, `padj_condition`).
#' @examples
#' ko  <- deTable(data.frame(gene = c("a", "b", "c"), lfc = c(-2, 1.5, -1),
#'                           padj = c(1e-5, 1e-3, 0.5)), lfc = "lfc")
#' age <- deTable(data.frame(gene = c("a", "b", "c"), lfc = c(-1, 1, 2),
#'                           padj = c(1e-3, 1e-4, 1e-6)), lfc = "lfc")
#' defineGroups(ko, age, conditionName = "age")
#' defineGroups(ko)                                   # perturbation only
#' defineGroups(condition = age, conditionName = "age")  # condition only
#' @export
defineGroups <- function(perturbation = NULL, condition = NULL, targets = NULL,
                         alpha = 0.05, minLFC = 0, perturbationType = c("loss", "gain"),
                         conditionName = "condition") {
    perturbationType <- match.arg(perturbationType)
    if (is.null(perturbation) && is.null(condition))
        stop("Give at least one contrast. To rank loci without expression data, ",
             "use buildCandidates(groups = NULL).", call. = FALSE)
    minLFC <- rep_len(minLFC, 2)
    sig <- function(d, m) d[d$padj <= alpha & abs(d$lfc) > m, , drop = FALSE]

    if (!is.null(perturbation)) {
        sigP <- sig(.asDE(perturbation), minLFC[1])
        if (!is.null(targets)) sigP <- sigP[sigP$gene %in% targets, ]
        names(sigP)[2:3] <- c("lfc_perturbation", "padj_perturbation")
        actSign <- if (perturbationType == "loss") -1 else 1
        sigP$regulation <- ifelse(sign(sigP$lfc_perturbation) == actSign,
                                  "activated", "repressed")
    }
    if (!is.null(condition)) {
        sigC <- sig(.asDE(condition), minLFC[2])
        if (!is.null(targets) && is.null(perturbation)) sigC <- sigC[sigC$gene %in% targets, ]
        names(sigC)[2:3] <- c("lfc_condition", "padj_condition")
        sigC$condition_dir <- ifelse(sigC$lfc_condition > 0, "up", "down")
    }

    if (!is.null(perturbation) && !is.null(condition)) {
        m <- merge(sigP, sigC, by = "gene")
        if (!nrow(m)) warning("No genes are significant in both contrasts.", call. = FALSE)
        m$group <- if (nrow(m)) paste(m$regulation, conditionName, m$condition_dir, sep = "_") else character()
        lev <- paste(rep(c("activated", "repressed"), each = 2), conditionName,
                     c("down", "up"), sep = "_")
        cols <- c("gene", "regulation", "condition_dir", "group", "lfc_perturbation",
                  "padj_perturbation", "lfc_condition", "padj_condition")
    } else if (!is.null(perturbation)) {
        m <- sigP
        m$group <- m$regulation
        lev <- c("activated", "repressed")
        cols <- c("gene", "regulation", "group", "lfc_perturbation", "padj_perturbation")
    } else {
        m <- sigC
        m$group <- if (nrow(m)) paste(conditionName, m$condition_dir, sep = "_") else character()
        lev <- paste(conditionName, c("down", "up"), sep = "_")
        cols <- c("gene", "condition_dir", "group", "lfc_condition", "padj_condition")
    }
    if (!nrow(m)) warning("No significant genes.", call. = FALSE)
    m$group <- factor(m$group, levels = lev)
    m <- m[order(m$group, m$gene), cols]
    rownames(m) <- NULL
    m
}

.asDE <- function(x) {
    if (is.data.frame(x) && identical(sort(intersect(names(x), c("gene", "lfc", "padj"))),
                                      c("gene", "lfc", "padj")))
        return(x)
    readDE(x)
}
