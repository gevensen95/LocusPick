#' Harmonize gene symbols to current official symbols
#'
#' Gene symbols change between annotation releases (for example Riken clones
#' that receive names). Mixing sources from different releases silently drops
#' genes from joins. This maps aliases and previous symbols to the current
#' official symbol using an `OrgDb`.
#'
#' Symbols that are already official, unknown, or ambiguous (an alias of more
#' than one gene) are returned unchanged.
#'
#' @param x Character vector of gene symbols.
#' @param orgDb An `OrgDb` object, e.g. `org.Mm.eg.db::org.Mm.eg.db`.
#' @return Character vector the same length as `x`, with attribute `changes`:
#'   a `data.frame` of the symbols that were renamed.
#' @examples
#' if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
#'     harmonizeSymbols(c("Hnf4a", "1600002H07Rik"), org.Mm.eg.db::org.Mm.eg.db)
#' }
#' @export
harmonizeSymbols <- function(x, orgDb) {
    .requirePkg("AnnotationDbi", "for symbol harmonization")
    x <- as.character(x)
    u <- unique(x[!is.na(x)])
    cur <- AnnotationDbi::keys(orgDb, keytype = "SYMBOL")
    need <- setdiff(u, cur)
    need <- intersect(need, AnnotationDbi::keys(orgDb, keytype = "ALIAS"))
    map <- stats::setNames(u, u)
    if (length(need)) {
        al <- suppressMessages(AnnotationDbi::mapIds(orgDb, keys = need, column = "SYMBOL",
                                                     keytype = "ALIAS", multiVals = "list"))
        one <- vapply(al, function(v) {
            v <- unique(v[!is.na(v)])
            if (length(v) == 1) v else NA_character_
        }, character(1))
        ok <- !is.na(one)
        map[names(one)[ok]] <- one[ok]
    }
    out <- unname(map[x])
    out[is.na(x)] <- NA_character_
    ch <- data.frame(original = names(map), harmonized = unname(map),
                     stringsAsFactors = FALSE)
    attr(out, "changes") <- ch[ch$original != ch$harmonized, , drop = FALSE]
    out
}
