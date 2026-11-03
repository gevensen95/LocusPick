#' Plot a locus
#'
#' Draws signal tracks, gene models and an optional expression inset for one
#' candidate locus.
#'
#' @param locus One row of a candidate table (needs `seqnames`, `start`, `end`,
#'   `gene`), or a `GRanges` of length 1 with a `gene` column.
#' @param signal Named list of signal sources (see [signalSources]).
#' @param genes Gene models from [readGeneModels()].
#' @param trackGroups Optional named list mapping a panel label to the track
#'   names drawn in it, e.g. `list(CUT&RUN = c("CNR_young", "CNR_old"))`.
#'   Tracks in one panel share a y-axis, so conditions are directly comparable;
#'   different panels have independent y-axes. Default: one panel per track.
#' @param colors Optional named vector of colors per track.
#' @param window Optional length-2 numeric giving the plotted region. If `NULL`,
#'   the window spans the peak and the gene TSS plus `flank`, or only the peak
#'   plus `flank` when that would exceed `maxWindow`.
#' @param flank,maxWindow See `window`.
#' @param expression Optional long `data.frame` with columns `gene`, `panel`,
#'   `group`, `value` (one row per sample) for the expression inset.
#' @param logExpression Plot `log2(value + 1)` in the inset.
#' @param yLimits Optional named list of `c(min, max)` per panel label, for
#'   identical axes across several loci in a figure.
#' @param annotate Add a title and subtitle with the candidate's statistics.
#'   Set `FALSE` for figure panels.
#' @param highlight Shade the peak.
#' @return A `patchwork` object.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
#' plotLocus(cand[1, ], ex$signal, ex$genes,
#'           trackGroups = list(Signal = c("A", "B")),
#'           expression = ex$expression)
#' @export
plotLocus <- function(locus, signal, genes, trackGroups = NULL, colors = NULL,
                      window = NULL, flank = 3000, maxWindow = 40000,
                      expression = NULL, logExpression = TRUE, yLimits = NULL,
                      annotate = TRUE, highlight = TRUE) {
    loc <- .asLocus(locus)
    gm <- .asGeneModels(genes)
    if (!is.list(signal) || is.null(names(signal)))
        stop("signal must be a named list of sources.", call. = FALSE)
    if (is.null(trackGroups))
        trackGroups <- stats::setNames(as.list(names(signal)), names(signal))
    miss <- setdiff(unlist(trackGroups), names(signal))
    if (length(miss)) stop("trackGroups name unknown track(s): ",
                           paste(miss, collapse = ", "), call. = FALSE)
    colors <- .trackColors(names(signal), colors)

    win <- if (is.null(window)) .locusWindow(loc, gm, flank, maxWindow) else window
    wgr <- GenomicRanges::GRanges(loc$seqnames, IRanges::IRanges(win[1], win[2]))

    xsc <- ggplot2::scale_x_continuous(limits = win, expand = c(0, 0),
                                       labels = function(b) format(b, big.mark = ",",
                                                                   scientific = FALSE))
    panels <- lapply(names(trackGroups), function(lab) {
        tr <- trackGroups[[lab]]
        d <- do.call(rbind, lapply(tr, function(t) {
            s <- fetchSignal(.asSource(signal[[t]]), wgr)
            if (!length(s)) return(NULL)
            data.frame(track = t,
                       xmin = pmax(GenomicRanges::start(s), win[1]),
                       xmax = pmin(GenomicRanges::end(s) + 1, win[2]),
                       score = s$score)
        }))
        .trackPanel(d, tr, lab, colors, loc, xsc, yLimits[[lab]], highlight)
    })
    gp <- .genePanel(gm, loc, win, xsc)
    tracks <- patchwork::wrap_plots(c(panels, list(gp)), ncol = 1,
                                    heights = c(rep(2, length(panels)), 1.2))

    p <- if (!is.null(expression)) {
        patchwork::wrap_elements(tracks) | .exprPanel(expression, loc$gene, logExpression)
    } else tracks
    if (!is.null(expression)) p <- p + patchwork::plot_layout(widths = c(3, 1))
    if (annotate) p <- p + patchwork::plot_annotation(
        title = .locusTitle(loc), subtitle = .locusSubtitle(loc),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                               plot.subtitle = ggplot2::element_text(size = 7)))
    p
}

.asLocus <- function(locus) {
    if (is(locus, "GRanges")) {
        if (length(locus) != 1) stop("locus must have length 1.", call. = FALSE)
        df <- as.data.frame(locus)
        df$seqnames <- as.character(df$seqnames)
        locus <- df
    }
    if (!is.data.frame(locus) || nrow(locus) != 1)
        stop("locus must be one row of a candidate table.", call. = FALSE)
    .stopIfMissingCols(locus, c("seqnames", "start", "end", "gene"), "locus")
    as.list(locus)
}

.locusWindow <- function(loc, gm, flank, maxWindow) {
    g <- gm$genes[gm$genes$gene_name == loc$gene &
                  as.character(GenomicRanges::seqnames(gm$genes)) == loc$seqnames]
    tss <- if (length(g)) g$tss[which.min(abs(g$tss - loc$start))] else NA
    w <- c(min(loc$start, tss, na.rm = TRUE), max(loc$end, tss, na.rm = TRUE)) +
        c(-flank, flank)
    if (diff(w) > maxWindow) w <- c(loc$start - flank, loc$end + flank)
    c(max(1, w[1]), w[2])
}

.trackColors <- function(tracks, colors) {
    base <- grDevices::hcl.colors(max(3, length(tracks)), "Dark 3")
    out <- stats::setNames(base[seq_along(tracks)], tracks)
    if (!is.null(colors)) out[names(colors)] <- colors
    out
}

.theme <- function() {
    ggplot2::theme_classic(base_size = 7) +
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
}

.trackPanel <- function(d, tracks, label, colors, loc, xsc, ylim, highlight) {
    p <- ggplot2::ggplot()
    if (highlight)
        p <- p + ggplot2::annotate("rect", xmin = loc$start, xmax = loc$end,
                                   ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.25)
    if (is.null(d) || !nrow(d)) {
        p <- p + ggplot2::annotate("text", x = mean(c(loc$start, loc$end)), y = 0.5,
                                   label = "no signal in window", size = 2.5) +
            ggplot2::scale_y_continuous(label, limits = c(0, 1))
    } else {
        d <- do.call(rbind, lapply(split(d, d$track), .stepCoords))
        d$track <- factor(d$track, levels = tracks)
        p <- p + ggplot2::geom_ribbon(data = d, ggplot2::aes(
                x = .data$x, ymin = 0, ymax = .data$y, fill = .data$track),
                colour = NA, outline.type = "full") +
            ggplot2::facet_grid(track ~ ., drop = FALSE) +
            ggplot2::scale_fill_manual(values = colors, guide = "none") +
            ggplot2::scale_y_continuous(label, limits = ylim,
                                        expand = ggplot2::expansion(mult = c(0, 0.05)))
    }
    p + xsc + .theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank())
}

# Convert signal intervals to a step outline (drops to 0 across gaps) so the
# track renders as one filled shape instead of many abutting rectangles.
.stepCoords <- function(d) {
    d <- d[order(d$xmin), ]
    n <- nrow(d)
    gap <- c(TRUE, d$xmin[-1] > d$xmax[-n])   # interval starts after a gap
    end <- c(gap[-1], TRUE)                    # interval ends before a gap
    x <- y <- vector("list", n)
    for (i in seq_len(n)) {
        xi <- c(d$xmin[i], d$xmax[i]); yi <- c(d$score[i], d$score[i])
        if (gap[i]) { xi <- c(d$xmin[i], xi); yi <- c(0, yi) }
        if (end[i]) { xi <- c(xi, d$xmax[i]); yi <- c(yi, 0) }
        x[[i]] <- xi; y[[i]] <- yi
    }
    data.frame(track = d$track[1], x = unlist(x), y = unlist(y))
}

.genePanel <- function(gm, loc, win, xsc) {
    chr <- loc$seqnames
    inWin <- function(gr) gr[as.character(GenomicRanges::seqnames(gr)) == chr &
                             GenomicRanges::end(gr) >= win[1] &
                             GenomicRanges::start(gr) <= win[2]]
    g <- .matchSeqStyle(gm$genes, GenomicRanges::GRanges(chr, IRanges::IRanges(1, 1)))
    e <- .matchSeqStyle(gm$exons, GenomicRanges::GRanges(chr, IRanges::IRanges(1, 1)))
    g <- inWin(g); e <- inWin(e)
    p <- ggplot2::ggplot()
    if (length(g)) {
        lab <- paste0(g$gene_name, " (", as.character(GenomicRanges::strand(g)), ")")
        y <- as.integer(factor(lab, levels = unique(lab)))
        gd <- data.frame(x = pmax(GenomicRanges::start(g), win[1]),
                         xend = pmin(GenomicRanges::end(g), win[2]), y = y,
                         target = g$gene_name == loc$gene)
        p <- p + ggplot2::geom_segment(data = gd, ggplot2::aes(
            x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y,
            colour = .data$target))
        if (length(e)) {
            ed <- data.frame(xmin = pmax(GenomicRanges::start(e), win[1]),
                             xmax = pmin(GenomicRanges::end(e), win[2]),
                             gene = e$gene_name)
            ed$y <- y[match(ed$gene, g$gene_name)]
            ed <- ed[!is.na(ed$y), ]
            ed$target <- ed$gene == loc$gene
            p <- p + ggplot2::geom_rect(data = ed, ggplot2::aes(
                xmin = .data$xmin, xmax = .data$xmax, ymin = .data$y - 0.3,
                ymax = .data$y + 0.3, fill = .data$target), colour = NA)
        }
        p <- p + ggplot2::scale_y_continuous(breaks = seq_along(unique(lab)),
                                             labels = unique(lab),
                                             expand = ggplot2::expansion(add = 0.6)) +
            ggplot2::scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "black"),
                                         guide = "none") +
            ggplot2::scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "black"),
                                       guide = "none")
    }
    p + xsc + ggplot2::labs(x = paste0(chr, " (bp)"), y = NULL) + .theme()
}

.exprPanel <- function(expression, gene, logExpression) {
    .stopIfMissingCols(expression, c("gene", "panel", "group", "value"), "expression")
    d <- expression[expression$gene == gene, ]
    if (!nrow(d)) return(patchwork::plot_spacer())
    d$y <- if (logExpression) log2(d$value + 1) else d$value
    d$group <- factor(d$group, levels = unique(expression$group))
    d$panel <- factor(d$panel, levels = unique(expression$panel))
    ggplot2::ggplot(d, ggplot2::aes(.data$group, .data$y)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90") +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1),
                            size = 0.8) +
        ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free") +
        ggplot2::labs(x = NULL, y = if (logExpression) "log2(expression + 1)" else "expression") +
        .theme()
}

.locusTitle <- function(loc) {
    parts <- c(loc$gene, if (!is.null(loc$group)) as.character(loc$group),
               if (!is.null(loc$feature)) as.character(loc$feature))
    paste(parts, collapse = "  |  ")
}

.locusSubtitle <- function(loc) {
    f <- function(nm, lab, fmt = "%.2f")
        if (!is.null(loc[[nm]]) && !is.na(loc[[nm]])) sprintf(paste0(lab, "=", fmt), loc[[nm]])
    s <- c(paste0(loc$seqnames, ":", loc$start, "-", loc$end),
           f("lfc_perturbation", "perturbation log2FC"),
           f("lfc_condition", "condition log2FC"),
           f("mappability", "mappability"),
           f("score", "score"))
    paste(s, collapse = " | ")
}

#' Write a review report of candidate loci
#'
#' Plots each locus on its own PDF page and, optionally, writes the full
#' candidate table. A locus that fails to plot is reported and skipped.
#'
#' @param loci Rows to plot (e.g. from [topLoci()]).
#' @param file Output PDF path.
#' @param ... Passed to [plotLocus()] (`signal`, `genes`, ...).
#' @param candidates Optional full candidate table to write as TSV.
#' @param tsv Path for the candidate TSV.
#' @param width,height PDF page size in inches.
#' @return Invisibly, a logical vector: which loci were plotted.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
#' top <- topLoci(rankLoci(cand, c(A_mean = 1)), n = 1)
#' f <- tempfile(fileext = ".pdf")
#' reportLoci(top, f, signal = ex$signal, genes = ex$genes)
#' @export
reportLoci <- function(loci, file, ..., candidates = NULL, tsv = NULL,
                       width = 9, height = 6) {
    if (!is.null(candidates)) {
        if (is.null(tsv)) tsv <- sub("\\.pdf$", "_candidates.tsv", file)
        utils::write.table(candidates, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    grDevices::pdf(file, width = width, height = height)
    on.exit(grDevices::dev.off())
    ok <- vapply(seq_len(nrow(loci)), function(i) {
        tryCatch({
            print(plotLocus(loci[i, , drop = FALSE], ...))
            TRUE
        }, error = function(e) {
            message("Skipped ", loci$gene[i], " (", loci$peak_id[i], "): ",
                    conditionMessage(e))
            FALSE
        })
    }, logical(1))
    message("Plotted ", sum(ok), "/", length(ok), " loci to ", file)
    invisible(ok)
}
