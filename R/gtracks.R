#' Produce genomic tracks
#'
#' This function extract the coverages of bam files and print them in track
#' format. It is also possible to add custom annotation in the form of a
#' \code{list} of \code{GRanges} or a \code{GRangesList}. A TxDb object can also
#' be used to add the names of the genes in the selected region.
#'
#' @param region A \code{GRanges} representing the region to display.
#' @param bam_files A \code{list} of bam filename. They must be indexed with
#'                  \code{samtools index} (i.e.: if you have a bam file named
#'                  \code{bam_file.bam}, you must also have a file named
#'                  \code{bam_file.bam.bai} in the same directory). The names of
#'                  the elements of the list will be used as names for the
#'                  displayed tracks.
#' @param ranges A \code{list} of {GRanges} or a \code{GRangesList} of custom
#'               annotations to add to the display. The names of the elements
#'               of the \code{list}/\code{GRangesList} will be used for the
#'                  displayed tracks.
#' @param annotation A \code{TxDb} object to use for annotation must be used
#'                   \code{org} param.
#' @param org The \code{OrgDb} object matching the \code{annotation}.
#'
#' @return A \code{ggbio::Tracks} object.
#'
#' @examples
#' region <- GRanges("chr1", IRanges(34840000, 34860000))
#' bam_files <- get_demo_bam_files()[1]
#' names(bam_files) <- "align1_rep1"
#' gtracks(region = region, bam_files = bam_files)
#'
#' @import ggbio
#' @import tools
#' @import metagene
#'
#' @export
gtracks <- function(region, bam_files, ranges = NULL, annotation = NULL, org = NULL) {
    # Validate params
    stopifnot(class(region) == "GRanges")
    stopifnot(length(region == 1))
    if (!is.null(ranges)) {
        stopifnot(class(ranges) == "list" | class(ranges == "GRangesList"))
        if (class(ranges) == "list") {
            stopifnot(all(lapply(ranges, function(x) class(x) == "GRanges")))
        }
    }
    if (!is.null(annotation)) {
        stopifnot(class(annotation) == "TxDb")
        stopifnot(class(org) == "OrgDb")
    }

    # Extract coverages
    mg <- metagene$new(regions = region, bam_files = bam_files)
    coverages <- mg$get_raw_coverages()
    names(coverages) <- names(bam_files)
    get_cov_df <- function(cov) {
        current_cov <- cov[gr][[1]]
        pos <- seq(start(gr), end(gr))
        data.frame(position = pos, coverage = as.numeric(current_cov))
    }
    df_coverages <- lapply(coverages, get_cov_df)

    # Crunch annotation
    gr_anno <- crunch(txdb, which = resize(region, 100000, fix = "center"))
    symbols <- select(org, keys = as.character(mcols(gr_anno)[["gene_id"]]),
                      column = "SYMBOL", keytype = "ENTREZID")[["SYMBOL"]]
    mcols(gr_anno)[["symbols"]] <- symbols
    i <- mcols(gr_anno)[["type"]] == "gap"
    gr_anno <- gr_anno[!i]
    levels(gr_anno) <- c("cds", "exon", "utr")
    grl_anno <- split(gr_anno, mcols(gr_anno)[["symbols"]])

    # Prepare tracks
    tracks_coverages <- lapply(df_coverages, function(x) {
                                   ggplot(x, aes(x = position, y = coverage)) +
                                       geom_line() + theme_bw()
                               })
    tracks_custom <- lapply(ranges, function(x) autoplot(x) + theme_bw())
    track_anno <- autoplot(grl_anno, aes(type = type)) + theme_bw() + xlim(region)
    tracks(c(tracks_coverages, tracks_custom, Annotation = track_anno)) + xlim(region)
}
