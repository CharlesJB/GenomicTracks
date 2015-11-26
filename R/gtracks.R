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
#' @param annotation A \code{TxDb} object to use for annotation.
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
gtracks <- function(region, bam_files, ranges = NULL, annotation = NULL) {

}
