#' StructuralVariantAnnotation: a package for SV annotation
#'
#' StructuralVariantAnnotation contains useful helper functions for reading
#' and interpreting structural variants calls. The packages contains functions
#' for parsing VCFs from a number of popular caller as well as functions for
#' dealing with breakpoints involving two separate genomic loci. The package
#' takes a `GRanges` based breakend-centric approach.
#'
#'    * Parse VCF objects with the `breakpointRanges()` and `breakendRanges()`functions.
#'    * Find breakpoint overlaps with the `findBreakpointOverlaps()`
#'   and `countBreakpointOverlaps()` functions.
#'    * Generate BEDPE files for circos plot with `breakpointgr2pairs()` function.
#'    * ...
#'
#' For more details on the features of StructuralVariantAnnotation, read the vignette:
#' `browseVignettes(package = "StructuralVariantAnnotation")`
#'
#' @docType package
#' @name StructuralVariantAnnotation
#' @import BiocGenerics
#' @import VariantAnnotation
#' @import rtracklayer
#' @importFrom Biostrings DNAStringSet reverseComplement getSeq nchar subseq
#' @importFrom pwalign pairwiseAlignment nindel insertion deletion nucleotideSubstitutionMatrix
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import IRanges
#' @import S4Vectors
#' @import SummarizedExperiment
#' @importFrom dplyr %>%
#' @importFrom methods as is setMethod setGeneric
#' @importFrom rlang .data
NULL
