#' RefSeq NLR domains
#'
#' A dataset containing the domain boundaries found in the NLRs extracted from the plant NCBI RefSeq proteomes using NLRtracker.
#'
#' @format A data frame with 815304 rows and 14 variables:
#' \describe{
#'   \item{RefSeq}{RefSeq assembly accession}
#'   \item{seqname}{Protein_ID}
#'   \item{Type}{CHAIN, MOTIF, REGION, DOMAIN}
#'   \item{description}{Domain name}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   \item{Species}{Species scientific name}
#'   \item{Order}{Phylogenetic order to which the organism belongs}
#'   \item{coded_by}{Transcript_ID}
#'   \item{Locs}{Locus_ID}
#'   \item{Status}{NLR, degenrate NLR, likely non-plant NLR, NLR-associated, MLKL}
#'   \item{Subclass_putative}{Putative NLR subclass: TIR-NLR, CCR-NLR, CCG10-NLR, CC-NLR}
#'   \item{Domain}{Domain architecture}
#'   \item{Domain_simplified}{Domain architecture simplified}
#'   ...
#' }
#' @source \url{https://zenodo.org/record/3936022/}
"RefSeq_NLR"
