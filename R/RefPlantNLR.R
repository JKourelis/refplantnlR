#' RefPlantNLR domains
#'
#' A dataset containing the domain boundaries found in the RefPlantNLR dataset.
#'
#' @format A data frame with 14720 rows and 13 variables:
#' \describe{
#'   \item{Status}{NLR, NLR-associated}
#'   \item{seqname}{sequence name}
#'   \item{Subclass}{NLR subclass: TIR-NLR, CCR-NLR, CCG10-NLR, CC-NLR}
#'   \item{Domain_simplified}{Domain architecture}
#'   \item{Species}{Scientific species name from which the NLR was cloned}
#'   \item{Order}{Phylogenetic order to which the host belongs}
#'   \item{Pathogen}{Pathogen against which NLR confers immunity}
#'   \item{Genus}{Genus}
#'   \item{Representative}{Does the NLR belong to the representative NLR set at 90% identity treshold: YES, NO}
#'   \item{Type}{CHAIN, MOTIF, REGION, DOMAIN}
#'   \item{description}{Domain name}
#'   \item{start}{start position}
#'   \item{end}{end position}
#'   ...
#' }
#' @source \url{https://zenodo.org/record/3936022/}
"RefPlantNLR"

