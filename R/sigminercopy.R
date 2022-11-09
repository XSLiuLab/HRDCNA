#' @title Calling Copy Number Features
#'
#' @description Detail information described in Sigminer.
#'
#' @param data Absolute copy number profile
#' @param genome_build Genome build version, should be 'hg19', 'hg38', 'mm9' or 'mm10'.
#'
#' @return The matrix of copy number features
#' @export
#'
#' @examples
#' cn_swgs <- readRDS("./data/test/cn_60_wgs10x.rds")
#' nmfcn_swgs <- sigminercopy(data = cn_swgs, "hg19")
#'
#' cn_wgs <- readRDS("./data/test/cn_60_wgs.rds")
#' nmfcn_wgs <- sigminercopy(data = cn_wgs, "hg19")
#'
#' cn_snp <- readRDS("./data/test/cn_60_snp.rds")
#' nmfcn_snp <- sigminercopy(data = cn_snp, "hg19")
#'
sigminercopy <-function(data, genome_build = c("hg19", "hg38", "mm10", "mm9")){
  genome_build <- match.arg(genome_build)
  copy <- sigminer::read_copynumber(data,
                                    seg_cols = c("chromosome", "start", "end", "segVal"),
                                    genome_build, complement = FALSE, verbose = TRUE)


  cnall_tally_W <- sigminer::sig_tally(copy, method = "W")
  nmf <- cnall_tally_W$nmf_matrix
  nmf <- as.data.frame(nmf)
  nmf$sample <- rownames(nmf)
  rownames(nmf) <- NULL
  return(nmf=nmf)
}
