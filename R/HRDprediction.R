#' @title HRD Prediction
#' @description a robust HRD predictor "HRDCNA" (Homologous Recombination
#'     Deficiency Prediction by Copy Number Feature) based on CNA features.
#'
#' @param data From function signimercopy
#'
#' @return HRDCNA probability scores of the samples.
#' @export
#'
#'
#' @examples
#' score_swgs <- HRDprediction(data = nmfcn_swgs)
#' score_wgs <- HRDprediction(data = nmfcn_wgs)
#' score_snp <- HRDprediction(data = nmfcn_snp)
#'
HRDprediction <- function(data){
  library(sigminer)
  library(gbm)
  data("modeldata")
  HRDCNAScore=predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
  samplescore <- as.data.frame(HRDCNAScore)
  samplescore$sample <- data$sample
  return(HRDCNAScore=samplescore)
}

