#' Differential expression(TPM/FPKM/RPKM) analysis for genes, the input genes should be log-transformed, Notice that it maybe not the optimal methods for differential expression analysis.
#'
#' @param input_anno Input annotation, first column is group
#' @param input_expr the expression matrix
#' @param contrast contrast group
#'
#' @return data.frame of differential expression table
#' @export
#'
#' @examples
fun_deg_tpm_limma = function(input_anno,input_expr,contrast = c("High","Low")){
  suppressMessages(library(limma))
  # data preprocessing ------------------
  message_paste = paste0("The input conatrst group is: ",paste(contrast,collapse = "-"),"\n",collapse = "")
  message_paste = paste0(message_paste,
                         paste0("Found ",nrow(input_anno)," samples;\n"),
                         collapse = "")
  input_anno[[1]] = as.character(input_anno[[1]])
  # 1) remove NA
  input_anno = na.omit(input_anno)
  message_paste = paste0(message_paste,
                         paste0("Left ",nrow(input_anno)," samples, when removing NA;\n"),
                         collapse = "")

  # 2) remove group not in contrast
  input_anno = input_anno[input_anno[[1]] %in% contrast,,drop  = F]

  message_paste = paste0(message_paste,
                         paste0("Left ",nrow(input_anno)," samples, when filter with input_contrast;\n"),
                         collapse = "")

  # 3) get common sampls id
  common_samplesID = intersect(rownames(input_anno),colnames(input_expr))
  input_anno = input_anno[common_samplesID,,drop = F]
  message_paste = paste0(message_paste,
                         paste0("Left ",nrow(input_anno)," samples, when filter with expression data;"),
                         collapse = "")
  message(message_paste)
  # get expression data
  input_expr = input_expr[,common_samplesID]

  # deg step -----------------
  # 1) rename annotation
  input_anno[[1]] = ifelse(input_anno[[1]] == contrast[1],"contrast_Group1", input_anno[[1]])
  input_anno[[1]] = ifelse(input_anno[[1]] == contrast[2],"contrast_Group2", input_anno[[1]])

  # 2) make design, contrast
  group = input_anno[[1]]
  design = model.matrix(~0+factor(input_anno[[1]]))
  rownames(design) = rownames(input_anno)
  colnames(design) = levels(factor(input_anno[[1]]))
  #print(dim(design))

  contrast.matrix = limma::makeContrasts(contrast_Group1-contrast_Group2,levels = design)
  fit <- limma::lmFit(input_expr,design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  tempOutput = limma::topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput)
  if(!"ID" %in% colnames(nrDEG)) nrDEG$ID = rownames(nrDEG)
  return(nrDEG)
}


#' Deg cutoff function
#'
#' @param input_df  data
#' @param pval pval
#' @param pval_cutoff pval_cutoff
#' @param logFC logFC
#' @param logFC_cutoff logFC_cutoff
#'
#' @return add p_sig column into input_df
#' @export
#'
#' @examples
fun_deg_cutoff = function(input_df,
                          pval = "adj.P.Val",
                          pval_cutoff=0.01,
                          logFC = "logFC",
                          logFC_cutoff= 1){
  input_df$p_sig = dplyr::case_when(
    input_df[[pval]] > pval_cutoff ~"N.S.",
    abs(input_df[[logFC]]) < logFC_cutoff ~"N.S.",
    input_df[[logFC]] > abs(logFC_cutoff)  ~"Up",
    input_df[[logFC]] <  -abs(logFC_cutoff) ~"Down",
    TRUE~NA_character_
  )
  input_df
}
