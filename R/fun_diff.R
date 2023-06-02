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





#' Differential Methylation Region(DMR) using ChAMP
#'
#' @param input_df pd
#' @param input_matrix beta
#' @param outfile outfile
#' @param input_group group
#' @param input_sample sample column
#' @param input_batch input_batch
#' @param filterXY T or F
#' @param arraytype arraytype
#' @param cores cores
#' @param method method
#' @param adjPVal adjPVal
#' @param compare.group compare.group
#'
#' @return list with differntial expressed genes
#' @export
#'
#' @examples
fun_dmr_champ <-  function(input_df,input_matrix,

                           # save the DMR results
                           outfile,

                           # clinical features params
                           input_group="Group",
                           input_sample = "sample",
                           input_batch = NULL,

                           # DMR params
                           filterXY = F,
                           arraytype="450K",
                           cores = 2,
                           method ="BMIQ",
                           adjPVal = 0.05,
                           compare.group = c("normal","cancer")
){

  # if not found the results
  if(!file.exists(outfile)){
    message(paste0(outfile," not found, performing DMR using ChAMP!"))
    # 1) get new clinical data
    samples_found = intersect(input_df[[input_sample]],colnames(input_matrix))
    input_df = input_df %>%
      dplyr::filter(.[[input_sample]] %in% samples_found) %>%
      dplyr::filter(.[[input_group]] %in% compare.group)
    print(table(input_df[[input_group]]))
    stopifnot("The comparison group should be equal to 2"=length(unique(input_df[[input_group]])) == 2)

    input_df[[input_group]] = factor(input_df[[input_group]],
                                     levels = compare.group)

    input_df =
      dplyr::arrange(input_df,!!as.name(input_group))


    # 2) filter_matix
    samples_found = input_df[[input_sample]]
    input_matrix = input_matrix[,samples_found]

    # 3) run
    # 3.1) load data
    myLoad =
      ChAMP::champ.filter(
        beta = input_matrix,
        pd = input_df,
        filterXY = filterXY,
        arraytype = arraytype
      )

    # 3.2) remove empty
    myLoad = ChAMP::champ.impute(
      beta = myLoad$beta
      ,pd =myLoad$pd
      ,ProbeCutoff = 0.2
      ,SampleCutoff = 1
    )

    # 3.3) Norm
    myNorm <- ChAMP::champ.norm(
      beta=myLoad$beta,
      arraytype=arraytype,
      cores=cores,
      method = method)

    # 04) DMP
    myDMP <- ChAMP::champ.DMP(
      beta = myNorm,
      adjPVal = adjPVal,
      pheno=myLoad$pd[[input_group]],
      compare.group = compare.group,
      arraytype = arraytype)

    # results
    res =
      list(data = myLoad,
           myNorm = myNorm,
           DMP = myDMP)
    res %>%
      saveRDS(outfile)
  }
  invisible(readRDS(outfile))
}
