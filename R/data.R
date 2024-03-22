#sig1
if(F){
  fun_sig <- list()
  fun_sig$data = list()

  fun_sig$data$pathological_response_signature_pmid34143979 <- list()
  fun_sig$data$pathological_response_signature_pmid34143979$PD1 = 'PDCD1'
  fun_sig$data$pathological_response_signature_pmid34143979$PDL1 = 'CD274'
  fun_sig$data$pathological_response_signature_pmid34143979$CD68 <- c('CD68')
  fun_sig$data$pathological_response_signature_pmid34143979$T_cells_signature =
    c('CD3D','CD3E','CD3G', 'CD6','SH2D1A','TRAT1') %>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$B_cells_signature =
    c('BLK', 'CD19','FCRL2','KIAA0125','MS4A1', 'PNOC','SPIB', 'TCL1A','TNFRSF17')%>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$Dendric_cells_signature <-
    c('CCL13','CD209','HSD11B1')%>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$Mast_cells_signature <-
    c('CPA3', 'HDC','MS4A2','TPSAB1','TPSB2')%>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$TIS_signature <-
    c('TIGIT', 'CD27','CD8A','PDCD1LG2','CXCR6', 'LAG3','CD274','CMKLR1','NKG7',
      'CCL5','PSMB10', 'IDO1','PPBP', 'HLA-DQA1', 'CD276','STAT1', 'HLA-DRB1', 'HLA-E')%>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$STAT1_signature <-
    c('TAP1', 'GBP1','IFIH1', 'PSMB9','CXCL9', 'IRF1','CXCL11','CXCL10')%>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$Mitotic_signature <-
    c('PLK1', 'CDK1','BUB1B', 'NEK2','TTK', 'MELK','PLK4', 'CHEK1','AURKA',
      'AURKB', 'BUB1','PBK' )%>%
    fun_gene2symbol() %>% pull(symbol)

  fun_sig$data$pathological_response_signature_pmid34143979$ESR1_PGR_ave <-
    c('ESR1','PGR')%>%
    fun_gene2symbol() %>% pull(symbol)

  fun_sig$data$pathological_response_signature_pmid34143979$TAMsurr_signature <-
    c('CXCL10','CXCL11', 'CCL8','LAMP3')%>%
    fun_gene2symbol() %>% pull(symbol)
  fun_sig$data$pathological_response_signature_pmid34143979$TcClassII_signature <-
    c('CD2','CD3G', 'CD8A','IFNG', 'TNF','GZMB', 'GZMH','PRF1', 'ZAP70',
      'HLA-DMA', 'HLA-DOA', 'HLA-DOB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1',
      'HLA-DQB2','HLA-DRA', 'HLA-DRB1', 'HLA-DRB2', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
      'HLA-DRB6', 'CIITA','CD74')%>%
    fun_gene2symbol() %>%
    pull(symbol)

  fun_sig$data$pathological_response_signature_pmid34143979$revisedPARPi7_signature <-
    c('BRCA1','CHEK1','MAPKAPK2','XRCC4','RAD17','POLB','CIRBP')%>%
    fun_gene2symbol() %>%
    pull(symbol)

  fun_sig$data$pathological_response_signature_pmid34143979$parpi_norm_genes <-
    c( 'RPL24', 'ABI2', 'GGA1', 'E2F4', 'IPO8', 'CXXC1', 'RPS10')%>%
    fun_gene2symbol() %>%
    pull(symbol)

  fun_sig$fun_pathological_response_signature_pmid34143979 <- function(
    input_expr,
    input_sigs = fun_sig$data$pathological_response_signature_pmid34143979){
    # 1) the input expression should not be z
    message("The input expression should be log transformed if RNA-seq data supplied!")
    message("The input expression should not be z-normalized!")
    if(!is.matrix(input_expr)){
      warning("Input is not matrix, convert it into matrix")
      input_expr <- as.matrix(input_expr)
    }

    samples_found <- colnames(input_expr)
    genes_found <- rownames(input_expr)
    results_df = data.frame(row.names = samples_found,sample = samples_found)

    # 2) each_genes analysis
    for(eachgene in c("PD1","PDL1","CD68")){
      gene_symbol <-  input_sigs[[eachgene]]
      if(gene_symbol %in% genes_found) {
        results_df[[eachgene]] <- biodf::fun_norm_scale(input_expr[gene_symbol,samples_found])
      }else{
        warning(eachgene,": not found!")
      }
      rm(eachgene)
      rm(gene_symbol)
    }

    # 3)
    for(each_sig in c('T_cells_signature','B_cells_signature','Dendric_cells_signature',
                      'Mast_cells_signature','TIS_signature','STAT1_signature',
                      'Mitotic_signature','ESR1_PGR_ave')){
      # current genes
      current_genes <- input_sigs[[each_sig]]
      #print(current_genes)

      # not found genes
      not_found_genes <- dplyr::setdiff(current_genes,genes_found)
      genes_found_inter <- dplyr::intersect(current_genes,genes_found)

      if(biodf::len(genes_found_inter) == 0){
        warning("All genes not found for: ",each_sig," The genes is: ",paste0(not_found_genes,collapse = ","))
        rm(current_genes)
        rm(not_found_genes)
        next
      }else if(biodf::len(not_found_genes)>0){
        warning("Some genes not found: ",paste0(not_found_genes,collapse = ","))
      }else{
        message("All genes found for: ", each_sig)
      }
      # current matrix
      current_matrix <- input_expr[genes_found_inter,samples_found]
      mean_value <- apply(current_matrix,1,mean,na.rm=T)
      matrix_sub <- apply(current_matrix,2,function(x) x - mean_value)
      sample_mean <- apply(matrix_sub,2,mean,na.rm=T)
      results_df[[each_sig]] <- biodf::fun_norm_scale(sample_mean[samples_found])

      # remove the data
      rm(current_matrix)
      rm(current_genes)
      rm(not_found_genes)
      rm(mean_value)
      rm(sample_mean)
    }
    # 3)
    sig1 <- dplyr::intersect(genes_found,input_sigs$TAMsurr_signature)
    sig2 <- dplyr::intersect(genes_found,input_sigs$TcClassII_signature)
    if(biodf::len(sig1) >0 & biodf::len(sig2) > 0){
      # sig1
      current_matrix <- input_expr[sig1,samples_found]
      mean_value <- apply(current_matrix,1,mean,na.rm=T)
      matrix_sub <- apply(current_matrix,2,function(x) x - mean_value)
      sample_mean <- apply(matrix_sub,2,mean,na.rm=T)
      sig1 <- biodf::fun_norm_scale(sample_mean[samples_found])
      rm(current_matrix)
      rm(mean_value)
      rm(sample_mean)

      # sig2
      current_matrix <- input_expr[sig2,samples_found]
      mean_value <- apply(current_matrix,1,mean,na.rm=T)
      matrix_sub <- apply(current_matrix,2,function(x) x - mean_value)
      sample_mean <- apply(matrix_sub,2,mean,na.rm=T)
      sig2 <- biodf::fun_norm_scale(sample_mean[samples_found])

      rm(current_matrix)
      rm(mean_value)
      rm(sample_mean)
      data.frame(row.names = samples_found,
                 sig1 = sig1,
                 sig2 = sig2) %>%
        dplyr::mutate(adj =ifelse(.$sig1 > .$sig2,abs(.$sig2)+1,abs(.$sig1)+1)) %>%
        mutate(final = log2((.$sig1 + adj)/(.$sig2+adj))) %>%
        mutate(final = fun_norm_scale(final)) -> final_df

      results_df[,'TAMsurr_TcClassII_ratio_signature'] <- final_df[rownames(results_df),'final']

    }else{
      warning("TAMsurr, TcClassII gene set found to be empty at least in one datasets")
    }

    # parpi
    parpi_sig_found <- intersect(genes_found,input_sigs$revisedPARPi7_signature)
    parpi_nor_found <- intersect(genes_found,input_sigs$parpi_norm_genes)
    parpi_nor_found_n = biodf::len(parpi_nor_found)

    sig_values = data.frame(row.names = input_sigs$revisedPARPi7_signature) %>%
      mutate(w = c(-0.5320, 0.5806, 0.0713, -0.1396, -0.1976, -0.3937, -0.2335),
             b = c(-0.0153, -0.006, 0.0031, -0.0044, 0.0014, -0.0165, -0.0126))

    if(biodf::len(parpi_sig_found) > 0  & biodf::len(parpi_nor_found) >0){
      # get normalized genes
      normal_expr <- input_expr[parpi_nor_found,] %>% apply(2,function(x) exp(mean(log(x))))

      #normal_expr <- input_expr[parpi_nor_found,] %>% apply(2,function(x) mean(x))
      parpi_expr <- input_expr[parpi_sig_found,]
      parpi_expr_normalized <- apply(parpi_expr,1,function(x) x/normal_expr)
      rownames(parpi_expr_normalized) <- colnames(parpi_expr) # samples * genes

      parpi_expr_normalized %<>% log2() %>% apply(2,function(x) x - median(x,na.rm=T)) # samples * genes
      parpi_expr_normalized <- t(parpi_expr_normalized) # genes * samples
      for(each_var in rownames(parpi_expr_normalized)){
        parpi_expr_normalized[each_var,] <-
          parpi_expr_normalized[each_var,] - sig_values[each_var,'b']
      }
      parpi_expr_normalized <- t(parpi_expr_normalized) # samples *genes

      res = parpi_expr_normalized %*% as.matrix(sig_values[colnames(parpi_expr_normalized),"w",drop=F]) %>%
        as.data.frame() %>%
        magrittr::set_colnames("final") %>%
        dplyr::mutate(final = final/sd(.$final,na.rm=T))
      #dplyr::mutate(final = fun_norm_scale(.$final))
      results_df[,"PARPi7_signature"] <- res[rownames(results_df),'final']

    }else{
      warning("parpi signature genes not found!")
    }
    return(results_df)
  }

  #emt
  fun_sig$fun_sig_emt <- function(input_mat,outfile){
      go_db <- biodata::db_go
      ecm_pos_neg = list()
      ecm_pos_neg$pos = go_db$gene[go_db$term == "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"]
      ecm_pos_neg$neg = go_db$gene[go_db$term == "GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"]

      ecm_df = fun_sig_ssGSGA(
        input_mat =input_mat,
        input_genes_list = ecm_pos_neg,
        outfile = outfile)

      ecm_df$pos46 = ecm_df$pos * 46
      ecm_df$neg28 = ecm_df$neg * 28

      ecm_df$score = ecm_df$pos46 - ecm_df$neg28
      return(ecm_df)

    }

  use_data(fun_sig,overwrite = T)

}

# sig2
if(F){
  fun_sig$data$consensus_signature_pmid34965943 <-
    fread('~/data/tmp/consensus_sig.txt') %>%
    as.list() %>%
    purrr::map(\(x) {
      x = x[x!= '']
      x = fun_gene2symbol(x) %>%
        pull(symbol) %>%
        unique
      x
      }
      )

  fun_sig$fun_consensus_signature_pmid34965943 <- function(
    input_mat,
    input_sig_list = biodf::fun_sig$data$consensus_signature_pmid34965943){
    message("The input_mat should be in (log scale)!")
    all_sig_names = names(input_sig_list)
    purrr::map(all_sig_names,\(each_sig){
      # common_names
      common_names = intersect(input_sig_list[[each_sig]],rownames(input_mat))
      current_mat = input_mat[common_names,]
      current_sig = apply(current_mat,2,mean,na.rm=T)
      data.frame(sig = current_sig) %>%
        set_colnames(each_sig)
    }) %>%
      purrr::reduce(cbind)
  }
  #test_mat = readRDS("~/data/project/db/tcga/rds/TCGA_ACC.RDS")
  fun_sig$fun_consensus_signature_pmid34965943(input_mat = test_mat$expr)
  fun_sig$data$TcellInf <- fread('~/data/tmp/tcellsig.txt') %>%
    as.data.frame() %>%
    set_rownames(.[[1]]) %>%
    select(-1) %>%
    set_rownames(gene_alias2symbol[rownames(.),"symbol"])

  fun_sig$data$TcellInf_houseKeeping <- c(
    'STK11IP','ZBTB34','TBC1D10B','OAZ1','POLR2A','G6PD','ABCF1',
    'C14orf102','UBB','TBP','SDHA') %>%
    {gene_alias2symbol[.,"symbol"]}
  fun_sig$fun_sig_TcellInf <- function(
    input_mat,
    normalized = T,
    model_coef = biodf::fun_sig$data$TcellInf,
    housekeeping_genes = biodf::fun_sig$data$TcellInf_houseKeeping,
    results_scale = T){
    # 1) get current matrix
    common_genes_found <- intersect(rownames(input_mat),rownames(model_coef))
    current_matrix <- t(input_mat[common_genes_found,]) # sample * genes
    # 2) whether normalized
    if(normalized){
      message('normalized towards housekeeping genes!')
      normalized_genes =  intersect(rownames(input_mat),housekeeping_genes)
      normalized_matrix = colMeans(input_mat[normalized_genes,],na.rm = T)
      current_matrix = current_matrix - normalized_matrix
    }
    # get results
    model_coef = model_coef[common_genes_found,,drop=F]
    res = current_matrix %*% as.matrix(model_coef) %>%
      as.data.frame() %>%
      set_colnames("TcellInf")
    # results scale
    if(results_scale) res$TcellInf = fun_norm_scale(res$TcellInf)

    res
  }
  use_data(fun_sig,overwrite = T)
}
