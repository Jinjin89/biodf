#' read data into list
#'
#' @param input_clin clinical file postion
#' @param input_expr expression file postion
#' @param input_maf maf file position
#'
#' @return a list contain data
#' @export
#'
fun_read_cbio = function(input_clin=NULL,input_expr=NULL,input_maf =NULL){
  res_list = list()
  # 1) get clinical data
  if(!is.null(input_clin)){
    clin = data.table::fread(input_clin)
    clin$sample = clin$`Sample ID`
    clin = janitor::clean_names(clin)
    clin = as.data.frame(clin)
    rownames(clin) = clin$sample
    res_list$clin =clin
  }


  # 2) get expression data
  if(!is.null(input_expr)){
    expr_tmp = data.table::fread(input_expr)
    expr_tmp = as.data.frame(expr_tmp)
    expr_tmp = expr_tmp[!duplicated(expr_tmp[[1]]),]
    rownames(expr_tmp) = expr_tmp[[1]]
    expr_tmp = expr_tmp[,-c(1,2)]
    expr_tmp = as.matrix(expr_tmp)
    if(sum(is.na(expr_tmp)) > 0 ){
      expr_tmp = impute::impute.knn(expr_tmp,rng.seed = 123)
      expr_tmp = expr_tmp$data
    }
    res_list$expr = expr_tmp
  }


  # 3) get mutation data

  if(!is.null(input_maf)){
    res_list$maf = data.table::fread(input_maf)
  }
  res_list
}



#' read single cell data downloaded from
#'
#' @param input_dir input directory
#' @param output_file the saved filed
#'
#' @return list of single cell
#' @export
#'
fun_read_tisch = function(input_dir,output_file)
{

  if(!file.exists(output_file)){

    # 1) get files
    all_files = list.files(input_dir,full.names = T)
    h5_file = all_files[str_detect(all_files,"expression.h5$")][1]
    meta_file = all_files[str_detect(all_files,"Metainfo_table.tsv$")][1]

    # 2) read sc data
    h5_read = Seurat::Read10X_h5(h5_file)
    sc_meta = data.table::fread(meta_file)
    sc_meta = as.data.frame(sc_meta)
    rownames(sc_meta) = sc_meta[[1]]
    sc_meta = dplyr::select(sc_meta,-1)

    # 3) process h5_file
    h5_read = MAESTRO::RNARunSeurat(h5_read)

    # 4) merge data
    meta_info = names(sc_meta)
    sc_rownames = rownames(h5_read$RNA@meta.data)
    for(i in seq_along(meta_info)){
      meta_tmp = meta_info[i]
      h5_read$RNA@meta.data[[meta_tmp]] = sc_meta[sc_rownames,meta_tmp]
    }

    # 5) rename the cell_major and cell minor
    h5_read$RNA$cell_major = h5_read$RNA$`Celltype (major-lineage)`
    h5_read$RNA$cell_minor = h5_read$RNA$`Celltype (minor-lineage)`
    h5_read$RNA$malignancy = h5_read$RNA$`Celltype (malignancy)`

    file.remove("MAESTRO.scRNA.Seurat_cluster.png")
    file.remove("MAESTRO.scRNA.Seurat_DiffGenes.tsv")
    file.remove("MAESTRO.scRNA.Seurat_PCElbowPlot.png")

    # 5) return
      saveRDS(h5_read,output_file)
  }
  return(readRDS(output_file))
}




#' Change the data.frame format into maf format
#'
#' @param input_df data.frame or file position of maf file
#' @param ... other fields for maf column
#' @param Hugo_Symbol Hugo_Symbol
#' @param Chromosome Chromosome
#' @param Start_Position Start_Position
#' @param End_Position End_Position
#' @param Reference_Allele Reference_Allele
#' @param Tumor_Seq_Allele2 Tumor_Seq_Allele2
#' @param Variant_Classification Variant_Classification
#' @param Variant_Type Variant_Type
#' @param Tumor_Sample_Barcode Tumor_Sample_Barcode
#'
#' @return a maf-like format
#' @export
#'
fun_read_df2maf = function(
    input_df,...,
    Hugo_Symbol = "Hugo_Symbol",
    Chromosome = "Chromosome",
    Start_Position = "Start_Position",
    End_Position = "End_Position",
    Reference_Allele = "Reference_Allele",
    Tumor_Seq_Allele2 ="Tumor_Seq_Allele2",
    Variant_Classification = "Variant_Type",
    Variant_Type ="Variant_Type",
    Tumor_Sample_Barcode = "Tumor_Sample_Barcode"
){
  #
  message("Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.")

  if(is.character(input_df)){
    input_df = data.table::fread(input_df)
  }
  # return maf format
  reutur_df =
    data.table(
      Hugo_Symbol = input_df[[Hugo_Symbol]],
      Chromosome = input_df[[Chromosome]],
      Start_Position = input_df[[Start_Position]],
      End_Position = input_df[[End_Position]],
      Reference_Allele = input_df[[Reference_Allele]],
      Tumor_Seq_Allele2 =input_df[[Tumor_Seq_Allele2]],
      Variant_Type =input_df[[Variant_Type]],
      Variant_Classification = input_df[[Variant_Classification]],
      Tumor_Sample_Barcode = input_df[[Tumor_Sample_Barcode]])

  input_args = list(...)
  if(length(input_args) != 0){
    message('extra args:')
    arg_names =names(input_args)
    for(i in seq_along(arg_names)){
      arg_name = arg_names[i]
      if(arg_name == ""){
        reutur_df[[input_args[[i]]]] = input_df[[input_args[[i]]]]
      }else{
        reutur_df[[arg_name]] = input_df[[input_args[[arg_name]]]]
      }
    }
  }
  return(reutur_df)
}
