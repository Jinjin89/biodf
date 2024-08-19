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
    expr_tmp = expr_tmp[!is.na(expr_tmp[[1]]),]
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


#' read methylation idat from dir, return the beta metarix
#'
#' @param input_dir input_dir contains the idat files(Grn,Red)
#' @param outfile the beta matrix file position
#' @param sample_anno_file if supplied write annotation file
#' @param input_df data.frame contains sample name, and idat file name
#' @param input_sample sample_colmn
#' @param input_filename idata filename
#' @param input_sample_regex extract sample form idat filename, if input_df is NULL
#' @param idat_suffix_regex extract idat form input_dir, if input_df is NULL
#'
#' @return data.frame
#' @export
#'

fun_read_meth_idat = function(
    input_dir,
    outfile,
    sample_anno_file = "sample_idat_file.csv",
    input_df =NULL,
    input_sample = "sample",
    input_filename = "filename",
    input_sample_regex = "GSM\\d+",
    idat_suffix_regex = "_Grn\\.idat$|_Grn\\.idat\\.gz$|_Red\\.idat$|_Red\\.idat\\.gz$"
){
  suppressMessages(require(minfi))

  # 1) check file exists, if not build it
  if(!file.exists(outfile)){
    # 1) input meta data not found, try to generate it
    if(is.null(input_df)){
      message("Gussing the results from input directory")
      all_files_found = list.files(input_dir,full.names = F)
      all_files_found_index = stringr::str_detect(
        all_files_found,
        idat_suffix_regex)
      # filter idat file
      all_files_found = all_files_found[all_files_found_index]
      all_files_found = stringr::str_remove_all(all_files_found,idat_suffix_regex)
      all_files_found = unique(all_files_found)
      # get sample names
      samples_names =stringr::str_extract(
        all_files_found,
        input_sample_regex)
      # get annotation
      sample_anno = data.frame(
        row.names = samples_names,
        sample = samples_names,
        filename =all_files_found)

      # write meta data into outfile
      if(!is.null(sample_anno_file)){
        write.csv(sample_anno,sample_anno_file,quote = F)
      }

    }else{
      sample_anno = data.frame(
        row.names = samples_names,
        sample = input_df[[input_sample]],
        filename = input_df[[input_filename]]
      )
    }

    # 2) read idat
    idat_list =
      purrr::map(seq_along(sample_anno$sample),function(i){
        # 1) get meta data
        sample_tmp = sample_anno$sample[i]
        filename_tmp = sample_anno$filename[i]
        message(paste0("reading from ", filename_tmp))
        idata_read = minfi::read.metharray(
          basenames = paste0(input_dir,"/",filename_tmp))
        idata_read =  minfi::getBeta(idata_read)
        colnames(idata_read) = sample_tmp
        probe_names = rownames(idata_read)
        idata_read = data.table::as.data.table(idata_read)
        idata_read$probe = probe_names
        idata_read
      }) %>%
      purrr::reduce(data.table::merge.data.table,by = "probe",all=T) %>%
      as.data.frame() %>%
      set_rownames(.$probe) %>%
      select(-probe)
    write.table(idat_list,file = outfile,sep = "\t",quote = F)
  }
  message("You may want to convert the output into matrix(check if it is numeric)")
  read.table(outfile,sep = "\t",row.names = 1)
}


#' read file into dataframe
#'
#' @param input_file the file postion of large text
#' @param index_col which column is index
#' @param ... params pass to fread
#'
#' @return data.frame objc
#' @export
#'

fun_read_df <- \(input_file,index_col ='sample',...){
  df = fread(input_file,...)
  if(len(index_col) == 0){
    return(df)
  }else if(len(index_col) == 1){
    df %>%
      as.data.frame() %>%
      magrittr::set_rownames(.[[index_col]]) %>%
      return
  }else{
    stop('index column must be length 1 of NULL')
  }

}
