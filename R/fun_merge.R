#' Merge varibles of matrix into dataframe
#'
#' @param input_df input_df
#' @param input_matrix input_matrix
#' @param input_variables variables to merge
#' @param key
#'
#' @return data.frame
#' @export
#'
fun_merge_matrix = function(input_df,input_matrix,input_variables,output_variables_names=input_variables,key = "sample"){
  # 1) get variables
  input_var_n = length(unique(input_variables))

  # 2) variables found
  input_variables = intersect(input_variables,
                              rownames(input_matrix))
  input_var_n_found = length(unique(input_variables))
  message(paste0("found ",input_var_n_found," from ",input_var_n))

  output_variables_names = fun_utils_broadcast(input_variables,output_variables_names)

  # 3) convert into matrix
  input_matrix = as.matrix(input_matrix)
  matrix_names = colnames(input_matrix)

  # 4) merge into df
  for(each_index in seq_along(input_variables)){
    each_var = input_variables[each_index]
    each_output = output_variables_names[each_index]
    df_tmp =
      data.frame(row.names = matrix_names,
                 var_tmp = input_matrix[each_var,matrix_names])
    input_df[[each_output]] = df_tmp[input_df[[key]],"var_tmp"]
  }
  return(input_df)
}



#' Merge the mutation status of gene into clinical data.frame
#'
#' @param input_df input_data.frame
#' @param input_maf the maf format data
#' @param input_variables genes to retrive the mutation data or list
#' @param key the key(sample,Tumor_Sample_Barcode)
#'
#' @return a data.frame with genes mutation status
#' @export
#'
#' @examples
fun_merge_maf = function(input_df,input_maf,input_variables,key = "sample"){
  # 1) extract data
  if(class(input_maf)[1] == "MAF"){
    input_maf = input_maf@data
  }else{
    message("MAF class data is preferred, as it filtered out some variants for us!")
  }

  # 2) get all the samples in maf
  all_samples_found = as.character(unique(input_maf$Tumor_Sample_Barcode))
  mut_df = data.frame(row.names = all_samples_found)

  # 3) get all the mutation status
  for(each_gene in input_variables){
    # 1) get all the mutated sample for this  gene
    mutated_samples_for_gene = as.character(input_maf$Tumor_Sample_Barcode[input_maf$Hugo_Symbol %in% each_gene])
    mutated_samples_for_gene = unique(mutated_samples_for_gene)

    # 2)
    each_gene = paste0(each_gene,collapse ="&")
    mut_df[[each_gene]] = "wt"
    mut_df[mutated_samples_for_gene,each_gene] = "mt"
    input_df[[paste0(each_gene,"_status")]] = mut_df[input_df[[key]],each_gene]

  }
  input_df
}



#' Merge the varaibles found in data.frame or file position, merge it into
#'
#' @param input_df input data.frame
#' @param input_df_merge the data.frame to merge
#' @param input_variables variables to merge
#' @param input_regex regex of variables to merge
#' @param key the input_df key, used to merge
#' @param df_merge_key the input_df_merge key, used to merge
#'
#' @return data.frame with merged varaibles
#' @export
#'
#' @examples
fun_merge_df = function(
    input_df,
    input_df_merge,
    input_variables = NULL,
    input_regex = "xcell|cibersort|mcpcounter",
    key = "sample",
    df_merge_key = "sample"){
  # 1) check data
  if(is.character(input_df_merge)){
    if(file.exists(input_df_merge)){
      input_df_merge = data.table::fread(input_df_merge)
    }else{
      stop("input_df_merge should be either data.frame or the file path of data")
    }
  }

  # 2) parse the merged data
  input_final = c()
  if(length(input_variables)>0){
    input_final = c(input_final,
                    intersect(input_variables,colnames(input_df_merge)))
  }

  if(length(input_regex) == 1){
    input_regex_found_var = colnames(input_df_merge)
    input_regex_found_var = input_regex_found_var[stringr::str_detect(input_regex_found_var,
                                                                      stringr::regex(input_regex,ignore_case = T))]
    input_final = unique(
      c(input_final,
        input_regex_found_var))
  }

  # 3) check the input_df_merge data
  stopifnot("df_merge_key not found in input_df_merge" = df_merge_key %in% colnames(input_df_merge))
  stopifnot("The merged variables not found in input_df_merge" = length(input_final)>0)

  # 4) get new data
  df_new = dplyr::select(input_df_merge,dplyr::all_of(input_final))
  df_new = as.data.frame(df_new)
  rownames(df_new) = input_df_merge[[df_merge_key]]

  # 5)
  for(each_var in input_final){
    input_df[[each_var]] = df_new[input_df[[key]],each_var]

  }
  input_df
}
