
fun_h52df <- function(input_h5_file,input_df_name,row_name = "_index",categroyCode=c("categories","codes")){

  ft_data_parse <- function(data_parse){
    if(is.numeric(data_parse)){
      return(data_parse)
    }else if(is.list(data_parse)){
      label = data_parse[[categroyCode[1]]]
      value = data_parse[[categroyCode[2]]] + 1
      return(label[value])
    }else{
      return(data_parse)
    }
  }
  # 1）read df list
  df_list = rhdf5::h5read(input_h5_file,input_df_name)
  df_name = names(df_list)
  message(paste0("found: ",length(df_name)))
  print(length(df_list))

  # 2）
  for(i in seq_along(df_name)){
    current_name = df_name[i]
    df_list[[current_name]] = ft_data_parse(df_list[[current_name]])
  }
  df = as.data.frame(df_list)
  if(length(row_name) == 1 && row_name %in% df_name){
    message("setting rownames")
    rownames(df) = df_list[[row_name]]
  }
  return(df)
}
