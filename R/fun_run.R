#' Run gistic2
#'
#' @param input_df meta-data, used to extract data
#' @param input_cnv_seg_file cnv segment file
#' @param output_dir The results folder
#' @param input_group_name Group
#' @param input_sample_name sample
#' @param loading whether reading the results using maftools::readGistic
#' @param gistic2_command gistic2 commaon
#' @param gistic2_refgene gistic2 refegene
#' @param LD_LIBRARY_PATH gistic2 LD_LIBRARY_PATH
#'
#' @return
#' @export
#'
#' @examples
fun_run_gistic2 <-function(
    input_df,
    input_cnv_seg_file,
    output_dir,
    input_group_name = 'Group',
    input_sample_name = 'sample',
    loading=T,

    gistic2_command = '/home/r/miniconda3/envs/gistic2/bin/gistic2',
    gistic2_refgene = '/home/r/miniconda3/envs/gistic2/share/gistic2-2.0.23-0/refgenefiles/hg19.mat',
    LD_LIBRARY_PATH="/home/r/miniconda3/envs/gistic2/share/mcr-8.3-0/v83/runtime/glnxa64:/home/r/miniconda3/envs/gistic2/share/mcr-8.3-0/v83/bin/glnxa64:/home/r/miniconda3/envs/gistic2/share/mcr-8.3-0/v83/sys/os/glnxa64:"
){
  require(biodf)
  require(maftools)
  message('* install the gistic2 with conda, and find the refergene accomplied by package in share folder, then get the LD_LIBRARY_PATH')

  # 1) check the configs
  stopifnot('output_dir must exist' = dir.exists(output_dir))
  stopifnot('input_group data could not be numeri' = !is.numeric(input_df[[input_group_name]]))

  # 2) getting data
  df_filter <- input_df %>% dplyr::filter(!is.na(!!as.name(input_group_name)))
  all_group = df_filter[[input_group_name]]
  unique_groups = unique(all_group)
  seg = fread(input_cnv_seg_file)

  gis_return_list = list()


  # 3)
  for(each_group in unique_groups){
    # 1) data_init
    group_dir = file.path(output_dir,each_group)
    seg_file = file.path(group_dir,'seg.txt')
    done_file = file.path(group_dir,'done')
    res_dir = file.path(group_dir,'results')

    if(!file.exists(done_file)){
      message('>>> results completion file not found, run...')
      # 1) getting the cnv data
      samples_for_group = df_filter %>%
        dplyr::filter(.[[input_group_name]] == each_group) %>%
        dplyr::pull(!!as.name(input_sample_name))
      seg_group <- seg %>% dplyr::filter(.[[1]] %in% samples_for_group)
      if(!dir.exists(group_dir)){
        dir.create(group_dir)
      }

      if(!dir.exists(res_dir)){
        dir.create(res_dir)
      }

      # 2) write the seg data into disk
      seg_group %>% data.table::fwrite(seg_file,sep = '\t',quote = F)

      # 3) running the gistic2 step
      gistic2_run = paste(
        paste0('LD_LIBRARY_PATH=',LD_LIBRARY_PATH),';',
        gistic2_command,'-b', res_dir,
        '-seg', seg_file,
        '-refgene',  gistic2_refgene
      )
      cat('**********************************\n')
      cat(gistic2_run,'\n')
      cat('**********************************\n')
      run_code = system(gistic2_run)
      if(run_code == 0){
        system(paste('touch', done_file,sep = ' '))
      }
      # 4) done
    }else{
      message('>>> results completion file found, skip!')
    }

    # loading results
    if(loading){
      message('>>> loading from: ',res_dir)
      gisticAllLesionsFile = list.files(res_dir,pattern = 'all_lesions\\.conf',full.names = T)
      gisticAmpGenesFile =  list.files(res_dir,pattern = 'amp_genes\\.conf',full.names = T)
      gisticDelGenesFile =  list.files(res_dir,pattern = 'del_genes\\.conf',full.names = T)
      gisticScoresFile =  list.files(res_dir,pattern = 'scores\\.gistic',full.names = T)

      gis_return_list[[each_group]] = maftools::readGistic(
        gisticAllLesionsFile = gisticAllLesionsFile,
        gisticAmpGenesFile = gisticAmpGenesFile,
        gisticDelGenesFile = gisticDelGenesFile,
        gisticScoresFile = gisticScoresFile)
    }
  }
  # loading gistic result

  if(loading){
    return(gis_return_list)
  }
}
