#' export current enviroment
#'
#' @param outfile the csv files for storing packages
#' @param lib.loc
#'
#' @return data.frame
#' @export
#'
fun_utils_export_enviroment <- \(outfile,lib.loc = NULL){
  parse_package <- \(pkg_name,lib.loc = NULL){
    # 1) get description
    pkg_description <- packageDescription(pkg_name)
    stopifnot(!is.na(pkg_description))

    # 2) collection info
    repos <- packageDescription(pkg_name,lib.loc = lib.loc,fields = 'Repository')
    remote_type <- packageDescription(pkg_name,lib.loc = lib.loc,fields = 'RemoteType')
    if(is.na(repos) ){ # no respos found
      if(is.na(remote_type)){ # no remote tpe found
        repos ='others'
      }else if(remote_type == 'github'){ # check github repos
        repos <- 'github'
      }else{
        repos = 'others'
      }
    }

    # 3) get version
    version <- packageDescription(pkg_name,lib.loc = lib.loc,fields = 'Version')

    # 3) pase git_url
    if(repos == 'github'){
      RemoteUsername <- packageDescription(pkg_name,lib.loc = lib.loc,fields = 'RemoteUsername')
      RemoteRepo <- packageDescription(pkg_name,lib.loc = lib.loc,fields = 'RemoteRepo')
      RemoteRef <- packageDescription(pkg_name,lib.loc = lib.loc,fields = 'RemoteRef')
      git_id <- file.path(RemoteUsername,RemoteRepo)

    }else{
      git_id <- ''
      RemoteRef <- ''
    }

    # df, names, source, install methods
    data.frame(
      pkg_name = pkg_name,
      repos = repos,
      version = version,
      github = git_id,
      github_ref = RemoteRef
    )
  }

  # 1) get installed packages
  pkg_installed <- installed.packages(lib.loc = lib.loc) |>
    as.data.frame()
  all_pkgs <- pkg_installed$Package
  # 2)
  return_df <- parse_package(all_pkgs[1],lib.loc = lib.loc)
  for(i in 2:length(all_pkgs)){
    current_pkg <- all_pkgs[i]
    return_df <- rbind(return_df,parse_package(current_pkg,lib.loc = lib.loc))
  }
  # saving files
  message('save files')
  return_df |> write.csv(outfile,row.names = F)
  return_df
}


#' install packages from a file exported by fun_utils_export_enviroment
#'
#' @param pkg_enviroment_file the outfile of fun_utils_export_enviroment
#' @param failed_files packges not saved
#' @param cran_version T or F, whether install version
#' @param cran_repos_url repos
#' @param proxy proxy used for biocmanager
#'
#' @return NULL
#' @export
#'
fun_utils_install_enviroment <- \(pkg_enviroment_file,failed_files,
                                  cran_version = F,
                                  cran_repos_url = 'https://mirrors.bfsu.edu.cn/CRAN/',
                                  proxy = 'socks5://127.0.0.1:1080'){
  df <- read.csv(pkg_enviroment_file)
  stopifnot(all(c("pkg_name",'repos','version','github','github_ref') %in% colnames(df)))
  rownames(df) <- df$pkg_name
  library(devtools)
  library(BiocManager)
  failed_install <- c()
  # current ins
  for(each_pkg in rownames(df)){
    message('**********************************************************************')
    pkg_found <- requireNamespace(each_pkg)
    repos <- df[each_pkg,'repos']
    version <- df[each_pkg,'version']
    github <- df[each_pkg,'github']
    github_ref <- df[each_pkg,'github_ref']

    if(pkg_found){
      message('*   ',each_pkg, ' found!')
    }else{
      if(repos == 'CRAN'){
        message('*   ',each_pkg, ' not found, will be installed from CRAN!')
        tryCatch(
          {
            Sys.setenv(HTTP_PROXY = "")
            Sys.setenv(HTTPS_PROXY = "")
            if(!cran_version) version = F;
            devtools::install_version(each_pkg,version = version,repos = cran_repos_url);
          }, # end of try
          error = \(e){NULL}# end of error function
        ) # end of tryCatch
      }else if(grepl(repos,"Bioconductor")){
        message('*   ',each_pkg, ' not found, will be installed from Bioconductor!')
        pkg_info <- tryCatch(
          {
            if(length(proxy) >0){
              Sys.setenv(HTTP_PROXY = "socks5://127.0.0.1:1080")
              Sys.setenv(HTTPS_PROXY = "socks5://127.0.0.1:1080")
            }
            BiocManager::install(each_pkg,ask = F)
          },# end of  try
          error = \(e){NULL} # end of error function
        )# end of tryCatch
      }else if(repos == 'github'){
        message('*   ',each_pkg, ' not found, will be installed from github!')
        Sys.setenv(HTTP_PROXY = "")
        Sys.setenv(HTTPS_PROXY = "")
        devtools::install_github(github,ref = github_ref)
      }else{
        message('*   ',each_pkg, ' not found, and try using Bioconductor!')
      }
      # try loading
      pkg_found <- requireNamespace(each_pkg)
      if(!pkg_found){
        failed_install <- c(failed_install,each_pkg)
        if(length(proxy) >0){
          Sys.setenv(HTTP_PROXY = "socks5://127.0.0.1:1080")
          Sys.setenv(HTTPS_PROXY = "socks5://127.0.0.1:1080")
        }
        BiocManager::install(each_pkg,ask = F)
      }
    }
  } # end of for install for loop
  df_failed <- df[failed_install,]
  df_failed |> write.csv(failed_files,row.names = F)

  message('**********************************************************************')
  message('      Total packages: ',nrow(df))
  message('Successful installed: ',nrow(df) - nrow(df_failed))
  message('    Failed installed: ',nrow(df_failed))
  message('**********************************************************************')
}
