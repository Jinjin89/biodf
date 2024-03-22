#' retrive avaiable outcomes
#'
#' @param ao_file output file
#' @param available_outcomes_args args to available_outcomes
#'
#' @return df
#' @export
#'

fun_mr_available_outcomes <- function(ao_file = '~/base/db/mr/available_outcomes.rds',
                                      available_outcomes_args = list()){
  if(!file.exists(ao_file)){
    ao <- do.call(TwoSampleMR::available_outcomes,available_outcomes_args)
    ao %>% saveRDS(ao_file)
    rm(ao)
  }
  readRDS(ao_file) %>%
    as.data.frame() %>%
    set_rownames(.[[1]]) %>%
    return
}

#' Get instruments SNPs
#'
#' @param exposure the exposure id
#' @param instrument_dir the output instrument-snps dir
#' @param extract_instruments_args the args
#'
#' @return df
#' @export
#'

fun_mr_extract_instrument <- function(
    exposure,
    instrument_dir = '~/base/db/mr',
    extract_instruments_args = list(p1 = 1e-10)){
  # 1) check meta data
  if(!dir.exists(instrument_dir)){
    dir.create(instrument_dir)
  }
  file_save = paste0(c(instrument_dir,paste0(exposure,'.rds')),collapse = '/')

  # 2) get ins
  if(!file.exists(file_save)){
    message('The exposure results not found!')
    extract_instruments_args$outcomes = exposure
    exposure_dat = do.call(TwoSampleMR::extract_instruments,extract_instruments_args)
    list(data = exposure_dat,
         params = extract_instruments_args) %>%
      saveRDS(file_save)
  }

  readRDS(file_save)

}

#' Mendelian randomization
#'
#' @param exposure exposure
#' @param outcome outcome
#' @param res_rds the results file
#' @param echo where print summary
#' @param mr_presso whether MR_presso
#' @param fun_mr_available_outcomes_args fun_mr_available_outcomes_args
#' @param fun_mr_extract_instrument_args fun_mr_extract_instrument_args
#' @param extract_outcome_data_args extract_outcome_data_args
#'
#' @return a list of MR results
#' @export
#'

fun_mr <- \(exposure,
            outcome,
            res_rds,
            echo=T,
            mr_presso=T,
            fun_mr_available_outcomes_args = list(),
            fun_mr_extract_instrument_args = list(
              instrument_dir ='~/base/db/mr',
              extract_instruments_args = list(p1 = 1e-10)),
            extract_outcome_data_args  =list()){

  if(!file.exists(res_rds)){
    message("===>",'results not found build it!')
    message("*    step 0: getting avaiable outcomes!")
    ao <- do.call(fun_mr_available_outcomes,fun_mr_available_outcomes_args)

    res_data <- list()
    message("*    step 1: setting up the instruments data(SNPs)!")
    # 2) get exposure snps
    res_data$exposure = exposure
    res_data$outcome = outcome

    fun_mr_extract_instrument_args$exposure = exposure
    res_data$exposure_dat = do.call(fun_mr_extract_instrument,fun_mr_extract_instrument_args)
    res_data$exposure_dat = res_data$exposure_dat$data # get data, drop the params

    message("**   step 1.1: how many snps found(exposure)!")
    message("**   step 1.1: ==>",nrow(res_data$exposure_dat))
    if(length(res_data$exposure_dat) == 0){
      return("instrument snps not found")
    }

    # 3) get outcome snps
    message("*    step 2: get the outcome snps")
    extract_outcome_data_args$snps = res_data$exposure_dat$SNP
    extract_outcome_data_args$outcome = outcome

    suppressMessages(
      {
        res_data$outcome_dat <- do.call(TwoSampleMR::extract_outcome_data, extract_outcome_data_args)
      }
    )
    message("**   step 2.1: how many snps found(outcome)!")
    message("**   step 2.1: ==>",nrow(res_data$outcome_dat))
    if(length(res_data$outcome_dat) == 0){
      return("outcome snps not found")
    }

    # 4) get harmonise data
    message("*    step 3: harmonise data")
    suppressMessages(
      {
        res_data$harmonised_dat = TwoSampleMR::harmonise_data(exposure_dat = res_data$exposure_dat,outcome_dat = res_data$outcome_dat)
      }
    )

    message("**   step 3.1: how many snps found(harmonise)!")
    message("**   step 3.1: ==>",nrow(res_data$harmonised_dat))
    if(length(res_data$harmonised_dat) == 0){
      return("harmonised_dat snps not found")
    }

    # 4) MR
    message("*    step 4: MR")
    res_data$MR = list()

    if(nrow(res_data$harmonised_dat) > 5){
      suppressMessages(
        {
          res_data$MR$MR <-  TwoSampleMR::mr(res_data$harmonised_dat)
          # perform presso
          if(mr_presso){
            res_data$MR$MR_presso <- TwoSampleMR::run_mr_presso(res_data$harmonised_dat)
          }

          # get results summary
          res_data$summary <- res_data$MR$MR %>%
            dplyr::select(method,nsnp,b,se,pval)

          if(mr_presso){
            res_data$summary <-
              res_data$summary %>%
              rbind(data.frame(
                method = "MR_presso",
                nsnp = NA,
                b = res_data$MR$MR_presso[[1]]$`Main MR results`[1,"Causal Estimate"],
                se = res_data$MR$MR_presso[[1]]$`Main MR results`[1,"Sd"],
                pval = res_data$MR$MR_presso[[1]]$`Main MR results`[1,"P-value"]))
            }
          }
      )
    }else{
      res_data$MR$MR <-  TwoSampleMR::mr(res_data$harmonised_dat)
      res_data$summary <- res_data$MR$MR %>%
        dplyr::select(method,nsnp,b,se,pval)
    }
    # merge the exposure and outcome id into summary data
    res_data$summary %<>%
      as.data.frame() %>%
      mutate(exposure_id = res_data$exposure,
             outcome_id = res_data$outcome) %>%
      mutate(exposure = ao[.$exposure_id,"trait"]) %>%
      mutate(outcome = ao[.$outcome_id,"trait"]) %>%
      mutate(exp_out = paste0(exposure, '->',outcome)) %>%
      mutate(ID = paste0(exposure_id,'->',outcome_id))

    message('Done<===')
    res_data %>% saveRDS(res_rds)
  }
  res = readRDS(res_rds)
  if(echo){
    print(res$summary)
  }
  return(res)
}

