#' Merge gene mutation status into clinical data
#'
#' @param input_df imvigor clinical data
#' @param input_variables genes
#'
#' @return data.frame, the merged df with genes mutation status
#' @export
#'

fun_extend_imvigor_mut <- function(
    input_df,input_variables){
  # 1)
  require(IMvigor210CoreBiologies)
  data("fmone",package = "IMvigor210CoreBiologies")

  # 2) get_genes_mutation-status
  known_short = fmone@assayData$known_short
  known_short = as.data.frame(known_short)
  known_short = as.data.frame(t(known_short))
  known_short$accession = rownames(known_short)
  known_short = tidyr::pivot_longer(known_short,-accession)

  # 3) get samples_patient
  pheno_data = fmone@phenoData@data
  known_short$ID = pheno_data[known_short$accession,"ANONPT_ID"]
  known_short$sample = rownames(input_df)[match(known_short$ID,input_df$ANONPT_ID)]
  #known_short = as.data.frame(known_short)

  known_short$value = ifelse(known_short$value == "",0,1)
  known_short %<>%
    dplyr::filter(!is.na(.$sample)) %>%
    dplyr::filter(value != 0) %>%
    dplyr::count(name,sample)
  all_samples_mut = known_short$sample %>% unique

  if(is.null(input_variables)){
    input_variables = unique(known_short$name)
  }else{
    input_variables = intersect(unique(known_short$name),input_variables)
  }

  for(g in input_variables){
    samples_found =
      known_short %>%
      dplyr::filter(.$name == g) %>%
      dplyr::pull(sample) %>%
      unique

    input_df[[paste0(g,"_status")]] =
      ifelse(input_df$sample %in% samples_found,"mt",
             ifelse(input_df$sample %in% all_samples_mut,"wt",NA))
  }
  input_df

}




fun_extend_pRRophetic <- function(
    input_mat,out_file,
    pRRopheticPredict_args = list()){
  if(!file.exists(out_file)){
    require(pRRophetic)
    # 1)get drugs
    drugs = stringr::str_split(
      "A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439",
      ","
    ) %>% unlist() %>% stringr::str_remove_all(" ")
    # get args
    default_args = list(
      tissueType = "allSolidTumors",
      batchCorrect = "eb",
      powerTransformPhenotype = TRUE,
      removeLowVaryingGenes = 0.2,
      minNumSamples = 10, selection =1,
      printOutput = TRUE,
      removeLowVaringGenesFrom = "homogenizeData",
      dataset = "cgp2014"
    )

    # get new args
    new_arg_names = names(pRRopheticPredict_args)
    for(i in seq_along(new_arg_names)){
      each_name = new_arg_names[i]
      if(each_name != ""){
        message("current params:",each_name, " is :",pRRopheticPredict_args[[each_name]])
        default_args[[each_name]] =pRRopheticPredict_args[[each_name]]
      }
    }
    default_args$testMatrix = input_mat


    # run the prediction

    purrr::map_dfc(drugs,\(each_drug){
      print(each_drug)
      current_args <-  default_args
      current_args$drug <-  each_drug
      res <- do.call(pRRophetic::pRRopheticPredict,current_args)
      res <-  data.frame(res)
      colnames(res) <-  each_drug
      return(res)
    }
    ) -> final_res
    final_res$sample <-  rownames(final_res)
    final_res %>% data.table::fwrite(out_file)
  }

  data.table::fread(out_file) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$sample)
}



#' Signature deconvolution: TIP
#'
#' @param input_mat input matrix
#' @param outfile results saving files
#'
#' @return data.frame
#' @export
#'

fun_extend_sig_tip <- function(input_mat,outfile){
  # 1)
  if(!file.exists(outfile)){

    sig_names = names(biodata::sig_list$step)
    sig_df <-
      purrr::map(sig_names,\(each_sig){
        sig_df = biodata::sig_list$step[[each_sig]]
        data.frame(gene = sig_df) %>%
          dplyr::mutate(term = each_sig)
      }) %>%
      purrr::list_rbind() %>%
      dplyr::count(term,gene)
    # 2) deconvolute
    deconv_df = fun_sig_deconv(
      input_mat = input_mat,
      input_sig = sig_df,
      outfile = outfile,
      statistic_used = 'ssgsea',
      input_args = list(method = 'ssgsea',minsize  = 1)
    ) %>% t %>%
      as.data.frame()
    deconv_df <-
      deconv_df%>%
      dplyr::mutate(Step2 = Step2.positive - Step2.negative) %>%
      dplyr::mutate(Step3 = Step3.positive - Step3.negative) %>%
      dplyr::mutate(Step5 = Step5.positive - Step5.negative) %>%
      dplyr::mutate(Step6 = Step6.positive - Step6.negative) %>%
      dplyr::mutate(Step7 = Step7.positive - Step7.negative) %>%
      dplyr::select(-matches('positive$|negative$')) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(),fun_norm_scale))
    print(deconv_df)
    # 3) soring df
    tip_names = colnames(deconv_df) %>% sort
    deconv_df[,tip_names] %>%
      mutate(sample = rownames(.)) %>%
      saveRDS(outfile)
  }
  readRDS(outfile)
}



