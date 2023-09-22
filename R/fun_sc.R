#' Scissor analysis
#'
#' @param ... params passing to Scissor::Scissor
#'
#' @return seurat obj with scissor feature
#' @export
#'
#' @examples
fun_sc_seurat_scissor =
  function(...){
    suppressMessages(require(Scissor))
    input_args = list(...)
    # 1) scissor
    if(!("scissor" %in% colnames(input_args$sc_dataset@meta.data))){
      infos1 = ft_Scissor(...)

      # 2 ) append scissor results into seurat
      Scissor_select <- rep(0, ncol(input_args$sc_dataset))
      names(Scissor_select) <- colnames(input_args$sc_dataset)
      Scissor_select[infos1$Scissor_pos] <- 1
      Scissor_select[infos1$Scissor_neg] <- 2
      input_args$sc_dataset <- AddMetaData(input_args$sc_dataset,metadata = Scissor_select, col.name = "scissor")
    }

    input_args$sc_dataset
  }

#' Scissor plot
#'
#' @param input_sc  seurat plot
#' @param pt.size pt.size
#' @param group.by group.by
#' @param cols colors
#'
#' @return ggplot obj
#' @export
#'
#' @examples
fun_sc_plot_scissor = function(input_sc,pt.size = 0.8,group.by = 'scissor',
                               cols = c('grey','indianred1','royalblue')){
  DimPlot(input_sc, reduction = 'umap',
          group.by = group.by,
          cols = cols, pt.size = pt.size, order = c(2,1))

}
#' Plot trajectory using monocle2
#'
#' @param input_sc input_seurat
#' @param input_genes the genes
#' @param input_stat the cell status
#' @param input_root the root cell
#'
#' @return monocle single cell data obj
#' @export
#'
#' @examples
fun_sc_seurat_trajectory <-
  function(input_sc,input_genes,input_stat = "cell_major",input_root =NULL){
    suppressMessages(
      {
        require(Seurat)
        require(monocle)
      }
    )
    # 1)
    input_exprs_data = as.matrix(Seurat::GetAssayData(input_sc)[input_genes,])

    pd <- new("AnnotatedDataFrame", data = input_sc@meta.data)
    fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(input_exprs_data),
                                                      gene_short_name =  rownames(input_exprs_data)))

    cds <- monocle::newCellDataSet(
      input_exprs_data,
      phenoData = pd,
      featureData = fd,
      lowerDetectionLimit = 1,
      expressionFamily = VGAM::negbinomial.size())

    #cds <- estimateSizeFactors(cds)
    pData(cds)$Size_Factor <- 1
    pData(cds)$Total_mRNAs <- colSums(Biobase::exprs(cds))

    cds <- monocle::detectGenes(cds, min_expr = 1)
    expressed_genes <- row.names(subset(fData(cds),
                                        num_cells_expressed >=
                                          5))
    ordering_genes <- expressed_genes
    cds <- monocle::setOrderingFilter(cds, ordering_genes)
    cds <- monocle::reduceDimension(cds, norm_method = "log", method = "DDRTree",
                                    pseudo_expr = 1)
    cds <- monocle::orderCells(cds)

    cds$State = cds[[input_stat]]

    if(length(input_root) >0){
      cds = monocle::orderCells(cds,root_state = factor(input_root,levels = c(0,1)))
    }

    cds
  }

#' Single cell preprocessing wrapper, params set to NULL to skip
#'
#' @param input_sc input_seurat project or input matirx
#' @param filter_feature_low the minimum feature
#' @param filter_feature_high the maximum feature
#' @param filter_mt the mitocondrial gens cutoff
#' @param scale_factor the scale factor
#' @param variables_method variables_method
#' @param variables_feature variables_feature
#' @param scale_data scale_data
#' @param nPCA nPCA
#' @param clustering_dims clustering_dims
#' @param clustering_resolution clustering_resolution
#' @param umap_dims umap_dims
#'
#' @return
#' @export
#'
#' @examples
fun_sc_seurat_preprocessing <-
  function(input_sc,
           filter_feature_low = NULL,
           filter_feature_high = NULL,
           filter_mt = NULL,
           scale_factor = NULL,
           variables_method = "vst",
           variables_feature = 2000,
           scale_data = T,
           nPCA=50,
           clustering_dims = 1:15,
           clustering_resolution = 0.5,
           umap_dims = 1:15
  ){
    # 1) get count data
    if(class(input_sc) == "Seurat"){
      message("Extracting counts from seurat obj...")
      input_count = input_sc@assays$RNA@counts
    }else{
      message("using the input matrix")
      input_count = input_sc
    }

    # 2) create seurat project
    sc_seurat = Seurat::CreateSeuratObject(
      counts = input_count
    )

    # 3) filtering features and mitochondria genes effect
    sc_seurat[["percent.mt"]] = Seurat::PercentageFeatureSet(sc_seurat,pattern = "^MT-")

    # 4) filtering step
    if(!is.null(filter_feature_low)){
      sc_seurat = subset(sc_seurat,subset= nFeature_RNA > filter_feature_low)
    }
    if(!is.null(filter_feature_high)){
      sc_seurat = subset(sc_seurat,subset= nFeature_RNA < filter_feature_high)
    }
    if(!is.null(filter_mt)){
      sc_seurat = subset(sc_seurat,subset= percent.mt < filter_mt)
    }

    # 5) Normalized step
    if(!is.null(scale_factor)){
      sc_seurat = Seurat::NormalizeData(sc_seurat,normalization.method = "LogNormalize",scale.factor =scale_factor )
    }

    # 6) find varaibles
    sc_seurat = Seurat::FindVariableFeatures(sc_seurat, selection.method = variables_method,nfeatures = variables_feature)

    # 7) scale data
    if(scale_data){
      sc_seurat = Seurat::ScaleData(sc_seurat,features = rownames(sc_seurat))
    }

    # 8) PCA
    sc_seurat = Seurat::RunPCA(sc_seurat,features = Seurat::VariableFeatures(sc_seurat),npcs = nPCA)


    # 9) clustering step
    sc_seurat = Seurat::FindNeighbors(sc_seurat,dims = clustering_dims)
    sc_seurat = Seurat::FindClusters(sc_seurat,resolution = clustering_resolution)

    # 10) UMAP
    sc_seurat = Seurat::RunUMAP(sc_seurat,dims = umap_dims)


    # 11) merge previous meta data
    if(class(input_sc) == "Seurat"){
      # 1) old df
      df = input_sc@meta.data
      # 2) rename
      colnames(df) = paste0(colnames(df),"_")

      # 3) merge old meta data into new
      for(each in colnames(df)){
        sc_seurat[[each]] = df[colnames(sc_seurat),each]
      }
    }
    sc_seurat
  }

#' fun_sc_read_mtx
#'
#' @param input_dir directory contains the mtx file
#' @param input_anno data.frame, first column is sample name, second column is the mtx file position
#' @param mtx_file_search key word to find mtx dir position
#' @param outfile the outfile
#'
#' @return
#' @export
#'
#' @examples
fun_sc_read_mtx <-  function(input_dir,outfile,input_anno=NULL,mtx_file_search = "barcodes.tsv.gz"){

  if(!is.null(outfile) && file.exists(outfile)){
    return(readRDS(outfile))
  }


  if(is.null(input_anno)){
    # 1) get mtx directly
    if(file.exists(paste0(input_dir,"/",mtx_file_search))){
      message(paste0(mtx_file_search, "found"))
      return(Seurat::Read10X(input_dir))
    }else{
      message("loop the directory to get mtx file list!")
      # 1) get all the directory
      list_dirs = list.dirs(input_dir,recursive = F,full.names = T)
      # 2) get_annotations_name
      sc_names = purrr::map_chr(list_dirs,basename)
      input_anno = data.frame(row.names = sc_names)
      input_anno$sc = sc_names

      # 3) get the files format
      sc_dirs =
        purrr::map_chr(list_dirs,function(x){
          #1) get new_dirs
          file_search =
            list.files(x,recursive = T,all.files = T,full.names = T,include.dirs = T)
          # 2) get the found files
          files_each = file_search[stringr::str_detect(file_search,mtx_file_search)]

          # 3) get new data
          if(length(files_each) ==0){
            return(NA)
          }else{
            dirname(files_each[1])
          }
        })
      input_anno$sc_dir = sc_dirs
      input_anno = input_anno[!is.na(input_anno[[2]]),]

    }
  }
  message("1) if: input_anno is supplied, -> first column: sample_name,2) second column: mtx file position\n2) else if:: found mtx, read it,\n3) else: loop the dirctory to get it ")

  # 2 get from dirs
  sc_seurat_list =
    purrr::map(seq_along(input_anno[[1]]),function(i){
      sc_10x_dir = input_anno[[2]][i]
      message(paste0(i,") Reading from ",sc_10x_dir))
      sc_count = Seurat::Read10X(sc_10x_dir)
      sc_seurat = Seurat::CreateSeuratObject(counts = sc_count)
      sc_seurat$sample =input_anno[[1]][i]
      sc_seurat
    })
  names(sc_seurat_list) = input_anno[[1]]

  if(length(sc_seurat_list) == 1){
    sc_final =sc_seurat_list[[1]]
  }else{
    sc_final = merge(sc_seurat_list[[1]],sc_seurat_list[2:length(sc_seurat_list)],names(sc_seurat_list))
  }

  saveRDS(sc_final,file = outfile)

  sc_final

}

#' Calculate the progeny score for Seurat obj using decouplerR
#'
#' @param input_seurat seurat obj
#' @param out_file the calculated score output
#' @param fun_sc which expression data in Seurat to use, the default is RNA data
#' @param top how many progeny genes for each pathway to use, default is 100
#'
#' @return
#' @export
#'
#' @examples
fun_sc_seurat_progeny <- function(input_seurat,
                                  out_file,
                                  fun_sc = function(x) x@assays$RNA@data,
                                  top = 100){
  if(!file.exists(out_file)){
    require(decoupleR)

    # 1) get progeny pathway
    net <- decoupleR::get_progeny(organism = 'human', top = top)


    # 2) get assay data
    input_sc_data = do.call(fun_sc,list(x = input_seurat))

    mat <- as.matrix(input_sc_data)
    message("The matrix values range from: ")
    cat(max(mat))
    cat("---->")
    cat(min(mat))
    cat("\n")

    acts <-  decoupleR::run_wmean(mat=mat, net=net, .source='source', .target='target',
                                  .mor='weight', times = 100, minsize = 5)

    acts %>%
      data.table::fwrite(out_file)

  }
  acts =
    data.table::fread(out_file)
  # 3) get progeny profile
  input_seurat[['pathwayswmean']] <- acts %>%
    dplyr::filter(statistic == 'norm_wmean') %>%
    tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                       values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)

  DefaultAssay(object = input_seurat) <- "pathwayswmean"

  input_seurat <- ScaleData(input_seurat)
  input_seurat@assays$pathwayswmean@data <- input_seurat@assays$pathwayswmean@scale.data


  input_seurat
}

#' cell_chat analysis
#'
#' @param input_sc Seurat obj
#' @param cellchat_rds the output file used to save results
#' @param group.by cellchat idents
#' @param cell_search cell db use
#' @param workers parrallel
#' @param plan multisession or ...
#'
#' @return cellchat obj
#' @export
#'
#' @examples
fun_sc_cellchat <- function(input_sc,
                            cellchat_rds,
                            group.by = "cell_type",
                            cell_search = "Secreted Signaling",
                            plan = "multisession",
                            workers = 4){
  if(!file.exists(cellchat_rds)){
    print("the cellchat obj not found, build it")
    print(cellchat_rds)
    # 1)
    require(CellChat)
    require(parallel)
    require(future)

    # 2) prepare cellchat data
    cellchat <- createCellChat(
      object = input_sc,
      meta = input_sc@meta.data,
      group.by = group.by)

    # 3) add meta data
    #cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = group.by)
    groupSize <- as.numeric(table(cellchat@idents))

    # 4) get db
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)
    if(length(cell_search) == 0 ){
      print("using default cellChatDB")
      CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    }else{
      print("Using subseted data:")
      print(cell_search)
      CellChatDB.use <- subsetDB(CellChatDB, search = cell_search)
    }
    cellchat@DB <- CellChatDB.use

    # 5)
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan(future::sequential) # reset the plan with sequentail
    future::plan(plan, workers = workers) # do parallel

    # 6) get overexpressed genes
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # 7)
    cellchat <- computeCommunProb(cellchat)

    # 8)
    cellchat <- computeCommunProbPathway(cellchat)

    # 9) aggregate
    cellchat <- aggregateNet(cellchat)

    # 10) save it
    cellchat %>%
      saveRDS(cellchat_rds)
  }
  readRDS(cellchat_rds)
}

#' plot cell chat
#'
#' @param input_cellchat cell chat
#' @param input_cell which cell to plot, NULL or cell types
#'
#' @return NULL
#' @export
#'
#' @examples
fun_sc_cellchat_plot <- function(
    input_cellchat,
    input_cell = NULL){

  # 1) get group size
  groupSize <- as.numeric(table(input_cellchat@idents))

  # 2)
  if(length(input_cell) == 0 ){
    netVisual_circle(input_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

  }else{
    mat <- input_cellchat@net$weight
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[input_cell, ] <- mat[input_cell, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = input_cell)
  }
}




#' Seurat obj signature deconvolution usingUCell
#'
#' @param input_seurat the Seurat object
#' @param input_list  signature list, must has name
#' @param ucell_params parmaters passing to AddModuleScore_UCell
#'
#' @return Seurat obj with signature deconvoluted
#' @export
#'
fun_sc_seurat_ucell <- function(
    input_seurat,
    input_list,
    ucell_params = list(name = "")){
  suppressMessages({
    library(UCell)
  })

  # 1) get signature name
  signature.names <- names(input_list)

  # 2) check if all in single cell data, not run
  check_index = purrr::map_lgl(signature.names,function(x) x %in% colnames(input_seurat@meta.data))
  if(all(check_index)){
    message("all found, skip the step")
  }else{
    input_list = input_list[!check_index]
    message(paste0("Found ",sum(check_index),", remain " ,length(input_list)," to be done!"))
    ucell_params$obj = input_seurat
    ucell_params$features = input_list
    input_seurat = do.call(AddModuleScore_UCell,ucell_params)
  }

  return(input_seurat)
}

