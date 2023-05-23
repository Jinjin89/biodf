fun_sc_seurat_scissor =
  function(input_sc){

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
