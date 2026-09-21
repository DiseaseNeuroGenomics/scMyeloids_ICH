# Author: Dimitrios Kyriakis



# ================================  Double Detection ==============================================
GetScrubletScores <- function(mat,sample_name, min.molecules.per.gene=10,
                              method=c("scrublet", "doubletDetection"),workdir) {
    tf.in <- tempfile()
    tf.out <- tempfile()
    hist_out <- paste0(workdir,"/2.scrublet_hist_",sample_name,".pdf")
    umap_out <- paste0(workdir,"/2.scrublet_umap_",sample_name,".pdf")
    dt <- mat[Matrix::rowSums(mat)>=min.molecules.per.gene,] %>% Matrix::t() %>% data.table::data.table()
    data.table::fwrite(dt, file=tf.in)
    cmd <- paste0("python", " -c 'import sys; import pandas; import scrublet;import matplotlib.pyplot as plt; ",
                    "df = pandas.read_csv(\"", tf.in, "\"); scrub = scrublet.Scrublet(df); ",
                    "doublet_scores, predicted_doublets = scrub.scrub_doublets();",
                    "pandas.DataFrame(dict(score=doublet_scores, is_doublet=predicted_doublets)).to_csv(\"",tf.out,"\");",
                    "scrub.plot_histogram();plt.savefig(\"",hist_out,"\");",
                    "scrub.set_embedding(\"UMAP\", scrublet.get_umap(scrub.manifold_obs_, 10, min_dist=0.3));",
                    "scrub.plot_embedding(\"UMAP\", order_points=True);plt.savefig(\"",umap_out,"\");'",
                    sep='')
    system(cmd, intern=F)
    x <- data.table::fread(tf.out,sep=',')[,2:3] %>% as.data.frame() %>% as.list() %>% lapply(`names<-`, colnames(mat))
    return(x)
}


# Create a function that compute all required parameters
outliers_madEst <- function(x,
                            b = 1.4826,
                            threshold = 2,
                            na.rm = TRUE){

    # If x is numeric or integer, applying the function.
    # Otherwise, stopping it.
    if(inherits(x,c("numeric","integer")) == FALSE)
    stop("x is neither numeric nor integer")

    if (na.rm == TRUE) {
    data <- na.omit(x)   # incomplete cases are removed
    } else {data <- x}

    # Calculate the MAD
    # computing the median of the data
    # median = the center of the CI defining acceptable values
    center <- median(data)
    # computing the median of data centered around the sample median
    MAD <- b*median(abs(data-center))
    # how many MAD from the median are the limits of the IC?
    half_CI <- threshold*MAD

    # Calculate the range of acceptable values
    LL_CI_MAD <- center-half_CI # lower limit of the median CI
    UL_CI_MAD <- center+half_CI # upper limit of the median CI

    # calculate the outliers
    outliers <- c(data[data < LL_CI_MAD],data[data > UL_CI_MAD])
    outliers_pos <- c(which(data < LL_CI_MAD),which(data > UL_CI_MAD))

    # Return results in list()
    invisible(list(Median = center,
                    MAD = MAD,
                    LL_CI_MAD = LL_CI_MAD,
                    UL_CI_MAD = UL_CI_MAD,
                    L_outliers = sort(data[data < LL_CI_MAD]),
                    U_outliers = sort(data[data > UL_CI_MAD]),
                    outliers = sort(c(data[data < LL_CI_MAD],
                                    data[data > UL_CI_MAD])),
                    outliers_pos=outliers_pos,
                    threshold=threshold))

}

outliers_mad<- function(x,b = 1.4826,threshold = 3,na.rm = TRUE){

    threshold <- as.numeric(threshold)
    na.rm <- as.logical(na.rm)

    out <- outliers_madEst(x,b,threshold,na.rm)
    out$call <- match.call()
    out$median <- out$Median
    out$MAD <- out$MAD
    out$limits <- as.vector(c(lower = out$LL_CI_MAD,upper = out$UL_CI_MAD))
    out$nb <- c("extremely low" = length(out$L_outliers),
                "extremely high" = length(out$U_outliers),
                total = length(out$outliers))

    class(out) <- "outliers_mad"
    out
}




read_10x_fun <- function(h5_file){
    # ===================== READ DATA ======================
    DEM <- as.data.frame(Seurat::Read10X_h5(h5_file))
    DEM <- DEM  %>% .[, colSums(.) > 10]  %>% .[,order(colSums(.),decreasing=T)]
    # cat(dim(DEM)[2])
    # ------------------------------------------------------
    return(DEM)
}

# ================================  Mitochondrial fraction ==============================================
mit_frac_calc <- function(DEM,workdir){
    rownames(DEM)[grep("^Mt-", rownames(DEM),ignore.case=T)]
    mt_pos <- grep("^Mt-", rownames(DEM),ignore.case=T)
    mit_frac_per_dataset <-  colSums(DEM[mt_pos, ]) / colSums(DEM)
    pdf(paste0(workdir,"/2.Mit_fraction.pdf"),width= 12,height=12)
    qplot(mit_frac_per_dataset[mit_frac_per_dataset < 1], xlab="Mit. fraction", ylab="#Cells", xlim=c(-0.01, 1), bins=30)
    dev.off()
    # table(mit_frac_per_dataset<0.1)
    DEM<- DEM[,mit_frac_per_dataset<0.1]
    # dim(DEM)
    return(DEM)
}
# --------------------------------------------------------------------------------------



# ================================  Scrublet ==============================================
scrublet_run <- function(DEM,sample,workdir,min.molecules.per.gene=100){
    # Careful the number of cores you use. Maybe an error in cl
    scrublet_doublet_info <- GetScrubletScores(mat=as.matrix(DEM),sample_name=sample,min.molecules.per.gene=min.molecules.per.gene, method="scrublet",workdir=workdir)
    # REMOVE PREDICTED DOUBLETS
    DEM <- DEM[,!scrublet_doublet_info$is_doublet]
    # dim(DEM)
    return(list("DEM"=DEM,"scr_info"=scrublet_doublet_info))
}
# --------------------------------------------------------------------------------------




# ============================= REMOVE MT/RBS PREPARE SEURAT OBJECT ===================================
remove_mt_rb_genes <- function(DEM,sample,workdir){
    # REMOVE rIBOSOMAL AND mIT GENES
    seur <- CreateSeuratObject(counts = DEM, project = sample, min.cells = 0, min.features = 0)
    seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
    seur[["percent.ribo"]] <- PercentageFeatureSet(seur, pattern = "^RPL|^RPS")
    percent.mt <- seur$percent.mt
    percent.ribo <- seur$percent.ribo
    DEM_r <- seur@assays$RNA@counts
    mt_pos <- grep("^Mt-|^RPL|^RPS", rownames(DEM_r),ignore.case=T)
    # rownames(DEM_r)[mt_pos]
    DEM_remain <- DEM_r[-mt_pos,]
    mt_pos <- grep("^Mt-|^RPL|^RPS", rownames(DEM_remain),ignore.case=T)
    # rownames(DEM_remain)[mt_pos]
    # dim(DEM)
    # dim(DEM_remain)
    # length(mt_pos)
    return(list("DEM"=DEM_remain,"percent.mt"=percent.mt,"percent.ribo"=percent.ribo))
}
# -----------------------------------------




# ============================= PREPARE SEURAT OBJECT ===================================
prepare_seura_object <- function(DEM,sample,percent.mt,percent.ribo,scrublet_doublet_info,workdir){
    seur <- CreateSeuratObject(counts = DEM, project = sample, min.cells = 10, min.features = 200)
    percent.mt_remain <- percent.mt[names(percent.mt) %in%names(seur$nCount_RNA)]
    percent.ribo_remain <- percent.ribo[names(percent.ribo) %in%names(seur$nCount_RNA)]
    seur$percent.mt <- percent.mt_remain
    seur$percent.ribo <- percent.ribo_remain

    # ====================== ADD INFO ========================
    seur$auto_scrub_is_doublet <- scrublet_doublet_info$is_doublet
    # --------------------------------------------------------
    p1 <- VlnPlot(seur,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo"),group.by="orig.ident")
    nFeature_Min <- outliers_mad(seur$nFeature_RNA,threshold = 3)$limits[1]
    nFeature_Max <- outliers_mad(seur$nFeature_RNA,threshold = 3)$limits[2]
    nCount_Min <- outliers_mad(seur$nCount_RNA,threshold = 3)$limits[1]
    nCount_Max<- outliers_mad(seur$nCount_RNA,threshold = 3)$limits[2]

    seur <- subset(seur, subset = nCount_RNA  > nCount_Min  & nCount_RNA < nCount_Max & nFeature_RNA > nFeature_Min & nFeature_RNA < nFeature_Max & percent.mt < 10 )
    seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
    seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seur)
    seur <- ScaleData(seur, features = all.genes,verbose = FALSE)
    # ====== Cell Cycle
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    seur <- CellCycleScoring(seur, s.features = s.genes, g2m.features = g2m.genes)
    # ------------------------------
    seur <- RunPCA(seur, features = VariableFeatures(object = seur))
    seur <- FindNeighbors(seur, dims = 1:20)
    seur <- FindClusters(seur, resolution = 0.5)
    seur <- RunUMAP(seur, dims = 1:20)
    # ----------------------------------------------------------------------------------------
    return(seur)
}


# Define a function to calculate MAD (Median Absolute Deviation)-based cutoff limits
calc_mad_cutoff <- function(num_vec, k = 2.5) {
    log_n_num_vec <- log1p(num_vec)  # Log-transform the input vector
    mad_log_n_num_vec <- mad(log_n_num_vec)  # Calculate MAD of the log-transformed vector
    up_mad_log_n_num_vec <- median(log_n_num_vec) + k * mad_log_n_num_vec  # Upper limit
    down_mad_log_n_num_vec <- median(log_n_num_vec) - k * mad_log_n_num_vec  # Lower limit
    return(c(exp(down_mad_log_n_num_vec), exp(up_mad_log_n_num_vec)))  # Return cutoff limits in original scale
}
