# Author: Dimitrios Kyriakis
# ========================= COLORS ============================
CellPop_colors <- c('#FDB462','#FB8072','#80B1D3','#B3DE69','#8DD3C7')
Colors_of_CL <- c('#E69F00','#56B4E9','#B3DE69','#1C91D4','#80B1D3','#666666','#FDB462','#8DD3C7','#AD7700','#D55E00','#CC79A7','#009E73','#0072B2')
names(Colors_of_CL) <- c('0','1','NEUT','8','MONO','4','B','T_NK','7','2','3','11','12')
Colors_of_Subtype <- c('#56B4E9','#0072B2','#005685','#B3DE69','#80B1D3','#D5C711','#CC79A7','#FDB462','#007756','#8DD3C7','#F0E442','#E69F00','#009E73','#FB8072','#D55E00','#1C91D4','#666666')
names(Colors_of_Subtype) <- c('MG_Adapt_CCL3','MG_Adapt_HSPA1A','MG_PVM_GPNMB','NEUT','MONO','MG_PVM_CD163','MG_Adapt_TMEM163','B Cells','MG_Prolif_MKI67','T Cells','MG_Adapt_HIST','MG_Adapt_AIF1','MG_Adapt_HIF1A','MG_Homeo_FRMD4A','MG_Adapt_IFI44L','MG_Homeo_PICALM','MG_exAM_ERN1')


colors <- c('Astrocytes'='#D9D9D9', # gray
            'B Cells'='#FDB462', # orange
            'Myeloid'='#FB8072', # red
            'MONO'='#80B1D3', # blue
            'Murel'='#BEBADA', # purple
            'NEUT'='#B3DE69', # green
            'Oligo'='#FCCDE5', # pink
            'T Cells'='#8DD3C7' # ciel
            )


Project_Colors <- c(
    # Cell populations
    'Astrocytes'='#D9D9D9', # gray
    'B Cells'='#FDB462', # orange
    'Myeloid'='#FB8072', # red
    'MONO'='#80B1D3', # blue
    'Murel'='#BEBADA', # purple
    'NEUT'='#B3DE69', # green
    'Oligo'='#FCCDE5', # pink
    'T Cells'='#8DD3C7', # ciel
    'MONO/NEUT'='#80B1D3',  # blue
    # Cell subtypes
    'MG_Homeo' = '#279e68', # dark green
    'MG_PVM' = '#ff7f0e', # dark orange
    'MG_Prolif' = '#d62728', # red
    'MG_Adapt' = '#1f77b4', # blue
    'MG_exAM' = '#ff9896', # light pink

    'MG_Homeo_FRMD4A' = '#279e68', # dark green
    'MG_Homeo_PICALM' = '#98df8a', # light green

    'MG_PVM_CD163' = '#ff7f0e', # dark orange
    'MG_PVM_GPNMB' = '#ffbb78', # light orange

    'MG_Prolif_MKI67' = '#d62728', # red

    'MG_Adapt_TMEM163' = '#1f77b4', # blue
    'MG_Adapt_AIF1' = '#c5b0d5', # light purple
    'MG_Adapt_IFI44L' = '#17becf', # cyan
    'MG_Adapt_HIF1A' = '#aec7e8', # light blue
    'MG_Adapt_CCL3' = '#aa40fc', # purple
    'MG_Adapt_HSPA1A' = '#9edae5', # light blue
    'MG_Adapt_HIST' = '#dbdb8d', # light yellow

    'MG_exAM_ERN1' = '#ff9896'
)
# ---------------------------------------------------------------



# ===================== Load All Libraries ======================
load_libraries <- function(){
    # ## We load the required packages
    library(dplyr)
    library(Seurat)
    library(SeuratDisk)
    library(SingleCellExperiment)
    library(harmony)


    library(WGCNA)
    library(hdWGCNA)
    library(dreamlet)

    library(readr)
    library(cowplot)
    library(patchwork)
    library(RColorBrewer)
    # co-expression network analysis packages:
    library(igraph)
    library(dittoSeq)
    library(scCustomize)
    library(Nebulosa)
    # library(UCell)
    '%notin%' <- Negate('%in%')
    colors <- brewer.pal(n = 9, name = 'Set3')
    colors[2] <- colors[8]
    colors[8] <- "azure3"
    cl_colors <- DiscretePalette_scCustomize(num_colors = 26,palette = "polychrome")
}




# ===================== SAVING Functions ========================
# saveCRDS() and readCRDS() - writing R objects with compression
saveCRDS <- function(object, filename, filter='zstd') {
  stopifnot(filter %in% c('zstd', 'lz4'))
  con = archive::file_write(file = filename, filter=filter)
  open(con)
  saveRDS(object, con)
  close(con)
}

readCRDS <- function(filename) {
  con <- archive::file_read(file = filename)
  res <- readRDS(con)
  close(con)
  res
}
# ---------------------------------------------------------------


# =========== cellspecificity return n genes =============
my_fun_return_x_max <- function(x,n){
    order(x,decreasing = T)[1:n]
}
# ---------------------------------------------------------------


# ===================== My heatmap for cellspecificity ========================
MyHeat <- function(x, genes = rownames(x), color="darkblue", assays=colnames(x)){
    # intersect preserving order from assays
	assays = intersect(assays, colnames(x))
	if( length(assays) == 0) stop("No valid assays selected")
	x = x[,assays,drop=FALSE]
	# subset based on specified genes
	x = x[rownames(x) %in% unique(genes),,drop=FALSE]
	# pass R CMD check
	value = variable = gene = NA
	# omit column totalCPM, if it exists
	i = which(colnames(x) == "totalCPM")
	if( length(i) > 0) x = x[,-1,drop=FALSE]
    df = data.frame(gene = rownames(x), x, check.names=FALSE)
	df_melt = reshape2::melt(df, id.vars="gene")
	df_melt$gene = factor(df_melt$gene, unique(genes))
	df_melt$variable = factor(df_melt$variable, assays)
	ratio = nlevels(df_melt$gene) / nlevels(df_melt$variable)
	# heatmap of cell type specificity
	ggplot(df_melt, aes(variable, gene, fill=value)) +
		geom_tile() +
		theme_classic() +
          scale_fill_gradient2(low = "#075AFF",
                           mid = "#FFFFCC",
                           high = "#FF0000")+
        xlab('') + ylab('') +
		ggtitle("Cell type specificity scores")+theme(text=element_text(size=21))+
        RotatedAxis()
}
# ---------------------------------------------------------------

my_VroomPlot<-function(x, ncol=3, alpha=.5, assays = names(x)){
	# Pass R CMD check
	y = NULL
	# intersect preserving order from assays
	assays = intersect(assays, names(x))
	if( length(assays) == 0) stop("No valid assays selected")

	# get common range across all plots
	###################################
	df_range = lapply( assays, function(id){

		if( !is.null(x[[id]]$voom.xy) ){
			res = with(x[[id]]$voom.xy, data.frame(range(x),range(y), id=id))
		}else{
			res = NULL
		}
		res
	})
	df_range = do.call(rbind, df_range)

	if( is.null(df_range) ){
		stop("Voom was not run on this object")
	}

	xlim = range(df_range$range.x.)
	ylim = c(0, max(df_range$range.y.))

	xlab = bquote(log[2](counts + 0.5))
	ylab = bquote(sqrt(standard~deviation))

	# only included assays were voom succeeded
	validAssays = droplevels(factor(unique(df_range)$id, assays))

	# make data.frame of points
	df.list = lapply( validAssays, function(id){
		with(x[[id]]$voom.xy, data.frame(id, x,y))
	})
	df_points = do.call(rbind, df.list)
	df_points$id = droplevels(factor(df_points$id, assays))
	df_points = df_points[order(df_points$id),]

	# make data.frame of curves
	df.list = lapply( validAssays, function(id){
		with(x[[id]]$voom.line, data.frame(id, x,y))
	})
	df_curve = do.call(rbind, df.list)
    unique_ids <- unique(df_curve$id)
    plotlist_vroom <- lapply(unique_ids,function(sub_id){
        df_point_sub <- df_points %>% filter(id == sub_id)
        df_curve_sub <- df_curve %>% filter(id == sub_id)

        ggplot(df_point_sub, aes(x,y)) +
            geom_point(size=0.1, alpha=alpha) +
            theme_classic() +
            theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5))+
            facet_wrap(~id, ncol=ncol) +
            xlab(xlab) +
            ylab(ylab) +
            xlim(xlim) +
            ylim(ylim) +
            geom_line(data = df_curve_sub, aes(x,y), color="red")
    })
    ggarrange(plotlist = plotlist_vroom,ncol = ncol,align ="hv")

}



gene_score_calc <- function(dataframe_obj){
    dataframe_obj <- dataframe_obj%>% mutate(
                enrichment.ratio = pct.1 / (pct.2 + .000001),
                diff.pct = pct.1 - pct.2,
                tstat= -log10(p_val_adj + 1e-320) * sign(avg_log2FC) ,
                gene.score = avg_log2FC * enrichment.ratio) %>%
                dplyr::select(gene,tstat, gene.score,p_val, p_val_adj, diff.pct, enrichment.ratio,
                    pct.1, pct.2, avg_log2FC,cluster)
    return(dataframe_obj)
}



binnary_gene_score_calc <- function(dataframe_obj){
    dataframe_obj <- dataframe_obj%>% mutate(
                cluster = ifelse(avg_log2FC >0 , "PD" ,"HC"),
                enrichment.ratio = ifelse(pct.1>pct.2, pct.1 / (pct.2 + .000001), pct.2 / (pct.1 + .000001)),
                diff.pct = pct.1 - pct.2,
                tstat= -log10(p_val_adj + 1e-320) * sign(avg_log2FC) ,
                gene.score = avg_log2FC * enrichment.ratio) %>%
                dplyr::select(gene,tstat, gene.score,p_val, p_val_adj, diff.pct, enrichment.ratio,
                    pct.1, pct.2, avg_log2FC,cluster)
    return(dataframe_obj)
}




my_Vlc_plot <- function(AllDF,x_lim,cellType,n_genes=10){
    DF <- AllDF[AllDF$assay==cellType,]
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    DF$diffexpressed[DF$logFC > 1.5 & DF$adj.P.Val < 0.05] <- "UP"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    DF$diffexpressed[DF$logFC < -1.5 & DF$adj.P.Val < 0.05] <- "DOWN"
    DF$delabel <- NA
    # add a column of NAs
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
    DF$diffexpressed[(DF$logFC > 1 & DF$adj.P.Val < 0.05) | (DF$logFC < -1 & DF$adj.P.Val < 0.05)] <- "p-value and log2FC >1"
    DF$diffexpressed[(DF$logFC > 1.5 & DF$adj.P.Val < 0.05) | (DF$logFC < -1.5 & DF$adj.P.Val < 0.05)] <-  "p-value and log2FC >1.5"

    DF <- DF[DF$assay==cellType,]
    pos_sub_volc <- DF %>% filter(logFC > 1.5 , adj.P.Val<0.05 )
    up_labs   <- pos_sub_volc$gene.name[1:n_genes]
    down_sub_volc <- DF %>% filter(logFC < -1.5 , adj.P.Val < 0.05)
    down_labs   <- down_sub_volc$gene.name[1:n_genes]
    DF$delabel <- NA
    DF$delabel[DF$gene.name %in% c(up_labs,down_labs)] <- DF$gene.name[DF$gene.name %in% c(up_labs,down_labs)]
    x_lim <- ceiling(max(abs(DF$logFC)))


    ggplot(data=DF, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
        geom_point(size=0.2) +
        theme_cowplot()  + theme_light()+
        geom_text_repel(color = 'black') +
        scale_shape_manual(values =c(19,1,1)) +
        scale_color_manual(values=c("black","blue","red")) +
        geom_vline(xintercept=c(-1.5, 1.5), col="black", linetype = "longdash") +
        geom_vline(xintercept=c(-1, 1), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")+
        ggtitle(cellType, subtitle = paste0("Total genes tested: ",dim(DF)[1] ))+

        xlim(-x_lim,x_lim)+
        theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
                plot.title = element_text(size = 21,face = "bold"),
                axis.text=element_text(size=12),
                axis.title=element_text(size=14),
                legend.text=element_text(size=12),legend.title=element_blank())
}


dh_meta_analysis = function(tabList) {
    if (is.null(names(tabList))) {
        names(tabList) = as.character(seq(length(tabList)))
    }
    for (key in names(tabList)) {
        tabList[[key]]$Dataset = key
    }
    df = do.call(rbind, tabList)
    df %>%
        as_tibble() %>%
        group_by(assay) %>%
        do(tidy(rma(yi = logFC, sei = logFC / t, data = ., method = "FE"))) %>%
        select(-term, -type) %>%
        mutate(FDR = p.adjust(p.value, "fdr")) %>%
        mutate('log10FDR' = -log10(FDR))
}

dh_plotTree = function(tree, low = "grey90", mid = "red", high = "darkred", xmax.scale = 1.5) {
    fig = ggtree(tree, branch.length = "none") +
        geom_tiplab(color = "black", size = 4, hjust = 0, offset = 0.2) +
        theme(legend.position = "top left", plot.title = element_text(hjust = 0.5))
    xmax = layer_scales(fig)$x$range$range[2]
    fig + xlim(0, xmax * xmax.scale)
}

dh_plotCoef = function(tab, coef, fig.tree, low = "grey90", mid = "red", high = "darkred", ylab) {
    tab$logFC = tab$estimate
    tab$celltype = factor(tab$assay, rev(ggtree::get_taxa_name(fig.tree)))
    tab$se = tab$std.error
    fig.es = ggplot(tab, aes(celltype, logFC)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
        geom_errorbar(aes(ymin = logFC - 1.96 * se, ymax = logFC + 1.96 * se), width = 0) +
        geom_point2(aes(color = pmin(4, -log10(FDR)), size = pmin(4, -log10(FDR)))) +
        scale_color_gradient2(name = bquote(-log[10]~FDR),
                              limits = c(0, 4), low = low, mid = mid, high = high, midpoint = -log10(0.01)) +
        scale_size_area(name = bquote(-log[10]~FDR), limits = c(0, 4)) +
        geom_text2(aes(label = '+', subset = FDR < 0.05), color = "white", size = 6, vjust = 0.4, hjust = 0.5) +
        theme_classic() +
        coord_flip() +
        xlab('') +
        ylab(ylab) +
        theme(axis.text.y = element_blank(),
              axis.text = element_text(size = 12),
              axis.ticks.y = element_blank(),
              text = element_text(size = 20)) +
        scale_y_continuous(breaks = scales::breaks_pretty(3))
    return(fig.es)
}


JaccardSets<- function(set1, set2){
        length(intersect(set1, set2))/length(unique(c(set1, set2)))
}


PairWiseJaccardSets<- function(Metadata,ident.1, ident.2){
    ident1_vec <- as.factor(as.vector(Metadata[[ident.1]]))
    names(ident1_vec) <- rownames(Metadata)
    ident2_vec <- as.factor(as.vector(Metadata[[ident.2]]))
    names(ident2_vec) <- rownames(Metadata)

    ident1.list<- split(names(ident1_vec), ident1_vec)
    ident2.list<- split(names(ident2_vec), ident2_vec)
    res<- matrix(nrow = length(ident1.list), ncol = length(ident2.list),
                    dimnames = list(names(ident1.list), names(ident2.list)))
    for (i in seq_along(ident1.list)){
            res[i, ]<- purrr::map_dbl(ident2.list, ~JaccardSets(ident1.list[[i]], .x))
    }
    return(res)
}


CalculateSilhouette<- function(object, dims = 1:50){
        if (length(dims) > ncol(object@reductions$pca@cell.embeddings)) {
                stop("please specify PCA dims smaller than calculated")
        }
        cell_distance<- dist(object@reductions$pca@cell.embeddings[, dims])
        cell_cluster<- as.numeric(as.character(Idents(object)))
        silhouette_score<- cluster::silhouette(cell_cluster, cell_distance)
        silhouette_score<- tibble::tibble(cluster = silhouette_score[,1],
                                          width = silhouette_score[,3],
                                          cell = colnames(object)) %>%
                dplyr::mutate(cluster = as.factor(cluster))
        return(silhouette_score)
}


generate_gradient_colors <- function(low, high, num_colors) {
  color_ramp <- colorRampPalette(c(low, high))
  colors <- color_ramp(num_colors)
  return(colors)
}
