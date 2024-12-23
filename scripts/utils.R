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
#         scale_x_discrete(breaks=seq(-x_lim,x_lim,1))+
        theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
                plot.title = element_text(size = 21,face = "bold"),
                axis.text=element_text(size=12),
                axis.title=element_text(size=14),
                legend.text=element_text(size=12),legend.title=element_blank())
}







dh_meta_analysis = function(tabList) {
    # Assign sequential names to the elements of tabList if they are not named
    if (is.null(names(tabList))) {
        names(tabList) = as.character(seq(length(tabList)))
    }

    # Add a new column 'Dataset' to each data frame in tabList indicating the dataset name
    for (key in names(tabList)) {
        tabList[[key]]$Dataset = key
    }

    # Combine all data frames in tabList into a single data frame
    df = do.call(rbind, tabList)

    # Perform meta-analysis for each matching gene and assay
    # Calculate standard error (se) from logFC and t-statistics using the formula t = logFC / se
    df %>%
        as_tibble() %>%
        group_by(assay) %>%
        do(tidy(rma(yi = logFC, sei = logFC / t, data = ., method = "FE"))) %>%  # Fixed-effects meta-analysis
        select(-term, -type) %>%  # Remove unnecessary columns
        mutate(FDR = p.adjust(p.value, "fdr")) %>%  # Adjust p-values for false discovery rate
        mutate('log10FDR' = -log10(FDR))  # Calculate -log10 of the FDR
}

dh_plotTree = function(tree, low = "grey90", mid = "red", high = "darkred", xmax.scale = 1.5) {
    # Plot a phylogenetic tree with tips labeled and no branch lengths
    fig = ggtree(tree, branch.length = "none") +
        geom_tiplab(color = "black", size = 4, hjust = 0, offset = 0.2) +  # Add tip labels
        theme(legend.position = "top left", plot.title = element_text(hjust = 0.5))  # Customize plot theme

    # Get the default maximum value of the x-axis
    xmax = layer_scales(fig)$x$range$range[2]

    # Adjust the x-axis scale to increase its width
    fig + xlim(0, xmax * xmax.scale)
}

dh_plotCoef = function(tab, coef, fig.tree, low = "grey90", mid = "red", high = "darkred", ylab) {
    # Prepare data for plotting
    tab$logFC = tab$estimate  # Rename 'estimate' column to 'logFC'
    tab$celltype = factor(tab$assay, rev(ggtree::get_taxa_name(fig.tree)))  # Create a factor for cell types
    tab$se = tab$std.error  # Rename 'std.error' column to 'se'

    # Create a coefficient plot
    fig.es = ggplot(tab, aes(celltype, logFC)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +  # Add a horizontal reference line
        geom_errorbar(aes(ymin = logFC - 1.96 * se, ymax = logFC + 1.96 * se), width = 0) +  # Add error bars
        geom_point2(aes(color = pmin(4, -log10(FDR)), size = pmin(4, -log10(FDR)))) +  # Add points with size and color based on FDR
        scale_color_gradient2(name = bquote(-log[10]~FDR),
                              limits = c(0, 4), low = low, mid = mid, high = high, midpoint = -log10(0.01)) +  # Set color gradient
        scale_size_area(name = bquote(-log[10]~FDR), limits = c(0, 4)) +  # Scale point size
        geom_text2(aes(label = '+', subset = FDR < 0.05), color = "white", size = 6, vjust = 0.4, hjust = 0.5) +  # Add labels for significant points
        theme_classic() +  # Use classic theme
        coord_flip() +  # Flip axes
        xlab('') +
        ylab(ylab) +  # Set y-axis label
        theme(axis.text.y = element_blank(),
              axis.text = element_text(size = 12),
              axis.ticks.y = element_blank(),
              text = element_text(size = 20)) +  # Customize axis and text styles
        scale_y_continuous(breaks = scales::breaks_pretty(3))  # Adjust y-axis breaks

    return(fig.es)  # Return the coefficient plot
}




my_DotPlot <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "red4"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA
    ) {
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
    scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
    )
    feature.groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
        X = 1:length(features),
        FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
        }
    ))
    if (any(is.na(x = feature.groups))) {
        warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
        )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
    }
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

    data.features <- FetchData(object = object, vars = features, cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
    } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
        if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
            return(mean(x = expm1(x = x)))
        }
        )
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
    )
    names(x = data.plot) <- unique(x = data.features$id)
    if (cluster.idents) {
    mat <- do.call(
        what = rbind,
        args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
    }
    data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    }
    )
    data.plot <- do.call(what = 'rbind', args = data.plot)
    if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    ngroup <- length(x = levels(x = data.plot$id))
    if (ngroup == 1) {
    scale <- FALSE
    warning(
        "Only one identity present, the expression values will be not scaled",
        call. = FALSE,
        immediate. = TRUE
    )
    } else if (ngroup < 5 & scale) {
    warning(
        "Scaling data with a low number of groups may produce misleading results",
        call. = FALSE,
        immediate. = TRUE
    )
    }
    avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
        data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
        if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
        } else {
        data.use <- log1p(x = data.use)
        }
        return(data.use)
    }
    )
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
    )
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (split.colors) {
    splits.use <- vapply(
        X = as.character(x = data.plot$id),
        FUN = gsub,
        FUN.VALUE = character(length = 1L),
        pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
        ),
        replacement = '',
        USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
        FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
        },
        color = cols[splits.use],
        value = avg.exp.scaled
    )
    }
    color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
    if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
    }
    if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
    }
    if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
        x = feature.groups[data.plot$features.plot],
        levels = unique(x = feature.groups)
    )
    }
    plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(shape = 21, colour = "black",
                mapping = aes_string(size = 'pct.exp', fill = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
        x = 'Features',
        y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) +
    theme_cowplot()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+RotatedAxis() +
        scale_fill_gradient(low = cols[1], high = cols[2])
    if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
        facets = ~feature.groups,
        scales = "free_x",
        space = "free_x",
        switch = "y"
    ) + theme(
        panel.spacing = unit(x = 1, units = "lines"),
        strip.background = element_blank()
    )
    }
    if (split.colors) {
    plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
    } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
    }
    return(plot)
}





plot_SANKEY<-function(Metadata,ident.1,ident.2,method='Jaccard'){
    # Libraries
    require(tidyverse)
    require(viridis)
    require(patchwork)
    require(hrbrthemes)
    require(circlize)
    # Load dataset from github
    # Package
    require(networkD3)
    # I need a long format
    if(method=='Jaccard'){
        data <- as.data.frame.matrix(PairWiseJaccardSets(Metadata,ident.1=ident.1, ident.2=ident.2))
    }else if(method=='prop'){
        data <- as.data.frame.matrix(prop.table(table(Metadata[[ident.1]],Metadata[[ident.2]])))
    }else{
        print('Please select: Jaccard or prop')
    }    

    data_long <- data %>%
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname) %>%
    filter(value > 0)
    colnames(data_long) <- c("source", "target", "value")
    data_long$target <- paste(data_long$target, " ", sep="")

    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
    
    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    data_long$IDsource=match(data_long$source, nodes$name)-1 
    data_long$IDtarget=match(data_long$target, nodes$name)-1

    # prepare colour scale
    ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

    # Make the Network
    sankeyNetwork(Links = data_long, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)
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
        # or as.integer
        cell_cluster<- as.numeric(as.character(Idents(object)))
        silhouette_score<- cluster::silhouette(cell_cluster, cell_distance)
        silhouette_score<- tibble::tibble(cluster = silhouette_score[,1],
                                          width = silhouette_score[,3],
                                          cell = colnames(object)) %>%
                dplyr::mutate(cluster = as.factor(cluster))
        return(silhouette_score)
}


# Function to generate gradient colors
generate_gradient_colors <- function(low, high, num_colors) {
  color_ramp <- colorRampPalette(c(low, high))
  colors <- color_ramp(num_colors)
  return(colors)
}



plot_JaccardResults<-function(jaccard_results,flip=F,
    cols=c("white","purple4"),num_colors=10){

    gradient_colors <- generate_gradient_colors(low=cols[1], 
            high=cols[2], 
            num_colors=num_colors)
    if(flip){
        jaccard_results <- t(jaccard_results)
    }
    pheatmap(jaccard_results,
        scale = FALSE,
        cellwidth = 20,
        cellheight = 20,
        color=gradient_colors,fontsize = 12)

}






Dreamlet_Enrichement <- function(Dreamlet_obj, GMT_file, 
    pvalue_thres=0.05,
    lfc_thres=0.1,
    coef="Diff_B_vs_G",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500
    ){

    require(clusterProfiler)
    require(parallel)
    
    # ============== Read All results ===================
    DEGs_file_all <- as.data.frame(topTable(Dreamlet_obj,
        p.value = 1,lfc = 0,number = 100000000,
        coef = "Diff_B_vs_G" ))
    background_genes <- unique(DEGs_file_all$ID)

    # ============== Read significant results ===================
    DEGs_file_signi <- as.data.frame(topTable(Dreamlet_obj,
        p.value = pvalue_thres,
        lfc = lfc_thres,
        number = 100000000,
        coef = "Diff_B_vs_G" ))

    DEGs_file_signi <- DEGs_file_signi[order(DEGs_file_signi$adj.P.Val),]
    colnames(DEGs_file_signi)[2]<- "gene.name"    
    DEGs_file_signi$cluster <- ifelse(DEGs_file_signi$logFC>0,'Bad','Good')
    DEGs_file_signi$cluster <- paste0(DEGs_file_signi$cluster,"_",DEGs_file_signi$assay)
    DEGs_file_signi$gene < DEGs_file_signi$gene.name


    # ============== Read GMT File ===================
    Loaded_GMT_file <- read.gmt(GMT_file)
    
    # ============== Perform Enrichement for all clusters  ===================
    Unique_Clusters <- unique(DEGs_file_signi$cluster)
    
    res_gmt_enrich_list <- mclapply(Unique_Clusters,
        function(cl_name,background_genes,DEGs_file,Loaded_GMT_file,
        qvalue_thres_enrich,pvalueCutoff,qvalueCutoff,minGSSize,maxGSSize){
        
        # ==== Subset specific cluster
        cl_degs <- DEGs_file[DEGs_file$cluster==cl_name,]
        # ==== Order LogFC
        ord_cl_degs <- cl_degs[order(cl_degs$logFC,decreasing=T),]
        gene_ids <- ord_cl_degs$gene.name
        
        # === Run ClusterProfiler
        cl_res <- clusterProfiler::enricher( gene_ids, 
            TERM2GENE= Loaded_GMT_file,
            pAdjustMethod="BH",
            universe=background_genes,
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize)
        
        df_cl_res <- as.data.frame(cl_res@result)
        df_cl_res$Cluster <- cl_name
        return(df_cl_res)
    },background_genes=background_genes,
        DEGs_file=DEGs_file_signi,
        Loaded_GMT_file=Loaded_GMT_file,
        qvalue_thres_enrich=qvalue_thres_enrich,
        pvalueCutoff=pvalueCutoff,
        qvalueCutoff=qvalueCutoff,
        minGSSize=minGSSize,
        maxGSSize=maxGSSize)

    
    # === Merge Cluster enrichement results
    df <- do.call("rbind", res_gmt_enrich_list)
    # === Add Columns  PathwaySize, GeneRatio_num, BgRatio_num, Comb_Ratios, Log10Qval
    df$PathwaySize <- unlist(lapply(df$GeneRatio,function(x){unlist(str_split(x,pattern='/'))[2]}))
    df$GeneRatio_num <- unlist(lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}))
    df$BgRatio_num <- unlist(lapply(df$BgRatio,function(x){as.numeric(eval(parse(text = x)))}))
    df$Comb_Ratios <- df$GeneRatio_num / df$BgRatio_num
    df$Log10Qval <- -log10(df$qvalue)
    df$Significant <- ""
    df$Significant[df$qvalue< 0.05] <- "*"
    return(df)
}    


plot_Dreamlet_Enrichement<-function(Dreamlet_Enrichement_result,
    GMT_name,
    n_terms=1,
    wt="Comb_Ratios",
    qvalue_thres=0.05,
    cols=c('white','red4')){
    
    Dreamlet_Enrichement_result$ID <- sub("^GOBP_|^GOMF_|^HALLMARK_|^KEGG_|^REACTOME_", "", Dreamlet_Enrichement_result$ID)
    Dreamlet_Enrichement_result$Description <- sub("^GOBP_|^GOMF_|^HALLMARK_|^KEGG_|^REACTOME_", "", Dreamlet_Enrichement_result$Description)

    df_sign <- as.data.frame(Dreamlet_Enrichement_result[Dreamlet_Enrichement_result$qvalue<0.05,])
    if(wt=="Comb_Ratios"){
        top10 <- df_sign %>% group_by(Cluster) %>% top_n(n = n_terms, wt = Comb_Ratios)
    }else{
        top10 <- df_sign %>% group_by(Cluster) %>% top_n(n = n_terms, wt = Log10Qval)
    }
    df_to_plot <- Dreamlet_Enrichement_result[Dreamlet_Enrichement_result$ID %in%unique(top10$Description),]
    df_to_plot$Description <- factor(df_to_plot$Description,levels=unique(top10$Description))
    
    
    p <- ggplot(df_to_plot,aes(x=Cluster,
                      y=Description,
                      fill = -log10(qvalue),
                      label=Significant))+
        geom_tile()+
        theme_cowplot() +
        geom_text(alpha=1,size=10)+
        RotatedAxis()+
        ggtitle(GMT_name)+
        theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18))+
        coord_flip()+
         theme(plot.title = element_text(hjust = 0.5))+
        scale_fill_gradient(low = cols[1], high = cols[2])+
        coord_equal()

    return(p)
}








# ================================= ENRICHEMENT FUNCTIONS =====================================
Fun_DEGs_GMT_Enrichment <- function(cl_name,DEGs_file,ngenes=50){
    cl_degs <- DEGs_file[DEGs_file$cluster==cl_name,]
    ord_cl_degs <- cl_degs[order(cl_degs$gene.score,decreasing=T),]
    ids <- ord_cl_degs$gene[1:ngenes]
    background<-rownames(DEGs_file)
    cl_res <- clusterProfiler::enricher(ids, TERM2GENE= loaded_gmt_file,
            pAdjustMethod="BH",
            universe=background,
            pvalueCutoff = 0.05,qvalueCutoff = 0.05,
            minGSSize = 5,maxGSSize = 500)
    df_cl_res <- as.data.frame(cl_res@result)
    df_cl_res$Cluster <- cl_name
    return(df_cl_res)
}

GMT_Enrichment_Dreamlet <- function(cl_name,DEGs_file,
                pval_thres,
                lfc_thres,
                ngenes=50){
    Vlcano_df <- as.data.frame(topTable(res.dl_all_contrast_bin,p.value = pval_thres,
                        lfc = fc_thres,
                        number = 100000000,
                        coef = "Diff_B_vs_G" ))
    Vlcano_df <- Vlcano_df[order(Vlcano_df$adj.P.Val),]
    colnames(Vlcano_df)[2]<- "gene.name"
    
    # add a column of CLUSTER
    Vlcano_df$cluster <- if(Vlcano_df$logFC>1) 'Bad' else 'Good'



    cl_degs <- DEGs_file[DEGs_file$cluster==cl_name,]
    ord_cl_degs <- cl_degs[order(cl_degs$gene.score,decreasing=T),]
    ids <- ord_cl_degs$gene[1:ngenes]
    background<-rownames(DEGs_file)
    cl_res <- clusterProfiler::enricher(ids, TERM2GENE= loaded_gmt_file,
            pAdjustMethod="BH",
            universe=background,
            pvalueCutoff = 0.05,qvalueCutoff = 0.05,
            minGSSize = 5,maxGSSize = 500)
    df_cl_res <- as.data.frame(cl_res@result)
    df_cl_res$Cluster <- cl_name
    return(df_cl_res)
}
# -----------------------------------------------------------------------



# ================================= MAP SYMBOLS TO ENTREZ =====================================
entrezMapper <- function(DEGenes, geneIDs, species="mmu") {  
    require(tidyverse)

    # Kegg convention
    if (species == "mmu") {
        OrgDb <- org.Mm.eg.db
        require(org.Mm.eg.db)
    } else if (species == "hsa") {
        OrgDb <- org.Hs.eg.db
        require(org.Hs.eg.db)
    }
    keytypes(OrgDb) # UNIPROT, ALIAS etc
    
    idsFromENSEMBL <- bitr(geneIDs$ensembl.id,
                            fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL", "GENENAME"),
                            OrgDb=OrgDb, drop = F)
    
    idsFromSYMBOL <- bitr(geneIDs$gene.name,
                            fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL", "GENENAME"),
                            OrgDb=OrgDb, drop = F)

    entrezID1 <- list()
    entrezID2 <- list()
    entrezID <- list()

    mismatched <- 0
    # ENSEMBL FOR
    for (symbol_g in DEGenes$gene.name ) {
        inds <- which(idsFromENSEMBL$ENSEMBL == DEGenes$ensembl.id[DEGenes$gene.name == symbol_g])
        entrezID1[[symbol_g]] <- idsFromENSEMBL$ENTREZID[inds]
        if (length(entrezID1[[symbol_g]]) == 0){
            entrezID1[[symbol_g]] <- NA
        }
        entrezID[[symbol_g]] <- entrezID1[[symbol_g]]
    }
    
    # SYMBOL FOR
    for (symbol_g in DEGenes$gene.name ) {
        inds <- which(idsFromSYMBOL$SYMBOL == DEGenes$gene.name[DEGenes$gene.name == symbol_g])
        entrezID2[[symbol_g]] <- idsFromSYMBOL$ENTREZID[inds]
        if (length(entrezID2[[symbol_g]]) == 0){
            entrezID2[[symbol_g]] <- NA
        }
        # print(symbol_g)
        # print(entrezID1[[symbol_g]])
        # print(entrezID2[[symbol_g]])
        if(length(entrezID2[[symbol_g]]) >1 ){
            entrezID2[[symbol_g]] <- entrezID2[[symbol_g]][1]
        }
        if(length(entrezID1[[symbol_g]]) >1 ){
            entrezID1[[symbol_g]] <- entrezID1[[symbol_g]][1]
        }

        if ( !is.na(entrezID1[[symbol_g]]) && !is.na(entrezID2[[symbol_g]]) ) {
            if (entrezID1[[symbol_g]] != entrezID2[[symbol_g]]) {
                # print("WARNING: entrezID mismatch!!!!!!!!")
                mismatched <- mismatched + 1
                entrezID[[symbol_g]] <- c(entrezID2[[symbol_g]])
            }
        }


        if (all(is.na(entrezID1[[symbol_g]]))) {
            entrezID[[symbol_g]] <- entrezID2[[symbol_g]]
        }
        if(length(entrezID[[symbol_g]])==2){
            entrezID[[symbol_g]] <- entrezID[[symbol_g]][1]
        }
    }

    
    print(table(is.na(entrezID)))
    DEGenes$entrezID <- as.character(entrezID)
    return(DEGenes)
}