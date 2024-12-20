library(GEOquery)

library(Seurat)
library(tibble)
library(magrittr)
library(dittoSeq)
library(patchwork)
library(harmony)

setwd('External_Data/')
save_Fig_dir <- 'SFigures/SFigure4/'



# Read in data
Metadata <- read.csv(file ="GSE189432_annotations.csv.gz")
Metadata <- Metadata[Metadata$sample %in% c('ctrl_cns','stroke_cns_24h_1','stroke_cns_24h_2','stroke_cns_72h'),]
unique(Metadata$sample)
dim(Metadata)


# ============================== READ CONTROL ============================
CTRL <- CreateSeuratObject(ReadMtx(
  mtx = "ctrl_cns/GSM5701739_ctrl_cns_matrix.mtx.gz/",
  features = "ctrl_cns/GSM5701739_ctrl_cns_features.tsv.gz",
  cells = "ctrl_cns/GSM5701739_ctrl_cns_barcodes.tsv.gz"
))
CTRL$orig.ident <- 'ctrl_cns'
my_vector <- rownames(CTRL@meta.data)
# Replace all occurrences of "-1" with a specific string
replacement_string <- "_ctrl_cns"
my_vector <- gsub("-1", replacement_string, my_vector)
CTRL <- RenameCells(CTRL, new.names =my_vector )
# -----------------------------------------------------------------------



# ============================== READ S24_1 ============================
S24_1 <- CreateSeuratObject(ReadMtx(
  mtx = "stroke_24_cns/GSM5701742_stroke_cns_24h_1_matrix.mtx.gz",
  features = "stroke_24_cns/GSM5701742_stroke_cns_24h_1_features.tsv.gz",
  cells = "stroke_24_cns/GSM5701742_stroke_cns_24h_1_barcodes.tsv.gz",
))
S24_1$orig.ident <- 'stroke_cns_24h_1'
my_vector <- rownames(S24_1@meta.data)
# Replace all occurrences of "-1" with a specific string
replacement_string <- "_stroke_cns_24h_1"
my_vector <- gsub("-1", replacement_string, my_vector)
S24_1 <- RenameCells(S24_1, new.names =my_vector )
# -----------------------------------------------------------------------


# ============================== READ S24_2 ============================
S24_2 <- CreateSeuratObject(ReadMtx(
  mtx = "stroke_24_cns_2/GSM5701743_stroke_cns_24h_2_matrix.mtx.gz",
  features = "stroke_24_cns_2/GSM5701743_stroke_cns_24h_2_features.tsv.gz",
  cells = "stroke_24_cns_2/GSM5701743_stroke_cns_24h_2_barcodes.tsv.gz"
))
S24_2$orig.ident <- 'stroke_cns_24h_2'
my_vector <- rownames(S24_2@meta.data)
# Replace all occurrences of "-1" with a specific string
replacement_string <- "_stroke_cns_24h_2"
my_vector <- gsub("-1", replacement_string, my_vector)
S24_2 <- RenameCells(S24_2, new.names =my_vector )
# -----------------------------------------------------------------------


# ============================== READ S72h ============================
S72 <- CreateSeuratObject(ReadMtx(
  mtx = "stroke_72_cns/GSM5701746_stroke_cns_72h_matrix.mtx.gz",
  features = "stroke_72_cns/GSM5701746_stroke_cns_72h_features.tsv.gz",
  cells = "stroke_72_cns/GSM5701746_stroke_cns_72h_barcodes.tsv.gz"
))
S72$orig.ident <- 'stroke_cns_72h'
my_vector <- rownames(S72@meta.data)
# Replace all occurrences of "-1" with a specific string
replacement_string <- "_stroke_cns_72h"
my_vector <- gsub("-1", replacement_string, my_vector)
S72 <- RenameCells(S72, new.names =my_vector )
# -----------------------------------------------------------------------



table(Metadata$sample)

#add.cell.ids = c("ctrl_cns", "stroke_cns_24h_1","stroke_cns_24h_2","stroke_cns_72h"),
#MICE <- merge(x = CTRL,y = c(S24_1,S24_2,S72,PIA,DURA,stroke_dura,stroke_pia), project = "Mice")
MICE <- merge(x = CTRL,y = c(S24_1,S24_2,S72), project = "Mice")

head(Metadata)
head(MICE@meta.data)
MICE$barcode <- rownames(MICE@meta.data)
MICE <- subset(MICE,subset= barcode %in% Metadata$barcode)
MICE$UMAP_1 <- Metadata$UMAP_1
MICE$UMAP_2 <- Metadata$UMAP_2
MICE$cluster <- Metadata$cluster
MICE$sample <- Metadata$sample
MICE[["RNA"]] <- JoinLayers(MICE[["RNA"]])
counts <- MICE@assays$RNA$counts
rownames(counts) <- toupper(rownames(MICE))
metadata <- MICE@meta.data



MICE <- CreateSeuratObject(counts =counts,meta.data = metadata )
umap <- Metadata[,c('UMAP_1','UMAP_2')]
rownames(umap) <- Metadata$barcode
MICE[["ref.umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(MICE))
all.genes <- rownames(MICE)
# MICE <- ScaleData(MICE, features = all.genes)

ifnb <- MICE
# ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$sample)
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
set.seed(120120224)
set.seed(24022012)

ifnb <- RunHarmony(ifnb, group.by.vars = "sample",dims.use = 1:40)
ifnb <- RunUMAP(ifnb, dims = 1:40, reduction = "harmony", min.dist = 0.3,seed.use =24022012 ,
                reduction.name = "umap.harmony", reduction.key = "Uh_", return.model = T)

ifnb
DimPlot(ifnb,group.by = 'cluster',reduction = 'umap.harmony',label = T,pt.size = 0.0001)

sort(unique(ifnb$cluster))
colors_paper <-c(
  "#FFC312",
  "#C4E538",
  "#12CBC4",
  "#FDA7DF",
  "#ED4C67",
  "#F79F1F",
  "#A3CB38",
  "#1289A7",
  "#D980FA",
  "#B53471",
  "#EE5A24",
  "#009432",
  "#0652DD",
  "#9980FA",
  "#833471",
  "#EA2027",
  "#006266",
  "#1B1464",
  "#5758BB",
  "#6F1E51",
  "#40407A"
)
names(colors_paper) <- c("Micro_1",
                         "Micro_2","Micro_3",
                         "stress_Micro","CAM_1","CAM_2",
                         "SAMC","Macro_1","Macro_2","stress_Myeloid","mDC1",'mDC2',
                         'Granulo_1','Granulo_2',
                         'Mast','prolif_cells',
                         'Bc','gdTc','ILC2','Tc','NK')
DimPlot(ifnb,group.by = 'cluster',reduction = 'umap.harmony',pt.size = 0.0001,cols=colors_paper)
dittoBarPlot(ifnb,group.by = 'sample',var = 'cluster',color.panel =colors_paper  )


# saveCRDS(ifnb, 'MCAO_brain_parenchymal.rds.ztsd')

source('utils.R')
setwd('External_Data/')
Metacells <- readCRDS('2024_03_28_Myeloid_Metacells_Subclass_ADAM.rds.ztsd')
Idents(Metacells) <- 'Clusters'


anchors <- FindTransferAnchors(
  reference = ifnb,
  query = Metacells,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:40
)


Metacells <- MapQuery(anchorset = anchors,
                      query = Metacells,
                      reference = ifnb,
                      refdata = list(celltype = "cluster"),
                      reference.reduction = "harmony",
                      reduction.model = "umap.harmony")

Metacells
DimPlot(Metacells, reduction = "ref.umap", group.by = "Clusters")
DimPlot(Metacells, reduction = "ref.umap", group.by = "predicted.celltype")+
  DimPlot(ifnb,group.by = 'cluster',reduction = 'umap.harmony',label = T)

DimPlot(Metacells, reduction = "umap", group.by = "predicted.celltype",cols = colors_paper)



DimPlot(Metacells, reduction = "umap", group.by = "predicted.celltype",cols = colors_paper)+
  dittoBarPlot(Metacells,group.by = 'Clusters',var = 'predicted.celltype',color.panel = colors_paper)+coord_flip()+
  plot_layout(widths=c(6,3))

p1 <- DimPlot(ifnb,group.by = 'cluster',reduction = 'umap.harmony',label.box = F,pt.size = 0.0001,cols=colors_paper,label=F)+xlim(-15,10)+ylim(-15,10)+ggtitle('Reference')
p2 <- DimPlot(Metacells, reduction = "umap", group.by = "predicted.celltype",cols = colors_paper)+NoLegend()


pdf(paste0(save_Fig_dir,'SupFig.4.pdf'),width = 18,height=6)
p1 + p2 + dittoBarPlot(Metacells,group.by = 'Clusters',var = 'predicted.celltype',color.panel = colors_paper)+coord_flip()+NoLegend()+
plot_layout(widths=c(6,6,3),guides = 'collect')
dev.off()

graphics.off()
png(paste0(save_Fig_dir,'SupFig.4.png'),width = 18,height=6)
p1 + p2 + dittoBarPlot(Metacells,group.by = 'Clusters',var = 'predicted.celltype',color.panel = colors_paper)+coord_flip()+NoLegend()+
plot_layout(widths=c(6,6,3),guides = 'collect')
dev.off()
