setwd("/hpcnfs/scratch/temporary/DIMA/")
library(Seurat,  lib.loc = "/usr/local/lib/R/site-library")
library(destiny)  
library(SingleCellExperiment)
library(ggplot2)
library(Signac)
library(SCENT)

data = readRDS("scmultiome_mpi/scmultiome_230629.rds")
data@assays$ATAC@links


time=read.csv('/hpcnfs/scratch/DIMA/piva/mpi_freiburg/timepoint_from_ingestion.csv')
data@meta.data['inferred_timepoints']=time$time
data@meta.data['barcodes']=row.names(data@meta.data)


# New fragments paths
oldpath='/tungstenfs/scratch/shared/gchao_ggiorget/data_freiburg/multiome/'
newpath='/hpcnfs/scratch/temporary/DIMA/atac_fragments_freiburg/'
for (e in 1:12) {
  data@assays$ATAC@fragments[[e]]@path=gsub(oldpath,newpath, data@assays$ATAC@fragments[[e]]@path )
  print(e)
  print(data@assays$ATAC@fragments[[e]])
  print(data@assays$ATAC@fragments[[e]]@path )
}


# =========== UMAP on ATAC =============

DefaultAssay(data) <- "ATAC"
# Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = "q0")

# Run singular value decomposition (LSI reduction, latent semantic indexing)
data <- RunSVD(data, n = 80)

data <- RunUMAP(data, reduction = "lsi", 
                dims = 2:30,  # We exclude the first dimension as this is typically correlated with sequencing depth
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")

DimPlot(data, reduction = "umap.atac")
DimPlot(data, reduction = "umap.atac", group.by='sample')+
  DimPlot(data, reduction = "umap.atac", group.by='dataset')+
    DimPlot(data, reduction = "umap.atac", group.by='annotation_leiden12')


library(scater)
sce=as.SingleCellExperiment(data)
c=c( "nCount_RNA","nFeature_RNA",'nCount_ATAC','nFeature_ATAC', 
     "percent.mt","dataset" ,"sample"  ,"annotation_leiden12" )
expl_pca=scater::getExplanatoryPCs(sce, dimred ='LSI',variables = c)
plotExplanatoryPCs(expl_pca/100, title='LSI')


data <- RunUMAP(data, reduction = "lsi", 
                dims = 2:30,  # We exclude the first dimension as this is typically correlated with sequencing depth
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sce=as.SingleCellExperiment(data)
c=c( "nCount_RNA","nFeature_RNA",'nCount_ATAC','nFeature_ATAC', 
     "percent.mt","dataset" ,"sample"  ,"annotation_leiden12" )
expl_pca=scater::getExplanatoryPCs(sce, dimred ='LSI',variables = c)
plotExplanatoryPCs(expl_pca/100, title='LSI')

DimPlot(data, group.by = 'dataset', pt.size = 0.1)
DimPlot(data, group.by = 'sample', pt.size = 0.3)
DimPlot(data, split.by  = 'sample', pt.size = 0.3)
DimPlot(data, group.by = 'annotation_leiden12', pt.size = 0.3)



# Retrieve samples color palette
require(scales)
my_color_palette <- hue_pal()(3)
DimPlot(object = data,group.by = "sample" , do.return = T) + 
  scale_color_manual(values = my_color_palette)
my_color_palette
#"#F8766D" "#00BA38" "#619CFF"

data@meta.data$sample=gsub("nej","CBP",data@meta.data$sample)

colori_samples=c("wt"="#00BA38" , "ez"= "#619CFF", "nej"="#F8766D")


colori_germs= c('yolk' ='#efd129',
'mesoderm'= '#e24041',
'anterioposterior_ectoderm' ='#182953',
'dorsal_ectoderm' ='#277ea6',
'neuroectoderm'= '#3ac4e7',
'ventral_ectoderm' ='#ab93c6',
'anterior_endoderm'= '#7ac143',
'posterior_endoderm'= '#40733e',
'undifferentiated_cells'= '#b22987',
'pluripotent_progenitors'='#ff9900')


DimPlot(data, reduction = "umap",group.by = "ident", cols =colori_germs,
        label = FALSE, pt.size=0.8)
DimPlot(data, reduction = "umap",group.by = "inferred_timepoints",
        label = FALSE, pt.size=1)


p1 <- DimPlot(data, reduction = "umap",group.by = "ident", cols =colori_germs,label = FALSE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(data, reduction = "umap.atac",group.by = "ident", cols =colori_germs,label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2


DimPlot(data, reduction = "umap",group.by = "ident",cols =colori_germs, # label = TRUE, 
        cols =colori_germs )  + ggtitle("RNA")
DimPlot(subset(data, sample=='wt'), reduction = "umap",group.by = "ident", cols =colori_germs, # label = TRUE, 
        cols =colori_germs )  + ggtitle("RNA") #+ NoLegend()
DimPlot(subset(data, sample=='wt'), reduction = "umap",group.by = "ident",# label = TRUE,
        cols =colori_germs ) + ggtitle("RNA")

p1wt <- DimPlot(subset(data, sample=='wt'), reduction = "umap",group.by = "ident", label = FALSE,cols =colori_germs) + NoLegend() + ggtitle("RNA")
p2wt <- DimPlot(subset(data, sample=='wt'), reduction = "umap.atac",group.by = "ident", label = FALSE, cols =colori_germs) + NoLegend() + ggtitle("ATAC")
p1wt|p2wt


## 




## integration RNA + ATAC 
# FindMultiModalNeighbors function constructs a weighted nearest neighbor (WNN) graph. 
# For each cell, we identify the nearest neighbors based on a weighted combination of two modalities, 
# one computed for each assay.
data <- FindMultiModalNeighbors(data, 
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:30, 2:30))
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

p3 <- DimPlot(data,reduction = "wnn.umap", label = FALSE,cols =colori_germs, repel = TRUE) + ggtitle("WNN (RNA+ATAC)")

p1 + p2 + p3

p4 <- DimPlot(data, reduction = "umap",group.by = "sample",cols=colori_samples, label = FALSE) + NoLegend() + ggtitle("RNA")
p5 <- DimPlot(data, reduction = "umap.atac",group.by = "sample",cols=colori_samples, label = FALSE) + NoLegend() + ggtitle("ATAC")
p6 <- DimPlot(data,reduction = "wnn.umap",group.by = "sample", cols=colori_samples,label = FALSE, repel = TRUE) + ggtitle("WNN (RNA+ATAC)")

p4 + p5 + p6

p3wt <- DimPlot(subset(data, sample=='wt'), reduction = "wnn.umap",group.by = "ident", label = FALSE,
                cols =colori_germs) + ggtitle("WNN (RNA+ATAC)")
p1wt|p2wt|p3wt


p1wtl <- DimPlot(subset(data, sample=='wt'), reduction = "umap",group.by = "leiden_res12", label = TRUE) + NoLegend() + ggtitle("RNA")
p2wtl <- DimPlot(subset(data, sample=='wt'), reduction = "umap.atac",group.by = "leiden_res12", label = TRUE) + NoLegend() + ggtitle("ATAC")
p3wtl <- DimPlot(subset(data, sample=='wt'), reduction = "wnn.umap",group.by = "leiden_res12", label = TRUE)
+ ggtitle("WNN (RNA+ATAC)")+ NoLegend() # guides(color = guide_legend(override.aes = list(size=4), ncol=3) )

p1wt+p2wt+p3wt+p1wtl+p2wtl+p3wtl

p1 + p2 + p3 + p4 + p5 + p6

str(p1$theme)

p7 <- DimPlot(data, reduction = "umap",group.by = "leiden_res12", 
              label = FALSE) + NoLegend() + ggtitle("RNA")
p8 <- DimPlot(data, reduction = "umap.atac",group.by = "leiden_res12", 
              label = FALSE) + NoLegend() + ggtitle("ATAC")
p9 <- DimPlot(data,reduction = "wnn.umap",group.by = "leiden_res12",
              label = FALSE, repel = TRUE) + ggtitle("WNN (RNA+ATAC)") +guides(color = guide_legend(override.aes = list(size=4), ncol=3) ) #+ NoLegend()
p7 + p8 + p9

p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9


DimPlot(data, reduction = "umap.atac", group.by='sample') + DimPlot(data, reduction = "umap.atac", group.by='annotation_leiden12')+
  DimPlot(data, reduction = "umap.atac", group.by='dataset')

## clustering of atac assay
DefaultAssay(data)='ATAC'

data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 1.2 )
data=DimPlot(data,reduction = "umap.atac",group.by = "leiden_res12")









