library(Seurat)
library(Signac)
require(rtracklayer)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(dplyr)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(ChIPseeker)
library(tidyr)

setwd('/hpcnfs/scratch/DIMA/piva/mpi_freiburg/scmultiome/') # path to repository


tss=import('./data/atac_peaks_associated_to_HVgenes_500KB_TSS.bed')

enh=read.csv('./data/hvgs2000_1st_linked_peak_enh_within20kb.csv' ) # only the first linked peak
head(enh)

# add peaks column
tss=as.data.frame(tss)
tss$peaks=paste(tss$seqnames, tss$start, tss$end, sep='-')


# =========== UMAP on ATAC =============

# data with linking only on WT
data=readRDS('./results/linked_wt_only.rds')
DefaultAssay(data) <- "peaks"

data@assays$peaks@links

# esclude yolk and undiff.cells
data=subset(data, annotation_leiden12 != 'yolk' )
data=subset(data, annotation_leiden12 != 'undifferentiated_cells') 


## select most accessible TSS
tss=tss[!duplicated(tss),]
avgacc= AverageExpression(data, assays = 'peaks', slot = 'data', features = tss$peaks, group.by = 'sample')
avgacc=as.data.frame(avgacc)

if(identical(tss$peaks,  row.names(avgacc))){
  tss$avgacc=avgacc$all
}


dim(tss)
head(tss)

mosta_tss=tss %>%
  arrange(name, desc(avgacc)) %>%
  group_by(name) %>%
  slice_head(n = 1)
head(mosta_tss)
dim(mosta_tss)


## ==== TSS ====

# to use most accessible tss:
tss=mosta_tss

tss_data=subset(data, features = tss$peaks)
tss_data
DefaultAssay(tss_data) <- "peaks"
tss_data <- RunTFIDF(tss_data)
tss_data <- FindTopFeatures(tss_data,
                            #min.cutoff = NULL
                            min.cutoff='q5', # to set the top 95% most common features as the VariableFeatures for the object.
                            ) 
tss_data@assays$peaks@meta.features
# Run singular value decomposition (LSI reduction, latent semantic indexing)
tss_data <- RunSVD(tss_data)
tss_data <- RunUMAP(tss_data, reduction = "lsi", 
                dims = 2:50,  # We exclude the first dimension as this is typically correlated with sequencing depth
                reduction.name = "umap.atac", reduction.key = "tss_atacUMAP_")

DimPlot(tss_data,reduction = "umap.atac", group.by='annotation_leiden12')+
  DimPlot(tss_data, reduction = "umap.atac", group.by='dataset')


## ==== Enhancers ====

dim(enh)
ehn_data=subset(data, features = enh$peaks)
ehn_data <- RunTFIDF(ehn_data)
ehn_data <- FindTopFeatures(ehn_data, 
                            min.cutoff='q5', # to set the top 95% most common features as the VariableFeatures for the object.
                          #  min.cutoff = NULL
                            )
# Run singular value decomposition (LSI reduction, latent semantic indexing)
ehn_data <- RunSVD(ehn_data)

ehn_data <- RunUMAP(ehn_data, reduction = "lsi", 
                dims = 2:50,  # We exclude the first dimension as this is typically correlated with sequencing depth
                reduction.name = "umap.atac", reduction.key = "enhancer_atacUMAP_")

DimPlot(ehn_data, reduction = "umap.atac", group.by='annotation_leiden12') +
  DimPlot(ehn_data, reduction = "umap.atac", group.by='dataset')

###

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

pdf('./figures/umap_tss_enh_inWT.pdf', paper='a4r')
DimPlot(ehn_data, reduction = "umap.atac", group.by='annotation_leiden12',cols =colori_germs, pt.size = 0.3, shuffle=TRUE)+
  DimPlot(tss_data, reduction = "umap.atac", group.by='annotation_leiden12', cols =colori_germs, pt.size = 0.3, shuffle=TRUE)+
  DimPlot(ehn_data, reduction = "umap.atac", group.by='dataset',pt.size = 0.3, shuffle=TRUE)+
  DimPlot(tss_data, reduction = "umap.atac", group.by='dataset', pt.size = 0.3, shuffle=TRUE)

dev.off()



pdf('./figures/umap_enh_inWT.pdf')
DimPlot(ehn_data, reduction = "umap.atac", group.by='annotation_leiden12',
        cols =colori_germs, pt.size = 0.3, shuffle=TRUE)
dev.off()

pdf('./figures/umap_tss_inWT.pdf')
DimPlot(tss_data, reduction = "umap.atac", group.by='annotation_leiden12',
        cols =colori_germs, pt.size = 0.3, shuffle=TRUE)
dev.off()


###
