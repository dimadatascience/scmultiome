library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(Seurat)
library(Signac)
require(rtracklayer)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(dplyr)
library(ChIPseeker)
library(tidyr)

setwd('/hpcnfs/scratch/DIMA/piva/mpi_freiburg/scmultiome/') # path to repository
data=readRDS('./results/linked_wt_only.rds')
DefaultAssay(data) <- "peaks"

# subset single cell data
data=subset(data, sample=='wt')
data=subset(data, annotation_leiden12 != 'yolk' )
data=subset(data, annotation_leiden12 != 'undifferentiated_cells') 

table(data$sample, data$annotation_leiden12)

custom_colors <- c(  'yolk' ='#efd129',
                     'mesoderm'= '#e24041',
                     'anterioposterior_ectoderm' ='#182953',
                     'dorsal_ectoderm' ='#277ea6',
                     'neuroectoderm'= '#3ac4e7',
                     'ventral_ectoderm' ='#ab93c6',
                     'anterior_endoderm'= '#7ac143',
                     'posterior_endoderm'= '#40733e',
                     'undifferentiated_cells'= '#b22987',
                     'pluripotent_progenitors'='#ff9900',
                     'grey')

# import peaks
file_path='./data/CRM_redfly.bed'
enh=import(file_path)

enh$peaks=GRangesToString(enh)


## _____ overlapping peaks _____

ov = findOverlaps(subject = GRanges(data@assays$peaks@ranges), 
                  query = enh)

overlap_peaks=GRanges(data@assays$peaks@ranges)[subjectHits(ov)]
overlap_peaks$name=enh[queryHits(ov)]$name


# remove TSSs from overlap peaks
all_tss=import('./data/allgenes_tss.bed')
all_tss=resize(all_tss, width = width(all_tss)+(400*2), fix = "center")

ov2 = findOverlaps(subject = overlap_peaks, 
                  query = all_tss)
overlap_peaks[subjectHits(ov2)] ## these are peaks overlapping with TSSs


overlap_peaks_no_tss=overlap_peaks[-subjectHits(ov2)]
overlap_peaks_no_tss


## UMAP
ehn_data=subset(data, features = GRangesToString(overlap_peaks_no_tss) )

ehn_data <- RunTFIDF(ehn_data)
ehn_data <- FindTopFeatures(ehn_data, 
                            min.cutoff=NULL, # NULL to use all the peaks in the assay
)

# Run singular value decomposition (LSI reduction, latent semantic indexing)
ehn_data <- RunSVD(ehn_data)

ehn_data <- RunUMAP(ehn_data, reduction = "lsi", 
                    dims = 2:50,  # exclude the first dimension as this is typically correlated with sequencing depth
                    reduction.name = "umap.atac", reduction.key = "enhancer_atacUMAP_")

DimPlot(ehn_data, reduction = "umap.atac", group.by='annotation_leiden12', cols = custom_colors, pt=0.8)
