setwd("/hpcscratch/ieo/DIMA/scmultiome_mpi/")

library(Seurat,  lib.loc = "/usr/local/lib/R/site-library")
library(Signac)
require(rtracklayer)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(dplyr)

# wt only data - linking done on wt only cells
datawt=readRDS('/hpcscratch/ieo/DIMA/scmultiome_mpi/08_other/wt_complete/linked_wt_only.rds')
ranges_wt=datawt@assays$peaks@ranges

data = readRDS("/hpcscratch/ieo/DIMA/scmultiome_mpi/scmultiome_230629.rds")
oldpath='/tungstenfs/scratch/shared/gchao_ggiorget/data_freiburg/multiome/'
newpath='/hpcnfs/scratch/temporary/piva/atac_fragments_freiburg/'
for (e in 1:12) {
  data@assays$ATAC@fragments[[e]]@path=gsub(oldpath,newpath, data@assays$ATAC@fragments[[e]]@path )
}

# mutant object
srat.sub <-
  data %>%
  subset(sample != 'wt')

table(srat.sub@meta.data$sample)

DefaultAssay(srat.sub) <- "ATAC"
macs2_counts <- FeatureMatrix(
  fragments = Fragments(srat.sub),
  features = ranges_wt, # peaks from wt only obj
  cells = Cells(srat.sub) 
)

## add annotation
annotations= datawt@assays$peaks@annotation
# Create a new assay using the MACS2 peak set
srat.sub[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
                                            annotation = annotations)

DefaultAssay(srat.sub) <- "peaks"

## combine wt and mutants
combined <- merge(
  x = datawt,
  y = srat.sub,
  add.cell.ids = c("wt", "mutants")
)

combined@assays$peaks
saveRDS(combined, '/hpcscratch/ieo/DIMA/scmultiome_mpi/wt_only_linked_with_mutants.rds')
max(combined@assays$peaks@data)
sce=as.SingleCellExperiment(combined)
max(sce@assays@data$counts)
library(zellkonverter)
writeH5AD( sce, file='/hpcnfs/scratch/DIMA/piva/mpi_freiburg/wt_only_linked_with_mutants_raw.h5ad',
           compression = "none")


assay=combined
DefaultAssay(assay)
# normalization, umap and clustering
assay <- RunTFIDF(assay) # default is method 1: log(TF×IDF), and scale factor is 10000
assay <- FindTopFeatures(assay, min.cutoff = 'q0')
assay <- RunSVD(object = assay)

assay <- RunUMAP(object = assay,reduction = 'lsi',dims = 2:30)
assay <- FindNeighbors(object = assay, reduction = 'lsi', dims = 2:30)
#assay <- FindClusters( object = assay, algorithm = 3, #  SLM algorithm;
#                       resolution = 0.4, verbose = FALSE )

DimPlot(assay, group.by = 'sample')+DimPlot(assay, group.by = 'annotation_leiden12')

sce=as.SingleCellExperiment(assay)
# substitute counts slot to save normalized data
sce@assays@data$counts=assay@assays$peaks@data 

max(sce@assays@data$counts)
library(zellkonverter)
writeH5AD( sce, file='/hpcnfs/scratch/DIMA/piva/mpi_freiburg/wt_only_linked_with_mutants_norm.h5ad',
           compression = "none")










