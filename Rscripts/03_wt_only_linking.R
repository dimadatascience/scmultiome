setwd("/hpcscratch/ieo/DIMA/scmultiome_mpi/08_other/wt_complete/")

library(Seurat,  lib.loc = "/usr/local/lib/R/site-library")
library(Signac)
require(rtracklayer)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(dplyr)

# Peak linking
# Import called peaks

data = readRDS("/hpcnfs/scratch/temporary/DIMA/scmultiome_mpi/scmultiome_230629.rds")
oldpath='/tungstenfs/scratch/shared/gchao_ggiorget/data_freiburg/multiome/'
newpath='/hpcnfs/scratch/temporary/DIMA/atac_fragments_freiburg/'
for (e in 1:12) {
  data@assays$ATAC@fragments[[e]]@path=gsub(oldpath,newpath, data@assays$ATAC@fragments[[e]]@path )
}

sample_name='wt'
  
basename=paste0(sample_name,'_only')
print(sprintf("Running %s", basename))

srat.sub <-
  data %>%
  subset(stage == 5) %>%
  subset(sample == sample_name)

peaks <- import(sprintf("%s_peaks.narrowPeak", basename))

select = c("chr2L", "chr2R", "chr3L", "chr3R", "chrX", "chrY")
peaks = peaks[as.vector(seqnames(peaks)) %in% select]

DefaultAssay(srat.sub) <- "ATAC"
macs2_counts <- FeatureMatrix(
  fragments = Fragments(srat.sub),
  features = peaks,
  cells = Cells(srat.sub)
)
annotations= srat.sub@assays$ATAC@annotation
# Create a new assay using the MACS2 peak set
srat.sub[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
                                            annotation = annotations)


DefaultAssay(srat.sub) <- "peaks"
# srat.sub[['peaks']]=RegionStats(srat.sub[['peaks']], genome = BSgenome.Dmelanogaster.UCSC.dm6)

srat.sub %>%
  RegionStats(genome = BSgenome.Dmelanogaster.UCSC.dm6,  sep = c("-", "-"),assay='peaks') %>%
  LinkPeaks(peak.assay = "peaks",
            expression.assay = "RNA",
            gene.id = TRUE) %>%
  saveRDS(sprintf("/hpcscratch/ieo/DIMA/scmultiome_mpi/08_other/wt_complete/linked_%s.rds", basename))

## /hpcscratch/ieo/DIMA/scmultiome_mpi/08_other/wt_complete/linked_wt_only.rds 

srat.sub



