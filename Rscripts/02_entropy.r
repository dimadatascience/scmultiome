library(Seurat)
library(destiny)  
library(SingleCellExperiment)
library(ggplot2)
library(Signac)

setwd('/hpcnfs/scratch/DIMA/piva/mpi_freiburg/scmultiome/') # path to repository
data = readRDS("./data/scmultiome_230629.rds")

## ==== SCENT: entropy ====
library(AnnotationDbi)
library(org.Hs.eg.db)
library(SCENT)

m=read.csv('./data/hd_fly_ortholog.csv')
gs=m$Human_gene
gs2entrez=AnnotationDbi::select(org.Hs.eg.db, keys=gs, columns='ENTREZID', keytype='SYMBOL')
head(gs2entrez)

joined_df <- merge(gs2entrez, m, by.x = "SYMBOL", how='inner',
                   by.y = "Human_gene", all.x = FALSE, all.y = FALSE)
joined_df=joined_df[!duplicated(joined_df$FBgn),]
joined_df=joined_df[!duplicated(joined_df$ENTREZID),]
data(net13Jun12)
net13Jun12.m[1:5,1:5]
joined_df$ENTREZID

gen=intersect(joined_df$ENTREZID, rownames(net13Jun12.m))

joined_df=joined_df[!is.na(joined_df$ENTREZID),]
rownames(joined_df)=joined_df$ENTREZID
joined_df=joined_df[gen,]
net13Jun12.m=net13Jun12.m[gen,gen]
dim(net13Jun12.m)
colnames(net13Jun12.m)=joined_df$FBgn
rownames(net13Jun12.m)=joined_df$FBgn

length(unique(joined_df$ENTREZID))

net13Jun12.m[1:8,1:8]

library(qlcMatrix)  
expr=data@assays$RNA@data
dim(expr)
expr=as.matrix(expr)
expr[1:8,1:8]

#expr=expr[colnames(net13Jun12.m),]
commonEID.v <- intersect(rownames(net13Jun12.m),rownames(expr))
k.v <- rowSums(net13Jun12.m[match(commonEID.v, rownames(net13Jun12.m)),])
ccat.v <- as.vector(corSparse(expr[match(commonEID.v,rownames(expr)),],Matrix(matrix(k.v,ncol=1))))

sum(rownames(data@meta.data)==colnames(expr))==length(colnames(expr))
data@meta.data['entropy_ccat']=ccat.v

pot.o <- InferPotencyStates(potest.v=data@meta.data['entropy_ccat'], 
                            pheno.v=data@meta.data['sample'],
                            type='CCAT', 
)
pot.o$het

FeaturePlot(data, features = 'entropy_ccat')
FeaturePlot(data, features = 'entropy_ccat', split.by='sample')
VlnPlot(data, features = 'entropy_ccat',sort='increasing',  pt.size =0.1)
VlnPlot(data, features = 'entropy_ccat', group.by='sample',sort='increasing',  pt.size =0.1)
VlnPlot(data, features = 'entropy_ccat', split.by='sample', sort='increasing',  pt.size =0.1)

write.csv(data@meta.data, './results/obs_entropy_ccat1.csv')
