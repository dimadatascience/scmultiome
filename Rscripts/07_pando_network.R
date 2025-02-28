library(Seurat)
library(Pando)
library(JASPAR2020)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(GenomicRanges)
library(Signac)
library(tidyr)
library(glmnet)
library(doParallel)
library(xgboost)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(grr)
library(motifmatchr)

setwd('/hpcnfs/scratch/DIMA/piva/mpi_freiburg/scmultiome/')
seurat_object_tot=readRDS('./data/scmultiome_230629.rds')

## ---- Pando Processing ----

# get motifs
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 7227,  # taxonomy id drosophila
                                                    all_versions = FALSE))
motif_tfs=as.data.frame(cbind(motifs=ID(pwm_set),tfs=name(pwm_set) )) 
pwm_set=pwm_set[motif_tfs$motifs]


# preparing single cell data for Pando processing
RenameGenesSeurat <- function(obj = seurat_object, newnames = newnames) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  ATAC <- obj@assays$ATAC
  if (nrow(ATAC) == length(newnames)) {
    if (length(ATAC@counts)) ATAC@counts@Dimnames[[1]]            <- newnames
    if (length(ATAC@data)) ATAC@data@Dimnames[[1]]                <- newnames
    if (length(ATAC@scale.data)) ATAC@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(ATAC) != nrow(newnames)"}
  obj@assays$ATAC <- ATAC
  return(obj)
}
# rename peaks "chrUn-xx-start-end", otherwise it will try to use xx as peak start
newnames=gsub('chrUn-', 'chrUn_' , rownames(GetAssay(seurat_object_tot, assay="ATAC")) )
seurat_object_tot=RenameGenesSeurat(obj = seurat_object_tot, newnames = newnames)


## set Params
method="glm"
filename='updown20k'

# Select only WT cells
seurat_object=subset(seurat_object_tot, subset = sample %in% c('wt'))

# Initiate GRN object and select candidate regions
seurat_object <- initiate_grn(seurat_object, 
                              peak_assay = "ATAC", 
                              rna_assay = "RNA", 
                              exclude_exons = TRUE #  Whether to consider exons for binding site inference
)

# Scan candidate regions for TF binding motifs
seurat_object <- find_motifs(
  seurat_object,
  pfm = pwm_set,
  motif_tfs = motif_tfs, #	df matching motifs with TFs. The first column: name of the motif, the second the name of the TF.
  genome = BSgenome.Dmelanogaster.UCSC.dm6
)
# Infer gene regulatory network
# gene_name column must be flybase
seurat_object@assays$ATAC@annotation$gene_symbol=seurat_object@assays$ATAC@annotation$gene_name
seurat_object@assays$ATAC@annotation$gene_name=seurat_object@assays$ATAC@annotation$gene_id


registerDoParallel(4)
seurat_object <- infer_grn(seurat_object,
                           peak_to_gene_method = 'Signac',
                           method =method, 
                           genes= seurat_object@assays$integrated@var.features , # high variable genes
                           parallel = TRUE, 
                           upstream=2e+04,
                           downstream=2e+04, 
                           extend=1e+05,
                           only_tss = TRUE
)

# Find gene and regulatory modules 
seurat_object <- find_modules(seurat_object, 
                              p_thresh = 0.05,  #  threshold on adjusted p-value
                              rsq_thresh = 0.1, #  R2 threshold on adjusted p-value
                              nvar_thresh = 10,  #  minimum number of variables in the model
                              min_genes_per_module = 5,  # minimum number of genes in a module
)

# Print modules
modules <- NetworkModules(seurat_object) 
# Plot goodness-of-fit metrics
plot_gof(seurat_object, point_size=3)
plot_module_metrics(seurat_object)

# Visualizing the GRN
seurat_object <- get_network_graph(seurat_object)
plot_network_graph(seurat_object, label_nodes=F,   node_size = c(2, 8))


save(seurat_object,  file='./results/networks_wt.rds')


# --- Plot graph with nodes colored by germlayer ---
# This function return the germlayer in which the gene is more expressed
from_gene_to_highestgerm = function(sobj, sample='wt'){
  network_node=as.data.frame((sobj@grn@networks$glm_network@graphs$module_graph %>% activate(nodes)))
  gene_conversion_tmp2 = gene_conversion$validated_id 
  names(gene_conversion_tmp2) = gene_conversion$current_symbol
  network_node$flybase= gene_conversion_tmp2[network_node$name]

  df_avg_gex=AverageExpression( subset(seurat_object_tot, sample==sample) ,
                                assays = 'RNA', group.by = 'annotation_leiden12', features =network_node$flybase  )$RNA
  
  max_expression_df <- data.frame(
    gene = rownames(df_avg_gex),
    max_cell_type = apply(df_avg_gex, 1, function(x) names(which.max(x))),
    max_expression = apply(df_avg_gex, 1, max)
  )
  max_expression_df$symbol=gene_conversion_tmp[max_expression_df$gene]
  max_expression=max_expression_df$max_cell_type
  names(max_expression)=max_expression_df$gene
  
  return(max_expression)
}


graph= sobj_wt@grn@networks$glm_network@graphs$module_graph
max_expression=from_gene_to_highestgerm(sobj_wt)
names = graph %>% activate(nodes) %>% pull(name)
graph = graph %>% activate(nodes) %>% mutate(highest_germ = max_expression[names(names)] )
sobj_wt@grn@networks$glm_network@graphs$module_graph=graph

pdf(sprintf('/hpcnfs/scratch/temporary/piva/plots/graphs_wt_only_%s.pdf', file_name), 
    width=8, height=8)
plot_network_graph(sobj_wt, label_nodes=TRUE,color_edges=TRUE) +
  geom_node_point(aes(color = factor(highest_germ)), size = 5) +
  scale_color_manual(values = custom_colors) 
dev.off()


# --- Plot graph with nodes colored by gene ID score ---
# The gene ID score represents the proportion of cells expressing the gene within the appropriate germ layer. 
# Specifically, it is calculated as the ratio of the number of cells expressing the gene in the proper germ layer 
# respect to the total number of cells expressing the gene.

graph=sobj_wt@grn@networks$glm_network@graphs$module_graph
max_expression=from_gene_to_highestgerm(sobj_wt)

# defined as "proper expression" the germ layer in which the gene is more expressed in wt condition
proper_expr=as.data.frame(max_expression)
colnames(proper_expr)='wt'
head(proper_expr)

# ---- EZ ----
# compute the germlayer in which the genes within the wt network are more expressed in ez-kd mutant
max_ez=as.data.frame(from_gene_to_highestgerm(sobj_wt,sample='ez' ))
if(identical(rownames(max_ez),rownames(proper_expr))){
  proper_expr$ez=max_ez$`from_gene_to_highestgerm(sobj_wt, sample = "ez")`
}

# compute the fraction of cell expressing the gene in the proper germlayer
subset_sobj=subset(seurat_object_tot, sample=='ez')
expression_matrix <- GetAssayData(subset_sobj, 
                                  slot = "counts", assay = 'RNA')  # or "data" for normalized
proper_expr$gene=rownames(proper_expr)

# Compute the metric
result <- proper_expr %>%
  rowwise() %>%
  mutate(
    num_expressing_in_cell_type_ez = sum(expression_matrix[gene, ] > 0 & subset_sobj$annotation_leiden12 == wt),
    total_num_expressing_ez = sum(expression_matrix[gene, ] > 0),
    metric_ez = num_expressing_in_cell_type_ez / total_num_expressing_ez
  ) 

ectopic=result$metric_ez
names(ectopic)=result$gene
# ---- 

# ---- WT ----
subset_sobj=subset(seurat_object_tot, sample=='wt')
expression_matrix <- GetAssayData(subset_sobj, 
                                  slot = "counts", assay = 'RNA')  # or "data" for normalized
# Compute the metric
result <- result %>%
  rowwise() %>%
  mutate(
    num_expressing_in_cell_type_wt = sum(expression_matrix[gene, ] > 0 & subset_sobj$annotation_leiden12 == wt),
    total_num_expressing_wt = sum(expression_matrix[gene, ] > 0),
    metric_wt = num_expressing_in_cell_type_wt / total_num_expressing_wt
  ) 

ectopic_wt=result$metric_wt
names(ectopic_wt)=result$gene

# ---- CBP ----
subset_sobj=subset(seurat_object_tot, sample=='nej')
expression_matrix <- GetAssayData(subset_sobj, 
                                  slot = "counts", assay = 'RNA')  # or "data" for normalized
# Compute the metric
result <- result %>%
  rowwise() %>%
  mutate(
    num_expressing_in_cell_type_cbp = sum(expression_matrix[gene, ] > 0 & subset_sobj$annotation_leiden12 == wt),
    total_num_expressing_cbp = sum(expression_matrix[gene, ] > 0),
    metric_cbp = num_expressing_in_cell_type_cbp / total_num_expressing_cbp
  ) 

ectopic_cbp=result$metric_cbp
names(ectopic_cbp)=result$gene

head(result)


# plot the network

graph=sobj_wt@grn@networks$glm_network@graphs$module_graph
names = graph %>% activate(nodes) %>% pull(name)
graph = graph %>% activate(nodes) %>% mutate(metric_ez = ectopic[names(names)] )
graph = graph %>% activate(nodes) %>% mutate(metric_wt = ectopic_wt[names(names)] )
graph = graph %>% activate(nodes) %>% mutate(metric_cbp = ectopic_cbp[names(names)] )


plot_network_graph(sobj_wt, label_nodes=TRUE,color_edges=TRUE) +
  geom_node_point(aes(color = factor(highest_germ)), size = 5) +
  scale_color_manual(values = custom_colors) + ggtitle("WT network") | 
  plot_network_graph(sobj_wt, label_nodes=TRUE, color_edges=FALSE ) +
  geom_node_point(aes(color = metric_wt, size=centrality), size = 5) +
  theme_void() + ggtitle("WT\n N.cell in proper germ / Total n. of cells expressing the gene") |
  plot_network_graph(sobj_wt, label_nodes=TRUE, color_edges=FALSE ) +
  geom_node_point(aes(color = metric_ez, size=centrality), size = 5) +
  theme_void() + ggtitle("Ez\nN.cell in proper germ / Total n. of cells expressing the gene") |
  plot_network_graph(sobj_wt, label_nodes=TRUE, color_edges=FALSE ) +
  geom_node_point(aes(color = metric_cbp, size=centrality), size = 5) +
  theme_void() + ggtitle("CBP\nN.cell in proper germ / Total n. of cells expressing the gene")

# ---



