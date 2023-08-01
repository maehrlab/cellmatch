## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("cellmatch")

knitr::opts_chunk$set(
  results = "hide", 
  fig.show='hide',
  warning = FALSE ,
  message = FALSE
)

## -----------------------------------------------------------------------------
# These packages help find and manipulate publicly available data for the demo.
library("magrittr")
library("Seurat")
devtools::install_github('satijalab/seurat-data')

# Obtain pancreas data for demo.
SeuratData::AvailableData()
SeuratData::InstallData("panc8")
library("SeuratData")
data("panc8")
panc8@meta.data$orig.ident %>% table

panc8 = Seurat::NormalizeData(panc8)
panc1 = subset(panc8, orig.ident %in% "human1")
panc2 = subset(panc8, orig.ident %in% "human2")
rm(panc8); gc()

## -----------------------------------------------------------------------------
# Aggregate by cluster
panc1$celltype %>% table
panc2$celltype %>% table
panc1_bulk = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(panc1, "data")), panc1$celltype, method = "average")
panc2_bulk = cellmatch::AggregateByCluster(expm1(Seurat::GetAssayData(panc2, "data")), panc2$celltype, method = "average")
panc1_bulk = as.matrix(panc1_bulk)
panc2_bulk = as.matrix(panc2_bulk)

# Convert counts per 10k to cpm
colSums(panc1_bulk)
colSums(panc2_bulk)
panc1_bulk = panc1_bulk*100
panc2_bulk = panc2_bulk*100

# Where do you want your output?
temp_dir = getwd() %>% file.path("demo")
dir.create(temp_dir)

## -----------------------------------------------------------------------------
matching_output = RunCellMatch(
  query = panc1_bulk,
  reference = panc2_bulk,
  results_path = temp_dir,
  K = 2000,
  num_init = 2
)


## -----------------------------------------------------------------------------

# gene selection
variable_genes = cellmatch::SelectInformativeGenes(panc1_bulk, K = 2000)

# Model selection
matchmodel = cellmatch::SelectModelThoroughly(panc1_bulk[variable_genes,],
                                              panc2_bulk[variable_genes,],
                                              verbose = T,
                                              num_init = 2,
                                              compute_penalty = "correlation_distance")
matchmodel$x

# Show plots justifying the model.
eval_results = cellmatch::EvaluateByGene(panc1_bulk[variable_genes,],
                                         panc2_bulk[variable_genes,],
                                         equivalents = matchmodel$x,
                                         do_heatmaps = T,
                                         results_path = temp_dir,
                                         compute_penalty = "correlation_distance" )

# Check neighboring models for goodness of fit.
neighbors = cellmatch::MutateModel(panc1_bulk[variable_genes,],
                                   panc2_bulk[variable_genes,],
                                   init = matchmodel$x)
cellmatch::PlotNeighboringModels(neighbors, results_path = temp_dir)


# Plot whatever genes you like
p = PairedHeatmap(
  query = panc1_bulk ,
  reference = panc2_bulk[,matchmodel$x],
  genes = c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10", "GCG", "PTPRC" , "INS")
)
p

# If you're working with human-mouse comparisons, check out the utilities for gene name conversion.
cellmatch::get_ortholog( head( rownames( panc2 ) ), from = "human", to = "mouse")
cellmatch::has_ortholog( head( rownames( panc2 ) ), from = "human", to = "mouse")
panc2_mousified = cellmatch::convert_species_rownames(head(panc2_bulk, 100), from = "human", to = "mouse")
dim(panc2_mousified)
rownames(panc2_mousified)

## -----------------------------------------------------------------------------
unlink(temp_dir, recursive = T)

