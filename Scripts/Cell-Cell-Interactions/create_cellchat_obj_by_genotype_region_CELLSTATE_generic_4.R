#!/usr/bin/env Rscript
#$ -cwd

library(CellChat)
library(patchwork)
library(Seurat)
library(SingleCellExperiment)
library(ggalluvial)
library(repr)
options(stringsAsFactors = FALSE)

options(future.globals.maxSize = 20000 * 1024 ^ 2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Require Genotype and region to be specified.", call.=FALSE)
} else if (length(args) > 2) {
  stop("Too many arguments provided. Require Genotypte and region to be specified.", call.=FALSE)
} else {
  GENE <- args[1]
  REGION_GROUP <- args[2]
  if (REGION_GROUP == "RV") {
    REGION <- "RV"
  } else if (REGION_GROUP == "LV") {
    REGION <- c("FW","AP","S")
  } else {
    stop("Invalid region specified, must be LV or RV", call.=FALSE)
  }
}

INDIR <- "/fast/AG_Huebner/huebner3/ANALYSES/20190926_gp_BadOyenhausen/scanpy/20210302/Cellchat/"
PREFIX <- "global_all_ANNOTATED_CELLSTATES_RAW_V6_"
SUFFIX <- ".rds"

OUTDIR <- "/fast/AG_Huebner/huebner3/ANALYSES/20190926_gp_BadOyenhausen/scanpy/20210302/Cellchat/"


INFILE <- paste(INDIR, PREFIX, GENE, SUFFIX, sep='')
OUTFILE <- paste(OUTDIR, "cellchat_global_V6_", GENE, "_", REGION_GROUP, "_CELLSTATES.rds", sep='')

if (!file.exists(INFILE)){
  ERROR.MSG <- paste(".h5seurat file does not exist: ", INFILE, sep='')
  stop(ERROR.MSG, call.=FALSE)
}

seurat_obj_orig <-  readRDS(INFILE)
seurat_obj_orig

seurat_obj <- subset(x = seurat_obj_orig, subset = validated == TRUE)
seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(object = seurat_obj) <- "Sample"

##########
### RV ###
##########

seurat_obj_RV <- subset(x = seurat_obj, subset = Region %in% REGION)
seurat_obj_RV
Idents(object = seurat_obj_RV) <- "cell_states"
cellchat <- createCellChat(object = seurat_obj_RV, group.by = "ident")
cellchat
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 16) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, file = OUTFILE)
