suppressPackageStartupMessages({
    library(splatter)
    library(Seurat)
})

batch.facLoc <- 0.4
de.facLoc <- 0.1
nBatches <- 4
nGroups <- 3
nGenes <- 5000
groupCells <- 300

sim <- splatSimulate(
    batchCells = rep(300 * nGroups, nBatches),
    batch.facLoc = batch.facLoc,
    group.prob = rep(1/nGroups, nGroups),
    nGenes = nGenes,
    batch.facScale = 0.1,
    method = "groups",
    verbose = FALSE,
    out.prob    = 0.001,
    de.prob     = 0.1,
    de.facLoc   = de.facLoc,
    de.facScale = 0.1,
    bcv.common  = 0.1
)
sim

se <- CreateSeuratObject(counts = counts(sim), meta.data = as.data.frame(colData(sim)))

save(se, file = "inst/extdata/seurat_splatter_sim.RData", version = 2)
