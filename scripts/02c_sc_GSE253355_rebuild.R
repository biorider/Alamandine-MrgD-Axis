# 02c: GSE253355 (Bandyopadhyay 2024). The authors' processed Seurat object lacks MRGPRD (gene-filtered),
# so per-cell counts are rebuilt from the RAW 10x matrices (33,538 features) and joined to the authors' labels by barcode.
# Run from repo root:  Rscript scripts/02c_sc_GSE253355_rebuild.R   (env var PHASE0_OPEN may point to the download dir)
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
open_dir <- Sys.getenv("PHASE0_OPEN", "data/open"); interim <- "data/interim"
rds <- file.path(open_dir, "GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz")
tmp <- tempfile(fileext = ".rds"); R.utils::gunzip(rds, destname = tmp, remove = FALSE, overwrite = TRUE)  # outer gzip layer
atl <- readRDS(tmp); unlink(tmp)                                                                         # inner layer read natively
cat("MRGPRD in authors' object:", "MRGPRD" %in% rownames(atl), "\n")
md <- atl@meta.data; md$cell <- rownames(md); rm(atl); invisible(gc())
md$compartment <- with(md, ifelse(cluster_anno_coarse %in% c("Mesenchymal","Endothelial"), "stroma",
                     ifelse(cluster_anno_l1 == "HSPC", "HSPC",
                     ifelse(cluster_anno_l2 %in% c("Monocyte","Macrophages"), "mono/macro",
                     ifelse(cluster_anno_l2 %in% c("CD4+ T-Cell","CD8+ T-Cell"), "T/NK",
                     ifelse(cluster_anno_coarse == "Muscle", "VSMC (perivascular)", paste0("other: ", cluster_anno_l1)))))))
rawdir <- file.path(interim, "GSE253355_raw"); dir.create(rawdir, showWarnings = FALSE)
untar(file.path(open_dir, "GSE253355_RAW.tar"), exdir = rawdir)
track <- c("MRGPRD","MAS1","CD34","MPO","CXCL12","ALB")
mats <- list.files(rawdir, pattern = "matrix.mtx.gz$"); cells <- list(); detnum <- NULL; detden <- c()
for (m in mats) {
  pre <- sub("_matrix.mtx.gz$", "", m); samp <- sub("^GSM[0-9]+_", "", pre)
  X <- as(readMM(file.path(rawdir, m)), "CsparseMatrix")
  ft <- read.delim(file.path(rawdir, paste0(pre, "_features.tsv.gz")), header = FALSE)
  bc <- read.delim(file.path(rawdir, paste0(pre, "_barcodes.tsv.gz")), header = FALSE)$V1
  colnames(X) <- paste0(samp, "_", bc); stopifnot(all(c("MRGPRD","MAS1") %in% ft$V2))
  keep <- intersect(colnames(X), md$cell); X <- X[, keep, drop = FALSE]; sm <- md[keep, ]
  df <- data.frame(cell = keep, sample = samp, patient = sub("_MACS$", "", samp), compartment = sm$compartment,
                   label_l2 = sm$cluster_anno_l2, total_umi = colSums(X), n_genes = colSums(X > 0))
  for (g in track) df[[g]] <- colSums(X[ft$V2 == g, , drop = FALSE])
  cells[[samp]] <- df; B <- X > 0
  for (k in unique(sm$compartment)) {
    idx <- which(sm$compartment == k); v <- rowSums(B[, idx, drop = FALSE])
    if (is.null(detnum)) detnum <- matrix(0, nrow(X), 0, dimnames = list(make.unique(ft$V2), NULL))
    if (!k %in% colnames(detnum)) { detnum <- cbind(detnum, 0); colnames(detnum)[ncol(detnum)] <- k; detden[k] <- 0 }
    detnum[, k] <- detnum[, k] + v; detden[k] <- detden[k] + length(idx)
  }
}
ct <- do.call(rbind, cells); write.csv(ct, file.path(interim, "sc_cells_GSE253355.csv"), row.names = FALSE)
write.csv(sweep(detnum, 2, detden[colnames(detnum)], "/"), file.path(interim, "sc_det_GSE253355.csv"))
cat("cells", nrow(ct), "of atlas", nrow(md), "\n"); print(table(ct$compartment))
