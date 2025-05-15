#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)  # à inclure dans l'env Conda
  library(rhdf5)
  library(Matrix)
  library(infercnv)
  library(IRanges)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# ─── Définir les options de ligne de commande ───────────────────────────
option_list <- list(
  make_option("--NormalCellFile", type="character", default = "data/RUN1_S1_hFF_WT.dna.h5"),
  make_option("--TumorCellFile",  type="character", default = "data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5"),
  make_option("--chrom",          type="character", default="10"),
  make_option("--max_cells",      type="integer",   default=5000),
  make_option("--min_reads",      type="integer",   default=100),
  make_option("--max_reads",      type="integer",   default=50000),
  make_option("--target_norm",    type="integer",   default=10000),
  make_option("--target_tum",     type="integer",   default=5000),
  make_option("--seed",           type="integer",   default=42),
  make_option("--threads",        type="integer",   default=4),
  make_option("--cutoff",         type="double",    default=0.1),
  make_option("--workdir",        type="character", default=getwd()),
  make_option("--HMM",            type="character", default="i6"),
  make_option("--out_dir",        type="character", default="./results/infercnv_out")
)

opt <- parse_args(OptionParser(option_list = option_list))
unlink(opt$out_dir, recursive = TRUE)  # nettoyage éventuel
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)


# ─── Initialisation globale ─────────────────────────────────────────────
set.seed(opt$seed)
options(expressions = 5e5)
setwd(opt$workdir)


# ─── Fonction de lecture d’un fichier .h5 spécifique à un chromosome ────
read_counts_h5 <- function(path, chrom, max_cells) {
  all_chr <- h5read(path, "/assays/dna_read_counts/ca/CHROM")
  idx_chr <- which(all_chr == chrom)
  stopifnot(length(idx_chr) > 0)
  
  all_bc <- h5read(path, "/assays/dna_read_counts/ra/barcode")
  sampled_cols <- if (length(all_bc) > max_cells) sample(seq_along(all_bc), max_cells) else seq_along(all_bc)
  bc_keep <- all_bc[sampled_cols]
  
  mat <- h5read(path, "/assays/dna_read_counts/layers/read_counts", index = list(idx_chr, sampled_cols))
  gene_ids <- h5read(path, "/assays/dna_read_counts/ca/id", index = list(idx_chr))
  starts   <- h5read(path, "/assays/dna_read_counts/ca/start_pos", index = list(idx_chr))
  ends     <- h5read(path, "/assays/dna_read_counts/ca/end_pos", index = list(idx_chr))
  
  rownames(mat) <- make.unique(gene_ids)
  colnames(mat) <- bc_keep
  
  gene_df <- data.frame(gene = rownames(mat), chr = chrom, start = as.numeric(starts), end = as.numeric(ends))
  list(mat = mat, gene_df = gene_df)
}

# ─── Lecture des matrices normal & tumor ─────────────────────────────────
cat(">>> Lecture des matrices de comptage H5...\n")
run1 <- read_counts_h5(opt$NormalCellFile, opt$chrom, opt$max_cells)
run2 <- read_counts_h5(opt$TumorCellFile,  opt$chrom, opt$max_cells)

# ─── Préfixes N_ et T_ pour les barcodes ────────────────────────────────
colnames(run1$mat) <- paste0("N_", colnames(run1$mat))
colnames(run2$mat) <- paste0("T_", colnames(run2$mat))

# ─── Fusion des matrices de comptage ─────────────────────────────────────
counts <- cbind(run1$mat, run2$mat)

# ─── Filtrage reads-per-cell ─────────────────────────────────────────────
reads <- Matrix::colSums(counts)
keep_cells <- which(reads >= opt$min_reads & reads <= opt$max_reads)
counts <- counts[, keep_cells]

# ─── Identification origine (normal / tumor) ─────────────────────────────
bc_all    <- colnames(counts)
is_normal <- startsWith(bc_all, "N_")
is_tumor  <- startsWith(bc_all, "T_")

cat(">>> Répartition après filtrage reads-per-cell :\n")
print(table(normal = is_normal, tumor = is_tumor))

# ─── Sous-échantillonnage équilibré ──────────────────────────────────────
n_norm_sel <- min(opt$target_norm, sum(is_normal))
n_tum_sel  <- min(opt$target_tum,  sum(is_tumor))

cat(sprintf(">>> Sélection de %d normales et %d tumorales\n", n_norm_sel, n_tum_sel))

sel_norm <- sample(which(is_normal), n_norm_sel)
sel_tum  <- sample(which(is_tumor),  n_tum_sel)

keep_balanced <- sort(c(sel_norm, sel_tum))
counts <- counts[, keep_balanced]

# ─── Génération des fichiers d’annotation & gene_order ───────────────────
bc_all <- colnames(counts)

annotation_df <- data.frame(
  cell_id   = bc_all,
  cell_type = ifelse(startsWith(bc_all, "T_"), "tumor", "normal")
)

write.table(annotation_df, file.path(opt$out_dir, "annotation.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(run1$gene_df, file.path(opt$out_dir, "gene_order.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ─── Création de l’objet InferCNV ────────────────────────────────────────

# ─── Création objet Seurat avec les mêmes données ────────────────────────
seu <- CreateSeuratObject(counts = counts)

annotations <- read.table(file.path(opt$out_dir, "annotation.txt"), sep = "\t", header = FALSE, col.names = c("cell", "type"))
seu$cell_type <- annotations$type[match(colnames(seu), annotations$cell)]

# (Optionnel pour la suite : UMAP / duplication)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)
seu <- RunUMAP(seu, dims = 1:20)

inf_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = file.path(opt$out_dir, "annotation.txt"),
  gene_order_file   = file.path(opt$out_dir, "gene_order.txt"),
  ref_group_names   = "normal"
)


cat(">>> Lancement d'inferCNV...\n")
inf_obj <- infercnv::run(
  infercnv_obj       = inf_obj,
  cutoff             = opt$cutoff,
  denoise            = TRUE,
  HMM                = TRUE,
  out_dir            = opt$out_dir,
  HMM_type           = opt$HMM,
  num_threads        = opt$threads,
  leiden_resolution  = 0.001
)


# ─── Ajouter résultats inferCNV au Seurat object ─────────────────────────
seu <- infercnv::add_to_seurat(seu, infercnv_output_path = opt$out_dir, top_n = 10)

# ─── Lecture des métadonnées annotées par inferCNV ───────────────────────
meta <- read.table(
  file.path(opt$out_dir, "map_metadata_from_infercnv.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

meta_tumor <- meta[grep("^tumor_", meta$subcluster), ]
gain_chr10_global_pct <- mean(meta_tumor$proportion_dupli_10, na.rm = TRUE) * 100
cat(sprintf("→ Gain global moyen chr10 (proportion_dupli_10) : %.2f %%\n", gain_chr10_global_pct))

hmm_file <- Sys.glob(file.path(opt$out_dir, "HMM_CNV_predictions.*.pred_cnv_genes.dat"))[1]
hmm <- fread(hmm_file)

# Harmoniser le nom de la colonne des cellules
cell_col <- grep("^cell", names(hmm), value = TRUE)[1]
setnames(hmm, cell_col, "cell_id")

# Ajouter type de cellule
hmm[, cell_type := ifelse(grepl("^tumor\\.", cell_id), "tumor", "normal")]

# Ajouter colonne « copies supplémentaires »
hmm[, copies_extra := dplyr::case_when(
  state == 4 ~ 1,
  state == 5 ~ 2,
  state == 6 ~ 3,
  TRUE       ~ 0
)]

hmm_roi <- hmm[chr == opt$chrom & start >= 35e6 & end <= 125e6]

mean_by_cell <- hmm_roi[, .(mean_extra = mean(copies_extra)), by = .(cell_id, cell_type)]
mean_tum <- mean(mean_by_cell[cell_type == "tumor"]$mean_extra)
mean_norm <- mean(mean_by_cell[cell_type == "normal"]$mean_extra)
gain_vs_norm <- mean_tum - mean_norm

cat(sprintf("\n→ Gain moyen dans 35–125 Mb (tumor – normal) : %.2f copies supplémentaires\n", gain_vs_norm))
cat(sprintf("   ≈ %.1f %% d’ADN en plus dans cette région (vs 2 copies attendues)\n", 100 * gain_vs_norm / 2))

résumer_gain_vs_norm <- function(hmm_tumor, hmm_norm, label, 
                                 start_mb = NULL, end_mb = NULL, 
                                 barcodes_cluster, barcodes_tumor) {
  if (!is.null(start_mb)) {
    hmm_tumor <- hmm_tumor[start >= start_mb * 1e6 & end <= end_mb * 1e6]
    hmm_norm  <- hmm_norm[start >= start_mb * 1e6 & end <= end_mb * 1e6]
  }
  
  gain_mean <- mean(hmm_tumor$copies_extra, na.rm = TRUE) - 
    mean(hmm_norm$copies_extra, na.rm = TRUE)
  
  pct <- 100 * length(barcodes_cluster) / length(barcodes_tumor)
  
  cat(sprintf("\n[%s] %d cellules (%.1f%% des tumorales)\n", label, length(barcodes_cluster), pct))
  if (!is.null(start_mb)) {
    cat(sprintf("[%s] Gain moyen (%.0f–%.0f Mb) : %.2f copies (≈ %.1f%% d’ADN en plus)\n", 
                label, start_mb, end_mb, gain_mean, 100 * gain_mean / 2))
  } else {
    cat(sprintf("[%s] Gain moyen chr10 : %.2f copies (≈ %.1f%% d’ADN en plus)\n", 
                label, gain_mean, 100 * gain_mean / 2))
  }
}
# ─── Trouver le sous-cluster avec le plus de gains 5 ou 6 ────────────────

# Extraire uniquement les lignes tumorales
hmm_tumor <- hmm[grepl("^tumor\\.", cell_id)]
hmm_tumor <- hmm_tumor[order(hmm_tumor$start,)]


# Extraire le nom du sous-cluster (ex: tumor.s9 → s9)
hmm_tumor[, subcluster := sub("^tumor\\.", "", cell_id)]

# Filtrer pour les états 5 et 6 (forts gains)
hmm_gain56 <- hmm_tumor[state %in% c(5, 6)]

# Calcul : pourcentage de bins "5 ou 6" par sous-cluster
gain_stats <- hmm_gain56[, .N, by = subcluster]
gain_stats <- gain_stats[order(-N)]

# Récupérer le top 1
best_subcluster <- gain_stats$subcluster[1]
cat(sprintf(">>> Sous-cluster avec le plus de gains 5/6 : %s\n", best_subcluster))


subclust <- paste0(best_subcluster)
start_int <- 35
end_int <- 127

cells_sclust <- rownames(meta)[meta$subcluster == subclust]
cells_tumor  <- rownames(meta)[startsWith(rownames(meta), "T_")]

hmm_sclust <- hmm[cell_id == paste0("tumor.", subclust)]
hmm_sclust <- hmm_sclust[order(hmm_sclust$start), ]

hmm_norm <- hmm[startsWith(cell_id, "normal.")]

résumer_gain_vs_norm(hmm_sclust, hmm_norm, subclust,
                     start_mb = start_int, end_mb = end_int,
                     barcodes_cluster = cells_sclust,
                     barcodes_tumor = cells_tumor)

résumer_gain_vs_norm(hmm_sclust, hmm_norm, subclust,
                     barcodes_cluster = cells_sclust,
                     barcodes_tumor = cells_tumor)

compter_cellules_avec_gain_region <- function(hmm, meta, start_mb, end_mb, chrom = 10) {
  start_bp <- start_mb * 1e6
  end_bp   <- end_mb * 1e6
  
  gain_bins <- hmm[
    chr == chrom & start >= start_bp & end <= end_bp &
      grepl("^tumor\\.", cell_id) & copies_extra > 0
  ]
  
  gain_subclusters <- unique(sub("^tumor\\.", "", gain_bins$cell_id))
  if (length(gain_subclusters) == 0) {
    cat(sprintf("\n⚠️ Aucun sous-cluster tumoral n’a de gain entre %d–%d Mb.\n", start_mb, end_mb))
    return(NULL)
  }
  
  tumor_barcodes <- rownames(meta)[startsWith(rownames(meta), "T_")]
  subclust_vec <- meta[tumor_barcodes, "subcluster"]
  cells_with_gain <- tumor_barcodes[subclust_vec %in% gain_subclusters]
  
  n_gain <- length(cells_with_gain)
  n_total <- length(tumor_barcodes)
  pct <- 100 * n_gain / n_total
  
  cat(sprintf("\n→ %d cellules tumorales (%.1f%% sur %d) présentent un gain entre %d–%d Mb\n",
              n_gain, pct, n_total, start_mb, end_mb))
  cat("   Sous-clusters impliqués :", paste(gain_subclusters, collapse = ", "), "\n")
  
  return(invisible(list(pct_gain = pct)))
}

compter_cellules_avec_gain_region(hmm, meta, start_int, end_int)
# Seuil sur proportion_dupli_10 (modifiable selon contexte)
dup_treshold <- 0.2

dupTumor <- WhichCells(
  seu,
  expression = cell_type == "tumor" & proportion_dupli_10 >= dup_treshold
)

cat(sprintf("\n>>> %d cellules tumorales ont une duplication ≥ %.1f\n", length(dupTumor), dup_treshold))

# UMAP avec surlignage des cellules dupliquées
p <- DimPlot(
  seu,
  reduction = "umap",
  cells.highlight = dupTumor,
  cols.highlight = "red",
  cols = "lightgrey"
) + ggtitle("Cellules tumorales avec un gain ≥ 0.2")

ggsave(
  filename = file.path(opt$out_dir, "duplication_chr10_umap.pdf"),
  plot = p, width = 8, height = 6
)

# Fonctions utilitaires
percent_extra <- function(x) round(x / 2 * 100, 2)  # copies supplémentaires → % ADN

# Compilation des valeurs
summary_df <- data.frame(
  métrique = c(
    "Gain global chr10 (tumor) [% ADN +]",
    "Gain moyen 35–125 Mb (tumor – normal) [% ADN +]",
    sprintf("Gain moyen chr10 – %s", subclust),
    sprintf("Gain moyen 35–127 Mb – %s", subclust),
    "Cellules tumorales avec gain 35–127 Mb [% cellules]"
  ),
  valeur = c(
    round(gain_chr10_global_pct, 2),
    percent_extra(gain_vs_norm),
    percent_extra(mean(hmm_sclust$copies_extra, na.rm = TRUE)),
    percent_extra(mean(hmm_sclust[start >= start_int*1e6 & end <= end_int*1e6]$copies_extra, na.rm = TRUE)),
    round(compter_cellules_avec_gain_region(hmm, meta, start_int, end_int)$pct_gain, 1)
  )
)

# Sauvegardes
write.csv(summary_df, file = file.path(opt$out_dir, "final_compare.csv"), row.names = FALSE)
write.csv(hmm_sclust, file = file.path(opt$out_dir, "résumé_CNV_tumor_subscluster_max_dupli.csv"), row.names = FALSE)
write.csv(hmm_tumor, file = file.path(opt$out_dir, "résumé_CNV_all_subscluster_tumor.csv"), row.names = FALSE)
file.create(file.path(opt$out_dir, ".done"))



