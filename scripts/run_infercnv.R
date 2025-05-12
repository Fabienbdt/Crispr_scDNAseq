# === LIBRAIRIES ========================================================


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("S4Arrays")

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("infercnv")

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("rhdf5")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)      # <-- ajoute au yml Conda
  library(rhdf5); library(Matrix); library(infercnv)
  library(GenomicRanges); library(ggplot2); library(dplyr)
  library(IRanges); library(Seurat)
})

## ─── Définition des options ────────────────────────────────────────────
option_list <- list(
  make_option("--NormalCellFile",      type="character", default = "RUN1_S1_hFF_WT.dna.h5"),
  make_option("--TumorCellFill",      type="character", default = "RUN2_S8_hFF_clone_6_KOfluo.dna.h5" ),
  make_option("--chrom",          type="character", default="10"),
  make_option("--max_cells",      type="integer",   default=5000),
  make_option("--min_reads",      type="integer",   default=100),
  make_option("--max_reads",      type="integer",   default=50000),
  make_option("--target_norm",    type="integer",   default=10000),
  make_option("--target_tum",     type="integer",   default=5000),
  make_option("--seed",           type="integer",   default=42),
  make_option("--threads",        type="integer",   default=4),
  make_option("--cutoff",        type="integer",   default=0.1),
  make_option("--workdir",        type="character", default=getwd()),
  make_option("--HMM",        type="character", default="i6"),
  make_option("--out_dir",        type="character", default="infercnv_out")
)
opt <- parse_args(OptionParser(option_list = option_list))

## ─── Appliquer les paramètres globaux ──────────────────────────────────
set.seed(opt$seed)
options(expressions = 5e5)
setwd(opt$workdir)


# === PARAMETRES OPTIONNELS ==============================================
set.seed(opt$seed)
options(expressions = 5e5)
setwd("/Volumes/LaCie")



######################### PARTIE 1 #######################################

# === FONCTION LECTURE .H5 CHR10 =========================================
read_counts_h5 <- function(path, chrom = opt$chrom, max_cells = opt$max_cells) {
  all_chr <- h5read(path, "/assays/dna_read_counts/ca/CHROM")
  idx_chr <- which(all_chr == chrom)
  stopifnot(length(idx_chr) > 0)
  
  all_bc <- h5read(path, "/assays/dna_read_counts/ra/barcode")
  sampled_cols <- if (length(all_bc) > max_cells) sample(seq_along(all_bc), max_cells) else seq_along(all_bc)
  bc_keep <- all_bc[sampled_cols]
  
  mat <- h5read(path, "/assays/dna_read_counts/layers/read_counts", index = list(idx_chr, sampled_cols))
  
  
  gene_ids <- h5read(path, "/assays/dna_read_counts/ca/id", index = list(idx_chr))
  starts   <- as.numeric(h5read(path, "/assays/dna_read_counts/ca/start_pos", index = list(idx_chr)))
  ends     <- as.numeric(h5read(path, "/assays/dna_read_counts/ca/end_pos", index = list(idx_chr)))
  
  rownames(mat) <- make.unique(gene_ids)
  colnames(mat) <- bc_keep
  
  gene_df <- data.frame(gene = rownames(mat), chr = chrom, start = starts, end = ends)
  list(mat = mat, gene_df = gene_df)
}



# === LECTURE & FILTRAGE ===================================================
run1_file <- opt$NormalCellFile   # normal
run2_file <- opt$TumorCellFill    # tumor
out_dir   <- "infercnv_chr10_out_filtered"

run1 <- read_counts_h5(run1_file, chrom = opt$chrom)
run2 <- read_counts_h5(run2_file, chrom = opt$chrom)

# Ajouter les préfixes avant de stocker les barcodes
colnames(run1$mat) <- paste0("N_", colnames(run1$mat))
colnames(run2$mat) <- paste0("T_", colnames(run2$mat))

# ➕ Définition des barcodes d’origine (avec préfixes)
all_bc_run1 <- colnames(run1$mat)
all_bc_run2 <- colnames(run2$mat)

# Trouver les barcodes communs et uniques (avec préfixes)
common_barcodes <- intersect(all_bc_run1, all_bc_run2)
unique_barcodes_run1 <- setdiff(all_bc_run1, all_bc_run2)
unique_barcodes_run2 <- setdiff(all_bc_run2, all_bc_run1)

# Concaténation des deux matrices
counts <- cbind(run1$mat, run2$mat)

# Filtrage reads-per-cell
reads <- Matrix::colSums(counts)
min_reads <- opt$min_reads
max_reads <- opt$max_reads
keep_cells <- which(reads >= min_reads & reads <= max_reads)
counts <- counts[, keep_cells]

# Vérifier provenance des barcodes filtrés
bc_all <- colnames(counts)
is_normal <- bc_all %in% all_bc_run1
is_tumor  <- bc_all %in% all_bc_run2

cat(">>> Répartition après filtrage reads-per-cell :\n")
print(table(normal = is_normal, tumor = is_tumor))



# === SOUS-ÉCHANTILLONNAGE DYNAMIQUE ======================
# Objectif : 5000 normales, 3000 tumorales (ou moins si pas assez disponibles)
target_norm <- opt$target_norm
target_tum  <- opt$target_tum

available_norm <- sum(is_normal)
available_tum  <- sum(is_tumor)

# Ajustement dynamique
n_norm_sel <- min(target_norm, available_norm)
n_tum_sel  <- min(target_tum,  available_tum)

cat(sprintf(">>> Sélection de %d normales (sur %d disponibles) et %d tumorales (sur %d disponibles)\n",
            n_norm_sel, available_norm, n_tum_sel, available_tum))

sel_norm <- sample(which(is_normal), n_norm_sel)
sel_tum  <- sample(which(is_tumor),  n_tum_sel)

keep_balanced <- sort(c(sel_norm, sel_tum))
counts <- counts[, keep_balanced]

bc_all <- colnames(counts)
is_normal <- bc_all %in% all_bc_run1



# === ANNOTATION & GENE ORDER =============================================
annotation_df <- data.frame(
  cell_id   = bc_all,
  cell_type = ifelse(startsWith(bc_all, "T_"), "tumor", "normal")
)

write.table(annotation_df, "annotation.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)


write.table(run1$gene_df, "gene_order.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)



# ======== CREATION OBJET SEURAT =========================================
# 1. Créer l'objet Seurat avec la matrice de comptage
seu <- CreateSeuratObject(counts = counts)

# 2. Ajouter l’annotation de type de cellule
annotations <- read.table("annotation.txt", sep = "\t", header = FALSE, col.names = c("cell", "type"))
seu$cell_type <- annotations$type[match(colnames(seu), annotations$cell)]



# === CREATION OBJET INFERCNV =============================================
unlink(opt$out_dir, recursive = TRUE)

inf_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = "annotation.txt",
  gene_order_file   = "gene_order.txt",
  ref_group_names   = "normal"
)

inf_obj <- infercnv::run(
  infercnv_obj = inf_obj,
  cutoff = opt$cutoff,
  denoise = TRUE,
  HMM = TRUE,
  out_dir = opt$out_dir,
  HMM_type = opt$HMM,
  num_threads = opt$threads,
  leiden_resolution = 0.001
)

# Ajout des résultats InferCNV à un objet Seurat
seu <- infercnv::add_to_seurat(seurat_obj = seu, infercnv_output_path = opt$out_dir, top_n = 10)




######################### PARTIE 2 ########################################

# ========== DETECTION DU GAIN ============================================
# 1) Lire la table d’annotations map_metadata_from_infercnv.txt
# Charger les données d'annotations
meta <- read.table(
  file = file.path(opt$out_dir, "map_metadata_from_infercnv.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extraire uniquement les cellules tumorales si besoin (ex: si elles commencent par "T_")
meta_tumor <- meta[grep("^tumor_", meta$subcluster), ]

# Calcul du total de gain observé (en %)
total_gain_percent <- mean(meta_tumor$proportion_dupli_10, na.rm = TRUE) * 100

# Affichage
print(total_gain_percent)


library(data.table)
library(GenomicRanges)
library(dplyr)

### A. Charger les états HMM par cellule×bin (fichier *genes.dat*)
genes_file <- Sys.glob(file.path(
  opt$out_dir, "HMM_CNV_predictions.*pred_cnv_genes.dat"))[1]
hmm <- fread(genes_file)

# ── Harmoniser les labels
cell_col  <- grep("^cell", names(hmm), value = TRUE)[1]
setnames(hmm, cell_col, "cell_id")

# ── Ajouter type de cellule (normal / tumor)
hmm_roi <- hmm[chr == 10 & start >= 35e6 & end <= 125e6]
hmm_roi[, cell_type := ifelse(grepl("^tumor\\.", cell_id), "tumor", "normal")]



### B.  Convertir l’état en « copies au-dessus du diploïde »
#  (i3: état 1=1 copie, 2=2 copies, 3=≥3 copies => gain = state-2)
hmm_roi[, copies_extra := case_when(
  state == 4 ~ 1,
  state == 5 ~ 2,
  state == 6 ~ 3,
  TRUE       ~ 0
)]
### C. Moyenne par cellule puis par groupe
mean_by_cell <- hmm_roi[, .(mean_extra = mean(copies_extra)),
                        by = .(cell_id, cell_type)]

mean_tum  <- mean(mean_by_cell[cell_type == "tumor"]$mean_extra)
mean_norm <- mean(mean_by_cell[cell_type == "normal"]$mean_extra)

gain_vs_norm <- mean_tum - mean_norm    # normalement ≈ mean_tum (norm ~0)

cat(sprintf(
  "\nGain moyen sur 35–125 Mb : %.2f copies supplémentaires\n", gain_vs_norm))
cat(sprintf(
  "→ Cela signifie qu'en moyenne, les cellules tumorales présentent %.1f %% d'ADN en plus dans cette région par rapport aux cellules normales (2 copies attendues).\n",
  100 * gain_vs_norm / 2))
###############################################################################
##  PARAMÈTRES
###############################################################################
chrom      <- 10               # chromosome d’intérêt (hg19)
start_int  <- 35               # borne Mb inférieure de la zone d’intérêt
end_int    <- 127              # borne Mb supérieure de la zone d’intérêt
subclust   <- "tumor_s7"       # sous-cluster à analyser

###############################################################################
## 1 • AJOUT DE LA COLONNE copies_extra DANS hmm ------------------------------
###############################################################################
hmm[, copies_extra := dplyr::case_when(
  state == 4 ~ 1,     # “duplication” (state-4 dans i6)
  state == 5 ~ 2,
  state == 6 ~ 3,
  TRUE       ~ 0
)]

###############################################################################
## 2 • CELLULES DU SOUS-CLUSTER tumor_s7 --------------------------------------
###############################################################################
meta <- read.table(
  file.path(opt$out_dir, "map_metadata_from_infercnv.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# barcodes de tumor_s7 (préfixe T_)
cells_s7 <- rownames(meta)[meta$subcluster == subclust]

# total de cellules tumorales (préfixe T_)
cells_tumor <- rownames(meta)[startsWith(rownames(meta), "T_")]

###############################################################################
## 3 • FILTRAGE DE hmm POUR LES LIGNES CONCERNANT CE SOUS-CLUSTER -------------
##     (Les prédictions HMM sont déjà agrégées par sous-cluster : une fois par
##      bin génomique – on utilise donc l’étiquette "tumor.<subcluster>").
###############################################################################
# Exemple si tu connais le cluster d’intérêt
subclust <- "tumor_s7"
hmm_s7 <- hmm[cell_id == paste0("tumor.", subclust)]
# lignes pour tumor_s7
hmm_norm <- hmm_roi[startsWith(cell_id, "normal.")]

###############################################################################
## 4 • FONCTION DE RÉSUMÉ : gain moyen (vs normales) + effectifs --------------
###############################################################################
résumer_gain_vs_norm <- function(dt_tumor, dt_norm, label,
                                 start_mb = NULL, end_mb = NULL,
                                 barcodes_cluster, barcodes_tumor) {
  if (!is.null(start_mb) && !is.null(end_mb)) {
    dt_tumor <- dt_tumor[start >= start_mb*1e6 & end <= end_mb*1e6]
    dt_norm  <- dt_norm [start >= start_mb*1e6 & end <= end_mb*1e6]
  }
  
  # moyenne par bin (déjà agrégée par sous-cluster) puis moyenne globale
  mean_tumor <- dt_tumor[, mean(copies_extra, na.rm = TRUE)]
  mean_norm  <- dt_norm [, mean(copies_extra, na.rm = TRUE)]
  gain_vs_norm <- mean_tumor - mean_norm
  
  # taille et pourcentage du sous-cluster (en nombre de cellules, issu de meta)
  n_cluster <- length(barcodes_cluster)
  pct       <- 100 * n_cluster / length(barcodes_tumor)
  
  cat(sprintf("\n[%s] %d cellules tumorales (%.1f %% du total tumor)\n",
              label, n_cluster, pct))
  cat(sprintf("[%s] Gain moyen%s : %.2f copies supplémentaires (vs normales)\n",
              label,
              ifelse(!is.null(start_mb),
                     sprintf(" sur %d–%d Mb", start_mb, end_mb),
                     " sur l'ensemble du chr10"),
              gain_vs_norm))
  cat(sprintf("→ Cela correspond à %.1f %% d'ADN en plus (vs 2 copies attendues)\n",
              100 * gain_vs_norm / 2))
}

###############################################################################
## 5 • RÉSULTATS
###############################################################################
cat("\n==== SOUS-CLUSTER", subclust, "====\n")
# Zone d’intérêt 35–127 Mb
résumer_gain_vs_norm(hmm_s7, hmm_norm, subclust,
                     start_mb = start_int, end_mb = end_int,
                     barcodes_cluster = cells_s7,
                     barcodes_tumor   = cells_tumor)

# Chromosome entier
résumer_gain_vs_norm(hmm_s7, hmm_norm, subclust,
                     barcodes_cluster = cells_s7,
                     barcodes_tumor   = cells_tumor)
compter_cellules_avec_gain_region <- function(hmm, meta, start_mb, end_mb) {
  # Conversion borne en bases
  start_bp <- start_mb * 1e6
  end_bp   <- end_mb * 1e6

  # Sous-ensemble des bins dans la région
  hmm_region <- hmm[start >= start_bp & end <= end_bp]

  # Binariser la présence d’un gain dans chaque bin
  hmm_region[, gain_bin := copies_extra > 0]

  # Regrouper par cellule et vérifier si AU MOINS UN bin a un gain
  cells_with_gain <- hmm_region[, .(gain_present = any(gain_bin, na.rm = TRUE)), by = cell_id]
  cells_with_gain <- cells_with_gain[gain_present == TRUE]

  # Récupérer les barcodes tumoraux uniquement
  tumor_barcodes <- rownames(meta)[startsWith(rownames(meta), "T_")]
  cells_with_gain_tumor <- cells_with_gain[cell_id %in% paste0("tumor.", meta[tumor_barcodes, "subcluster"])]

  # Résumé
  n_gain  <- nrow(cells_with_gain_tumor)
  n_total <- length(tumor_barcodes)
  pct     <- 100 * n_gain / n_total

  cat(sprintf("\n→ %d cellules tumorales (%0.1f %% du total) présentent un gain dans la région %d–%d Mb\n",
              n_gain, pct, start_mb, end_mb))
}
#############
compter_cellules_avec_gain_region <- function(hmm, meta,
                                              start_mb, end_mb,
                                              chrom = 10) {
  # 1 • Région en pb ----------------------------------------------------------
  start_bp <- start_mb * 1e6
  end_bp   <- end_mb   * 1e6
  
  # 2 • Sous-ensemble des bins tumoraux avec gain dans la fenêtre -------------
  gain_bins <- hmm[
    chr == chrom &
      start >= start_bp & end <= end_bp &                # fenêtre
      grepl("^tumor\\.", cell_id) &                      # clusters tumoraux
      copies_extra > 0                                   # gain
  ]
  
  # 3 • Liste des sous-clusters tumoraux impliqués ----------------------------
  #    (on enlève le préfixe "tumor." pour récupérer le nom du subcluster seul)
  gain_subclusters <- unique(sub("^tumor\\.", "", gain_bins$cell_id))
  if (length(gain_subclusters) == 0) {
    cat(sprintf(
      "\n⚠️  Aucun sous-cluster tumoral ne montre de gain entre %d et %d Mb.\n",
      start_mb, end_mb))
    return(invisible(NULL))
  }
  
  # 4 • Barcodes tumoraux et comptage -----------------------------------------
  tumor_barcodes <- rownames(meta)[startsWith(rownames(meta), "T_")]
  subclust_vec   <- meta[tumor_barcodes, "subcluster"]
  
  # cellules dont le subcluster fait partie de la liste avec gain
  cells_with_gain <- tumor_barcodes[subclust_vec %in% gain_subclusters]
  
  n_gain  <- length(cells_with_gain)
  n_total <- length(tumor_barcodes)
  pct     <- 100 * n_gain / n_total
  
  # 5 • Affichage -------------------------------------------------------------
  cat(sprintf(
    "\n→ %d cellules tumorales (%.1f %% sur %d) présentent un gain dans %d–%d Mb\n",
    n_gain, pct, n_total, start_mb, end_mb))
  cat("   Sous-clusters impliqués :", paste(gain_subclusters, collapse = ", "), "\n")
  
  invisible(list(
    n_gain          = n_gain,
    n_total_tumor   = n_total,
    pct_gain        = pct,
    subclusters_hit = gain_subclusters,
    cells_with_gain = cells_with_gain      # vecteur de barcodes utiles si besoin
  ))
}



# Fenêtre 35–127 Mb sur le chromosome 10
compter_cellules_avec_gain_region(hmm, meta, start_int, end_int)


# === DETECTION DES POSITIONS ============================================
# Normalisation & Identification des features variables
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)
seu <- RunUMAP(seu, dims = 1:20)

dup_treshold <- 0.2                         # ajustez si besoin
dupTumor <- WhichCells(
  seu,
  expression = cell_type == "tumor" &
    proportion_dupli_10 >= dup_treshold
)
length(dupTumor)      # combien de cellules retenues

DimPlot(
  seu,
  reduction = "umap",           # ou "pca", "tsne"
  cells.highlight = dupTumor,
  cols.highlight = "red",
  cols = "lightgrey"
) + ggtitle("Tumor cells with chr10 duplication ≥ 0.2")

go10 <- read.table("gene_order.txt", header = FALSE, sep = "\t")
colnames(go10) <- c("gene", "chr", "start", "end")
regions_file <- Sys.glob(file.path(opt$out_dir, "HMM_CNV_predictions.*.pred_cnv_regions.dat"))[1]
hmm10 <- read.table(regions_file, header = TRUE, sep = "\t")

# 4) Matrice logique : gain si state > 4 
gains_chr10 <- subset(hmm10, state >= 4)

genes_gain <- go10$gene[
  sapply(1:nrow(gains_chr10), function(i) {
    any(
      go10$start <= gains_chr10$end[i] &
        go10$end >= gains_chr10$start[i]
    )
  })
]


# Sauvegarde  
ggsave(
  file.path(opt$out_dir, "duplications_chr10.pdf"),
  width = 8, height = 4
)

###############################################################################
##  Résumé écrit + figure
###############################################################################
###############################################################################
##  RÉSUMÉ EN % + FIGURE
###############################################################################
library(data.table)
library(ggplot2)

###  A • Fonctions utilitaires ----------------------------------------------
percent_extra <- function(x) round(x / 2 * 100, 2)   # copies_extra → %

###  B • Table des métriques -------------------------------------------------
summary_df <- data.frame(
  métrique = c("Gain global chr10 (tumeur) [% ADN +]",
               "Gain moyen 35–125 Mb (tumeur – normal) [% ADN +]",
               paste0("Gain moyen chr10 – ", subclust, " [% ADN +]"),
               paste0("Gain moyen 35–127 Mb – ", subclust, " [% ADN +]"),
               "Cellules tumorales avec gain 35–127 Mb [% cellules]"),
  valeur = c(
    round(total_gain_percent, 2),                                        # déjà en %
    percent_extra(gain_vs_norm),                                         # copies → %
    percent_extra(mean(hmm_s7$copies_extra, na.rm = TRUE)),              # chr10
    percent_extra(mean(hmm_s7[start >= start_int*1e6 &
                                end   <= end_int*1e6]$copies_extra,
                       na.rm = TRUE)),                                   # fenêtre
    round(compter_cellules_avec_gain_region(hmm, meta,
                                            start_int, end_int)$pct_gain, 1)
  )
)

summary_df

file.create(file.path(opt$out_dir, ".done"))




