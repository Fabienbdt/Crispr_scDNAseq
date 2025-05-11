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
  make_option("--cutoff",        type="integer",   default=1),
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
  num_threads = opt$threads
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

###############################################################################
##  Pré-requis library(data.table)
library(GenomicRanges)
library(dplyr)
library(data.table)


### 1. charger les sous-groupes tumoraux ------------------------------
meta <- read.table(file.path(opt$out_dir, "map_metadata_from_infercnv.txt"))
tumor_groups <- meta$subcluster[grepl("^tumor_", meta$subcluster)]

### 2. lire les régions CNV ------------------------------------------
regions_file <- Sys.glob(file.path(
  opt$out_dir, "HMM_CNV_predictions.*.pred_cnv_regions.dat"))[1]
cnv <- fread(regions_file)

## ➜ créer un label sans le préfixe “tumor.” / “normal.”
cnv[, subcluster := sub("^.+\\.", "", cell_group_name)]   # garde “tumor_s16”

### 3. garder les gains sur chr10 dans les groupes tumoraux ----------
gains_chr10 <- cnv[state >= 4 & chr == 10 &
                     subcluster %in% tumor_groups]

### 4. définir la fenêtre 35–125 Mb et calculer la fraction ----------
roi      <- GRanges("10", IRanges(35e6, 127e6))
roi_len  <- width(roi)          # 90 000 001 bp

gain_gr  <- GRanges(seqnames = gains_chr10$chr,
                    ranges   = IRanges(gains_chr10$start, gains_chr10$end),
                    grp      = gains_chr10$subcluster)

ov        <- findOverlaps(gain_gr, roi)
ov_len    <- pmin(end(gain_gr)[queryHits(ov)],  end(roi)) -
  pmax(start(gain_gr)[queryHits(ov)], start(roi)) + 1

prop_by_grp <- data.frame(
  subcluster = mcols(gain_gr)$grp[queryHits(ov)],
  overlap_bp = ov_len
) %>% 
  group_by(subcluster) %>% 
  summarise(prop_gain_roi = sum(overlap_bp) / roi_len)   # fraction 0-1

print(prop_by_grp)

### 5. moyenne simple ou pondérée -----------------------------
mean_gain_pct <- mean(prop_by_grp$prop_gain_roi) * 100
cat(sprintf("\nGain moyen (35–125 Mb) = %.1f %%\n", mean_gain_pct))

# moyenne pondérée par le nombre de cellules, si dispo
if ("n_cells" %in% names(meta)) {
  prop_by_grp <- left_join(prop_by_grp,
                           meta[, c("subcluster", "n_cells")],
                           by = "subcluster")
  w_gain_pct <- weighted.mean(prop_by_grp$prop_gain_roi,
                              w = prop_by_grp$n_cells) * 100
  cat(sprintf("Gain moyen pondéré (par n_cells) = %.1f %%\n", w_gain_pct))
}

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

# ── Garder chr10 et fenêtre 35-125 Mb
hmm_roi <- hmm[chr == 10 & start >= 35e6 & end <= 125e6]

# ── Ajouter type de cellule (normal / tumor)
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

file.create("results/infercnv/.done")






