# === LIBRAIRIES ========================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("S4Arrays")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)      # <-- ajoute au yml Conda
  library(rhdf5); library(Matrix); library(infercnv)
  library(GenomicRanges); library(ggplot2); library(dplyr)
  library(IRanges); library(Seurat)
})

## ─── Définition des options ────────────────────────────────────────────
option_list <- list(
  make_option("--run1_file",      type="character"),
  make_option("--run2_file",      type="character"),
  make_option("--chrom",          type="character", default="10"),
  make_option("--max_cells",      type="integer",   default=5000),
  make_option("--min_reads",      type="integer",   default=100),
  make_option("--max_reads",      type="integer",   default=50000),
  make_option("--target_norm",    type="integer",   default=5000),
  make_option("--target_tum",     type="integer",   default=5000),
  make_option("--seed",           type="integer",   default=42),
  make_option("--threads",        type="integer",   default=4),
  make_option("--cutoff",        type="integer",   default=0),
  make_option("--workdir",        type="character", default=getwd()),
  make_option("--HMM",        type="character", default="i3"),
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
run1_file <- "RUN1_S1_hFF_WT.dna.h5"   # normal
run2_file <- "RUN2_S8_hFF_clone_6_KOfluo.dna.h5"        # tumor
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


# === DETECTION DES POSITIONS ============================================
go10 <- read.table("gene_order.txt", header = FALSE, sep = "\t")
colnames(go10) <- c("gene", "chr", "start", "end")
regions_file <- file.path(opt$out_dir,"HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
hmm10 <- read.table(regions_file, header = TRUE, sep = "\t")

# 4) Matrice logique : gain si state > 2  
gains_chr10 <- subset(hmm10, state == 3)

genes_gain <- go10$gene[
  sapply(1:nrow(gains_chr10), function(i) {
    any(
      go10$start <= gains_chr10$end[i] &
        go10$end >= gains_chr10$start[i]
    )
  })
]

go_gain <- go10[go10$gene %in% genes_gain, ]

"""# 1. Sélection des cellules dupliquées
cells_dupli <- rownames(meta_tumor)[meta_tumor$top_dupli_1 == 1]

# 2. Restreindre la matrice aux gènes de go_gain et cellules dupliquées
counts_gain_genes <- counts[rownames(counts) %in% go_gain$gene, colnames(counts) %in% cells_dupli]

# 3. Calcul du total de duplication (somme des comptages) pour chaque gène dans toutes les cellules tumorales dupliquées
total_duplication <- rowSums(counts_gain_genes)

# 4. Ajout dans le tableau go_gain
go_gain$total_duplication <- total_duplication"""













