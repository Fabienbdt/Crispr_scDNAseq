#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(karyotapR)
  library(rhdf5)
  library(SummarizedExperiment)
})

# --- Définir les options ---
option_list <- list(
  make_option("--run1_file", type = "character", default = "data/RUN1_S1_hFF_WT.dna.h5"),
  make_option("--run2_file", type = "character", default = "data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5"),
  make_option("--panel",     type = "character", default = "CO261"),
  make_option("--design",    type = "character", default = "data/6969-design-summary.csv"),
  make_option("--out_dir",   type = "character", default = "results/karyotapr")
)

opt <- parse_args(OptionParser(option_list = option_list))
unlink(opt$out_dir, recursive = TRUE)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Charger les fichiers HDF5 directement ---
te_run1 <- createTapestriExperiment(opt$run1_file, panel.id = opt$panel)
te_run2 <- createTapestriExperiment(opt$run2_file, panel.id = opt$panel)

# --- QC & Normalisation ---
qc_filter <- function(te, min_reads = 500, max_mt = 5) {
  cts <- assay(te, "counts")
  total_reads <- colData(te)$total.reads
  mt_idx <- which(rowData(te)$chr %in% c("MT", "chrM"))
  pct_mito <- colSums(cts[mt_idx, , drop = FALSE]) / total_reads * 100
  keep_cells <- which(total_reads >= min_reads & pct_mito <= max_mt)
  return(te[, keep_cells])
}

te_run1_norm <- calcNormCounts(qc_filter(te_run1))
te_run2_norm <- calcNormCounts(qc_filter(te_run2))

# --- Charger le design ---
design_df <- read.csv(opt$design, stringsAsFactors = FALSE, check.names = FALSE)

# --- Log2FC ---
med1 <- pmax(apply(assay(te_run1_norm), 1, median), 1e-3)
med2 <- pmax(apply(assay(te_run2_norm), 1, median), 1e-3)
log2fc1 <- log2((assay(te_run1_norm) + 1e-3) / med1)
log2fc2 <- log2((assay(te_run2_norm) + 1e-3) / med2)

# --- Info des sondes & filtrage chr10 ---
probe_info <- as.data.frame(rowData(te_run1_norm))
probe_info$probe <- rownames(probe_info)
probe_info <- merge(probe_info, design_df[, c("AmpID", "chr", "amplicon_start", "amplicon_end")],
                    by.x = "probe", by.y = "AmpID")
is_chr10 <- probe_info$chr.y %in% c("chr10", "10")
centro <- is_chr10 & probe_info$amplicon_start >= 38e6 & probe_info$amplicon_end <= 48e6
telo   <- is_chr10 & probe_info$amplicon_start >= 120e6 & probe_info$amplicon_end <= 135e6

# --- Métriques ---
region_metrics <- function(log2fc, region) {
  avg <- (2^mean(colMeans(log2fc[region, , drop = FALSE], na.rm = TRUE)) - 1) * 100
  pct <- sum(colMeans(log2fc[region, , drop = FALSE], na.rm = TRUE) > 1) / ncol(log2fc) * 100
  return(c(avg, pct))
}
run1 <- cbind(centro = region_metrics(log2fc1, centro),
              telo   = region_metrics(log2fc1, telo))
run2 <- cbind(centro = region_metrics(log2fc2, centro),
              telo   = region_metrics(log2fc2, telo))

# --- Résumé & sauvegarde ---
summary_df <- data.frame(
  Region = c("Centromérique", "Télomérique"),
  Run1_Gain_Moyen_pct = run1[1, ],
  Run2_Gain_Moyen_pct = run2[1, ],
  Run1_Cellules_pct   = run1[2, ],
  Run2_Cellules_pct   = run2[2, ]
)


write.csv(summary_df, file.path(opt$out_dir, "karyotapR_results.csv"), row.names = FALSE, quote = FALSE)

# --- Figure ---
png(file.path(opt$out_dir, "gain_percentage_chr10_barplot.png"), width = 800, height = 600)
mat <- t(as.matrix(summary_df[, c("Run1_Cellules_pct", "Run2_Cellules_pct")]))
colnames(mat) <- summary_df$Region
barplot(mat, beside = TRUE, col = c("steelblue", "tomato"),
        legend.text = rownames(mat), args.legend = list(x = "topright", bty = "n"),
        ylab = "Pourcentage de cellules avec gain (%)",
        main = "Comparaison gain chr10 centromérique / télomérique")
dev.off()
