# Script R annoté en français pour upload sur GitLab

# Charger les bibliothèques nécessaires
library(karyotapR)        # Pour créer et manipuler des objets Tapestri
library(rhdf5)            # Pour lire les fichiers h5
library(SummarizedExperiment)  # Pour manipuler les objets SummarizedExperiment

# Définir les répertoires des données pour chaque RUN
run1_dir  <- file.path(getwd(), "RUN1_S1_hFF_WT")    # Répertoire de la première expérience
run2_dir  <- file.path(getwd(), "RUN2_S8_hFF_clone_6_KOfluo")  # Répertoire de la deuxième expérience

# Localiser le premier fichier .dna.h5 dans chaque répertoire
run1_h5 <- list.files(run1_dir, "\\.dna\\.h5$", full.names=TRUE)[1]
run2_h5 <- list.files(run2_dir, "\\.dna\\.h5$", full.names=TRUE)[1]

# Créer les objets TapestriExperiment pour chaque RUN
te_run1 <- createTapestriExperiment(run1_h5, panel.id = "CO261")  # Chargement RUN1
te_run2 <- createTapestriExperiment(run2_h5, panel.id = "CO261")  # Chargement RUN2

# Afficher le nombre de cellules et de sondes pour chaque RUN
message("RUN1 : ", dim(te_run1)[2], " cellules × ", dim(te_run1)[1], " sondes")
message("RUN2 : ", dim(te_run2)[2], " cellules × ", dim(te_run2)[1], " sondes")

# Définir une fonction de filtrage QC générique
qc_filter <- function(te, min_reads=500, max_mt=5){
  # Extraction des comptes de lectures et calcul du pourcentage mitochondrial
  cts <- assay(te, "counts")
  total_reads <- colData(te)$total.reads  # Nombre total de lectures par cellule
  mt_idx <- which(rowData(te)$chr %in% c("MT","chrM"))  # Indices des sondes mitochondriales
  pct_mito <- colSums(cts[mt_idx,,drop=FALSE]) / total_reads * 100  # Pourcentage de lectures mitochondriales
  
  # Sélectionner les cellules satisfaisant aux critères QC
  keep_cells <- which(total_reads >= min_reads & pct_mito <= max_mt)
  message("Nombre de cellules avant/après QC : ", ncol(te), " → ", length(keep_cells))
  return(te[,keep_cells])  # Retourner l'objet filtré
}

# Appliquer la QC à chaque RUN
te_run1_qc <- qc_filter(te_run1)  # QC pour RUN1
te_run2_qc <- qc_filter(te_run2)  # QC pour RUN2

# Normaliser les données via calcNormCounts
te_run1_norm <- calcNormCounts(te_run1_qc)  # Normalisation RUN1
te_run2_norm <- calcNormCounts(te_run2_qc)  # Normalisation RUN2

# Charger le fichier de conception du panel (design)
design_csv <- file.path(getwd(), "6969-design-summary.csv")  # Chemin vers le CSV de design

design_df <- read.csv(design_csv, stringsAsFactors = FALSE, check.names = FALSE)  # Lecture du design
head(design_df)  # Vérification rapide des premières lignes

# Calculer la médiane des comptes par probe pour chaque RUN
probe_med_run1 <- apply(assay(te_run1_norm), 1, median)  # Médiane RUN1
probe_med_run2 <- apply(assay(te_run2_norm), 1, median)  # Médiane RUN2

# Éviter les valeurs nulles pour le calcul du log2FC
probe_med_run1[probe_med_run1 < 1e-3] <- 1e-3
probe_med_run2[probe_med_run2 < 1e-3] <- 1e-3

# Calcul du log2 fold-change pour chaque cellule et chaque probe
log2fc_run1 <- log2((assay(te_run1_norm) + 1e-3) / probe_med_run1)
log2fc_run2 <- log2((assay(te_run2_norm) + 1e-3) / probe_med_run2)

# Extraire les informations de position des sondes à partir du design
probe_info <- as.data.frame(rowData(te_run1_norm))  # Informations initiales des sondes
probe_info$probe <- rownames(probe_info)  # Ajouter le nom de la sonde en colonne
probe_info <- merge(probe_info,
                    design_df[, c("AmpID","chr","amplicon_start","amplicon_end")],
                    by.x="probe", by.y="AmpID")  # Fusion avec le design

# Définir les régions centromérique et télomérique sur le chromosome 10
centro_region <- probe_info$chr.y %in% c("chr10","10") & 
  probe_info$amplicon_start >= 38000000 & probe_info$amplicon_end <= 48000000  # Région centromérique

telo_region <- probe_info$chr.y %in% c("chr10","10") & 
  probe_info$amplicon_start >= 120000000 & probe_info$amplicon_end <= 135000000  # Région télomérique

# Vérifier la définition des régions
table(centro_region)  # Nombre de sondes centromériques
table(telo_region)    # Nombre de sondes télomériques

# Calcul du log2FC moyen par cellule dans chaque région pour RUN1
run1_centro_log2fc <- colMeans(log2fc_run1[centro_region, , drop=FALSE], na.rm=TRUE)
run1_telo_log2fc   <- colMeans(log2fc_run1[telo_region, , drop=FALSE], na.rm=TRUE)

# Conversion en pourcentage de gain d'ADN moyen
run1_centro_gain_pct <- (2^mean(run1_centro_log2fc) - 1) * 100  # Gain moyen centromérique
run1_telo_gain_pct   <- (2^mean(run1_telo_log2fc) - 1) * 100    # Gain moyen télomérique

# Proportion de cellules avec gain significatif (log2FC > 1)
run1_centro_pct_cells <- sum(run1_centro_log2fc > 1) / length(run1_centro_log2fc) * 100
run1_telo_pct_cells   <- sum(run1_telo_log2fc > 1) / length(run1_telo_log2fc) * 100

# Répéter les calculs pour RUN2
run2_centro_log2fc <- colMeans(log2fc_run2[centro_region, , drop=FALSE], na.rm=TRUE)
run2_telo_log2fc   <- colMeans(log2fc_run2[telo_region, , drop=FALSE], na.rm=TRUE)

run2_centro_gain_pct <- (2^mean(run2_centro_log2fc) - 1) * 100
run2_telo_gain_pct   <- (2^mean(run2_telo_log2fc) - 1) * 100

run2_centro_pct_cells <- sum(run2_centro_log2fc > 1) / length(run2_centro_log2fc) * 100
run2_telo_pct_cells   <- sum(run2_telo_log2fc > 1) / length(run2_telo_log2fc) * 100

# Préparer le tableau récapitulatif des résultats pour chaque RUN
summary_runs <- data.frame(
  Region = c("Centromérique","Télomérique"),
  Run1_Gain_Moyen_pct = c(run1_centro_gain_pct, run1_telo_gain_pct),
  Run2_Gain_Moyen_pct = c(run2_centro_gain_pct, run2_telo_gain_pct),
  Run1_Cellules_pct = c(run1_centro_pct_cells, run1_telo_pct_cells),
  Run2_Cellules_pct = c(run2_centro_pct_cells, run2_telo_pct_cells)
)

# ----- Ajout de la section pour sauvegarder les résultats et générer des figures -----

# Créer le répertoire de sortie si nécessaire
out_dir <- file.path(getwd(), "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Générer et sauvegarder un barplot comparatif du pourcentage de cellules avec gain
png(file.path(out_dir, "gain_percentage_chr10_barplot.png"),
    width = 800, height = 600)

# Préparer les données pour le barplot
mat <- t(as.matrix(summary_runs[, c("Run1_Cellules_pct", "Run2_Cellules_pct")]))
colnames(mat) <- summary_runs$Region
rownames(mat) <- c("Run1", "Run2")

# Tracer le barplot
barplot(mat,
        beside = TRUE,
        col = c("steelblue", "tomato"),       # Couleurs pour distinguer les RUNs
        legend.text = rownames(mat),
        args.legend = list(x = "topright", bty = "n"),
        ylab = "Pourcentage de cellules avec gain (%)",
        main = "Comparaison du pourcentage de cellules avec gain\n(chr10 centromérique vs télomérique)",
        ylim = c(0, max(mat, na.rm = TRUE) * 1.1))

dev.off()
message("✅ Figure enregistrée : gain_percentage_chr10_barplot.png")

# 2. Sauvegarder le tableau récapitulatif en CSV
write.csv(summary_runs,
          file.path(out_dir, "comparaison_run1_run2_chr10.csv"),
          row.names = FALSE, quote = FALSE)
message("✅ CSV enregistré : comparaison_run1_run2_chr10.csv")