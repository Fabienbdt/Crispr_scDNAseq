# =============================================
# SCRIPT EXPÉRIMENTAL – ANALYSE CNV CHR10 AVEC MOSAIC SDK (SI VERSION COMPLÈTE)
# =============================================

# === 1. Import des bibliothèques nécessaires ===
import missionbio.mosaic as ms              # SDK officiel Mosaic
import pandas as pd                         # Manipulation des matrices
import numpy as np                          # Outils mathématiques
import seaborn as sns                       # Visualisation
import matplotlib.pyplot as plt             # Affichage des graphiques

# === 2. Chargement des fichiers .h5 ===
# On suppose ici qu’on a accès à la version complète du SDK Mosaic
# On désactive tous les filtres automatiques (variants, cellules, etc.)

wt_path = "/net/cremi/redery/Projet/RUN1_S1_hFF_WT.dna.h5"
crispr_path = "/net/cremi/redery/Projet/RUN2_S8_hFF_clone_6_KOfluo.dna.h5"

sample_wt = ms.load(wt_path, raw=True, filter_variants=False, filter_cells=False, whitelist=[], single=True)
sample_crispr = ms.load(crispr_path, raw=True, filter_variants=False, filter_cells=False, whitelist=[], single=True)

# === 3. Extraction des données CNV ===
# On récupère les barcodes et les amplicons sous forme de DataFrame
# Ces objets sont normalement accessibles dans le SDK complet

df_wt = sample_wt.dna.counts               # Matrice de lecture : cellules x amplicons
df_crispr = sample_crispr.dna.counts

amplicons_wt = df_wt.columns
amplicons_crispr = df_crispr.columns

# === 4. Filtrage du chromosome 10 ===
# On isole les amplicons du chromosome 10 en utilisant leur nom

chr10_amplicons = [a for a in amplicons_wt if a.startswith("chr10") and a in amplicons_crispr]

df_wt_chr10 = df_wt[chr10_amplicons]
df_crispr_chr10 = df_crispr[chr10_amplicons]

# === 5. Normalisation des lectures par amplicon ===
# Chaque valeur est divisée par la moyenne de son amplicon → permet de corriger les biais de couverture

norm_wt = df_wt_chr10.div(df_wt_chr10.mean(axis=0), axis=1)
norm_crispr = df_crispr_chr10.div(df_crispr_chr10.mean(axis=0), axis=1)

# === 6. Calcul de la moyenne CNV par cellule ===

mean_cnv_wt = norm_wt.mean(axis=1)
mean_cnv_crispr = norm_crispr.mean(axis=1)

# === 7. Détection des cellules avec gain ===
# On considère une cellule comme "à gain" si sa moyenne CNV dépasse 1.5

threshold = 1.5
gain_cells_wt = mean_cnv_wt[mean_cnv_wt > threshold]
gain_cells_crispr = mean_cnv_crispr[mean_cnv_crispr > threshold]

# === 8. Calcul du pourcentage de cellules avec gain ===

pct_wt = 100 * len(gain_cells_wt) / len(norm_wt)
pct_crispr = 100 * len(gain_cells_crispr) / len(norm_crispr)

print("Résultats CNV chr10 (avec SDK Mosaic complet simulé)")
print(f"WT     : {len(gain_cells_wt)} cellules à gain / {len(norm_wt)} total → {pct_wt:.2f}%")
print(f"CRISPR : {len(gain_cells_crispr)} cellules à gain / {len(norm_crispr)} total → {pct_crispr:.2f}%")
print(f"Différence absolue : {pct_crispr - pct_wt:.2f}%")

# === 9. Identification des amplicons significativement affectés ===
# On regarde les amplicons où la moyenne dépasse le seuil chez les cellules à gain

gain_amplicons_crispr = norm_crispr.loc[gain_cells_crispr.index].mean(axis=0)
gain_amplicons_crispr = gain_amplicons_crispr[gain_amplicons_crispr > threshold]

print(f"{len(gain_amplicons_crispr)} amplicons gagnés détectés dans CRISPR")

# === 10. Heatmap des cellules avec gain (optionnelle) ===

plt.figure(figsize=(12, 6))
sns.heatmap(norm_crispr.loc[gain_cells_crispr.index], cmap="coolwarm", center=1.0)
plt.title("Heatmap CNV normalisé – cellules à gain (chr10 – CRISPR)")
plt.xlabel("Amplicons chr10")
plt.ylabel("Cellules CRISPR avec gain")
plt.tight_layout()
plt.show()
