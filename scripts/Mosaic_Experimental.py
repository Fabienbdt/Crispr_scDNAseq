# =============================================
# SCRIPT EXPÉRIMENTAL – ANALYSE CNV CHR10 AVEC MOSAIC SDK (SI VERSION COMPLÈTE)
# =============================================

import argparse
import missionbio.mosaic as ms
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# === Argument parsing ===
parser = argparse.ArgumentParser(description="Analyse CNV chr10 via Mosaic SDK complet (expérimental)")
parser.add_argument("--wt", default="data/RUN1_S1_hFF_WT.dna.h5", help="Fichier .h5 WT")
parser.add_argument("--crispr", default="data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5", help="Fichier .h5 CRISPR")
parser.add_argument("--bed", default="data/6969-amplicon.bed", help="Fichier .bed des amplicons")
parser.add_argument("--threshold", type=float, default=1.5)
parser.add_argument("--distance", type=int, default=50000)
parser.add_argument("--output", default="mosaic_experimental", help="Fichier résumé .txt")
args = parser.parse_args()

# === Chargement des fichiers .h5 (avec SDK complet) ===
sample_wt = ms.load(args.wt, raw=True, filter_variants=False, filter_cells=False, whitelist=[], single=True)
sample_crispr = ms.load(args.crispr, raw=True, filter_variants=False, filter_cells=False, whitelist=[], single=True)

# === Extraction des données CNV ===
df_wt = sample_wt.dna.counts
df_crispr = sample_crispr.dna.counts

# === Filtrage chr10 via .bed ===
bed = pd.read_csv(args.bed, sep="\t", header=None, names=["chr", "start", "end", "amplicon"])
amplicons_chr10 = bed[bed["chr"] == "chr10"]["amplicon"].tolist()
amplicons_common = [a for a in amplicons_chr10 if a in df_wt.columns and a in df_crispr.columns]

df_wt_chr10 = df_wt[amplicons_common]
df_crispr_chr10 = df_crispr[amplicons_common]

# === Analyse CNV ===
def compute_gain(df, threshold):
    norm = df.div(df.mean(axis=0), axis=1)
    mean = norm.mean(axis=1)
    gain_cells = mean[mean > threshold]
    return norm, gain_cells.index.tolist(), 100 * len(gain_cells) / len(df)

norm_wt, ids_wt, pct_wt = compute_gain(df_wt_chr10, args.threshold)
norm_crispr, ids_crispr, pct_crispr = compute_gain(df_crispr_chr10, args.threshold)

# === Moyenne CNV par amplicon (CRISPR) ===
mean_amplicon = norm_crispr.loc[ids_crispr].mean(axis=0)
amplicons_gain = mean_amplicon[mean_amplicon > args.threshold].index.tolist()
bed_gain = bed[bed["amplicon"].isin(amplicons_gain)].sort_values(by="start").reset_index(drop=True)

# === Regroupement en régions continues ===
regions = []
if not bed_gain.empty:
    start = bed_gain.loc[0, "start"]
    end = bed_gain.loc[0, "end"]
    count = 1
    for i in range(1, len(bed_gain)):
        s, e = bed_gain.loc[i, "start"], bed_gain.loc[i, "end"]
        if s - end < args.distance:
            end = max(end, e)
            count += 1
        else:
            regions.append((int(start), int(end), count))
            start, end, count = s, e, 1
    regions.append((int(start), int(end), count))
# === Export résultats texte ===
os.makedirs(os.path.dirname(args.output), exist_ok=True)
with open(args.output, "w") as f:
    f.write("Résultats Mosaic CNV (expérimental)\n")
    f.write(f"Amplicons analysés : {len(amplicons_common)}\n")
    f.write(f"WT     : {len(ids_wt)} / {len(df_wt_chr10)} → {pct_wt:.2f}%\n")
    f.write(f"CRISPR : {len(ids_crispr)} / {len(df_crispr_chr10)} → {pct_crispr:.2f}%\n")
    f.write(f"Différence : {pct_crispr - pct_wt:.2f}%\n\n")
    f.write("Régions détectées :\n")
    for s, e, n in regions:
        f.write(f"chr10:{s}-{e} ({n} amplicons)\n")
    if not bed_gain.empty:
        f.write(f"\nPremière position : chr10:{int(bed_gain['start'].min())}\n")
        f.write(f"Dernière position : chr10:{int(bed_gain['end'].max())}\n")

# === Export CSV pour intégration finale ===
summary_df = pd.DataFrame({
    "Metric": [
        "Amplicons analysés",
        "Cellules WT avec gain (%)",
        "Cellules CRISPR avec gain (%)",
        "Différence (%)",
        "Nombre de régions détectées"
    ],
    "Value": [
        len(amplicons_common),
        round(pct_wt, 2),
        round(pct_crispr, 2),
        round(pct_crispr - pct_wt, 2),
        len(regions)
    ]
})

csv_output_path = os.path.join(os.path.dirname(args.output), "final_compare.csv")
summary_df.to_csv(csv_output_path, index=False)
print(f"✅ CSV enregistré : {csv_output_path}")

done_path = os.path.join(os.path.dirname(args.output), ".done")
with open(done_path, "w") as f:
    f.write("done\n")

