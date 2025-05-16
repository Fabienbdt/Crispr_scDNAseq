# =============================================
# SCRIPT FONCTIONNEL – Analyse CNV chr10 avec h5py
# =============================================

import h5py
import pandas as pd
import numpy as np
import argparse
import os

# === Argument parsing ===
parser = argparse.ArgumentParser(description="Analyse CNV sur chr10 depuis fichiers .h5")
parser.add_argument("--wt", default="data/RUN1_S1_hFF_WT.dna.h5", help="Fichier .h5 de l'échantillon WT")
parser.add_argument("--crispr", default="data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5", help="Fichier .h5 de l'échantillon CRISPR")
parser.add_argument("--bed", default="6969-amplicon.bed", help="Fichier BED des amplicons")
parser.add_argument("--threshold", type=float, default=1.5, help="Seuil de gain CNV")
parser.add_argument("--distance", type=int, default=50000, help="Distance max entre amplicons pour les regrouper")
parser.add_argument("--output", default="Analyse_manuelle_par h5py", help="Fichier texte de sortie")
args = parser.parse_args()

# === Chargement du fichier BED ===
bed = pd.read_csv(args.bed, sep="\t", header=None, names=["chr", "start", "end", "amplicon"])
amplicons_chr10 = bed[bed["chr"] == "chr10"]["amplicon"].tolist()

# === Fonction pour lire un fichier .h5 ===
def load_counts(h5_path):
    with h5py.File(h5_path, 'r') as f:
        counts = f['assays/dna_read_counts/layers/read_counts'][:]
        barcodes = [b.decode() for b in f['assays/dna_read_counts/ra/barcode'][:]]
        amplicons = [a.decode() for a in f['assays/dna_read_counts/ca/id'][:]]
    df = pd.DataFrame(counts, index=barcodes, columns=amplicons)
    return df[[a for a in amplicons_chr10 if a in df.columns]]

# === Chargement des données ===
df_wt = load_counts(args.wt)
df_crispr = load_counts(args.crispr)

# === Fonction d'analyse CNV ===
def compute_gain_stats(df, threshold):
    norm = df.div(df.mean(axis=0), axis=1)
    means = norm.mean(axis=1)
    gain_cells = means[means > threshold]
    return norm, gain_cells.index.tolist(), 100 * len(gain_cells) / len(df)

norm_wt, gain_ids_wt, pct_wt = compute_gain_stats(df_wt, args.threshold)
norm_crispr, gain_ids_crispr, pct_crispr = compute_gain_stats(df_crispr, args.threshold)

# === Identification des amplicons gagnés (CRISPR) ===
amplicon_mean = norm_crispr.loc[gain_ids_crispr].mean(axis=0)
amplicons_gain = amplicon_mean[amplicon_mean > args.threshold].index.tolist()
bed_gain = bed[bed["amplicon"].isin(amplicons_gain)].sort_values(by="start").reset_index(drop=True)

# === Regroupement des régions gagnées ===
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

# === Sauvegarde des résultats texte ===
os.makedirs(os.path.dirname(args.output), exist_ok=True)
with open(args.output, "w") as f:
    f.write("Résultats CNV chr10 (pipeline fonctionnel)\n")
    f.write(f"Nombre d'amplicons analysés : {df_crispr.shape[1]}\n")
    f.write(f"WT     : {len(gain_ids_wt)} / {len(df_wt)} → {pct_wt:.2f}%\n")
    f.write(f"CRISPR : {len(gain_ids_crispr)} / {len(df_crispr)} → {pct_crispr:.2f}%\n")
    f.write(f"Différence : {pct_crispr - pct_wt:.2f}%\n\n")

    f.write("Régions de gains détectées (CRISPR) :\n")
    for s, e, n in regions:
        f.write(f"chr10:{s}-{e} ({n} amplicons)\n")

    if not bed_gain.empty:
        f.write(f"\nPremière position : chr10:{int(bed_gain['start'].min())}\n")
        f.write(f"Dernière position : chr10:{int(bed_gain['end'].max())}\n")

# === Export CSV pour comparaison finale ===
df_summary = pd.DataFrame({
    "Metric": [
        "Amplicons analysés",
        "Cellules WT avec gain (%)",
        "Cellules CRISPR avec gain (%)",
        "Différence (%)",
        "Nombre de régions détectées"
    ],
    "Value": [
        df_crispr.shape[1],
        round(pct_wt, 2),
        round(pct_crispr, 2),
        round(pct_crispr - pct_wt, 2),
        len(regions)
    ]
})

csv_output_path = os.path.join(os.path.dirname(args.output), "final_compare.csv")
df_summary.to_csv(csv_output_path, index=False)
print(f"✅ CSV enregistré : {csv_output_path}")


done_path = os.path.join(os.path.dirname(args.output), ".done")
with open(done_path, "w") as f:
    f.write("done\n")
