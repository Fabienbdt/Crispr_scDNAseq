# =============================================
# SCRIPT FONCTIONNEL – Analyse CNV chr10 avec h5py
# =============================================

import h5py
import pandas as pd
import numpy as np

# === Chemins des fichiers ===
run1_path = "/net/cremi/redery/Projet/RUN1_S1_hFF_WT.dna.h5"
run2_path = "/net/cremi/redery/Projet/RUN2_S8_hFF_clone_6_KOfluo.dna.h5"
bed_path = "/net/cremi/redery/Projet/6969-amplicon.bed"

# === Chargement du fichier BED pour identifier les amplicons chr10 ===
bed = pd.read_csv(bed_path, sep="\t", header=None, names=["chr", "start", "end", "amplicon"])
amplicons_chr10 = bed[bed["chr"] == "chr10"]["amplicon"].tolist()

# === Fonction de chargement des données à partir du fichier .h5 ===
def load_counts(h5_path):
    with h5py.File(h5_path, 'r') as f:
        counts = f['assays/dna_read_counts/layers/read_counts'][:]
        barcodes = [b.decode() for b in f['assays/dna_read_counts/ra/barcode'][:]]
        amplicons = [a.decode() for a in f['assays/dna_read_counts/ca/id'][:]]
    df = pd.DataFrame(counts, index=barcodes, columns=amplicons)
    return df[[a for a in amplicons_chr10 if a in df.columns]]

# === Chargement des deux runs ===
df_wt = load_counts(run1_path)
df_crispr = load_counts(run2_path)

print(f"Nombre d'amplicons analysés sur chr10 : {df_wt.shape[1]}")

# === Fonction d'analyse : normalisation + détection de gain ===
def compute_gain_stats(df_chr10, threshold=1.5):
    norm = df_chr10.div(df_chr10.mean(axis=0), axis=1)
    means = norm.mean(axis=1)
    gain_cells = means[means > threshold]
    return norm, gain_cells.index.tolist(), 100 * len(gain_cells) / len(df_chr10)

# === Analyse CNV sur chr10 pour WT et CRISPR ===
norm_wt, gain_wt_ids, percent_wt = compute_gain_stats(df_wt)
norm_crispr, gain_crispr_ids, percent_crispr = compute_gain_stats(df_crispr)

# === Identification des amplicons avec gain dans la condition CRISPR ===
amplicon_gains = norm_crispr.loc[gain_crispr_ids].mean(axis=0)
amplicons_gain = amplicon_gains[amplicon_gains > 1.5].index.tolist()
bed_gain = bed[bed["amplicon"].isin(amplicons_gain)].sort_values(by="start").reset_index(drop=True)

# === Fusion des amplicons proches en régions continues (< 50kb) ===
regions = []
start = bed_gain.loc[0, "start"]
end = bed_gain.loc[0, "end"]
count = 1

for i in range(1, len(bed_gain)):
    s, e = bed_gain.loc[i, "start"], bed_gain.loc[i, "end"]
    if s - end < 50000:
        end = max(end, e)
        count += 1
    else:
        regions.append((int(start), int(end), count))
        start, end, count = s, e, 1
regions.append((int(start), int(end), count))

# === Affichage des résultats ===
print("\n--- Résultats CNV chr10 (lecture brute .h5) ---")
print(f"WT     : {len(gain_wt_ids)} cellules à gain / {len(df_wt)} total → {percent_wt:.2f}%")
print(f"CRISPR : {len(gain_crispr_ids)} cellules à gain / {len(df_crispr)} total → {percent_crispr:.2f}%")
print(f"Différence absolue : {percent_crispr - percent_wt:.2f}%")

print("\n➤ Régions de gain détectées sur chr10 :")
for s, e, n in regions:
    print(f" - chr10:{s}-{e}  ({n} amplicons)")

print(f"\n➤ Première position détectée : chr10:{int(bed_gain['start'].min())}")
print(f"➤ Dernière position détectée  : chr10:{int(bed_gain['end'].max())}")
