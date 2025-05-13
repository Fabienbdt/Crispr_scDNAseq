import os
import glob
import pandas as pd

# Chemin racine des résultats
base_dir = "results"

# Dictionnaire pour stocker les résultats par outil
summaries = {}

# Rechercher tous les fichiers CSV dans les sous-dossiers
for csv_path in glob.glob(os.path.join(base_dir, "*", "*.csv")):
    tool_name = os.path.basename(os.path.dirname(csv_path))  # ex: infercnv
    try:
        df = pd.read_csv(csv_path)
        summaries[tool_name] = df
    except Exception as e:
        print(f"⚠️ Erreur lors de la lecture de {csv_path}: {e}")

# Fusionner les résultats si possible
with open("results/comparison/summary.txt", "w") as f:
    for tool, df in summaries.items():
        f.write(f"=== Résultats pour {tool.upper()} ===\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")
