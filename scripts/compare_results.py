import os
import glob
import pandas as pd

# Chemin racine des résultats
base_dir = "results"
output_path = os.path.join(base_dir, "comparison", "summary.txt")
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# Dictionnaire pour stocker les résultats
summaries = {}

# Recherche uniquement les fichiers nommés final_compare.csv
for csv_path in glob.glob(os.path.join(base_dir, "*", "final_compare.csv")):
    tool_name = os.path.basename(os.path.dirname(csv_path))  # ex: infercnv_out
    try:
        df = pd.read_csv(csv_path)
        summaries[tool_name] = df
    except Exception as e:
        print(f"⚠️ Erreur lors de la lecture de {csv_path}: {e}")

# Écrire la synthèse dans un fichier texte
with open(output_path, "w") as f:
    if summaries:
        for tool, df in summaries.items():
            f.write(f"=== Résultats pour {tool.upper()} ===\n")
            f.write(df.to_string(index=False))
            f.write("\n\n")
    else:
        f.write("Aucun fichier final_compare.csv trouvé dans les sous-dossiers de 'results/'.\n")
        print("⚠️ Aucun fichier final_compare.csv trouvé.")

done_path = os.path.join(os.path.dirname(args.output), ".done")
with open(done_path, "w") as f:
    f.write("done\n")

