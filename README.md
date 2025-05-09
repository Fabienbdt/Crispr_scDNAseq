# Snake — Pipeline Snakemake pour exécuter et comparer des scripts R



## Ce projet fournit un pipeline reproductible basé sur **Snakemake**. Il exécute automatiquement plusieurs scripts R présents dans `scripts/`, accepte des paramètres personnalisables via `config.yaml`, puis compare leurs résultats

## Structure du projet

crispr/snake/
├── config.yaml # Configuration des scripts à exécuter
├── envs/ # Environnement Conda requis pour R
│ └── r.yml
├── scripts/ # Scripts R (et Python) à exécuter
│ ├── script1.R
│ ├── script2.R
│ ├── script3.R
│ └── compare_results.py
├── results/ # Résultats générés automatiquement
├── snakefile # Fichier Snakefile principal
└── README.md # Ce fichier


## Prérequis

* Python ≥ 3.8  
* [Conda / Mamba](https://docs.conda.io/en/latest/)  
* Snakemake ≥ 7  
* Accès SSH au dépôt GitLab (si privé)

Installation rapide de Snakemake :

```bash
conda create -n snake_env snakemake -c bioconda -c conda-forge
conda activate snake_env
## Collaborate with your team
```

##  Étape 1 — Configurer config.yaml
Modifie le fichier config.yaml pour lister tes scripts et leurs arguments :

```yaml
scripts:
  A: "scripts/script1.R"
  B: "scripts/script2.R"
  C: "scripts/script3.R"

params:
  A: "--input dataA.csv"
  B: "--threshold 0.05"
  C: "--mode full"
```

scripts: : clé = label (dossier de sortie), valeur = chemin vers le script R.

params: (optionnel) : arguments à passer au script correspondant.

## Étape 2 Lancer le pipeline

Depuis le dossier contenant le snakefile, exécute :
```bash
snakemake --use-conda --cores 4
```

## Étape 3 — Consulter les résultats
Chaque script est exécuté dans results/<label>/ (ex. results/A/).

Le résumé de comparaison est généré dans :
results/comparison/summary.txt

## Personnalisation
Ajoute ou remplace des scripts R dans scripts/.

Ajuste compare_results.py pour des analyses/graphes spécifiques.

Modifie envs/r.yml si tes scripts nécessitent d’autres packages R.

## Crédits & Contact
Développé par Fabien Bidet ; Raphael Edery ; Ziyi Zhao ; Tom Bourrachot — Université de Bordeaux
✉️ fabien.bidet@etu.u-bordeaux.fr


