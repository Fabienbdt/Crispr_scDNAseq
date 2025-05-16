![Last Commit](https://img.shields.io/github/last-commit/Fabienbdt/Crispr_scDNAseq)

🔬 Crispr_scDNAseq — Pipeline Snakemake pour l’analyse CNV sur données scDNA-seq
Ce dépôt propose un pipeline modulaire et reproductible basé sur Snakemake, combinant plusieurs outils R et Python pour l’analyse des variations du nombre de copies (CNV) à partir de données Mission Bio Tapestri.

Il prend en charge :

* infercnv (via Docker),

* karyotapR,

* Deux versions de Mosaic (expérimentale et fonctionnelle),

* Une phase de comparaison automatique des résultats.
---

##  Structure du projet

```bash
Crispr_scDNAseq/
├── Snakefile
├── config.yaml
├── data/
│   ├── RUN1_S1_hFF_WT.dna.h5
│   ├── RUN2_S8_hFF_clone_6_KOfluo.dna.h5
│   ├── 6969-amplicon.bed
│   └── 6969-design-summary.csv
├── envs/
│   └── r.yml
├── scripts/
│   ├── run_infercnv.R
│   ├── run_karyotapR.R
│   ├── Analyse_CNV_Manuelle.py
│   ├── Mosaic_Experimental.py
│   └── compare_results.py
├── results/
└── Dockerfile  # utilisé pour exécuter infercnv

```
## Prérequis

* Python ≥ 3.8  
* [Conda / Mamba](https://docs.conda.io/en/latest/)  
* Snakemake ≥ 7  
* Accès SSH au dépôt GitHub

## 1. Cloner le dépôt

```bash
git clone https://github.com/Fabienbdt/Crispr_scDNAseq.git
cd Crispr_scDNAseq
```


### 2.Ajouter vos fichiers .h5 (non versionnés)
Créez le dossier data/ et ajoutez vos fichiers .h5 :

```bash
mkdir -p data/
cp /chemin/vers/tes/fichiers/*.h5 data/
```
Le dossier `data/` est ignoré par Git grâce à `.gitignore`, pour éviter de versionner des fichiers volumineux.  
Assurez-vous que les chemins spécifiés dans `config.yaml` pointent bien vers vos fichiers `.h5` dans ce dossier.




### 3. Installation de Snakemake

```bash
conda create -n snake_env snakemake -c bioconda -c conda-forge
conda activate snake_env
```

### 4. Environnements spécifiques

```bash

conda env create -f envs/r.yml
conda env create -f envs/python_h5.yml
conda env create -f envs/python_mosaic.yml

```

### 5. Docker pour infercnv (obligatoire)
🐳 Construction du conteneur Docker pour infercnv
```bash

docker build -t crispr_infercnv .

```
L’image inclut R 4.3, infercnv, rhdf5 et toutes ses dépendances. Elle sera utilisée automatiquement par Snakemake pour cette tâche.



## Configurer config.yaml
Modifie le fichier config.yaml pour lister tes scripts et leurs arguments :

```yaml
scripts:
  infercnv: "scripts/run_infercnv.R"
  karyotapr: "scripts/run_karyotapR.R"
  mosaic_functional: "scripts/Analyse_CNV_Manuelle.py"
  mosaic_experimental: "scripts/Mosaic_Experimental.py"
  compare: "scripts/compare_results.py"

params:
  infercnv:
    run1_file: "data/RUN1_S1_hFF_WT.dna.h5"
    run2_file: "data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5"
    chrom: "10"
    out_dir: "results/infercnv"
  karyotapr:
    run1_file: "..."
    design: "data/6969-design-summary.csv"
  ...

```

##  Lancer le pipeline

Depuis le dossier contenant le snakefile, exécute :
```bash
snakemake --use-conda --cores 4
```
⚠️ --use-singularity fonctionne aussi avec Docker sur les systèmes disposant de Docker Desktop.

### TEST D'UN SCRIPTS UNIQUE  : EX avec run_infercnv.R:

```bash

docker run --rm -v $(pwd):/work -w /work crispr_infercnv \
  Rscript scripts/run_infercnv.R \
    --NormalCellFile data/RUN1_S1_hFF_WT.dna.h5 \
    --TumorCellFile  data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5 \
    --out_dir results/infercnv


```

## Étape 3 — Consulter les résultats
Chaque outil écrit ses résultats dans results/<outil>/.

S’il réussit, un fichier final_compare.csv est généré.

La commande :
```bash
scripts/compare_results.py
```

génère automatiquement :

```bash
results/comparison/summary.txt
```

## résultats et comparaison

Chaque outil écrit ses résultats dans results/<outil>/

Un fichier standardisé final_compare.csv est généré pour chaque outil (si tout se passe bien)

Un script Python (compare_results.py) centralise automatiquement tous les final_compare.csv disponibles et génère :
```bash
results/comparison/summary.txt
```

##  Personnalisation
Ajoute tes propres scripts dans le dossier scripts/

Adapte les paramètres dans config.yaml

Le Dockerfile peut être étendu si des paquets R supplémentaires sont nécessaires


## Crédits & Contact
👥 Crédits
Développé par :
📍 Fabien Bidet
📍 Raphael Edery
📍 Ziyi Zhao
📍 Tom Bourrachot
Université de Bordeaux

✉️ Contact : fabien.bidet@etu.u-bordeaux.fr

