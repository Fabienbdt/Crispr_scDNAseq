![Last Commit](https://img.shields.io/github/last-commit/Fabienbdt/Crispr_scDNAseq)

# Snake — Pipeline Snakemake pour exécuter et comparer des scripts R



## Ce projet fournit un pipeline reproductible basé sur **Snakemake**. Il exécute automatiquement plusieurs scripts R présents dans `scripts/`, accepte des paramètres personnalisables via `config.yaml`, puis compare leurs résultats

## Structure du projet
```arduino
Snake/
├── Snakefile
├── config.yaml
├── envs/
│   └── r.yml
├── scripts/
│   ├── run_infercnv.R
│   ├── karyotapR.R
│   ├── mosaic.py
│   └── compare_results.py
└── results/
```

## Prérequis

* Python ≥ 3.8  
* [Conda / Mamba](https://docs.conda.io/en/latest/)  
* Snakemake ≥ 7  
* Accès SSH au dépôt GitHub

### Cloner le dépôt

```bash
git clone https://github.com/Fabienbdt/Crispr_scDNAseq.git
cd Crispr_scDNAseq
```


### Ajouter vos fichiers .h5 (non versionnés)
Créez le dossier data/ et ajoutez vos fichiers .h5 :

```bash
mkdir -p data/
cp /chemin/vers/tes/fichiers/*.h5 data/
```
Le dossier `data/` est ignoré par Git grâce à `.gitignore`, pour éviter de versionner des fichiers volumineux.  
Assurez-vous que les chemins spécifiés dans `config.yaml` pointent bien vers vos fichiers `.h5` dans ce dossier.




Installation rapide de Snakemake (si besoin):

```bash
conda create -n snake_env snakemake -c bioconda -c conda-forge
conda activate snake_env
```

Puis installe l’environnement R :

```bash

conda env create -f envs/r.yml
conda activate snake_env
```

### Pour l'utilisation de inferCNV sur différents système d'exploitation : 

🐳 Construction du conteneur Docker pour infercnv

Installez docker si message d'erreur du type : zsh: command not found: docker
```bash

docker build -t crispr_infercnv .

```
L’image inclut R 4.3, infercnv, rhdf5 et toutes ses dépendances. Elle sera utilisée automatiquement par Snakemake pour cette tâche.



##  Étape 1 — Configurer config.yaml
Modifie le fichier config.yaml pour lister tes scripts et leurs arguments :

```yaml
scripts:
  infercnv:  "scripts/run_infercnv.R"
  mosaic:    "scripts/run_mosaic.py"
  karyotapr: "scripts/run_karyotapR.R"


params:
  infercnv: "--input dataA.csv"
  mosaic: "--threshold 0.05"
  karyotapr: "--mode full"
```

scripts: : clé = label (dossier de sortie), valeur = chemin vers le script R.

params: (optionnel) : arguments à passer au script correspondant.

## Étape 2 Lancer le pipeline

Depuis le dossier contenant le snakefile, exécute :
```bash
snakemake --use-singularity --cores 4
```
⚠️ --use-singularity fonctionne aussi avec Docker sur les systèmes disposant de Docker Desktop.

## Étape 3 — Consulter les résultats
Chaque script est exécuté dans results/<label>/ (ex. results/A/).

Le résumé de comparaison est généré dans :
results/comparison/summary.txt

### TEST D'UN SCRIPTS UNIQUE  : EX avec run_infercnv.R:



```bash

docker run --rm -v $(pwd):/work -w /work crispr_infercnv \
  Rscript scripts/run_infercnv.R \
    --NormalCellFile data/RUN1_S1_hFF_WT.dna.h5 \
    --TumorCellFile  data/RUN2_S8_hFF_clone_6_KOfluo.dna.h5 \
    --out_dir results/infercnv


```
## Personnalisation
Ajoute ou remplace des scripts R dans scripts/.

Ajuste compare_results.py pour des analyses/graphes spécifiques.

Adaptez le Dockerfile si vos scripts requièrent d’autres packages


## Crédits & Contact
Développé par Fabien Bidet ; Raphael Edery ; Ziyi Zhao ; Tom Bourrachot — Université de Bordeaux
✉️ fabien.bidet@etu.u-bordeaux.fr


