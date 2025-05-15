################################################################################
## Snakefile – exécute InferCNV (via Docker), Mosaic, KaryotapR, puis compare ##
################################################################################

configfile: "config.yaml"
SCRIPTS = config["scripts"]

def build_args(tool):
    return " ".join(f"--{k} {v}" for k, v in config["params"].get(tool, {}).items())

rule all:
    input: "results/comparison/summary.txt"

# ──────────────── Rule pour InferCNV (via Docker) ────────────────────────
rule infercnv:
    input:
        script = SCRIPTS["infercnv"]
    output:
        "results/infercnv/.done"
    params:
        args = build_args("infercnv")
    container:
        "crispr_infercnv:latest"
    shell:
        """
        mkdir -p results/infercnv
        Rscript {input.script} {params.args} && touch {output}
        """

# ──────────────── Rule générique pour Mosaic & KaryotapR ────────────────
rule run_tool:
    input:
        script = lambda wc: SCRIPTS[wc.label]
    output:
        "results/{label}/.done"
    params:
        args = lambda wc: build_args(wc.label)
    conda:
        lambda wc: "envs/r.yml"
    shell:
        """
        mkdir -p results/{wildcards.label}
        {{ 'Rscript' if input.script.endswith('.R') else 'python' }} {input.script} {params.args} && touch {output}
        """

# ──────────────── Comparaison finale, même si certains résultats manquent ────────────────
rule compare:
    input:
        script = SCRIPTS["compare"]
    output:
        "results/comparison/summary.txt"
    conda:
        "envs/r.yml"
    shell:
        """
        mkdir -p results/comparison
        python {input.script} > {output} || echo '[⚠️ compare_results.py terminé avec avertissements]'
        """
