configfile: "config.yaml"
SCRIPTS = config["scripts"]

def build_args(tool):
    return " ".join(f"--{k} {v}" for k, v in config["params"].get(tool, {}).items())

rule all:
    input:
        "results/comparison/summary.txt"

# ────────────────────────────
# INFERCNV (via Docker)
# ────────────────────────────
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
        Rscript {input.script} {params.args}
        touch {output}
        """

# ────────────────────────────
# Analyse_CNV_Manuelle (via python_h5)
# ────────────────────────────
rule python_h5:
    input:
        script = config["scripts"]["python_h5"]
    output:
        "results/python_h5/.done"
    params:
        args = build_args("python_h5")
    conda:
        "envs/mosaic.yml"
    shell:
        """
        mkdir -p results/python_h5
        python {input.script} {params.args}
        touch {output}
        """


# ────────────────────────────
# MOSAIC EXPERIMENTAL (via python_mosaic)
# ────────────────────────────
rule manual_analysis:
    input:
        script = SCRIPTS["mosaic_experimental"]
    output:
        "results/mosaic_experimental/.done"
    params:
        args = build_args("mosaic_experimental")
    conda:
        "envs/python_mosaic.yml"
    shell:
        """
        mkdir -p results/mosaic_experimental
        python {input.script} {params.args}
        touch {output}
        """

# ────────────────────────────
# KARYOTAPR (via R)
# ────────────────────────────
rule karyotapr:
    input:
        script = SCRIPTS["karyotapr"]
    output:
        "results/karyotapr/.done"
    params:
        args = build_args("karyotapr")
    conda:
        "envs/r.yml"
    shell:
        """
        mkdir -p results/karyotapr
        Rscript {input.script} {params.args}
        touch {output}
        """

# ────────────────────────────
# COMPARAISON (via R ou Python)
# ────────────────────────────
rule compare:
    input:
        expand("results/{label}/.done", label=[
            "infercnv", 
            "karyotapr", 
            "mosaic_functional", 
            "mosaic_experimental"
        ])
    output:
        "results/comparison/summary.txt"
    conda:
        "envs/python_h5.yml"
    script:
        SCRIPTS["compare"]
