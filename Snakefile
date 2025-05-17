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
rule h5:
    input:
        script = config["scripts"]["h5"]
    output:
        "results/h5/.done"
    params:
        args = build_args("h5")
    conda:
        "envs/h5.yml"
    shell:
        """
        mkdir -p results/h5
        python {input.script} {params.args}
        touch {output}
        """


# ────────────────────────────
# MOSAIC EXPERIMENTAL (via python_mosaic)
# ────────────────────────────
rule mosaic:
    input:
        script = SCRIPTS["mosaic"]
    output:
        "results/mosaic/.done"
    params:
        args = build_args("mosaic")
    conda:
        "envs/mosaic.yml"
    shell:
        """
        mkdir -p results/mosaic
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
        infercnv = "results/infercnv/.done",
        karyotapr = "results/karyotapr/.done",
        h5 = "results/h5/.done",
        mosaic = temp("results/mosaic/.done")  # marqué comme temporaire (facultatif)
    output:
        "results/comparison/summary.txt"
    conda:
        "envs/h5.yml"
    script:
        SCRIPTS["compare"]

