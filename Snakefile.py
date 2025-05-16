configfile: "config.yaml"
SCRIPTS = config["scripts"]

def build_args(tool):
    return " ".join(f"--{k} {v}" for k, v in config["params"].get(tool, {}).items())

rule all:
    input: "results/comparison/summary.txt"

# ──────────── RULE INFERCNV : via DOCKER ─────────────
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

# ───────────── RULE MOSAIC : envs/python_mosaic.yml ─────────────
rule mosaic:
    input:
        script = SCRIPTS["mosaic"]
    output:
        "results/mosaic/.done"
    params:
        args = build_args("mosaic")
    conda:
        "envs/python_mosaic.yml"
    shell:
        """
        mkdir -p results/mosaic
        python {input.script} {params.args}
        touch {output}
        """

# ───────────── RULE KARYOTAPR : envs/r.yml ─────────────
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

rule compare:
    input:
        expand("results/{label}_out/final_compare.csv", label=["infercnv", "mosaic", "karyotapr"])
    output:
        "results/comparison/summary.txt"
    conda:
        "envs/r.yml"
    shell:
        "python {SCRIPTS['compare']}"
