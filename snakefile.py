################################################################################
## Snakefile – exécute InferCNV, Mosaic, KaryotapR, puis compare les résultats ##
################################################################################
configfile: "config.yaml"
SCRIPTS = config["scripts"]

def build_args(tool):
    return " ".join(f"--{k} {v}" for k, v in config["params"].get(tool, {}).items())

rule all:
    input: "results/comparison/summary.txt"

rule run_tool:
    input:
        script=lambda wc: SCRIPTS[wc.label]
    output:
        touch("results/{label}/.done")
    params:
        args=lambda wc: build_args(wc.label)
    conda:
        lambda wc: "envs/r.yml"       # un seul environnement R commun (au début)
    shell:
        """
        mkdir -p results/{wildcards.label}
        {{ 'Rscript' if input.script.endswith('.R') else 'python' }} {input.script} {params.args}
        """

# Crée les cibles .done pour chaque outil :
expand("results/{label}/.done", label=["infercnv", "mosaic", "karyotapr"])

rule compare:
    input:
        expand("results/{label}/.done", label=["infercnv", "mosaic", "karyotapr"])
    output:
        "results/comparison/summary.txt"
    conda:
        "envs/r.yml"
    shell:
        "python {SCRIPTS['compare']} > {output}"
