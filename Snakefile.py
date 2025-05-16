configfile: "config.yaml"
SCRIPTS = config["scripts"]

def build_args(tool):
    return " ".join(f"--{k} {v}" for k, v in config["params"].get(tool, {}).items())

rule all:
    input: "results/comparison/summary.txt"

rule infercnv:
    input:
        script = SCRIPTS["infercnv"]
    output:
        "results/infercnv_out/final_compare.csv"
    params:
        args = build_args("infercnv")
    container:
        "docker://crispr_infercnv:latest"
    shell:
        """
        mkdir -p results/infercnv_out
        Rscript {input.script} {params.args} || echo "infercnv failed" > {output}
        """

rule run_tool:
    input:
        script = lambda wc: SCRIPTS[wc.label]
    output:
        "results/{label}_out/final_compare.csv"
    params:
        args = lambda wc: build_args(wc.label)
    conda:
        "envs/r.yml"
    shell:
        """
        mkdir -p results/{wildcards.label}_out
        {{ 'Rscript' if input.script.endswith('.R') else 'python' }} {input.script} {params.args} || echo "{wildcards.label} failed" > {output}
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
