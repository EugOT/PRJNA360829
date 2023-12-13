"""individual"""
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from pandas import read_table
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.5.1")

##### load config and sample sheets #####
configfile: "config.yaml"

samples = read_table(config["samples"]).set_index(["alias"], drop=False)

proj_dir_path=abspath("../")
refdata_dir_path=abspath("/mnt/pub/GENCODE/M27/")
whitelists_dir_path=abspath("/mnt/pub/whitelists/")

##### target rules #####

shell.executable("/bin/bash")
shell.prefix("source /home/etretiakov/.bashrc; export LC_ALL=en_US.utf-8 && export LANG=en_US.utf-8 &&")

rule all:
    input:
        expand(["kb_{sample}/counts_filtered/adata.h5ad"],
                sample=samples["alias"]),
        # expand(["bus_output_{sample}/output.bus"],
        #         sample=samples["alias"]),
        # expand(["bus_output_{sample}/umicorr.bus"],
        #         sample=samples["alias"]),
        # expand(["count_{sample}/output.mtx"],
        #         sample=samples["alias"]),
        # expand(["count_ds_pred_{sample}/output.mtx"],
        #         sample=samples["alias"])

# localrules: all, kallisto_g
ruleorder: kb_gc #> kallisto_gc > bustools_sort_gc > bustools_umicorrect_gc > bustools_count_gc > bustools_downsample_gc > bustools_predict_gc

##### load rules #####

rule kb_gc:
    input:
        r1=join("fastq_gz","{sample}_R1.fastq.gz"),
        r2=join("fastq_gz","{sample}_R2.fastq.gz")
    output:
        "kb_{sample}/counts_filtered/adata.h5ad"
    log:
        "logs/kallisto-bus_{sample}.log"
    params:
        idx=join(refdata_dir_path, "transcriptome.idx"),
        t2g=join(refdata_dir_path, "transcripts_to_genes.txt"),
        tech=lambda wildcards: str(samples["technology"][wildcards.sample]).upper(),
        outdir="kb_{sample}"
    conda: "kallisto-bus.yaml"
    threads: 10
    shell:
        ("kb count \
            -i {params.idx} \
            -g {params.t2g} \
            -o {params.outdir} \
            -x {params.tech} \
            -t {threads} \
            --h5ad \
            --filter bustools \
            {input.r1} {input.r2}")

rule kallisto_gc:
    input:
        r1=join("fastq_gz","{sample}_R1.fastq.gz"),
        r2=join("fastq_gz","{sample}_R2.fastq.gz")
    output:
        "bus_output_{sample}/output.bus",
        "bus_output_{sample}/matrix.ec",
        "bus_output_{sample}/transcripts.txt"
    log:
        "logs/kallisto-bus_{sample}.log"
    params:
        idx=join(refdata_dir_path, "transcriptome.idx"),
        tech=lambda wildcards: samples["technology"][wildcards.sample],
        outdir="bus_output_{sample}"
    conda: "kallisto-bus.yaml"
    threads: 10
    shell:
        ("kallisto bus \
            -i {params.idx} \
            -o {params.outdir} \
            -x {params.tech} \
            -t {threads} \
            {input.r1} {input.r2}")

rule bustools_sort_gc:
    input:
        "bus_output_{sample}/output.bus"
    output:
        "bus_output_{sample}/sort.bus"
    log:
        "logs/bustools-sort_{sample}.log"
    conda: "kallisto-bus.yaml"
    threads: 4
    shell:
        ("bustools sort -T /tmp/ -t {threads} -o {output} {input} ")

rule bustools_collapse_gc:
    input:
       bus="bus_output_{sample}/sort.bus",
       ec="bus_output_{sample}/matrix.ec",
       t="bus_output_{sample}/transcripts.txt"
    output:
        "bus_output_{sample}/coll.bus"
    log:
        "logs/bustools-umicorr_{sample}.log"
    params:
        t2g=join(refdata_dir_path, "transcripts_to_genes.txt"),
        outfile="bus_output_{sample}/coll"
    conda: "kallisto-bus.yaml"
    shell:
        ("bustools collapse -e {input.ec} -g {params.t2g} -t {input.t} -o {params.outfile} {input.bus} ")

rule bustools_umicorrect_gc:
    input:
       bus="bus_output_{sample}/sort.bus",
       ec="bus_output_{sample}/matrix.ec",
       t="bus_output_{sample}/transcripts.txt"
    output:
        "bus_output_{sample}/umicorr.bus"
    log:
        "logs/bustools-umicorr_{sample}.log"
    params:
        t2g=join(refdata_dir_path, "transcripts_to_genes.txt")
    conda: "kallisto-bus.yaml"
    shell:
        ("bustools umicorrect  -e {input.ec} -g {params.t2g} -t {input.t} -o {output} {input.bus} ")

rule bustools_count_gc:
    input:
        bus="bus_output_{sample}/umicorr.bus",
        ec="bus_output_{sample}/matrix.ec",
        t="bus_output_{sample}/transcripts.txt"
    output:
        "count_{sample}/output.mtx"
    log:
        "logs/bustools-count_{sample}.log"
    params:
       t2g=join(refdata_dir_path, "transcripts_to_genes.txt"),
       outdir="count_{sample}/"
    conda: "kallisto-bus.yaml"
    shell:
        ("bustools count --hist -m --genecounts -e {input.ec} -g {params.t2g} -t {input.t} -o {params.outdir} {input.bus} ")

rule bustools_downsample_gc:
    input:
        bus="bus_output_{sample}/umicorr.bus",
        ec="bus_output_{sample}/matrix.ec",
        t="bus_output_{sample}/transcripts.txt"
    output:
        "count_ds_{sample}/output.mtx"
    log:
        "logs/bustools-count-ds_{sample}.log"
    params:
       t2g=join(refdata_dir_path, "transcripts_to_genes.txt"),
       ds=0.1,
       outdir="count_ds_{sample}/"
    conda: "kallisto-bus.yaml"
    shell:
        ("bustools count --downsample {params.ds} --hist -m --genecounts -e {input.ec} -g {params.t2g} -t {input.t} -o {params.outdir} {input.bus} ")

rule bustools_predict_gc:
    input:
        "count_ds_{sample}/output.mtx"
    output:
        "count_ds_pred_{sample}/output.mtx"
    log:
        "logs/bustools-predict_{sample}.log"
    params:
        t=10,
        indir="count_ds_{sample}/",
        outdir="count_ds_pred_{sample}/"
    conda: "kallisto-bus.yaml"
    shell:
        ("bustools predict -t {params.t} -o {params.outdir} {params.indir} 2>{log}")
