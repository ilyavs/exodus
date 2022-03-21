import pandas as pd
import re
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()


wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in config['samples'].keys()]),
    ref_id = '|'.join([re.escape(x) for x in config['references'].keys()]),

rule all:
    input:
        expand("{sample}/0.reads/{R}_fastqc.html", sample=config['samples'].keys(), R=['R1', 'R2']),
        expand('Summary/Combined_relative_abundances.csv')

rule cp_reads:
    input:
        r1_path = lambda wildcards: S3.remote(config['samples'][wildcards.sample]['r1_path']) if config['samples'][wildcards.sample]['r1_path'].startswith('s3') else config['samples'][wildcards.sample]['r1_path'],
        r2_path = lambda wildcards: S3.remote(config['samples'][wildcards.sample]['r2_path']) if config['samples'][wildcards.sample]['r2_path'].startswith('s3') else config['samples'][wildcards.sample]['r2_path']
    output:
        r1_path = temporary("{sample}/0.reads/R1.fq"),
        r2_path = temporary("{sample}/0.reads/R2.fq")
    shell:
        "if [[ {input.r1_path} == *.zst ]]; then zstd -d {input.r1_path} -o {output.r1_path}; elif [[ {input.r1_path} == *.gz ]]; then zcat {input.r1_path} > {output.r1_path}; else cp {input.r1_path} {output.r1_path}; fi;"
        "if [[ {input.r2_path} == *.zst ]]; then zstd -d {input.r2_path} -o {output.r2_path}; elif [[ {input.r2_path} == *.gz ]]; then zcat {input.r2_path} > {output.r2_path}; else cp {input.r2_path} {output.r2_path}; fi;"

rule cp_refs:
    input: lambda wildcards: S3.remote(config['references'][wildcards.ref_id]) if config['references'][wildcards.ref_id].startswith('s3') else config['references'][wildcards.ref_id]
    output: temporary("references/{ref_id}.fasta")
    shell:
        "cat '{input}' | sed -e 's/[[:space:]]\+/_/g' | sed -e 's/[=|,\/]\+/_/g' | sed 's/^>/>{wildcards.ref_id}_/g' > {output}"

rule fastqc_reads:
    input:
        "{sample}/0.reads/R1.fq",
        "{sample}/0.reads/R2.fq"
    output:
        "{sample}/0.reads/R1_fastqc.html",
        "{sample}/0.reads/R2_fastqc.html"
    threads: 2
    params: out_dir = "{sample}/0.reads"
    benchmark:
        '{sample}/benchmarks/0.fastqc_reads.tsv'
    shell:
        "fastqc -q --extract -t {threads} -o {params.out_dir} {input};"
        "rm {params.out_dir}/*.zip"

rule bwa_index_refs:
    input:
        "references/{ref_id}.fasta"
    output:
        temporary("references/{ref_id}.fasta.amb"),
        temporary("references/{ref_id}.fasta.ann"),
        temporary("references/{ref_id}.fasta.bwt"),
        temporary("references/{ref_id}.fasta.pac"),
        temporary("references/{ref_id}.fasta.sa"),
    shell:
        "bwa index {input} &> /dev/null"

rule samtools_index_reference:
    input:
        "references/{ref_id}.fasta"
    output:
        "references/{ref_id}.fasta.fai"
    shell:
        "samtools faidx {input}"

rule bwa_map_reads_to_ref:
    input:
        r1 = "{sample}/0.reads/R1.fq",
        r2 = "{sample}/0.reads/R2.fq",
        index = ["references/{ref_id}.fasta.amb",
                 "references/{ref_id}.fasta.ann",
                 "references/{ref_id}.fasta.bwt",
                 "references/{ref_id}.fasta.pac",
                 "references/{ref_id}.fasta.sa"]
    output:
        temporary("{sample}/1.align/{sample}_{ref_id}.sam")
    log: "{sample}/logs/1.align_{sample}_{ref_id}.log"
    benchmark: "{sample}/benchmarks/1.align_{sample}_{ref_id}.tsv"
    threads: 10
    params:
        prefix = "references/{ref_id}.fasta"
    shell:
        'bwa mem -t {threads} -o {output[0]} {params.prefix} {input[0]} {input[1]} 2> {log}'

rule sort_sam:
    input:
        '{sample}/1.align/{sample}_{ref_id}.sam'
    output:
        temporary('{sample}/1.align/{sample}_{ref_id}.sorted.bam')
    shell:
        'samtools sort -o {output} -O BAM {input}'

rule filter_bam:
    input:
        '{sample}/1.align/{sample}_{ref_id}.sorted.bam'
    output:
        bam = temporary('{sample}/1.align/{sample}_{ref_id}.sorted.filtered.bam'),
        bai = temporary('{sample}/1.align/{sample}_{ref_id}.sorted.filtered.bam.bai')
    params:
        minq = 0
    shell:
         'samtools view -u -q {params.minq} {input[0]} > {output.bam};'
         'samtools index {output.bam}'

rule merge_bam:
    input:
        expand('{{sample}}/1.align/{{sample}}_{ref_id}.sorted.filtered.bam', ref_id=config['references'].keys())
    output:
        bam = temporary('{sample}/1.align/{sample}.merged.bam'),
        sorted_bam = temporary('{sample}/1.align/{sample}.merged.sorted.bam'),
    threads: 5
    shell:
        'samtools merge -@ $(expr {threads} - 1) {output.bam} {input};'
        'samtools sort -@ $(expr {threads} - 1) -n -o {output.sorted_bam} -O BAM {output.bam}'

rule calculate_relative_abundances:
    input:
        '{sample}/1.align/{sample}.merged.sorted.bam'
    output:
        relative_abundances_path = temporary('{sample}/2.relative_abundances/{sample}_relative_abundances.csv'),
    benchmark: "{sample}/benchmarks/2.relative_abundances_{sample}.tsv"
    script:
        'Scripts/ProcessMergedBam.py'

rule combine_relative_abundances:
    input:
        expand('{sample}/2.relative_abundances/{sample}_relative_abundances.csv', sample=config['samples'].keys())
    output:
        'Summary/Combined_relative_abundances.csv'
    run:
        import pandas as pd
        combined = []
        for i in input:
            combined.append(pd.read_csv(i))
        combined = pd.concat(combined)
        combined.to_csv(output[0], index=False)