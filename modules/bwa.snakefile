rule bwa:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        sam="analysis/bwa/{sample}/{sample}.sam"
    wildcard_constraints:
        sample="[^/]+",
    threads: 8
    message: "BWA using bwa mem on {wildcards.sample}"
    container: "docker://docker.io/biocontainers/bwa:v0.7.17_cv1"
    log:
        stdout="logs/bwa/{sample}.stdout",
        stderr="logs/bwa/{sample}.stderr"
    resources:
        tasks=2,
        mem_mb=1,
        mpi="srun",
        runtime=2880,
        tasks_per_node=1
    params:
        tmp="analysis/bwa/{sample}/temp",
        rg="\'@RG\\tID:{sample}\\tLB:1\\tSM:{sample}\\tPL:ILLUMINA'",
    shell:
        "mkdir -p {params.tmp} && "
        "bwa mem -M"
        " -R {params.rg}"
        " -t {threads}"
        " {config[ref_bwaidx]}"
        " {input.fq1} {input.fq2} > {output.sam}"

rule sortbam:
    input:
        sam="analysis/bwa/{sample}/{sample}.sam"
    output:
        bam="analysis/bwa/{sample}/{sample}.s.bam"
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/sortbam.{sample}.txt"
    threads: 4
    resources:
        tasks=2,
        mem_mb=1,
        mpi="srun",
        runtime=2880,   
        tasks_per_node=1   
    container: "docker://docker.io/biocontainers/samtools:v1.9-4-deb_cv1"
    params:
        tmp="analysis/bwa/{sample}/temp",
        ref_genome=config['ref_genome']
    shell:
        "mkdir -p {params.tmp} && "
        "samtools view -Shu"
        " -@ {threads}"
        " {input.sam} |"
        " samtools sort"
        " -@ {threads} - > {output.bam}"

rule index_sbam:
    input:
        bam="analysis/bwa/{sample}/{sample}.s.bam"
    output:
        bai="analysis/bwa/{sample}/{sample}.s.bam.bai"
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/index_sbam.{sample}.txt"
    resources:
        tasks=2,
        mem_mb=1,
        mpi="srun",
        runtime=2880,   
        tasks_per_node=1
    threads: 1
    container: "docker://docker.io/biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "samtools index {input.bam}"

rule md_sbam:
    input:
        bam="analysis/bwa/{sample}/{sample}.s.bam"
    output:
        bam="analysis/bwa/{sample}/{sample}.s.md.bam",
        txt="analysis/bwa/{sample}/{sample}.metrics.txt"
    wildcard_constraints:
        sample="[^/]+",
    benchmark: "benchmarks/md_sbam.{sample}.txt"
    threads: 1
    container: "docker://docker.io/biocontainers/picard:v2.3.0_cv3"
    resources:
        tasks=2,
        mem_mb=1,
        mpi="srun",
        runtime=2880,   
        tasks_per_node=1
    params:
        tmp="analysis/bwa/{sample}/temp",
        java_opts="-Xms8g -Xmx30g"
    shell:
        "mkdir -p {params.tmp} && "
        "picard"
        " '{params.java_opts}'"
        " -XX:ConcGCThreads={threads}"
        " MarkDuplicates"
        " TMP_DIR={params.tmp}"
        " I={input.bam}"
        " O={output.bam}"
        " M={output.txt}"
        " VALIDATION_STRINGENCY=LENIENT"
        " REMOVE_DUPLICATES=false"
        " REMOVE_SEQUENCING_DUPLICATES=false"

rule index_mdbam:
    input:
        bam="analysis/bwa/{sample}/{sample}.s.md.bam"
    output:
        bai="analysis/bwa/{sample}/{sample}.s.md.bam.bai"
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/index_mdbam.{sample}.txt"
    threads: 1
    resources:
        tasks=2,
        mem_mb=1,
        mpi="srun",
        runtime=2880,   
        tasks_per_node=1
    container: "docker://docker.io/biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "samtools index {input.bam}"

