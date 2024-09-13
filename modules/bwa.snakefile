rule bwa:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        sam="analysis/bwa/{sample}/{sample}.sam"
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/bwa.{sample}.txt"
    threads: 3
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
    threads: 3
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
    threads: 1
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
    shell:
        "samtools index {input.bam}"

