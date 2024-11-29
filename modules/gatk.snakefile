rule haplotypecaller:
    input:
        bam="analysis/gatk/{sample}/{sample}.s.md.ir.bam"
    output:
        vcf="analysis/gatk/{sample}/{sample}.vcf"
    wildcard_constraints:
        sample="[^/]+",
    benchmark: "benchmarks/haplotypecaller.{sample}.txt"
    threads: 16
    resources:          
        tasks=2,        
        mem_mb=1,       
        mpi="srun",     
        runtime=2880,   
        tasks_per_node=1
    container: "docker://docker.io/broadinstitute/gatk:4.0.5.1"
    params:
        tmp="analysis/gatk/{sample}/temp",
        java_opts="-Xmx4g"
    shell:
        "gatk"
        " --java-options '{params.java_opts}'"
        " HaplotypeCaller"
        " -R {config[ref_genome]}"
        " -I {input.bam}"
        " -O {output.vcf}"
        " --native-pair-hmm-threads {threads}"

rule realigner_tc:
    input:
        bam="analysis/bwa/{sample}/{sample}.s.md.bam",
        bai="analysis/bwa/{sample}/{sample}.s.md.bam.bai"
    output:
        intervals="analysis/gatk/{sample}/{sample}.s.md.bam.intervals"
    wildcard_constraints:
        sample="[^/]+",
    benchmark: "benchmarks/realigner_tc.{sample}.txt"
    threads: 3
    resources:          
        tasks=2,        
        mem_mb=1,       
        mpi="srun",     
        runtime=2880,   
        tasks_per_node=1
    params:
        tmp="analysis/haplotypecaller/{sample}/temp",
        java_opts="-Xms8g -Xmx30g"
    shell:
        "java"
        " {params.java_opts}"
        " -jar {config[ref_gatk35]}"
        " -T RealignerTargetCreator"
        " -R {config[ref_genome]}"
        " -I {input.bam}"
        " -o {output.intervals}"
        " -nt {threads}"

rule indelrealigner:
    input:
        bam="analysis/bwa/{sample}/{sample}.s.md.bam",
        intervals="analysis/gatk/{sample}/{sample}.s.md.bam.intervals"
    output:
        bam="analysis/gatk/{sample}/{sample}.s.md.ir.bam"
    wildcard_constraints:                            
        sample="[^/]+",                              
    benchmark: "benchmarks/realigner_tc.{sample}.txt"
    resources:          
        tasks=2,        
        mem_mb=1,       
        mpi="srun",     
        runtime=2880,   
        tasks_per_node=1
    params:
        tmp="analysis/haplotypecaller/{sample}/temp",
        java_opts="-Xms8g -Xmx30g"
    shell:
        "java"
        " {params.java_opts}"
        " -jar {config[ref_gatk35]}"
        " -T IndelRealigner"
        " -R {config[ref_genome]}"
        " -I {input.bam}"
        " -targetIntervals {input.intervals}"
        " -o {output.bam}"

rule index_irbam:
    input:
        bam="analysis/gatk/{sample}/{sample}.s.md.ir.bam"
    output:
        bai="analysis/gatk/{sample}/{sample}.s.md.ir.bam.bai"
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/index_irbam.{sample}.txt"
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

rule zip_vcf:
    input:
        vcf="analysis/gatk/{sample}/{sample}.vcf"
    output:
        vcf="analysis/gatk/{sample}/{sample}.vcf.gz"
    resources:          
        tasks=2,        
        mem_mb=1,       
        mpi="srun",     
        runtime=2880,   
        tasks_per_node=1
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/index_irbam.{sample}.txt"
    shell:
        "bgzip -c {input.vcf} > {output.vcf}"

rule index_vcf:
    input:
        vcf="analysis/gatk/{sample}/{sample}.vcf.gz"
    output:
        tbi="analysis/gatk/{sample}/{sample}.vcf.gz.tbi"
    wildcard_constraints:
        sample="[^/]+",
    #benchmark: "benchmarks/index_irbam.{sample}.txt"
    resources:          
        tasks=2,        
        mem_mb=1,       
        mpi="srun",     
        runtime=2880,   
        tasks_per_node=1
    shell:
        "tabix -p vcf {input.vcf}"
