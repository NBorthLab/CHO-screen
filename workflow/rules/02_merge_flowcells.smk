### merge samples from different flowcells
rule merge_flowcells:
    input:
        flowcell1="resources/unzipped/flowcell_1_S{sample}_val_{strand}.fastq",
        flowcell2="resources/unzipped/flowcell_2_S{sample}_val_{strand}.fastq"
    output:
        ("resources/flowcell_merged/S{sample}_val_{strand}_merged.fq")
    threads:
        1
    resources:
        mem_mb_per_cpu= 10000
    shell:
        "cat {input.flowcell1} {input.flowcell2} > {output}"
