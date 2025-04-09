### merge overlapping reads: merge paired end reads into an overlapping fragment, and also increase quality score
rule merge_overlapping_reads:
    input:
        forward_read="resources/flowcell_merged/S{sample}_val_1_merged.fq",
        reverse_read="resources/flowcell_merged/S{sample}_val_2_merged.fq"
    output:
        temp("resources/merged_overlapping_reads/S{sample}.assembled.fastq"),
        temp("resources/merged_overlapping_reads/S{sample}.discarded.fastq"),
        temp("resources/merged_overlapping_reads/S{sample}.unassembled.forward.fastq"),
        temp("resources/merged_overlapping_reads/S{sample}.unassembled.reverse.fastq")
    threads: 4
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    conda: "../envs/03_merge_overlapping.yaml"
    shadow: "shallow"
    shell:
        "pear -f {input.forward_read} -r {input.reverse_read} -y {resources.mem_mb_per_cpu}  -o ./resources/merged_overlapping_reads/S{wildcards.sample} -j {threads}"
        #-f, --forward-fastq: Forward paired-end FASTQ file.
        #-r, --reverse-fastq: Reverse paired-end FASTQ file.
        #-o, --output: Output filename.
        #-j Number of threads to use
