### beginning with tct are in forward direction, separating merged between forward and reverse direction
rule fastq_grep:
    input:
        "resources/merged_overlapping_reads/S{sample}.assembled.fastq"
    output:
        forward_amplicon=temp("resources/fastq_greped/S{sample}_forward_amplicon.fastq"),
        reverse_amplicon=temp("resources/fastq_greped/S{sample}_reverse_amplicon.fastq")
    threads:
        1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu=1000
    conda: "../envs/04_sort_directions.yaml"
    shadow: "shallow"
    shell:
        "fastq-grep ^TCT {input} > {output.forward_amplicon}; "
        "fastq-grep -v ^TCT {input} > {output.reverse_amplicon}"
        #^: beginn of a line
        #-v: --invert match; all lines NOT matching
