## trim sequenze of everything except gRNA
## trimming down to guideRNA
rule trim_all:
    input:
        "resources/fastq_greped/S{sample}_amplicon.fastq"
    output:
        ("resources/trim_all/S{sample}_amplicon_trim.fastq")
    threads:
        8
    resources:
        mem_mb_per_cpu= 1000
    params:
        "-m 18 -M 22 "
        "-g GGAAAGGACGAAACACCG...GTTTCAGAGCTATGCTGGA "
        # -m: minimum length
        # -M: maximum length
    conda: "../envs/07_trim.yaml"
    shadow: "shallow"
    shell:
        "cutadapt {params} -j {threads} -o {output} {input}"
        # -o: name of outputfile
        # -j: number of CPU cores to use
