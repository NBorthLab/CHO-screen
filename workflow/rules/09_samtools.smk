### sam to bam file with header
rule sam_to_bam:
    input:
        "resources/aligned/S{sample}.sam"
    output:
        "resources/aligned/S{sample}.bam"
    params:
        "-bh"
        # -b: output in bam-format
        # -h --with-header: Include the header in the output
        # -o --output: Output to FILE [stdout]
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 1000
    conda: "../envs/09_samtools.yaml"
    shadow: "shallow"
    shell:
        "samtools view {params} {input} -o {output}"
