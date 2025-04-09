### combine forward and reverse amplicon to one .fastq file 
rule combine_amplicon:
    input:
        forward_amplicon="resources/fastq_greped/S{sample}_forward_amplicon.fastq",
        reverse_amplicon="resources/fastq_greped/S{sample}_reverse_amplicon_revcomp.fastq"
    output:
        temp("resources/fastq_greped/S{sample}_amplicon.fastq")
    threads:
        1
    resources:
        mem_mb_per_cpu= 1000
    shadow: "shallow"
    shell:
        "cat {input.forward_amplicon} {input.reverse_amplicon} > {output}"

