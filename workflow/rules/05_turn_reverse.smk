##turning the merged, who are facing the reversed direction around, 
##so all of them are in the same direction. 
REV_AMP_CMD = r"""{print "@"($name);print revcomp($seq);print "+";print reverse($qual)}"""

rule bioawk:
    input:
        reverse_amplicon="resources/fastq_greped/S{sample}_reverse_amplicon.fastq"
    output:
        temp("resources/fastq_greped/S{sample}_reverse_amplicon_revcomp.fastq")
    params:
        "-c fastx",
    threads:
        1
    resources:
        mem_mb_per_cpu= 1000
    conda: "../envs/05_turn_reverse.yaml"
    shadow: "shallow"
    shell:
        "bioawk {params} {REV_AMP_CMD:q} {input.reverse_amplicon} > {output}"

