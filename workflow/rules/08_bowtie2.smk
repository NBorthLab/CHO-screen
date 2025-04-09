### creating index for alignment
REFERENCE = config["add_gene_name_path_to_gene_table"]

rule create_bowtie_index:
    input:
        config["add_gene_name_path_to_gene_table"]
    output:
        #output0= REFERENCE,
        output1= REFERENCE + ".1.bt2",
        output2= REFERENCE + ".2.bt2",
        output3= REFERENCE + ".3.bt2",
        output4= REFERENCE + ".4.bt2",
        outputrev1= REFERENCE + ".rev.1.bt2",
        outputrev2= REFERENCE + ".rev.2.bt2"
    params:
        basename= config["add_gene_name_path_to_gene_table"]
    conda: "../envs/08_bowtie.yaml"
    shadow: "shallow"
    shell:
        "bowtie2-build {input} {params.basename}"



### align sequenze to reference
rule align_reads:
    input:
        sequence_file="resources/trim_all/S{sample}_amplicon_trim.fastq",
        reference_file= config["add_gene_name_path_to_gene_table"],
        for_rule_calling= REFERENCE + ".1.bt2"
    output:
        log=("resources/aligned/S{sample}_amplicon_trim.fastq_2"),
        sam=temp("resources/aligned/S{sample}.sam"),
    params:
        "--norc -D30 -R6 --score-min L,0,0" 
        #--norc: no reverse compliment
        # -D: attempts to find a new best or second best alignment, bevor moving on
        # -R: number of re-seeds
    threads:
        16
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    conda: "../envs/08_bowtie.yaml"
    shadow: "shallow"
    shell:
        "bowtie2 -x {input.reference_file} {params} -p {threads} "
        "-U {input.sequence_file} -S {output.sam} 2> {output.log}"
        #-U: list of files containing unpaired reads to be aligned
        #-p: Number of threads to use

 
