MAGECK_LIBRARY = config["mageck_library_config"]

### methode to get thecounttable from .bam files
rule mageck_count_bam:
    input:
        library = MAGECK_LIBRARY,
        in1 = "resources/aligned/S11.bam",
        in2 = "resources/aligned/S12.bam",
        in3 = "resources/aligned/S13.bam",
        inrem = "resources/aligned/S14.bam",
        in4 = "resources/aligned/S15.bam",
        in5 = "resources/aligned/S16.bam",
        in6 = "resources/aligned/S21.bam",
        in7 = "resources/aligned/S23.bam",
        inrem2 = "resources/aligned/S24.bam",
        in8 = "resources/aligned/S25.bam",
        in9 = "resources/aligned/S26.bam",
    output:
        count_table_norm = "resources/table_bowtieM/" + RUN_NAME + ".count_normalized.txt",
        count_table = "resources/table_bowtieM/" + RUN_NAME + ".count.txt",
        count_summary = "resources/table_bowtieM/" + RUN_NAME + ".countsummary.txt",
        log = "resources/table_bowtieM/" + RUN_NAME + ".log",
        report_rmd = "resources/table_bowtieM/" + RUN_NAME + ".count_report.Rmd",
        summary_rnw = "resources/table_bowtieM/" + RUN_NAME + "_countsummary.Rnw",
        summary_r = "resources/table_bowtieM/" + RUN_NAME + "_countsummary.R",
    params:
        "-n resources/table_bowtieM/" + RUN_NAME + " --sample-label "
        "treat_1,treat_2,treat_3,treat_rem,treat_4,treat_5,ctrl_1,ctrl_2,ctrl_rem,ctrl_3,ctrl_4"
    conda: "../envs/10_mageck.yaml"
    shadow: "shallow"
    threads:
        8
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 30000
    shell:
        "mageck count -l {input.library} {params} --fastq {input.in1} {input.in2} "
        "{input.in3} {input.inrem} {input.in4} "
        "{input.in5} {input.in6} {input.in7} {input.inrem2} {input.in8} {input.in9}"

