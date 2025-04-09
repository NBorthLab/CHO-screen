### rank gRNA based on the read count tables provided
rule mageck_rra:
    input:
        mageck_count_table = "resources/table_{version}/results/" + RUN_NAME + ".count_table_{durchgang}.txt",
    output:
        something = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME + ".R",
        gene_summary = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME + ".gene_summary.txt",
        sgrna_summary = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME + ".sgrna_summary.txt",
        log_summary = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME + ".log",
        report_Rmd = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME + ".report.Rmd",
        summary_summary = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME + "_summary.Rnw",
    params:
        mageck_rra_control = config["mageck_test_mageck_control"],
        mageck_rra_treated = config["mageck_test_mageck_treated"],
        output_folder = "resources/table_{version}/results/{durchgang}_rra/" + RUN_NAME
    conda: "../envs/10_mageck.yaml"
    shadow: "shallow"
    threads: 8
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 50000
    shell:
        "mageck test -k {input.mageck_count_table} "
        "-c {params.mageck_rra_control} "
        "-t {params.mageck_rra_treated} "
        "-n {params.output_folder} --keep-tmp "

