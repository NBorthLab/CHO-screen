## descriptive statistics
rule descriptive_statistics_plots:
    input:
        path_to_input_file_counttable = "resources/table_{version}/" + RUN_NAME + ".count_nozero.txt",
        path_to_input_file_summary = "resources/table_{version}/" + RUN_NAME + ".countsummary.txt",
        path_to_input_file_unfiltered = "resources/table_{version}/" + RUN_NAME + ".count.txt",
        path_to_config_file = CONFIGFILE,
    output:
       output_ = temp("results/table_{version}/" + RUN_NAME + "fakefile.txt"),
       plot_ = "results/table_{version}/" + RUN_NAME + "_correlation_all.png"
    params:
        path_to_output_folder = "results/table_{version}/",
        set_for_mageck = "TRUE"
    conda: "../envs/13_descriptive_statistics.yaml"
    shadow: "shallow"
    threads: 1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    shell:
        "Rscript ./workflow/scripts/13_16_analyse_counttable.R " 
        "{input.path_to_input_file_counttable} "
        "{input.path_to_input_file_summary} "
        "{input.path_to_config_file} "
        "{params.path_to_output_folder} {params.set_for_mageck} "
        "{input.path_to_input_file_unfiltered}"
