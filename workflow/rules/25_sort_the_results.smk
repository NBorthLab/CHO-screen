# restructure directory to results_overview
rule symlinks_for_sorted_results:
    input:
        count_table_mageck = "resources/table_mageck/" + RUN_NAME + ".count_nozero.txt", 
        count_table_bowtie = "resources/table_bowtieM/" + RUN_NAME + ".count_nozero.txt", 
        correlation_plot_of_samples_bowtie = "results/table_bowtieM/" + RUN_NAME + "_correlation_all.png", 
        correlation_plot_of_samples_mageck = "results/table_mageck/" + RUN_NAME + "_correlation_all.png", 
        sgRNA_summary_bowtie = "resources/table_bowtieM/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt", 
        sgRNA_summary_mageck = "resources/table_mageck/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt", 
    output:
        count_table_mageck = "results_overview/table_mageck/" + RUN_NAME + ".count_nozero.txt", 
        count_table_bowtie = "results_overview/table_bowtieM/" + RUN_NAME + ".count_nozero.txt", 
        correlation_plot_of_samples_bowtie = "results_overview/table_bowtieM/" + RUN_NAME + "_correlation_all.png", 
        correlation_plot_of_samples_mageck = "results_overview/table_mageck/" + RUN_NAME + "_correlation_all.png", 
        sgRNA_summary_bowtie = "results_overview/table_bowtieM/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt",
        sgRNA_summary_mageck = "results_overview/table_mageck/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt",
    threads:
        1
    resources:
        mem_mb_per_cpu= 500,
    run:
        import os
        for input_file, output_file in zip(input, output):
                os.symlink(os.path.abspath(input_file), output_file)

