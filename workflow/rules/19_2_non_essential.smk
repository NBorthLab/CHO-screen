## categorise non-essential
rule cat_non_essential:
    input:
        input_file_counttable = "resources/table_" + NON_ESS_BASIS + "/results/" + RUN_NAME + ".count_table_region.txt",
        input_file_sgrna_summary = "resources/table_" + NON_ESS_BASIS + "/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt",
        path_to_config_file = CONFIGFILE,
        input_file_gene_target = "results/gene_all.txt",
        input_file_gene_all = "results/gene_summary_all.txt",
    output:
        ("results/selected_least_essential.txt")
    params:
        input_folder = "resources/table_" + RUN_NAME + "/results/region_rra/table_bowtieM/results/",
        path_to_output_folder = "results/",
        input_folder_essential = "resources/",
        outputfolder = "results/",
        project_run_name = "find_least_essential",
    conda:
        "../envs/19_cat_non_essential.yaml"
    threads: 1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    shell:
        """
        Rscript ./workflow/scripts/19_cat_non_essential.R \
        {params.input_folder} \
        {input.input_file_counttable} \
        {input.input_file_sgrna_summary} \
        {params.path_to_output_folder} \
        {input.path_to_config_file} \
        {params.input_folder_essential} \
        {input.input_file_gene_target} \
        {input.input_file_gene_all} \
        {params.outputfolder} \
        {params.project_run_name}
        """
