MAGECK_LIBRARY = config["mageck_library_config"]

##combine bowtie table and mageck table
rule combine_bowtie_mageck:
    input:
        gene_mageck = "resources/table_mageck/results/region_rra/" + RUN_NAME + ".gene_summary.txt",
        gene_bowtie = "resources/table_bowtieM/results/region_rra/" + RUN_NAME + ".gene_summary.txt",
        annotation = "resources/raw/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff",
        guide_info = MAGECK_LIBRARY ,
    output:
        table_long = "results/genes_in_windows.txt",
        all_regions_summary = "results/gene_summary_all.txt",
        summary = "results/summary_all_tables.txt",
        neg_results = "results/negative_results_table.txt",
        essential_all = "results/gene_all.txt",
    params:
        path_to_output_folder = "results/",
        run_name = RUN_NAME,
    conda: "../envs/18_combine_tables.yaml"
    threads: 1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    shell:
        "Rscript ./workflow/scripts/18_combine_tables.R " 
        "{input.gene_mageck} {input.gene_bowtie} {input.annotation} {input.guide_info} "
        "{params.run_name} {params.path_to_output_folder}"

