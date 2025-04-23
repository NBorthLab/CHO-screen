MAGECK_LIBRARY = config["mageck_library_config"]

# create an overview table
rule summary_table_overview:
    input:
        results_table = "results/summary_table.txt",
        genes_chrom_states = "results/genes_chrom_states.txt",
        sgrna_summary_bowtie = "resources/table_bowtieM/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt",
        sgrna_summary_mageck = "resources/table_mageck/results/region_rra/" + RUN_NAME + ".sgrna_summary.txt",
        guides_design = MAGECK_LIBRARY
    output:
        compare_bowtie_mageck = "results_overview/compare_bowtie_mageck.png",
        summary_table_overview = "results_overview/summary_table_overview.txt",
        results_table_conclusion = "results_overview/results_all_regions.txt",
        sgrna_summary_conclusion = "results_overview/sgRNA_summary.txt",
    threads:
        1
    resources:
        mem_mb_per_cpu= 5000,
    conda:
        "../envs/26_summary_table.yaml"
    shell:
        """
        Rscript ./workflow/scripts/26_create_summary_table.R \
        {input.results_table} \
        {input.genes_chrom_states} \
        {input.sgrna_summary_bowtie} \
        {input.sgrna_summary_mageck} \
        {input.guides_design} \
        {output.compare_bowtie_mageck} \
        {output.summary_table_overview} \
        {output.results_table_conclusion} \
        {output.sgrna_summary_conclusion}
        """
