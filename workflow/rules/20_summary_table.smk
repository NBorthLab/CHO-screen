# summary table - combine tables and put information together
rule summary_tables:
    input:
        summary_all_tables = "results/summary_all_tables.txt",
##        essential_table = "results/selected_least_essential.txt",
        target_table = "results/gene_all.txt",
        chrom_states = "resources/raw/53h_dense.bed",
    output:
        summary = "results/summary_table.txt",
        chrom = "results/genes_chrom_states.txt"
    params:
        outputfolder = "results/",
        folder_chrom_states = "resources/",
#    benchmark:
#        "../benchmarks/summary_table.benchmark"
    conda:
        "../envs/20_summary_table.yaml"
    threads: 1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    shell:
        """
        Rscript ./workflow/scripts/20_summary_table.R \
        {input.summary_all_tables} \
        {input.target_table} \
        {input.chrom_states} \
        {params.folder_chrom_states} \
        {params.outputfolder}
        """


        ##   {input.essential_table} \
    
