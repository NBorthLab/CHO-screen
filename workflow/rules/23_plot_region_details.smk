# plot all interesting regions incl annotation, chrom_state, transcripts
rule plot_region_details:
    input:
        expressed_genes = "resources/raw/expressed-genes.txt",
        not_expressed_genes = "resources/raw/not-expressed-genes.txt",
        chromatin_states =  "resources/raw/53h_dense.bed",
        transcripts = "resources/normalised_transcripts_{selection}.txt",
        essentiality = "results/genes_chrom_states.txt",
        gtf = "resources/raw/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff",
        selection_table = "results/{selection}.txt",
        #table with regions in column called "id"
    output:
        "results_overview/plot_region_details_{selection}.pdf"
    params:
        txdb = "resources/raw/GCF_003668045.3_CriGri-PICRH-1.0_genomic.sqlite",
    conda:
        "../envs/23_plot_region_details.yaml"
    shell:
        """
        Rscript ./workflow/scripts/23_plot_region_details.R \
        {input.expressed_genes} \
        {input.not_expressed_genes} \
        {input.chromatin_states} \
        {input.transcripts} \
        {input.essentiality} \
        {input.gtf} \
        {input.selection_table} \
        {output} \
        {params.txdb}
        """



