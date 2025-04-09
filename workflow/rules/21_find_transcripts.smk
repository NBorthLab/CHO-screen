## compare transcripts to position within some datasets
rule find_transcripts:
    input:
        path_to_bam = "resources/raw/",
        region_table = "results/summary_table.txt",
        selection_table = "results/{selection}.txt",
        # path "resources/raw/" hardcoded in Rscript... :(
        # also hardcode : "/aligned_sorted_INH04" ,x ,".bam")
    output:
        temp("resources/transcripts_{selection}.txt")
    params:
        folder_for_selection_file = "results/",
        outputfolder = "resources/"
    conda:
        "../envs/21_find_transcripts.yaml"
    shell:
        """
        Rscript ./workflow/scripts/21_find_transcripts.R \
        {input.region_table} \
        {input.selection_table} \
        {input.path_to_bam} \
        {params.outputfolder} \
        {params.folder_for_selection_file}
        """
