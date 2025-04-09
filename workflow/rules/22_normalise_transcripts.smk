# normalise to 
rule normalise_transcripts:
    input:
        table_transcripts = "resources/transcripts_{selection}.txt",
        ## hardcoded samplenames inside
    output:
        "resources/normalised_transcripts_{selection}.txt"
    params:
        outputfolder = "resources"
    conda:
        "../envs/22_normalise_transcripts.yaml"
    shell:
        """
        Rscript ./workflow/scripts/22_normalise_transcripts.R \
        {input.table_transcripts} \
        {output}
        """
