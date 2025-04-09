### removes all regions which counts in the control samples does not match a cut off

rule remove_zeros:
    input:
        table = "resources/table_{version}/" + RUN_NAME + ".count.txt",
        config_file = CONFIGFILE
    output:
        nozero = "resources/table_{version}/" + RUN_NAME + ".count_nozero.txt",
        allzero = "resources/table_{version}/count_table_all_zero.txt",
        # outputfolder hardcoded in Rscript
    params:
        run_name = RUN_NAME,
        outputfolder = "resources/",
    conda: "../envs/12_remove_zeros.yaml"
    shadow: "shallow"
    threads: 1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    shell:
        """
        Rscript ./workflow/scripts/12_remove_zeros.R \
        {input.table} \
        {input.config_file} \
        {wildcards.version} \
        {params.run_name} \
        {params.outputfolder}
        """
