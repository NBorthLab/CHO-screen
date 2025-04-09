# create the sample table for mageck:
rule prepare_for_mageck:
    input:
        counttable = "resources/table_{version}/" + RUN_NAME + ".count_nozero.txt", 
        configfile = CONFIGFILE
    output:
        mageck_count_table_region = "resources/table_{version}/results/" + RUN_NAME + ".count_table_region.txt"
    params:
        path_to_output_folder = "resources/table_{version}/results/"
    conda: "../envs/14_prepare_mageck_test.yaml"
    shadow: "shallow"
    threads: 1
    resources:
        cores = lambda wc, threads: threads,
        mem_mb_per_cpu= 10000
    shell:
        "Rscript ./workflow/scripts/14_prepare_files_for_mageck.R {input.counttable} {input.configfile} "
        "{params.path_to_output_folder}"
