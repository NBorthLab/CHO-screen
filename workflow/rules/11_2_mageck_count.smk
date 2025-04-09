MAGECK_LIBRARY = config["mageck_library_config"]

### methode to get counttable from fastq
rule mageck_count:
    input:
        library = MAGECK_LIBRARY,
        in1= "resources/raw/flowcell_1/215716_S1_L004_R1_001.fastq.gz", bn1 = "resources/raw/flowcell_2/215716_S1_L001_R1_001.fastq.gz",
        in2= "resources/raw/flowcell_1/215717_S2_L004_R1_001.fastq.gz", bn2 = "resources/raw/flowcell_2/215717_S2_L001_R1_001.fastq.gz",
        in3= "resources/raw/flowcell_1/215718_S3_L004_R1_001.fastq.gz", bn3 = "resources/raw/flowcell_2/215718_S3_L001_R1_001.fastq.gz",
        in4= "resources/raw/flowcell_1/215719_S4_L004_R1_001.fastq.gz", bn4 = "resources/raw/flowcell_2/215719_S4_L001_R1_001.fastq.gz",
        in5= "resources/raw/flowcell_1/215720_S5_L004_R1_001.fastq.gz", bn5 = "resources/raw/flowcell_2/215720_S5_L001_R1_001.fastq.gz",
        in6= "resources/raw/flowcell_1/215721_S6_L004_R1_001.fastq.gz", bn6 = "resources/raw/flowcell_2/215721_S6_L001_R1_001.fastq.gz",
        in7= "resources/raw/flowcell_1/215722_S7_L004_R1_001.fastq.gz", bn7 = "resources/raw/flowcell_2/215722_S7_L001_R1_001.fastq.gz",
        in8= "resources/raw/flowcell_1/215724_S8_L004_R1_001.fastq.gz", bn8 = "resources/raw/flowcell_2/215724_S8_L001_R1_001.fastq.gz",
        in9= "resources/raw/flowcell_1/215725_S9_L004_R1_001.fastq.gz", bn9 = "resources/raw/flowcell_2/215725_S9_L001_R1_001.fastq.gz",
        in10= "resources/raw/flowcell_1/215726_S10_L004_R1_001.fastq.gz", bn10 = "resources/raw/flowcell_2/215726_S10_L001_R1_001.fastq.gz",
        in11= "resources/raw/flowcell_1/215727_S11_L004_R1_001.fastq.gz", bn11 = "resources/raw/flowcell_2/215727_S11_L001_R1_001.fastq.gz",

        rin1= "resources/raw/flowcell_1/215716_S1_L004_R2_001.fastq.gz", rbn1 = "resources/raw/flowcell_2/215716_S1_L001_R2_001.fastq.gz",
        rin2= "resources/raw/flowcell_1/215717_S2_L004_R2_001.fastq.gz", rbn2 = "resources/raw/flowcell_2/215717_S2_L001_R2_001.fastq.gz",
        rin3= "resources/raw/flowcell_1/215718_S3_L004_R2_001.fastq.gz", rbn3 = "resources/raw/flowcell_2/215718_S3_L001_R2_001.fastq.gz",
        rin4= "resources/raw/flowcell_1/215719_S4_L004_R2_001.fastq.gz", rbn4 = "resources/raw/flowcell_2/215719_S4_L001_R2_001.fastq.gz",
        rin5= "resources/raw/flowcell_1/215720_S5_L004_R2_001.fastq.gz", rbn5 = "resources/raw/flowcell_2/215720_S5_L001_R2_001.fastq.gz",
        rin6= "resources/raw/flowcell_1/215721_S6_L004_R2_001.fastq.gz", rbn6 = "resources/raw/flowcell_2/215721_S6_L001_R2_001.fastq.gz",
        rin7= "resources/raw/flowcell_1/215722_S7_L004_R2_001.fastq.gz", rbn7 = "resources/raw/flowcell_2/215722_S7_L001_R2_001.fastq.gz",
        rin8= "resources/raw/flowcell_1/215724_S8_L004_R2_001.fastq.gz", rbn8 = "resources/raw/flowcell_2/215724_S8_L001_R2_001.fastq.gz",
        rin9= "resources/raw/flowcell_1/215725_S9_L004_R2_001.fastq.gz", rbn9 = "resources/raw/flowcell_2/215725_S9_L001_R2_001.fastq.gz",
        rin10= "resources/raw/flowcell_1/215726_S10_L004_R2_001.fastq.gz", rbn10 = "resources/raw/flowcell_2/215726_S10_L001_R2_001.fastq.gz",
        rin11= "resources/raw/flowcell_1/215727_S11_L004_R2_001.fastq.gz", rbn11 = "resources/raw/flowcell_2/215727_S11_L001_R2_001.fastq.gz",
    output:
        count_table_norm = "resources/table_mageck/" + RUN_NAME + ".count_normalized.txt",
        count_table = "resources/table_mageck/" + RUN_NAME + ".count.txt",
        count_summary = "resources/table_mageck/" + RUN_NAME + ".countsummary.txt",
        log = "resources/table_mageck/" + RUN_NAME + ".log",
        report_rmd = "resources/table_mageck/" + RUN_NAME + ".count_report.Rmd",
        summary_rnw = "resources/table_mageck/" + RUN_NAME + "_countsummary.Rnw",
        summary_r = "resources/table_mageck/" + RUN_NAME + "_countsummary.R",
    params:
        "-n resources/table_mageck/" + RUN_NAME + " --sample-label "
        "treat_1,treat_2,treat_3,treat_rem,treat_4,treat_5,ctrl_1,ctrl_2,ctrl_rem,ctrl_3,ctrl_4"
    threads:
        8
    resources:
        mem_mb_per_cpu= 30000
    conda: "../envs/10_mageck.yaml"
    shadow: "shallow"
    shell:
        "mageck count -l {input.library} {params} --fastq {input.in1},{input.bn1} {input.in2},{input.bn2} {input.in3},{input.bn3} {input.in4},{input.bn4} " 
        "{input.in5},{input.bn5} {input.in6},{input.bn6} {input.in7},{input.bn7} {input.in8},{input.bn8} {input.in9},{input.bn9} "
        "{input.in10},{input.bn10} {input.in11},{input.bn11} "
        "--fastq-2 {input.rin1},{input.rbn1} {input.rin2},{input.rbn2} {input.rin3},{input.rbn3} {input.rin4},{input.rbn4} "
        "{input.rin5},{input.rbn5} {input.rin6},{input.rbn6} {input.rin7},{input.rbn7} {input.rin8},{input.rbn8} {input.rin9},{input.rbn9} " 
        "{input.rin10},{input.rbn10} {input.rin11},{input.rbn11}"
