rule rename_files_0:
    input:          
        s1= "resources/raw/flowcell_1/215716_S1_L004_R{strand}_001.fastq.gz",
        s2= "resources/raw/flowcell_1/215717_S2_L004_R{strand}_001.fastq.gz",
        s3= "resources/raw/flowcell_1/215718_S3_L004_R{strand}_001.fastq.gz",
        s4= "resources/raw/flowcell_1/215719_S4_L004_R{strand}_001.fastq.gz",
        s5= "resources/raw/flowcell_1/215720_S5_L004_R{strand}_001.fastq.gz",
        s6= "resources/raw/flowcell_1/215721_S6_L004_R{strand}_001.fastq.gz",
        s7= "resources/raw/flowcell_1/215722_S7_L004_R{strand}_001.fastq.gz",
        s8= "resources/raw/flowcell_1/215724_S8_L004_R{strand}_001.fastq.gz",
        s9= "resources/raw/flowcell_1/215725_S9_L004_R{strand}_001.fastq.gz", 
        s10= "resources/raw/flowcell_1/215726_S10_L004_R{strand}_001.fastq.gz",
        s11= "resources/raw/flowcell_1/215727_S11_L004_R{strand}_001.fastq.gz",
      # #s12= "resources/raw/flowcell_1/216037_S12_L004_R{strand}_001.fastq.gz",

        s21= "resources/raw/flowcell_2/215716_S1_L001_R{strand}_001.fastq.gz",
        s22= "resources/raw/flowcell_2/215717_S2_L001_R{strand}_001.fastq.gz",
        s23= "resources/raw/flowcell_2/215718_S3_L001_R{strand}_001.fastq.gz",
        s24= "resources/raw/flowcell_2/215719_S4_L001_R{strand}_001.fastq.gz",
        s25= "resources/raw/flowcell_2/215720_S5_L001_R{strand}_001.fastq.gz",
        s26= "resources/raw/flowcell_2/215721_S6_L001_R{strand}_001.fastq.gz",
        s27= "resources/raw/flowcell_2/215722_S7_L001_R{strand}_001.fastq.gz",
        s28= "resources/raw/flowcell_2/215724_S8_L001_R{strand}_001.fastq.gz",
        s29= "resources/raw/flowcell_2/215725_S9_L001_R{strand}_001.fastq.gz",
        s210= "resources/raw/flowcell_2/215726_S10_L001_R{strand}_001.fastq.gz",
        s211= "resources/raw/flowcell_2/215727_S11_L001_R{strand}_001.fastq.gz",
      #  #s212= "resources/raw/flowcell_2/216037_S12_L001_R{strand}_001.fastq.gz",
    output:
        out1= temp("resources/unzipped/flowcell_1_S11_val_{strand}.fastq"),
        out2= temp("resources/unzipped/flowcell_1_S12_val_{strand}.fastq"),
        out3= temp("resources/unzipped/flowcell_1_S13_val_{strand}.fastq"),
        out4= temp("resources/unzipped/flowcell_1_S14_val_{strand}.fastq"),  
        out5= temp("resources/unzipped/flowcell_1_S15_val_{strand}.fastq"),
        out6= temp("resources/unzipped/flowcell_1_S16_val_{strand}.fastq"),
        out7= temp("resources/unzipped/flowcell_1_S21_val_{strand}.fastq"),
        out8= temp("resources/unzipped/flowcell_1_S23_val_{strand}.fastq"),
        out9= temp("resources/unzipped/flowcell_1_S24_val_{strand}.fastq"),  
        out10= temp("resources/unzipped/flowcell_1_S25_val_{strand}.fastq"),
        out11= temp("resources/unzipped/flowcell_1_S26_val_{strand}.fastq"),
       # #out12= temp("resources/raw/flowcell_1/215711_S11_L001_R{strand}_001.fastq.gz"),
        
        out21= temp("resources/unzipped/flowcell_2_S11_val_{strand}.fastq"),
        out22= temp("resources/unzipped/flowcell_2_S12_val_{strand}.fastq"),
        out23= temp("resources/unzipped/flowcell_2_S13_val_{strand}.fastq"),
        out24= temp("resources/unzipped/flowcell_2_S14_val_{strand}.fastq"),  
        out25= temp("resources/unzipped/flowcell_2_S15_val_{strand}.fastq"),
        out26= temp("resources/unzipped/flowcell_2_S16_val_{strand}.fastq"),
        out27= temp("resources/unzipped/flowcell_2_S21_val_{strand}.fastq"),
        out28= temp("resources/unzipped/flowcell_2_S23_val_{strand}.fastq"),
        out29= temp("resources/unzipped/flowcell_2_S24_val_{strand}.fastq"),  
        out210= temp("resources/unzipped/flowcell_2_S25_val_{strand}.fastq"),
        out211= temp("resources/unzipped/flowcell_2_S26_val_{strand}.fastq"),
       # #out212= temp("resources/raw/flowcell_2/215711_S11_L001_R{strand}_001.fastq.gz"),
    threads:
        32
    resources:
        mem_mb_per_cpu= 500
    conda:
        "../envs/03_merge_overlapping.yaml"
    shell:
        "pigz -cd -p {threads} {input.s1} > {output.out1}; "
        "pigz -cd -p {threads} {input.s2} > {output.out2}; "
        "pigz -cd -p {threads} {input.s3} > {output.out3}; "
        "pigz -cd -p {threads} {input.s4} > {output.out4}; "
        "pigz -cd -p {threads} {input.s5} > {output.out5}; "
        "pigz -cd -p {threads} {input.s6} > {output.out6}; "
        "pigz -cd -p {threads} {input.s7} > {output.out7}; "
        "pigz -cd -p {threads} {input.s8} > {output.out8}; "
        "pigz -cd -p {threads} {input.s9} > {output.out9}; "
        "pigz -cd -p {threads} {input.s10} > {output.out10}; "
        "pigz -cd -p {threads} {input.s11} > {output.out11}; "
        "pigz -cd -p {threads} {input.s21} > {output.out21}; "
        "pigz -cd -p {threads} {input.s22} > {output.out22}; "
        "pigz -cd -p {threads} {input.s23} > {output.out23}; "
        "pigz -cd -p {threads} {input.s24} > {output.out24}; "
        "pigz -cd -p {threads} {input.s25} > {output.out25}; "
        "pigz -cd -p {threads} {input.s26} > {output.out26}; "
        "pigz -cd -p {threads} {input.s27} > {output.out27}; "
        "pigz -cd -p {threads} {input.s28} > {output.out28}; "
        "pigz -cd -p {threads} {input.s29} > {output.out29}; "
        "pigz -cd -p {threads} {input.s210} > {output.out210}; "
        "pigz -cd -p {threads} {input.s211} > {output.out211} "

