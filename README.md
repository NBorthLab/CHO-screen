# CHO-screen

This repository contains the analysis pipeline and relevant scripts described in the [original publication title] [doi]. 
The work aims to identify using a genome-scale CRISPR screen approach.

[abstract]


```
## Filestructure
./
├── config/                                 # (provide) configurations
├── resources/                              # Raw-data and other resources
│   ├── aligned/                            # Intermediate results
│   ├── flowcell_merged/                    # Intermediate results
│   ├── raw/                                # (provide) Raw-data
│   │  ├── 53h_dense.bed                    # (provide) Chromatin states****
│   │  ├── aligned_sorted_INH040x.bam       # (provide) Transcripomics data
│   │  ├── reactive-genes.txt               # (provide) List of reactive Genes
│   │  ├── non-reactive-genes.txt           # (provide) List of non reactive Genes
│   │  ├── library_for_mageck_count.txt     # (provide) pgRNA library used for MAGeCK count
│   │  ├── pgRNA-gene-id-sequence-list.fa   # (provide) pgRNA library used for Bowtie2
|   |  ├── flowcell_1/                      # (provide) .fasta.gz from sequencing
|   |  └── flowcell_2/                      # (provide) .fasta.gz from sequencing
|   ├── table_bowtieM/                      # Intermediate results 
|   |  └── results/                         # Intermediate results
|   └── table_mageck/                       # Intermediate results
|      └── results/                         # Intermediate results
├── results/                                # Intermediate results
│   ├── table_bowtieM/                      # Intermediate results
|   └── table_mageck/                       # Intermediate results 
├── results_overview/                       # Results produced in the workflow
│   └── table_bowtieM/                      # Results produced in the workflow
|      ├── results/                         # Results produced in the workflow
|      └── table_mageck/                    # Results produced in the workflow
├── workflow/                               # Workflow definitions
|   ├── profiles/                           # Snakemake profiles
│   ├── envs/                               # Conda environments
│   ├── rules/                              # Snakemake rules
│   ├── scripts/                            # Scripts
│   └── Snakefile                           # Main Snakefile
└── README.md
```

### preprocess for workflow:
- update parameters within ./config/config.yaml
- provide ./raw/53h_dense.bed in structure
```
$ head 53h_dense.bed 
#track name="Tp4" description="Tp4 (Emission ordered)" visibility=1 itemRgb="On"
chr10   0       21000   Repressed heterochromatin       0       .       0       21000   #3399ff
chr10   21000   27400   Quiescent/low   0       .       21000   27400   #d9d9d9
chr10   27400   47200   Repressed heterochromatin       0       .       27400   47200   #3399ff
chr10   47200   77800   Quiescent/low   0       .       47200   77800   #d9d9d9
chr10   77800   79200   Repressed heterochromatin       0       .       77800   79200   #3399ff
chr10   79200   146600  Quiescent/low   0       .       79200   146600  #d9d9d9
chr10   146600  155400  Repressed heterochromatin       0       .       146600  155400  #3399ff
chr10   155400  171000  Quiescent/low   0       .       155400  171000  #d9d9d9
chr10   171000  173000  Repressed heterochromatin       0       .       171000  173000  #3399ff
```
- provide ./raw/aligned_sorted_INH040x.bam
```
$ samtools view aligned_sorted_INH0401.bam | head
7001253F:661:CD4LWANXX:1:2313:13165:87523#ACATCCGACTGCGGAT      16      NW_023276806.1  7549    255     20M1671N33M7289N22M     *       0       0       ACTAACCCTAACCCTAACCCTACCCCTCTAACCCTAAC
CCTATCCCTAACCCTCTAACCCTAACTCTAACCCTCT   FFFFFFFFFFFFFFFFFFFFBBFFFFFFFFFFFFFFFFFFFFFFFFFFBFFF<FFFFFFFFF<FFFFFFFBBBBB     NH:i:1  HI:i:1  AS:i:48 nM:i:8
7001253F:661:CD4LWANXX:1:2313:13165:87523#ACATCCGACTGCGGAT      16      NW_023276806.1  7549    255     20M1671N33M7289N22M     *       0       0       ACTAACCCTAACCCTAACCCTACCCCTCTAACCCTAAC
CCTATCCCTAACCCTCTAACCCTAACTCTAACCCTCT   FFFFFFFFFFFFFFFFFFFFBBFFFFFFFFFFFFFFFFFFFFFFFFFFBFFF<FFFFFFFFF<FFFFFFFBBBBB     NH:i:1  HI:i:1  AS:i:48 nM:i:8
7001253F:661:CD4LWANXX:2:2205:17783:33201#ACATCCGACTGCGGAT      0       NW_023276806.1  29163   255     100M    *       0       0       AGTGGACTGAAGATTTTAAAATAGCGGAGTGGATTAATGACACTATGACAACAG
AGAAATGGTAAATATATGAGCTGAGAGATCTGACAACCATTCACAT  BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NH:i:1  HI:i:1  AS:i:98 nM:i:0
7001253F:661:CD4LWANXX:2:2205:17783:33201#ACATCCGACTGCGGAT      0       NW_023276806.1  29163   255     100M    *       0       0       AGTGGACTGAAGATTTTAAAATAGCGGAGTGGATTAATGACACTATGACAACAG
AGAAATGGTAAATATATGAGCTGAGAGATCTGACAACCATTCACAT  BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NH:i:1  HI:i:1  AS:i:98 nM:i:0
7001253F:661:CD4LWANXX:1:1215:19053:36813#ACATCCGACTGCGGAT      16      NW_023276806.1  29960   255     100M    *       0       0       AGTTAATTATTTAAATGTTAATCATCAACATTCTCATTCTCTGAAGACTGTGAT
ACATTTATTTTATATTTCTTATTTAAAATTATTTATTTAATGCATG  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB    NH:i:1  HI:i:1  AS:i:98 nM:i:0
7001253F:661:CD4LWANXX:1:1215:19053:36813#ACATCCGACTGCGGAT      16      NW_023276806.1  29960   255     100M    *       0       0       AGTTAATTATTTAAATGTTAATCATCAACATTCTCATTCTCTGAAGACTGTGAT
ACATTTATTTTATATTTCTTATTTAAAATTATTTATTTAATGCATG  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB    NH:i:1  HI:i:1  AS:i:98 nM:i:0
7001253F:661:CD4LWANXX:1:2205:14024:84479#ACATCCGACTGCGGAT      256     NW_023276806.1  32720   1       99M1S   *       0       0       GCCTTAATTTGTCACCACCATGGGATGGGTTAGTAGAAGCCTTAACTCTAGTCT
TGATTATTCTTCTTTTGTTCATATTGGTCTTATATTGTGCTCATCA  BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFB    NH:i:3  HI:i:3  AS:i:97 nM:i:0
7001253F:661:CD4LWANXX:1:2205:14024:84479#ACATCCGACTGCGGAT      256     NW_023276806.1  32720   1       99M1S   *       0       0       GCCTTAATTTGTCACCACCATGGGATGGGTTAGTAGAAGCCTTAACTCTAGTCT
TGATTATTCTTCTTTTGTTCATATTGGTCTTATATTGTGCTCATCA  BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFB    NH:i:3  HI:i:3  AS:i:97 nM:i:0
7001253F:661:CD4LWANXX:1:1101:1941:75923#ACATCCGACTGCGGAT       272     NW_023276806.1  34912   3       3S96M   *       0       0       TCTTGGCAAGACACCGGCCAGGGTCCACGCTAGTTATGGGAATCCAGCAAGAAA
ACAGACAGATCAGAGAATTGCAACAGGAGAACAAAGAACTGCGAA   7///FFFB///B<//</<<///<///7<///<FBFFFFF<//F</F/F<//<B/F/FFFFB<BF<FFFB/<</FF/F/FFB</<FFFF<FFF/FBBB/B     NH:i:2  HI:i:2  AS:i:90 nM:i:2
7001253F:661:CD4LWANXX:1:1101:1941:75923#ACATCCGACTGCGGAT       272     NW_023276806.1  34912   3       3S96M   *       0       0       TCTTGGCAAGACACCGGCCAGGGTCCACGCTAGTTATGGGAATCCAGCAAGAAA
ACAGACAGATCAGAGAATTGCAACAGGAGAACAAAGAACTGCGAA   7///FFFB///B<//</<<///<///7<///<FBFFFFF<//F</F/F<//<B/F/FFFFB<BF<FFFB/<</FF/F/FFB</<FFFF<FFF/FBBB/B     NH:i:2  HI:i:2  AS:i:90 nM:i:2
```
- provide ./raw/reactive-genes.txt & ./raw/non-reactive-genes.txt
```
$ head non-reactive-genes.txt 
Rpe65
Nexn
Ifi44l
Ttll7
Dnai3
LOC100762860
LOC100762574
LOC100752806
Gbp5
LOC100750831
```
- provide ./raw/library_for_mageck_count.txt
```
$ head library_for_mageck_count.txt 
>595-w10_NC_048595_1_1499500-1499519_NC_048595_1_1650461-1650480        GAAAACCTGATTGACACATG    595-w10
>595-w10_NC_048595_1_1499143-1499162_NC_048595_1_1650498-1650517        CCACCACTAGAAACAGGGAG    595-w10
>595-w10_NC_048595_1_1499243-1499262_NC_048595_1_1650893-1650912        CAGATGGTGGGACTATTGGG    595-w10
>595-w10_NC_048595_1_1499516-1499535_NC_048595_1_1650140-1650159        CATGTGGTGACCTTGGAACG    595-w10
>595-w10_NC_048595_1_1499812-1499831_NC_048595_1_1650218-1650237        CACAGGGAAAGAGCCTGGTG    595-w10
>595-w10_NC_048595_1_1499796-1499815_NC_048595_1_1650041-1650060        CTGTAGACTCTACAATCACA    595-w10
>595-w10_NC_048595_1_1499810-1499829_NC_048595_1_1650385-1650404        ATCACAGGGAAAGAGCCTGG    595-w10
>595-w10_NC_048595_1_1499368-1499387_NC_048595_1_1650079-1650098        CCACTGTGACACTCTGCCAT    595-w10
>595-w100_NC_048595_1_14999285-14999304_NC_048595_1_15150001-15150020   TGTGTGTGCACATACCCCTG    595-w100
>595-w100_NC_048595_1_14999286-14999305_NC_048595_1_15150169-15150188   GTGTGTGCACATACCCCTGA    595-w100
```
- provide ./raw/library_for_mageck_count.txt
```
$ head pgRNA_gene_id_sequence_list.fa
>>595-w10_NC_048595_1_1499500-1499519_NC_048595_1_1650461-1650480
GAAAACCTGATTGACACATG
>>595-w10_NC_048595_1_1499143-1499162_NC_048595_1_1650498-1650517
CCACCACTAGAAACAGGGAG
>>595-w10_NC_048595_1_1499243-1499262_NC_048595_1_1650893-1650912
CAGATGGTGGGACTATTGGG
>>595-w10_NC_048595_1_1499516-1499535_NC_048595_1_1650140-1650159
CATGTGGTGACCTTGGAACG
>>595-w10_NC_048595_1_1499812-1499831_NC_048595_1_1650218-1650237
CACAGGGAAAGAGCCTGGTG
```
### general workflow:
1. Rename and unzip files using 'pigz'. ([01_prep_data.smk](workflow/rules/01_prep_data.smk)
2. Merge both flowcells using 'cat'. ([02_merge_flowcells.smk](workflow/rules/02_merge_flowcells.smk)
3. Merge paired end reads into an overlapping fragment using ''. ([03_merge_overlapping.smk](workflow/rules/03_merge_overlapping.smk)
4. Sort for read directions using 'awk'. ([04_sort_directions.smk](workflow/rules/04_sort_directions.smk)
5. Turn reversed reads using ''. ([05_turn_reverse.smk](workflow/rules/05_turn_reverse.smk)
6. Combine amplicons using 'cat'. ([06_combine_amp.smk](workflow/rules/06_combine_amp.smk)
7. Trimming reads using ''. ([07_trim.smk](workflow/rules/07_trim.smk)
8. Alignment using 'Bowtie2' ([08_bowtei2.smk](workflow/rules/08_bowtie2.smk)
9. Change .sam to .bam using 'samtools'. ([09_samtools.smk](workflow/rules/09_samtools.smk)
10. Generate counttable out of .bam-files using 'mageck' ([10_2_mageck_count_bam.smk](workflow/rules/10_2_mageck_count_bam.smk)
11. Generate counttable out of .fasta-files using 'mageck' ([11_2_mageck_count.smk](workflow/rules/11_2_mageck_count.smk)
12. Remove low-count guides ([12_2_remove_zeros.smk](workflow/rules/12_2_remove_zeros.smk)
13. Descriptive statistics ([13_2_descriptive_statistics.smk](workflow/rules/13_2_descriptive_statistics.smk)
14. Prepare for ranking ([14_2_prepare_mageck_test.smk](workflow/rules/14_2_prepare_mageck_test.smk)
15. Move files to different folder ([15_2_move_files.smk](workflow/rules/15_2_move_files.smk)
16. Guide ranking ([17_2_test_mageck_rra.smk](workflow/rules/17_2_test_mageck_rra.smk)
17. Combine results from "counttable .bam" and "counttable .fasta" ([18_2_combine_tables.smk](workflow/rules/18_2_combine_tables.smk)
18. Classify non-essential regions ([19_2_non_essential.smk](workflow/rules/19_2_non_essential.smk)
19. combine informatino in a table ([20_summary_table.smk](workflow/rules/20_summary_table.smk)
20. Find transcripts in essential regions using 'Rsubread' ([21_find_transcripts.smk](workflow/rules/21_find_transcripts.smk)
21. Normalisation of transcript-data ([22_normailse_transcripts.smk](workflow/rules/22_normailse_transcripts.smk)
22. Plot essential regions ([23_plot_region_details.smk](workflow/rules/23_plot_region_details.smk)
23. Sort results ([24_sort_results.smk](workflow/rules/24_sort_results.smk)
24. Create a summary table ([25_create_summary_table.smk](workflow/rules/25_create_summary_table.smk)

