folder_args <- commandArgs(trailingOnly = TRUE)

path_to_input_summary_all_tables <- folder_args[1]
path_to_input_file_least_essential <- folder_args[2]
target_gene_all <- folder_args[3]
chromatin_states_filename <- folder_args[4]
input_folder_chromatin_states <- folder_args[5]
outputfolder <- folder_args[6]

library(dplyr)
library(data.table)
library(stringr)


read_new_file <- function(path_to_file) {spaltenbeschriftung_newfile <-
  scan(file = path_to_file, what = character(),
       nlines = 1, sep = "\t")
newfile <- matrix(scan( file = path_to_file, skip = 1,
                        what = character(), sep = "\t"),
                  ncol = length(spaltenbeschriftung_newfile), byrow = TRUE)
newfile <- newfile <- as.data.frame(newfile, stringsAsFactors = FALSE)
colnames(newfile) <- spaltenbeschriftung_newfile
return(newfile)
}

save_table_to_csv <-  function(filename_table_input, outputfolder, c_row_names = FALSE, c_col_names = FALSE) {
  write.table(as.matrix(filename_table_input),
              sep = "\t", file = paste0(outputfolder, substitute(filename_table_input), ".txt"),
              quote = FALSE, row.names = c_row_names, col.names = c_col_names
  )
}



## update averything summary_tabla with classification:

everything_for_fede <- read_new_file(path_to_input_summary_all_tables)
selected_least_essential <- read_new_file(path_to_input_file_least_essential)
target_gene_all <- read_new_file(target_gene_all)

everything_for_fede$classification <- rep("undefined", length(everything_for_fede$id))
everything_for_fede$classification <- ifelse(everything_for_fede$id %in% target_gene_all$id, "essential",
                                             ifelse(everything_for_fede$id %in% selected_least_essential$id, "non_essential", "undefined"))

everything_for_fede <- select(everything_for_fede, id, chromosome, windowbegin, windowend, LFC, FDR, classification, Genes, everything())
colnames(everything_for_fede) <- c("ID", "Chromosome", "Region_Start", "Region_End", "LFC", "FDR", "classification", "Genes",
                                     "num.bowtie", "neg_score.bowtie", "neg_p_value.bowtie", "neg_fdr.bowtie",
                                     "neg_rank.bowtie", "neg_goodsgrna.bowtie", "neg_lfc.bowtie", "pos_score.bowtie",
                                     "pos_p_value.bowtie", "pos_fdr.bowtie", "pos_rank.bowtie", "pos_goodsgrna.bowtie",
                                     "pos_lfc.bowtie", "table.bowtie", "num.mageck", "neg_score.mageck", "neg_p_value.mageck",
                                     "neg_fdr.mageck", "neg_rank.mageck", "neg_goodsgrna.mageck", "neg_lfc.mageck",
                                     "pos_score.mageck", "pos_p_value.mageck", "pos_fdr.mageck", "pos_rank.mageck",
                                     "pos_goodsgrna.mageck", "pos_lfc.mageck", "table.mageck", "LFC.bowtie", "FDR.bowtie",
                                     "LFC.mageck", "FDR.mageck")          

summary_table <- everything_for_fede


# Subset the data
filtered_data <- summary_table %>%
  filter(FDR <= 0.25, LFC <= log2(0.85))


chromatin_states_list <- lapply(chromatin_states_filename, function(path_to_file) {
  spaltenbeschriftung_newfile <- c("chromosome", "start", "end", "name", "something_1", "some_dot", "start_2", "end_2", "colorcode" )
  newfile <- matrix(scan(file = path_to_file, skip = 1, what = character(), sep = "\t"),
                    ncol = length(spaltenbeschriftung_newfile), byrow = TRUE)
  newfile <- as.data.frame(newfile, stringsAsFactors = FALSE)
  colnames(newfile) <- spaltenbeschriftung_newfile
  return(newfile)
})
names(chromatin_states_list) <- str_sub(chromatin_states_filename, 42, -5)



chrrrs <- c( "chr1_0", "chr1_1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
             "chr8", "chr9", "chr10", "chrX", "unplaced_scaffold_1")
chromosomessss <- c("NW_023276806.1", "NW_023276807.1", "NC_048595.1",
                    "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1",
                    "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1",
                    "NC_048604.1", "NW_023276808.1")

chrrrs <- setNames(chrrrs, chromosomessss)
chromosomessss <- setNames(chromosomessss, chrrrs)


chrom_state_for_cytobands <- as.data.table(chromatin_states_list[[1]])

chrom_state_for_cytobands$chromosome <- chromosomessss[chrom_state_for_cytobands$chromosome]

chrom_state_for_cytobands$chromosome[is.na(chrom_state_for_cytobands$chromosome)] <- "unplaced_scaff_xxx"
chrom_state_for_cytobands$gieStain <- chrom_state_for_cytobands$name


all_gene_info<- summary_table
all_gene_info$Region_Start <- as.numeric(all_gene_info$Region_Start)
all_gene_info$Region_End <- as.numeric(all_gene_info$Region_End)
dt_all_gene_info <- data.table(all_gene_info[,1:8])
dput(head(all_gene_info))

dt_chrom_states <- data.table(chrom_state_for_cytobands[,1:4])
colnames(dt_chrom_states) <- c("Chromosome", "Region_Start", "Region_End", "status")
dt_chrom_states$Region_Start <- as.numeric(dt_chrom_states$Region_Start)
dt_chrom_states$Region_End <- as.numeric(dt_chrom_states$Region_End)

dt_all_gene_info_new <- dt_all_gene_info

chromatin_states <- c("ID")

empty_df <- data.frame(matrix(ncol = length(chromatin_states), nrow = 0))
colnames(empty_df) <- chromatin_states
empty_df$ID <- as.character(empty_df$ID)

region_table <- empty_df
sel_row <- 1


for (sel_row in 1:nrow(all_gene_info)) { 
  
ID <- as.character(dt_all_gene_info[sel_row,"ID"])
bereich_start <- as.numeric(dt_all_gene_info[sel_row,"Region_Start"])
bereich_ende <- as.numeric(dt_all_gene_info[sel_row,"Region_End"])
bereich_laenge <- length(bereich_start:bereich_ende)
chromosome <- as.character(dt_all_gene_info[sel_row,"Chromosome"])

bereich_df <- data.table(Gene = rep(ID, bereich_laenge),
                         Region_Start = c(bereich_start:bereich_ende),
                         Region_End = c(bereich_start:bereich_ende),
                         Chromosome = rep(chromosome, bereich_laenge))

setkey(dt_chrom_states, Chromosome, Region_Start, Region_End)
states_in_region <- foverlaps(bereich_df, dt_chrom_states, type="any")
states <- summary(as.factor(states_in_region$status))
states_perc <- states / sum(states) * 100

region_table <- full_join(region_table,
                      as.data.frame(cbind(ID, t(states_perc))))
print(sel_row)
}

dt_all_gene_info_new <- full_join(dt_all_gene_info, region_table)

genes_chrom_states <- dt_all_gene_info_new

colnames(genes_chrom_states) <- c("ID",  "Chromosome",  "Region_Start",  "Region_End",  "LFC",  "FDR",  "classification",  "Genes",
                                    "Active Promoter",  "Flanking TSS downstream",
                                   "Polycomb repressed regions",  "Quiescent/low",  "Strong Transcription",  "Weak enhancer",
                                  "Repressed heterochromatin", "Weak genic enhancer", "Active Enhancer 1",  "Flanking TSS upstream",
                                   "Active Enhancer 2",  "unknown")



save_table_to_csv(summary_table, outputfolder, FALSE, TRUE)
save_table_to_csv(genes_chrom_states, outputfolder, FALSE, TRUE)
