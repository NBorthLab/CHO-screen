folder_args <- commandArgs(trailingOnly = TRUE)

inputfile <- folder_args[1]
outputfolder <- folder_args[2]

library(Rsubread)
library(data.table)
library(dplyr)
library(stringr)
library(edgeR)
library(TTR)


convert_df_to_numeric <- function(count_table, skip_columns = c(1,2)) {
  
  if (!require("dplyr", character.only=T, quietly=T)) {
    install.packages("dplyr", repos="http://cran.at.r-project.org")
    library("dplyr", character.only=T)
  }
  library("dplyr", character.only=T)
  
  count_table <- dplyr::select(count_table, all_of(skip_columns), everything())
  count_table_colnames <- colnames(count_table)
  gene_names <- as.data.frame(count_table[,1:length(skip_columns)])
  colnames(gene_names) <- count_table_colnames[1:length(skip_columns)]
  s_ <- sapply(count_table[,!(colnames(count_table) %in% colnames(gene_names))],as.integer)
  data_data_frame <- data.frame(gene_names, s_)
  colnames(data_data_frame) <- count_table_colnames
  return(data_data_frame)
}

read_new_file <- function(path_to_file) {
  spaltenbeschriftung_newfile <-
    scan(
      file = path_to_file,
      what = character(),
      nlines = 1,
      sep = "\t"
    )
  newfile <-
    matrix(
      scan(
        file = path_to_file ,
        skip = 1,
        what = character(),
        sep = "\t"
      ),
      ncol = length(spaltenbeschriftung_newfile),
      byrow = TRUE
    )
  newfile <- as.data.frame(newfile)
  colnames(newfile) <- spaltenbeschriftung_newfile
  return(newfile)
}

transcripts_numbers <- read_new_file(inputfile)
transcripts_numbers <- convert_df_to_numeric(transcripts_numbers, c(1:5))

colnames(transcripts_numbers) <- 
  ifelse(grepl("^aligned", colnames(transcripts_numbers)), paste0("s_", gsub("\\D","", colnames(transcripts_numbers))),colnames(transcripts_numbers))


tmm_factor <- calcNormFactors(transcripts_numbers[,grepl("^s_", colnames(transcripts_numbers))]) 
tmm_norm_table <- t(t(transcripts_numbers[,grepl("^s_", colnames(transcripts_numbers))]) * tmm_factor)

# Matrix of normalized logCPM values
normalized_logcpm <- cpm(tmm_norm_table, log = TRUE)
rm(tmm_norm_table)
data_logcpm <- cbind(transcripts_numbers[,!grepl("^s_", colnames(transcripts_numbers))],normalized_logcpm)
rm(normalized_logcpm)
data_logcpm$geo_mean <- rowMeans((data_logcpm[, grep("^s_", names(data_logcpm))]))


trans_numbers <- data_logcpm


# Calculate the rolling mean with a window size of 3
rolling_mean <- runMin(trans_numbers$geo_mean, n = 9, cumulative = FALSE)

rolling_mean[is.na(rolling_mean)] <- -1


trans_numbers$min <- rolling_mean
# Print the result
setDT(trans_numbers)

trans_numbers_ha <- trans_numbers[, ':=' (ID = str_split_i(GeneID, "_", 1),
                                          subset = str_split_i(GeneID, "_",2)),
                                 ]

transcripted <- trans_numbers_ha[,.(any_transcripts = any(min > 0)), by = ID]

genome_data <- merge(trans_numbers_ha, transcripted, all=TRUE)

haha_new <- genome_data

write.table(
  as.matrix(haha_new),
  sep = "\t",
  file = paste0(outputfolder),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

