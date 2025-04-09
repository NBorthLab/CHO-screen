
folder_args <- commandArgs(trailingOnly = TRUE)

path_to_input_file_counttable <- folder_args[1]
path_to_config_file <- folder_args[2]
table_version <- folder_args[3]
run_name <- folder_args[4]
path_outputfolder <- folder_args[5]

outputfolder <- paste0(path_outputfolder, "table_", table_version, "/")

library(dplyr)
library(stringr)


read_new_file <- 
function(path_to_file) {
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

convert_df_to_numeric <- 
function(count_table, skip_columns = c(1,2)) {
  count_table_colnames <- colnames(count_table)
  gene_names <- as.data.frame(count_table[,skip_columns])
  s_ <- sapply(count_table[,-skip_columns],as.numeric)
  data_data_frame <- data.frame(gene_names, s_)
  colnames(data_data_frame) <- count_table_colnames
  return(data_data_frame)
}

scan_config_file <- 
function(path_to_config_file, sep = ";") {
  scan(path_to_config_file, sep = ";", what = character())
}

read_from_config_table <- 
function(config_content, conf_file = config_table){
  x <- (sub(config_content, "", conf_file[substr(conf_file, 1, length(unlist(strsplit(config_content,NULL))))
                                          == config_content]))
  x <- unlist(strsplit(x, "[,]"))
  return(x)
}

save_table_to_txt <- 
function(filename_table_input, outputfolder, c_row_names = FALSE, c_col_names = FALSE) {
  write.table(as.matrix(filename_table_input),
              sep = "\t", file = paste0(outputfolder, substitute(filename_table_input), ".txt"),
              quote = FALSE, row.names = c_row_names, col.names = c_col_names
  )
}


####data:

count_table <- read_new_file(path_to_input_file_counttable) 
str(count_table)
count_table <- convert_df_to_numeric(count_table, c(1,2))

config_table <- scan_config_file(path_to_config_file)
count_cut_off <- as.numeric(read_from_config_table("count_cut_off_per_guide: "))
samplenumbers_control <- read_from_config_table("mageck_test_d0_control: ")
samplenumbers_treated <- read_from_config_table("mageck_test_d0_treated: ")
samplenumbers_day_0_as_control <- read_from_config_table("mageck_test_day_0_as_control: ") 

#remove all rows all zeros:
if (ncol(count_table) > 3) {
    count_table_no_zero <- subset(count_table[rowSums(count_table[,3:ncol(count_table)]) > 0,])
  count_table_all_zero <- subset(count_table[!rowSums(count_table[,3:ncol(count_table)]) > 0,])
} else {
  count_table_no_zero <- subset(count_table[(count_table[,3]) > 0,])
  count_table_all_zero <- subset(count_table[!(count_table[,3]) > 0,])
}

save_table_to_txt(count_table_all_zero, outputfolder, FALSE, TRUE)

#remove all rows which countsums in controls are below count_cut_off:
samples_to_analyse <- c(samplenumbers_day_0_as_control, samplenumbers_control)
samplenames <- colnames(count_table)
samplenames_to_select <- ifelse(samplenames %in% samples_to_analyse, TRUE, FALSE)

## save above and below cutoff, if it is provided, elso return to "no_zero":
if (length(count_cut_off) != 0) {
    if (ncol(count_table) > 3) {      
  count_table_above_cutoff <- subset(count_table_no_zero[rowSums(count_table_no_zero[,samplenames_to_select]) > count_cut_off,])
  count_table_below_cutoff <- subset(count_table_no_zero[!rowSums(count_table_no_zero[,samplenames_to_select]) > count_cut_off,])
    } else {
  count_table_above_cutoff <- subset(count_table_no_zero[(count_table_no_zero[,3]) > count_cut_off,])
  count_table_below_cutoff <- subset(count_table_no_zero[!(count_table_no_zero[,3]) > 0,])
    }

save_table_to_txt(count_table_below_cutoff, outputfolder, FALSE, TRUE)
    
write.table(as.matrix(count_table_above_cutoff),
              sep = "\t", file = paste0(outputfolder, run_name, ".count_nozero.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}

if (!exists("count_table_above_cutoff")) {
  write.table(as.matrix(count_table_no_zero),
              sep = "\t", file = paste0(outputfolder, run_name, ".count_nozero.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  }




