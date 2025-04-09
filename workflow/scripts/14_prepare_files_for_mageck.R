folder_args <- commandArgs(trailingOnly = TRUE)

path_to_input_file_counttable <- folder_args[1]
path_to_config_file <- folder_args[2]
path_to_output_folder <- folder_args[3]

library(dplyr)
library(stringr)


# functions:
rownames_from_txt <- function(name_of_file) {
              return(scan(
                file = name_of_file,
                what = character(),
                flush = TRUE,
                sep = "\t"
              ))
            }


read_new_file <- function(path_to_file) {spaltenbeschriftung_newfile <-
           scan(file = path_to_file, what = character(),
           nlines = 1, sep = "\t")
            newfile <- matrix(scan( file = path_to_file, skip = 1,
                            what = character(), sep = "\t"),
                            ncol = length(spaltenbeschriftung_newfile), byrow = TRUE)
          newfile <- as.data.frame(newfile)
          colnames(newfile) <- spaltenbeschriftung_newfile
          return(newfile)
          }


read_from_config_file <- function(config_file_id, conf_file = config_file){
              x <- (sub(config_file_id, "", conf_file[substr(conf_file, 1, length(unlist(strsplit(config_file_id,NULL))))
                                                      == config_file_id]))
              x <- unlist(strsplit(x, "[,]"))
              return(x)
            }


convert_df_to_numeric    <- function(count_table, skip_columns = c(1,2)) {
        count_table_colnames <- colnames(count_table)
        gene_names <- as.data.frame(count_table[,skip_columns])
        s_ <- sapply(count_table[,-skip_columns],as.numeric)
        data_data_frame <- data.frame(gene_names, s_)
        colnames(data_data_frame) <- count_table_colnames
        return(data_data_frame)
      }


save_table_to_csv <-  function(filename_table_input, outputfolder, c_row_names = FALSE, c_col_names = FALSE) {
            write.table(as.matrix(filename_table_input),
              sep = "\t", file = paste0(outputfolder, substitute(filename_table_input), ".txt"),
              quote = FALSE, row.names = c_row_names, col.names = c_col_names
            )
}


####data analyses####
count_table <- read_new_file(path_to_input_file_counttable)
count_table <- convert_df_to_numeric(count_table,c(1,2))

config_file <- scan(path_to_config_file, sep = ";", what = character())

project_run_name <- read_from_config_file("run_name: ")

samplenumbers_control <- (read_from_config_file("mageck_test_mageck_control: "))
samplenumbers_treated <- (read_from_config_file("mageck_test_mageck_treated: "))
samplenumbers_day_0_as_control <- (read_from_config_file
                                                     ("mageck_test_mageck_day_0_as_control: "))

### count table for analyses
samples_to_analyse <- c("Gene", "sgRNA",samplenumbers_day_0_as_control, samplenumbers_treated, samplenumbers_control)
samplenames <- colnames(count_table)
samplenames_to_select <- ifelse(samplenames %in% samples_to_analyse, TRUE, FALSE)
mageck_count_table_all <- subset(count_table, select =  samplenames_to_select)
rm(samples_to_analyse, samplenames, samplenames_to_select)

to_select_region <- grepl("-w", mageck_count_table_all$Gene) 
mageck_count_table_region <- subset(mageck_count_table_all, subset = to_select_region)
write.table(as.matrix(mageck_count_table_region),
              sep = "\t", file = paste0(path_to_output_folder, project_run_name, ".count_table_region.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE
            )
