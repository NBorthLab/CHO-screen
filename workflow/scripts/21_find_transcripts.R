folder_args <- commandArgs(trailingOnly = TRUE)

input_region_table <- folder_args[1]
input_selection_table <- folder_args[2]
input_path_to_bam <- folder_args[3]
input_path_to_output_folder <- folder_args[4]
input_folder_for_selection_file <- folder_args[5]


library(Rsubread)
library(data.table)
library(dplyr)
library(stringr)

for_selection <- TRUE

region_table <- input_region_table
selection_table <- input_selection_table
output_folder <- input_path_to_output_folder


convert_df_to_numeric <- function(count_table, skip_columns = c(1,2)) {
  
  if (!require("dplyr", character.only=T, quietly=T)) {
    install.packages("dplyr", repos="http://cran.at.r-project.org")
    library("dplyr", character.only=T)
  }
  library("dplyr", character.only=T)
  
  count_table <- select(count_table, all_of(skip_columns), everything())
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
  newfile <- as.data.frame(newfile, stringsAsFactors = FALSE)
  colnames(newfile) <- spaltenbeschriftung_newfile
  return(newfile)
}


bamfiles <- character()

for (i in 1:2){
  x <- str_pad(i, 2, side = "left", pad = "0")
  bamf <- paste0(input_path_to_bam, "/aligned_sorted_INH04" ,x ,".bam")
  bamfiles <- append(bamfiles, bamf)
}


###filter for interesting regions to reduce file-size
daten <- read_new_file(region_table)
print("read region_table")
if (for_selection == TRUE) {
    selection_sel <- read_new_file(selection_table)
    rownames( selection_sel) <-  selection_sel$id
    daten <- daten[daten$ID %in% rownames( selection_sel),]
} else {
  daten <- subset(daten, classification == "essential" | classification == "_" | is.na(Genes),
                     select = c(ID, Chromosome, Region_Start, Region_End))
}


empty_df <- data.frame(matrix(ncol = length(daten), nrow = 0))
colnames(empty_df) <- colnames(daten)
empty_df$ID <- as.character(empty_df$ID)
empty_df$Chromosome <- as.character(empty_df$Chromosome)

region_table <- empty_df

for (i in 1:nrow(daten)) {
  ID <- as.character(daten[i,"ID"]) 
  bereich_start <- as.integer(daten[i,"Region_Start"])
  bereich_ende <- as.integer(daten[i,"Region_End"])
  bereich_laenge <- length(bereich_start:bereich_ende)
  bereich_seq <- seq(bereich_start, bereich_ende, by =  10)    ##set region size here
  chromosome <- as.character(daten[i,"Chromosome"])
  
  id_vec <- as.character()
  for(kb in 1:(length(bereich_seq)-1)) {
   id_vec <- c(id_vec, as.character(paste0(ID, "_", kb)))
  }

  bereich_df <- data.table(Gene = ID,
                           id = id_vec,
                           Region_Start = bereich_seq[1:(length(bereich_seq)-1)],
                           Region_End = bereich_seq[2:length(bereich_seq)],
                           Chromosome = chromosome)
  
  region_table <- full_join(region_table, bereich_df)
  
  print(i)
}


region_table <- data.frame(
  GeneID=region_table$id,
  Gene=region_table$Gene,
  Chr=region_table$Chromosome,
  Start=as.integer(region_table$Region_Start),
  End=as.integer(region_table$Region_End),
  Strand="+",
  stringsAsFactors=FALSE)


transcripts_fc <- featureCounts(bamfiles,annot.ext=region_table, allowMultiOverlap = T)

transcripts_table <- data.frame(do.call(cbind, transcripts_fc$annotation)
                                , transcripts_fc$counts, stringsAsFactors = FALSE)

write.table(
  as.matrix(transcripts_table),
  sep = "\t",
 file = paste0(output_folder,"transcripts_", str_replace(input_selection_table, input_folder_for_selection_file, "")),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

