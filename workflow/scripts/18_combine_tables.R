
folder_args <- commandArgs(trailingOnly = TRUE)

gene_mageck <- folder_args[1]
gene_bowtie <- folder_args[2]
input_annotation_file <- folder_args[3]
path_to_input_file_genetable_with_names <- folder_args[4]
run_name <- folder_args[5]
outputfolder <- folder_args[6]


## cutoffs:
lfc_cutoff <- log2(0.85)
fdr_cutoff <- 0.25

library(ggplot2)
library(venn)
library(genomeIntervals)
library(dplyr)
library(readxl)
library(stringi)
library(stringr)
library(data.table)
library(tidyr)

colors_venn <- c('#d9f0d3',"#CCEBC5" )


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

plot_to_png <- function(plot_titel, outputfolder){
  png(file = paste0(outputfolder, run_name, "_", deparse(substitute(plot_titel)),".png"),width = 1200, height = 740)
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




#########data:
###bowtie:

gene_summary_bowtie <- read_new_file(gene_bowtie)
gene_summary_bowtie <- select(gene_summary_bowtie, 1, everything())
gene_summary_bowtie <- convert_df_to_numeric(gene_summary_bowtie, c(1))
gene_summary_bowtie$table <- rep("bowtie", length(gene_summary_bowtie$`neg|rank`))
colnames(gene_summary_bowtie) <- gsub("\\|" , "_", colnames(gene_summary_bowtie))
colnames(gene_summary_bowtie) <- gsub("\\|" , "_", colnames(gene_summary_bowtie))
colnames(gene_summary_bowtie) <- gsub("-" , "_", colnames(gene_summary_bowtie))


###mageck:

gene_summary_mageck <- read_new_file(gene_mageck)
gene_summary_mageck <- select(gene_summary_mageck, 1, everything())
gene_summary_mageck <- convert_df_to_numeric(gene_summary_mageck, c(1))
gene_summary_mageck$table <- rep("mageck", length(gene_summary_mageck$`neg|rank`))
colnames(gene_summary_mageck) <- gsub("\\|" , "_", colnames(gene_summary_mageck))
colnames(gene_summary_mageck) <- gsub("\\|" , "_", colnames(gene_summary_mageck))
colnames(gene_summary_mageck) <- gsub("-" , "_", colnames(gene_summary_mageck))


###combine both tables:
gene_summary_all <- full_join(gene_summary_bowtie, gene_summary_mageck
, by = c("id"), suffix = c(".bowtie", ".mageck"))



###find fdr and lfc for bowtie and mageck:
gene_summary_all$LFC.bowtie <- gene_summary_all$pos_lfc.bowtie
gene_summary_all$FDR.bowtie <- gene_summary_all$pos_fdr.bowtie

gene_summary_all$LFC.bowtie[abs(gene_summary_all$neg_lfc.bowtie) > gene_summary_all$pos_lfc.bowtie] <-
  gene_summary_all$neg_lfc.bowtie[abs(gene_summary_all$neg_lfc.bowtie) > gene_summary_all$pos_lfc.bowtie]
gene_summary_all$FDR.bowtie[abs(gene_summary_all$neg_lfc.bowtie) > gene_summary_all$pos_lfc.bowtie] <- 
  gene_summary_all$neg_fdr.bowtie[abs(gene_summary_all$neg_lfc.bowtie) > gene_summary_all$pos_lfc.bowtie]

gene_summary_all$LFC.mageck <- gene_summary_all$pos_lfc.mageck
gene_summary_all$FDR.mageck <- gene_summary_all$pos_fdr.mageck

gene_summary_all$LFC.mageck[abs(gene_summary_all$neg_lfc.mageck) > gene_summary_all$pos_lfc.mageck] <-
  gene_summary_all$neg_lfc.mageck[abs(gene_summary_all$neg_lfc.mageck) > gene_summary_all$pos_lfc.mageck]
gene_summary_all$FDR.mageck[abs(gene_summary_all$neg_lfc.mageck) > gene_summary_all$pos_lfc.mageck] <- 
  gene_summary_all$neg_fdr.mageck[abs(gene_summary_all$neg_lfc.mageck) > gene_summary_all$pos_lfc.mageck]


###find the more favourable fdr and it´s corresponding lfc for each region
 gene_summary_all$FDR = gene_summary_all$FDR.mageck
 gene_summary_all$FDR[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck] = 
  gene_summary_all$FDR.bowtie[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck]


 gene_summary_all$LFC = gene_summary_all$LFC.mageck
 gene_summary_all$LFC[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck] = 
 gene_summary_all$LFC.bowtie[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck]

print("subsetting")
J_gene_all <- subset(gene_summary_all, gene_summary_all$FDR <= fdr_cutoff & gene_summary_all$LFC <= lfc_cutoff, 
                              select = !(grepl("^pos", colnames(gene_summary_all))))


colors_venn <- c('#d9f0d3',"#CCEBC5" )


 venn_plot <- function(venn_diagramm_data, venn_diagram_title, colors_for_venn){
   venn(venn_diagramm_data, zcolor = colors_for_venn,  opacity = 0.6, ilcs = 3.5, sncs = 2.5,
      box = FALSE, ggplot = TRUE, ellipse = TRUE) +
   labs(title = venn_diagram_title) +
   theme(plot.title = element_text(size=34)) +
   theme(legend.title = element_blank()) +
   theme(plot.title = element_text(hjust = 0.5))
 }

venn_diagramm_data <- list(MAGeCK = J_gene_all$id[J_gene_all$FDR.mageck < fdr_cutoff],
                           Bowtie2 = J_gene_all$id[J_gene_all$FDR.bowtie < fdr_cutoff]
                           )

venn_diagram_title <- paste0("Regions lead to depletion (", length(J_gene_all$id), ")")

venn_category_names <- paste0(names(venn_diagramm_data)[1],"(", length(venn_diagramm_data[[1]]), ")")
for(i in 2:length(names(venn_diagramm_data))) {
  venn_category_names <- append(venn_category_names, 
                                paste0(names(venn_diagramm_data)[i],"(", length(venn_diagramm_data[[i]]), ")"))
}

names(venn_diagramm_data) <- venn_category_names

plot_to_png(venn_region_depletion, outputfolder)
venn_plot(venn_diagramm_data, venn_diagram_title, colors_venn) #[1:length(venn_diagramm_data)])
dev.off()


############ read annotation file 
cho_annotation <- readGff3(input_annotation_file, quiet = TRUE)
cho_annotation_table <- cbind(annotation(cho_annotation),
                              cho_annotation,
                              getGffAttribute(cho_annotation, 
                                              c("ID","Dbxref", "Name", "gbkey", "gene")))

####df mit fenster ...
sgRNA_pairs <- read_new_file(path_to_input_file_genetable_with_names) 
colnames(sgRNA_pairs) <- c("all", "seq", "genewindow")
sgRNA_pairs$start_g1 <- str_split_i(sgRNA_pairs$all, "_", 5)
sgRNA_pairs$start_g1 <- as.numeric(str_split_i(sgRNA_pairs$start_g1, "-", 1))
sgRNA_pairs$end_g2 <- str_split_i(sgRNA_pairs$all, "_", 9)
sgRNA_pairs$end_g2 <- as.numeric(str_split_i(sgRNA_pairs$end_g2, "-", 2))
sgRNA_pairs <- drop_na(sgRNA_pairs)
chromosome_names <- as.data.frame(str_split(sgRNA_pairs$all, "_", simplify = TRUE)[,6:7])
sgRNA_pairs$chromosome <- str_c(chromosome_names[,1], chromosome_names[,2], sep = "_")
sgRNA_pairs$chromosome <- str_c(sgRNA_pairs$chromosome, ".1")


###window for overlap
genewindows <- unique(sgRNA_pairs$genewindow)
windowbeginn <- 1
windowende <- 1

for (gw in genewindows) {
  windowbeginn <- append(windowbeginn, min(sgRNA_pairs$start_g1[sgRNA_pairs$genewindow == gw]))
  windowende <- append(windowende, max(sgRNA_pairs$end_g2[sgRNA_pairs$genewindow == gw]))
} 

genewindows <- data.frame(genewindow = genewindows, 
                          windowbeginn = windowbeginn[-1],
                          windowende = windowende[-1])


for (gw in genewindows$genewindow){
  genewindows$chromosome[genewindows$genewindow == gw] <- 
    unique(sgRNA_pairs$chromosome[sgRNA_pairs$genewindow == gw])
}

genewindows <- select(genewindows, chromosome, windowbeginn, windowende, genewindow)

###gen lookup for windows:
to_select <- cho_annotation_table$gbkey %in% grep("Gene", cho_annotation_table$gbkey, value =TRUE)
DF1 <- subset(cho_annotation_table, to_select, select = c('seq_name','1', '2','ID'))
colnames(DF1) <- c("ID","start", "ende", "name") 
DF1$ID <-  as.character(DF1$ID)
DF1$name <- gsub("^gene-","",DF1$name)


DF2 <- genewindows
colnames(DF2) <- c("ID", "start", "ende", "annot")


DTgene <- data.table(DF1) 
DTwindow <- data.table(DF2) 
setkey(DTgene, ID, start, ende)
genes_in_windows <- foverlaps(DTwindow, DTgene, type="any")
genes_in_windows$windowsize <- genes_in_windows$i.ende - genes_in_windows$i.start
genes_in_windows$essential <- rep("Gene", length(genes_in_windows$ID))


genes_in_windows_wide <- genes_in_windows[,c("name", "annot", "essential")]
genes_in_windows_wide <- setDT(genes_in_windows_wide)
genes_in_windows_wide <- nest(genes_in_windows_wide, data = c(name))
genes_in_windows_wide$data <- sapply(genes_in_windows_wide$data, function(x) 
  paste(unlist(x, use.names = FALSE), collapse = ", ") )
genes_in_windows_wide <- dcast(genes_in_windows_wide, annot ~ essential, value.var = "data")
colnames(genes_in_windows_wide) <- c("id", "Genes" )


gene_summary_all <- full_join(gene_summary_all, genes_in_windows_wide)
J_gene_all <- left_join(J_gene_all, genes_in_windows_wide)
colnames(genewindows) <- c("chromosome", "windowbegin", "windowend", "id")

J_results_table <- left_join(J_gene_all, genewindows)
J_results_table <- subset(J_results_table, select = c(chromosome, windowbegin, windowend, Genes) )

everything_for_fede <- left_join(gene_summary_all, genewindows)

summary_all_tables <- everything_for_fede

print(paste0(outputfolder, substitute(summary_all_tables), ".txt"))

results_table <- J_results_table

save_table_to_csv(genes_in_windows, outputfolder, FALSE, TRUE)
save_table_to_csv(gene_summary_all, outputfolder, FALSE, TRUE)
save_table_to_csv(summary_all_tables, outputfolder, FALSE, TRUE)
save_table_to_csv(results_table, outputfolder, FALSE, TRUE)

gene_all <- J_gene_all
save_table_to_csv(gene_all, outputfolder, FALSE, TRUE)

