
folder_args <- commandArgs(trailingOnly = TRUE)

gene_mageck <- folder_args[1]
gene_bowtie <- folder_args[2]
input_annotation_file <- folder_args[3]
path_to_input_file_genetable_with_names <- folder_args[4]
run_name <- folder_args[5]
outputfolder <- folder_args[6]
#path_to_input_file_least_essential <- folder_args[7]


# # ### debuggingg:
# setwd("C:/Users/akm/Documents/")
#input_folder_mageck <- "GitHub/lean_CHO/screen/table_mageck/results/windows_rra/"
## input_folder_bowtie <- "GitHub/lean_CHO/screen/table_bowtieM/results/windows_rra/"
# input_file_gene_summary <- "gws_test.gene_summary.txt"
# gene_mageck <- paste0(input_folder_mageck, input_file_gene_summary)
# gene_bowtie <- paste0(input_folder_bowtie, input_file_gene_summary)
# input_file_sgrna_summary <- "gws_test.sgrna_summary.txt"
# outputfolder <- "Github/lean_CHO/screen/"
# project_run_name <- "gws_test"
# path_to_input_file_genetable_with_names <- "GitHub/lean_CHO/final_output_LEAN_GW.xlsx"
# input_annotation_file <- "GitHub/lean_CHO/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff"
# path_to_input_file_least_essential <- "GitHub/lean_CHO/screen/selected_least_essential.txt"



#gene_mageck <- "/data/borth/amolin/leeaan/resources/table_mageck/results/region_rra/screen.gene_summary.txt"
#gene_bowtie <- "/data/borth/amolin/leeaan/resources/table_bowtieM/results/region_rra/screen.gene_summary.txt"
#input_annotation_file <- "/data/borth/amolin/leeaan/resources/raw/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff"
#path_to_input_file_genetable_with_names <- "/data/borth/amolin/leeaan/resources/raw/library_for_mageck_count_original_setting.txt"
#run_name <- "screen"
#outputfolder <- "results/"


#setwd("/data/borth/amolin/lean_cho_screen/")

#lfc_cutoff

## cutoffs:
lfc_cutoff <- log2(0.85)
fdr_cutoff <- 0.25

#library(MAGeCKFlute)
# library(clusterProfiler)
library(ggplot2)
library(venn)
#library("ggplot2")
library(genomeIntervals)
library(dplyr)
library(readxl)
library(stringi)
library(stringr)
#library("GenomicRanges")
library(data.table)
library(tidyr)
#library(ggrepel)


# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("genomeIntervals")



# control_isch <- "#1C6DAB"
# treated_isch <- "#E41A1C"
# time0_isch <- "#FFD92F"

# result_isch_light <- "#d9f0d3"
# result_isch_middle <- "#7fbf7b"
# result_isch_dark <- "#1b7837"

# grey_isch <- "gray80"
# gblack_isch <- "#999999"

colors_venn <- c('#d9f0d3',"#CCEBC5" )

# colors=c( "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
#           "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
#           "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
#           "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#762a83","#af8dc3",'#e7d4e8',
#           '#f7f7f7','#d9f0d3','#7fbf7b','#1b7837')






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

plot_to_png <- function(plot_titel, outputfolder){#, run_name = project_run_name){
  png(file = paste0(outputfolder, run_name, "_", deparse(substitute(plot_titel)),".png"),width = 1200, height = 740)
}


convert_df_to_numeric    <- function(count_table, skip_columns = c(1,2)) {
  ####achtung reihenfolge relevant für beschriftung...noch ändern
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




#########datazeugs:

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



###find the more favourable lfc and it´s corresponding fdr for each region
##gene_summary_all$LFC <- gene_summary_all$LFC.bowtie
##gene_summary_all$FDR <- gene_summary_all$FDR.bowtie

##gene_summary_all$LFC[abs(gene_summary_all$LFC.mageck) > abs(gene_summary_all$LFC.bowtie)] <- 
##    gene_summary_all$LFC.mageck[abs(gene_summary_all$LFC.mageck) > abs(gene_summary_all$LFC.bowtie)]
##gene_summary_all$FDR[abs(gene_summary_all$LFC.mageck) > abs(gene_summary_all$LFC.bowtie)] <-
##    gene_summary_all$FDR.mageck[abs(gene_summary_all$LFC.mageck) > abs(gene_summary_all$LFC.bowtie)]





#####################################################################the new old right one:
## ###find the more favourable fdr and it´s corresponding lfc for each region

 gene_summary_all$FDR = gene_summary_all$FDR.mageck
 gene_summary_all$FDR[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck] = 
  gene_summary_all$FDR.bowtie[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck]


 gene_summary_all$LFC = gene_summary_all$LFC.mageck
 gene_summary_all$LFC[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck] = 
 gene_summary_all$LFC.bowtie[gene_summary_all$FDR.bowtie < gene_summary_all$FDR.mageck]



negative_gene_all <- subset(gene_summary_all, gene_summary_all$FDR <= fdr_cutoff & gene_summary_all$LFC <= lfc_cutoff, 
                              select = !(grepl("^pos", colnames(gene_summary_all))))

str(negative_gene_all)

#str(compare)

#compare <- read_new_file("/data/borth/amolin/leeaan/targets_col.txt")

#summary(as.factor(negative_gene_all$id %in% compare$id))


#library(dplyr)

# Subset the data
#filtered_data <- gene_summary_all %>%
#  filter(FDR <= 0.25, LFC <= log(0.85))

# View the subset
#head(filtered_data)
##awk {if($6 < 0.25 && $5 < log(0.85)) print $1} results/summary_table.txt | wc -l

#str(filtered_data)

##subset(gene_summary_all, gene_summary_all



#str(negative_gene_all)








#head(negative_gene_all)
#str(negative_gene_all)

#summary(as.factor(negative_gene_all_old$id %in% negative_gene_all$id))

#gene_summary_all_old <- read_new_file("/data/borth/amolin/lean_cho_screen/results_original_settings/gene_summary_all.txt")
#negative_gene_all_old <- read_new_file("/data/borth/amolin/lean_cho_screen/results_original_settings/negative_gene_all.txt")

#str(gene_summary_all_old)
#str(negative_gene_all_old)

#fuck_old <- (gene_summary_all_old[
#  as.numeric(gene_summary_all_old$neg_fdr.bowtie) <= 0.25,])
#str(fuck_old)

#fuck_new <- gene_summary_bowtie[gene_summary_bowtie$neg_fdr <= 0.25,]
#str(fuck_new)
#summary(as.factor(fuck_new$id %in% fuck_old$id))









#summary(as.factor(gene_summary_all_old$id 
#%in% gene_summary_all$id))


#plot(gene_summary_bowtie$neg_fdr, gene_summary_all_old$neg_fdr.bowtie)

#colnames(gene_summary_mageck)
#colnames(gene_summary_all_old)

#str(neue_tabelle)


#summary(as.factor(negative_gene_all$id %in% neue_tabelle$id))



# positive_gene_all <- subset(gene_summary_all, gene_summary_all$LFC > lfc_cutoff *-1 & gene_summary_all$FDR <= fdr_cutoff,
#                               select = !(grepl("^neg", colnames(gene_summary_all))))

# topnames <- c("601-w290", "599-w370", "600-w83",
#               "599-w968", "595-w1249", "807-w1644")
#               #positive_gene_all$id

# plot_to_png(volcano_plot, outputfolder)
# p12 = VolcanoView_2(gene_summary_all, x = "LFC", y = "FDR", Label = "id",
#                   x_cutoff = lfc_cutoff *-1, y_cutoff = fdr_cutoff, top = 0,
#                   topnames = topnames,
#                   mycolour = c(grey_isch, result_isch_dark, result_isch_middle),
#                   alpha = 0.6) +
#   labs(title = "Negative and positive selection", xlab = "log2FC", ylab = "-log10(FDR)") +
#   theme(plot.title = element_text(size=26)) +
#   theme(plot.title = element_text(hjust = 0.5))

# plot(p12)
# dev.off()



# ############venn


# control_isch <- "#1C6DAB"
# treated_isch <- "#E41A1C"
# time0_isch <- "#FFD92F"

# result_isch_light <- "#d9f0d3"
# result_isch_middle <- "#7fbf7b"
# result_isch_dark <- "#1b7837"

# grey_isch <- "gray80"
# gblack_isch <- "#999999"

colors_venn <- c('#d9f0d3',"#CCEBC5" )

# colors_venn=c( "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
#           "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
#           "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
#           "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#762a83","#af8dc3",'#e7d4e8',
#           '#f7f7f7','#d9f0d3','#7fbf7b','#1b7837')

 venn_plot <- function(venn_diagramm_data, venn_diagram_title, colors_for_venn){
   venn(venn_diagramm_data, zcolor = colors_for_venn,  opacity = 0.6, ilcs = 3.5, sncs = 2.5,
      box = FALSE, ggplot = TRUE, ellipse = TRUE) +
   labs(title = venn_diagram_title) +
   theme(plot.title = element_text(size=34)) +
   theme(legend.title = element_blank()) +
   theme(plot.title = element_text(hjust = 0.5))
 }

venn_diagramm_data <- list(MAGeCK = negative_gene_all$id[negative_gene_all$FDR.mageck < fdr_cutoff],
                           Bowtie2 = negative_gene_all$id[negative_gene_all$FDR.bowtie < fdr_cutoff]
                           )

venn_diagram_title <- paste0("Regions lead to depletion (", length(negative_gene_all$id), ")")

venn_category_names <- paste0(names(venn_diagramm_data)[1],"(", length(venn_diagramm_data[[1]]), ")")
for(i in 2:length(names(venn_diagramm_data))) {
  venn_category_names <- append(venn_category_names, 
                                paste0(names(venn_diagramm_data)[i],"(", length(venn_diagramm_data[[i]]), ")"))
}

names(venn_diagramm_data) <- venn_category_names

plot_to_png(venn_region_depletion, outputfolder)
venn_plot(venn_diagramm_data, venn_diagram_title, colors_venn) #[1:length(venn_diagramm_data)])
dev.off()



# ###positive gene   positive enhancement
# venn_diagramm_data <- list(MAGeCK = positive_gene_all$id[positive_gene_all$FDR.mageck < fdr_cutoff],
#                            Bowtie2 = positive_gene_all$id[positive_gene_all$FDR.bowtie < fdr_cutoff])

# venn_diagram_title <- paste0("Regions lead to enrichment (", length(positive_gene_all$id), ")")

# venn_category_names <- paste0(names(venn_diagramm_data)[1],"(", length(venn_diagramm_data[[1]]), ")")
# for(i in 2:length(names(venn_diagramm_data))) {
#   venn_category_names <- append(venn_category_names, 
#                                 paste0(names(venn_diagramm_data)[i],"(", length(venn_diagramm_data[[i]]), ")"))
# }

# names(venn_diagramm_data) <- venn_category_names

# plot_to_png(venn_region_enrichment, outputfolder)
# venn_plot(venn_diagramm_data, venn_diagram_title, colors[1:length(venn_diagramm_data)])
# dev.off()



############ read annotatien file 

cho_annotation <- readGff3(input_annotation_file, quiet = TRUE)
cho_annotation_table <- cbind(annotation(cho_annotation),
                              cho_annotation,
                              getGffAttribute(cho_annotation, 
                                              c("ID","Dbxref", "Name", "gbkey", "gene")))

####df mit fenster ...
# sgRNA_pairs <- read_excel(path_to_input_file_genetable_with_names, 
#                           sheet = "final_output_perfect_pairs_gRNA")
# sgRNA_pairs$guide_id <- paste0(sgRNA_pairs$seqnames...2,":",sgRNA_pairs$start...3,"-",sgRNA_pairs$end...4,"_",
#                                sgRNA_pairs$seqnames...9,":",sgRNA_pairs$start...10,"-",sgRNA_pairs$end...11)

# sgRNA_pairs <- subset(sgRNA_pairs, select=c(seqnames...2,start...3, end...4, Window...5, start...10, end...11, guide_id))

# colnames(sgRNA_pairs) <- c("chromosome","start_g1","end_g1","window","start_g2","end_g2","guide_id")

# sgRNA_pairs$window <- gsub("Window", "w", sgRNA_pairs$window)

# summary(as.factor(sgRNA_pairs$chromosome))
# sgRNA_pairs$genewindow <- gsub("\\.1", "", sgRNA_pairs$chromosome)
# sgRNA_pairs$genewindow <- gsub("^NC_048", "", sgRNA_pairs$genewindow)
# sgRNA_pairs$genewindow <- gsub("^NW_023276", "", sgRNA_pairs$genewindow)
# sgRNA_pairs$genewindow <- paste0(sgRNA_pairs$genewindow,"-",sgRNA_pairs$window)

#neu ab da...................0601
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

str(sgRNA_pairs)


#head(chrom)
#is.matrix(chrom)
#str_split_i(string, pattern, i)
head(sgRNA_pairs)
#summary(is.na(sgRNA_pairs$start_g2))

# sgRNA_pairs$guide_id <- paste0(sgRNA_pairs$seqnames...2,":",sgRNA_pairs$start...3,"-",sgRNA_pairs$end...4,"_",
#                                sgRNA_pairs$seqnames...9,":",sgRNA_pairs$start...10,"-",sgRNA_pairs$end...11)
# sgRNA_pairs <- subset(sgRNA_pairs, select=c(seqnames...2,start...3, end...4, Window...5, start...10, end...11, guide_id))
# colnames(sgRNA_pairs) <- c("chromosome","start_g1","end_g1","window","start_g2","end_g2","guide_id")
# sgRNA_pairs$window <- gsub("Window", "w", sgRNA_pairs$window)
# summary(as.factor(sgRNA_pairs$chromosome))
# sgRNA_pairs$genewindow <- gsub("\\.1", "", sgRNA_pairs$chromosome)
# sgRNA_pairs$genewindow <- gsub("^NC_048", "", sgRNA_pairs$genewindow)
# sgRNA_pairs$genewindow <- gsub("^NW_023276", "", sgRNA_pairs$genewindow)
# sgRNA_pairs$genewindow <- paste0(sgRNA_pairs$genewindow,"-",sgRNA_pairs$window)




###frueester beginn bis spaetestes ende als window with overlap
genewindows <- unique(sgRNA_pairs$genewindow)
windowbeginn <- 1
windowende <- 1
#gw <- "807-w987"

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

head(genewindows)

genewindows <- select(genewindows, chromosome, windowbeginn, windowende, genewindow)

#rm(sgRNA_pairs, cho_annotation)

###checken welches gen in welchem window ist...

to_select <- cho_annotation_table$gbkey %in% grep("Gene", cho_annotation_table$gbkey, value =TRUE)
DF1 <- subset(cho_annotation_table, to_select, select = c('seq_name','1', '2','ID'))
colnames(DF1) <- c("ID","start", "ende", "name") 
DF1$ID <-  as.character(DF1$ID)
DF1$name <- gsub("^gene-","",DF1$name)


DF2 <- genewindows
colnames(DF2) <- c("ID", "start", "ende", "annot")

str(DF2)
str(DF1)

DTgene <- data.table(DF1) #y
DTwindow <- data.table(DF2) #x
setkey(DTgene, ID, start, ende)
genes_in_windows <- foverlaps(DTwindow, DTgene, type="any")#, which=TRUE)
genes_in_windows$windowsize <- genes_in_windows$i.ende - genes_in_windows$i.start
genes_in_windows$essential <- rep("Gene", length(genes_in_windows$ID))

#ifelse(tolower(genes_in_windows$name) %in% tolower(essential_genes$`Gene Name`), "essential", 
#                                      ifelse(tolower(genes_in_windows$name) %in% tolower(non_essential_genes$`Gene Name`),
#                                             "non_essential", 
#                                             ifelse(tolower(genes_in_windows$name) %in% tolower(silenced_genes$`Gene Name`), "silenced", 
#                                                    "unknown")))


genes_in_windows_wide <- genes_in_windows[,c("name", "annot", "essential")] #4,5,9)]
genes_in_windows_wide <- setDT(genes_in_windows_wide)
genes_in_windows_wide <- nest(genes_in_windows_wide, data = c(name))
genes_in_windows_wide$data <- sapply(genes_in_windows_wide$data, function(x) 
  paste(unlist(x, use.names = FALSE), collapse = ", ") )
genes_in_windows_wide <- dcast(genes_in_windows_wide, annot ~ essential, value.var = "data")
colnames(genes_in_windows_wide) <- c("id", "Genes" )

str(genes_in_windows_wide)

head(genes_in_windows_wide)

gene_summary_all <- full_join(gene_summary_all, genes_in_windows_wide)#, by = join_by(c("id","annot")))
negative_gene_all <- left_join(negative_gene_all, genes_in_windows_wide)
# positive_gene_all <- left_join(positive_gene_all, genes_in_windows_wide)

head(genewindows)
#colnames(genewindows) <- c("id", "windowbegin", "windowend", "chromosome")
colnames(genewindows) <- c("chromosome", "windowbegin", "windowend", "id")

# positive_results_table <- left_join(positive_gene_all, genewindows)
# positive_results_table <- subset(positive_results_table, select = c(chromosome, windowbegin, windowend, Genes) )

str(negative_gene_all)

negative_results_table <- left_join(negative_gene_all, genewindows)
negative_results_table <- subset(negative_results_table, select = c(chromosome, windowbegin, windowend, Genes) )

str(negative_gene_all)

everything_for_fede <- left_join(gene_summary_all, genewindows)

#str(gene_summary_all)

# Subset the data
#filtered_data <- everything_for_fede %>%
#  filter(FDR <= 0.25, LFC <= log2(0.85))

# View the subset
#head(filtered_data)
#awk '{if($6 < 0.25 && $5 < log(0.85)) print $1}' results/summary_table.txt | wc -l

#str(filtered_data)

##subset(gene_summary_all, gene_summary_all







#str(genewindows)


# ### update averything for fede with classification:

# selected_least_essential <- read_new_file(path_to_input_file_least_essential)

# everything_for_fede$classification <- rep("undefined", length(everything_for_fede$id))
# everything_for_fede$classification <- ifelse(everything_for_fede$id %in% negative_gene_all$id, "essential",
#                                              ifelse(everything_for_fede$id %in% selected_least_essential$id, "non_essential", "undefined"))

# everything_for_fede$classification[everything_for_fede$id %in% positive_gene_all$id] <- "growth_enhancing"

# summary(as.factor(everything_for_fede$classification))

# everything_for_fede <- select(everything_for_fede, id, chromosome, windowbegin, windowend, LFC, FDR, classification, Genes, everything())
# colnames(everything_for_fede) <- c("ID", "Chromosome", "Region_Start", "Region_End", "LFC", "FDR", "classification", "Genes",
#                                      "num.bowtie", "neg_score.bowtie", "neg_p_value.bowtie", "neg_fdr.bowtie",
#                                      "neg_rank.bowtie", "neg_goodsgrna.bowtie", "neg_lfc.bowtie", "pos_score.bowtie",
#                                      "pos_p_value.bowtie", "pos_fdr.bowtie", "pos_rank.bowtie", "pos_goodsgrna.bowtie",
#                                      "pos_lfc.bowtie", "table.bowtie", "num.mageck", "neg_score.mageck", "neg_p_value.mageck",
#                                      "neg_fdr.mageck", "neg_rank.mageck", "neg_goodsgrna.mageck", "neg_lfc.mageck",
#                                      "pos_score.mageck", "pos_p_value.mageck", "pos_fdr.mageck", "pos_rank.mageck",
#                                      "pos_goodsgrna.mageck", "pos_lfc.mageck", "table.mageck", "LFC.bowtie", "FDR.bowtie",
#                                      "LFC.mageck", "FDR.mageck")          

summary_all_tables <- everything_for_fede

head(summary_all_tables)

print(paste0(outputfolder, substitute(summary_all_tables), ".txt"))


# save_table_to_csv(positive_results_table, outputfolder, FALSE, TRUE)


# save_table_to_csv(positive_gene_all, outputfolder, FALSE, TRUE)





save_table_to_csv(genes_in_windows, outputfolder, FALSE, TRUE)
save_table_to_csv(gene_summary_all, outputfolder, FALSE, TRUE)
save_table_to_csv(summary_all_tables, outputfolder, FALSE, TRUE)
save_table_to_csv(negative_results_table, outputfolder, FALSE, TRUE)


gene_all <- negative_gene_all
save_table_to_csv(gene_all, outputfolder, FALSE, TRUE)







