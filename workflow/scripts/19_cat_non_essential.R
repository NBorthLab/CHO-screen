folder_args <- commandArgs(trailingOnly = TRUE)

input_folder <- folder_args[1]
input_file_counttable <- folder_args[2]
input_file_sgrna_summary <- folder_args[3]
path_to_output_folder <- folder_args[4]
path_to_config_file <- folder_args[5]
input_folder_essential <- folder_args[6]
input_file_gene_target <- folder_args[7]
input_file_gene_all <- folder_args[8]
outputfolder <- folder_args[9]
project_run_name <- folder_args[10]

target_gene <- input_file_gene_target
all_gene <-input_file_gene_all

## for GENES:
cutoff_fdr <- 0.25

## for GUIDES:
cutoff_padj <- 0.01                     
cutoff_baseMean <- 1                    
cutoff_log2FoldChange <- -log2(85/100)
cutoff_lfcSE <- 0.5              
nr_lost_guides <- 1              

suppressPackageStartupMessages({
library(venn)
library(ggplot2)
library(genomeIntervals)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(DESeq2)
library(plyr)
library(ggrepel)
})

control_isch <- "burlywood3"  #"#1C6DAB"
treated_isch <- "aquamarine3" #"#E41A1C" yesssss
time0_isch <- "#FFD92F"

result_isch_light <- "lightcyan3"  #"#d9f0d3"
result_isch_middle <- "seashell"  #"#7fbf7b"
result_isch_dark <- "cyan4"  #"#1b7837"

black_isch <- "#999999"

name_for_treated <- "treated"
name_for_ctrl <- "control"
 

# functions:
print("start functions")

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
newfile <- as.data.frame(newfile, stringsAsFactors = FALSE)
colnames(newfile) <- spaltenbeschriftung_newfile
return(newfile)
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


plot_to_png <- function(plot_titel, outputfolder, run_name = project_run_name){
  png(file = paste0(outputfolder, run_name, "_", deparse(substitute(plot_titel)),".png"),width = 1200, height = 740)
}


convert_df_to_numeric    <- function(count_table, skip_columns = c(1,2)) {
  count_table_colnames <- colnames(count_table)
  gene_names <- as.data.frame(count_table[,skip_columns], stringsAsFactors = FALSE)
  s_ <- sapply(count_table[,-skip_columns],as.numeric)
  data_data_frame <- data.frame(gene_names, s_, stringsAsFactors = FALSE)
  colnames(data_data_frame) <- count_table_colnames
  return(data_data_frame)
}

convert_df_to_integer    <- function(count_table, skip_columns = c(1,2)) {
  count_table_colnames <- colnames(count_table)
  gene_names <- as.data.frame(count_table[,skip_columns], stringsAsFactors = FALSE)
  s_ <- sapply(count_table[,-skip_columns],as.integer)
  data_data_frame <- data.frame(gene_names, s_, stringsAsFactors = FALSE)
  colnames(data_data_frame) <- count_table_colnames
  return(data_data_frame)
}

save_table_to_csv <-  function(filename_table_input, outputfolder, c_row_names = FALSE, c_col_names = FALSE) {
  write.table(as.matrix(filename_table_input),
              sep = "\t", file = paste0(outputfolder, substitute(filename_table_input), ".txt"),
              quote = FALSE, row.names = c_row_names, col.names = c_col_names
  )
}

find_least_essential <- function(least_essential, cutoff_padj, cutoff_baseMean, cutoff_log2FoldChange, cutoff_lfcSE) {
  
  if (cutoff_padj != 0) {
    least_essential <- least_essential[padj > cutoff_padj,,]  
  }
  
  if (cutoff_baseMean != 0) {
    least_essential <- least_essential[baseMean > cutoff_baseMean,,]  
  }
  
  if (cutoff_log2FoldChange != 0) {
    least_essential <- least_essential[log2FoldChange > -cutoff_log2FoldChange,,]  
  }
  
  if (cutoff_lfcSE != 0) {
    least_essential <- least_essential[lfcSE < cutoff_lfcSE,,]  
  }
  
  least_essential <-  left_join(least_essential, least_essential[, .N, by = Gene])
  
  setnames(least_essential, "N", "nr_ne_guides")
  
  return(least_essential)
}


######### volcano plot to show regions:
VolcanoView_3 <- function (df, x = "logFC", y = "adj.P.Val", id = "id", selectvector = "",  Label = NA, top = 5, 
                           topnames = NULL, x_cutoff = log2(1.5), y_cutoff = 0.05, 
                           mycolour = c("grey90", "grey90", "aquamarine3", "cyan4"),
                           #colours: "no", "up", "down", "select"
                           alpha = 0.8, force = 0.1, main = NULL,
                           xlab = bquote(~log[2] ~ "(FC)"), ylab = bquote(~ - ~ log[10] ~ "(FDR)"),
                           filename = NULL, width = 4, 
                           height = 2.5, ...) {
  requireNamespace("ggplot2", quietly = TRUE) || stop("need ggplot2 package")
  requireNamespace("ggrepel", quietly = TRUE) || stop("need ggrepel package")
  gg = df[, c(x, y, id)]
  gg$group = "no"
  gg$group[gg[, x] > x_cutoff & gg[, y] < y_cutoff] = "up"
  gg$group[gg[, x] < -x_cutoff & gg[, y] < y_cutoff] = "down"
  gg$group[gg[, id] %in% selectvector ] = "select"
  gg[, y] = -log10(gg[, y])
  if (!(top == 0 & is.null(topnames))) {
    gg$Label = rownames(gg)
    if (!is.na(Label)) 
      gg$Label = df[, Label]
    gg = gg[order(gg[, y], abs(gg[, x]), decreasing = TRUE), 
    ]
    idx1 = idx2 = c()
    if (top > 0) {
      idx1 = which(gg$group == "up")[1:min(top, sum(gg$group == 
                                                      "up"))]
      idx2 = which(gg$group == "down")[1:min(top, sum(gg$group == 
                                                        "down"))]
    }
    idx = unique(c(idx1, idx2, which(gg$Label %in% topnames)))
    gg$Label = as.character(gg$Label)
    gg$Label[setdiff(1:nrow(gg), idx)] = ""
  }
  gg$color = gg$group
  gg$color[gg$Label != ""] = "gray90"
  mycolour = c(mycolour)
  names(mycolour) = c("no", "up", "down", "select")
  
  levels(gg$group) <- c( "no", "up" , "down", "select" )
  p = ggplot(gg, aes_string(x = x, y = y, label = "Label"))
  p = p + geom_point(aes_string(fill = "group", colour = "color"), shape = 21, 
                     alpha = alpha, size = 4)
  p = p + geom_point(data = subset(gg, group == 'select'), colour = mycolour["select"],
                     alpha = (alpha - 0.1), size = 4)
  p = p + scale_color_manual(values = mycolour)
  p = p + scale_fill_manual(values = mycolour)
  p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
  p = p + geom_vline(xintercept = c(-x_cutoff), linetype = "dotted")
  p = p + xlim(min(gg[, x]) - 0.001, max(gg[, x]) + 0.001)
  p = p + labs(x = xlab, y = ylab, title = main)
  p = p + theme_bw(base_size = 22)
  p = p + theme(
                plot.title = element_text(hjust = 0.5),
                legend.position = "none",
                panel.grid = element_blank()  #removes gridlines
  )
  if (!(top == 0 & is.null(topnames))) {
    p = p + ggrepel::geom_text_repel(size = 6, segment.color = "grey90", segment.size = 0.3)
  }
  if (!is.null(filename)) {
    ggsave(plot = p, filename = filename, width = width, 
           height = height, units = "in", ...)
  }
  return(p)
}




#########data:
print("start_reading")
count_table <- read_new_file(input_file_counttable)
config_table <- scan_config_file(path_to_config_file)

samplenumbers_control <- read_from_config_table("mageck_test_d0_control: ")
samplenumbers_treated <- read_from_config_table("mageck_test_d0_treated: ")

count_table <- convert_df_to_integer(count_table, c(1,2))

# get the table of read counts
readcounts <- count_table[,-2]
row.names(readcounts) <- readcounts$sgRNA
readcounts <- readcounts[,-1]
readcounts <- select(readcounts, ctrl_1, ctrl_2, ctrl_3, ctrl_4, everything())

readmeta <- data.frame(project = rep("gws", length(names(readcounts))),
                       samplenumber = gsub(" ", "", names(readcounts)),
                       row.names = names(readcounts))

head(readmeta)

readmeta[,"treatment"] <- ifelse(readmeta$samplenumber %in% samplenumbers_treated, "treated", "control")
  
head(readmeta)

for (i in c("project", "treatment")) { 
  readmeta[,i] <- as.factor(readmeta[,i])  
}


####################settings###################

Auswertungsname <- "least_essential"
deseq_design <- ~treatment
first_order <- "control"


####generate_deseq_subset_object

rc <- readcounts 
rm <- readmeta

deseq.ds <- DESeqDataSetFromMatrix(countData = rc,
                                   colData = rm,
                                   design = deseq_design)

# setting the first - level - factor
colData(deseq.ds)$treatment <- relevel(colData(deseq.ds)$treatment, first_order)


#running dge-analysis
deseq.ds <- DESeq(deseq.ds)
resultsNames(deseq.ds)
res <- results(deseq.ds, contrast = c("treatment", "treated", "control"))

res_df <- as.data.frame(res, stringsAsFactors = FALSE)
res_df$sgRNA <- rownames(res_df)

count_table <- left_join(count_table, res_df)

ct <- as.data.table(count_table, stringsAsFactors = FALSE)
ct <- left_join(ct, ct[, .N, by = Gene])

print("set new names")
setnames(ct, "N", "nr_total_guides")


plot_to_png(histogramm_lfcSE_least_essential, outputfolder,project_run_name)
hist(ct$lfcSE[ct$log2FoldChange > cutoff_log2FoldChange & ct$padj > cutoff_padj ], breaks = 100,
     xlab = "lfc-SE" , ylab = " Frequency " , main = "Frequencies of standard errors")
dev.off()


selected_least_essential <- find_least_essential(ct, cutoff_padj, cutoff_baseMean, cutoff_log2FoldChange, cutoff_lfcSE)

setnames(selected_least_essential, "Gene", "id")

selected_least_essential_2 <-  selected_least_essential[nr_ne_guides >= nr_total_guides - 1,,]


#### comparison least essential and most essential

target_gene <- read_new_file(target_gene)
target_gene <- select(target_gene, id, Genes , everything())
target_gene <- convert_df_to_numeric(target_gene, c(1:2))

any(unique(selected_least_essential[nr_ne_guides >= nr_total_guides - nr_lost_guides , id,]) %in% unique(target_gene$id))

all_gene_summary <- read_new_file(all_gene)
all_gene_summary <- select(all_gene_summary, id, Genes, everything())
all_gene_summary <- convert_df_to_numeric(all_gene_summary, c(1:2))

res <- results(deseq.ds, contrast = c("treatment", "treated", "control"))

plot_to_png(log2fc_cutoff_least_essential, outputfolder, project_run_name)
plotMA(res , alpha = 0.05, main = "treatment_treated_vs_control" ,
       ylim = c ( -2, 2)) +
  abline(h = -cutoff_log2FoldChange, col = "red") 
dev.off()


## nonessential
region_data <- all_gene_summary
region_data <- select(region_data, id, Genes, everything())
region_data <- convert_df_to_numeric(region_data, c(1:2))
region_data_t <- as.data.table(region_data)


##combine both
testvector_a <- region_data_t[LFC > cutoff_log2FoldChange , id,]
testvector_b <- selected_least_essential[nr_ne_guides >= nr_total_guides - nr_lost_guides , id,]

testvector <- unique(c(testvector_a, testvector_b))  
png(file = paste0(outputfolder, project_run_name, "_region_data_.png"),width = 1200, height = 740)
p12 = VolcanoView_3(region_data, x = "LFC", y = "FDR", id = "id", testvector, Label = "id",
                  x_cutoff = cutoff_log2FoldChange, y_cutoff = cutoff_fdr, top = 0,
                  alpha = 0.6) +
  labs(title = "selection") +
  theme(plot.title = element_text(size = 26)) +
  theme(plot.title = element_text(hjust = 0.5))

plot(p12)
dev.off()

select_regions <- region_data$id %in% testvector

selected_least_essential <- region_data[select_regions,,]

save_table_to_csv(selected_least_essential,outputfolder, FALSE, TRUE)
