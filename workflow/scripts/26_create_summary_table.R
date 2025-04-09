
folder_args <- commandArgs(trailingOnly = TRUE)

input_results_table <- folder_args[1]
input_genes_chrom_states <- folder_args[2]

input_sgrna_summary_bowtie <- folder_args[3]
input_sgrna_summary_mageck <- folder_args[4]

input_guides_design <- folder_args[5]

output_compare_bowtie_mageck <- folder_args[6]
output_summary_table_overview <- folder_args[7]
output_results_table_conclusion <- folder_args[8]
output_sgrna_summary_conclusion <- folder_args[9]


library(data.table)
library(venn)
library(ggplot2)
library(ggpolypath)

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
    newfile <- as.data.frame(newfile, stringsAsFactors= FALSE)
    colnames(newfile) <- spaltenbeschriftung_newfile
    return(newfile)
  }

convert_df_to_numeric <- 
  function(count_table, skip_columns = c(1,2)) {
    count_table_colnames <- colnames(count_table)
    gene_names <- as.data.frame(count_table[,skip_columns], stringsAsFactors= FALSE)
    s_ <- sapply(count_table[,-skip_columns],as.numeric)
    data_data_frame <- data.frame(gene_names, s_, stringsAsFactors= FALSE)
    colnames(data_data_frame) <- count_table_colnames
    return(data_data_frame)
  }

save_table_to_txt <- 
  function(filename_table_input, outputfolder, c_row_names = FALSE, c_col_names = FALSE) {
    write.table(as.matrix(filename_table_input),
                sep = "\t", file = paste0(outputfolder, substitute(filename_table_input), ".txt"),
                quote = FALSE, row.names = c_row_names, col.names = c_col_names
    )
  }


gene_summary_all_df <- read_new_file( input_results_table )
genes_chrom_states_df <- read_new_file( input_genes_chrom_states )
 
setDT(gene_summary_all_df)
colnames(gene_summary_all_df)
setnames(gene_summary_all_df, "classification", "Classification")
 
setDT(genes_chrom_states_df)
colnames(genes_chrom_states_df)
setnames(genes_chrom_states_df,
         c("classification"),
         c("Classification")
         )

new_allmighty_table <- merge(gene_summary_all_df, genes_chrom_states_df,
                            all = TRUE)
str(new_allmighty_table)

new_allmighty_table[,unknown:=NULL]
columns_numeric <- c("LFC", "FDR", 
                     "Region_Start", "Region_End",
                     "num.bowtie", "neg_score.bowtie",
                     "neg_p_value.bowtie", "neg_fdr.bowtie",
                     "neg_rank.bowtie", "neg_goodsgrna.bowtie",
                     "neg_lfc.bowtie", "pos_score.bowtie",
                     "pos_p_value.bowtie", "pos_fdr.bowtie",
                     "pos_rank.bowtie", "pos_goodsgrna.bowtie",
                     "pos_lfc.bowtie", "num.mageck", "neg_score.mageck",
                     "neg_p_value.mageck", "neg_fdr.mageck",
                     "neg_rank.mageck", "neg_goodsgrna.mageck",
                     "neg_lfc.mageck", "pos_score.mageck",
                     "pos_p_value.mageck", "pos_fdr.mageck",
                     "pos_rank.mageck", "pos_goodsgrna.mageck",
                     "pos_lfc.mageck", "LFC.bowtie", "FDR.bowtie",
                     "LFC.mageck", "FDR.mageck",
                     "Active Promoter", "Flanking TSS downstream",
                     "Polycomb repressed regions", "Quiescent/low",
                     "Strong Transcription", "Weak enhancer",
                     "Repressed heterochromatin", "Weak genic enhancer",
                     "Active Enhancer 1", "Flanking TSS upstream",
                     "Active Enhancer 2"         
                     )

for(col in columns_numeric)
  set(new_allmighty_table, j = col, value = as.numeric(new_allmighty_table[[col]]))

sapply(new_allmighty_table, class)
str(new_allmighty_table)
colnames(new_allmighty_table)

new_allmighty_table <- as.data.frame(new_allmighty_table)


###find goodsgrna and score for bowtie and mageck:
new_allmighty_table$goodsgrna.bowtie <- new_allmighty_table$pos_goodsgrna.bowtie
new_allmighty_table$score.bowtie <- new_allmighty_table$pos_score.bowtie
new_allmighty_table$goodsgrna.bowtie[abs(new_allmighty_table$neg_lfc.bowtie) > new_allmighty_table$pos_lfc.bowtie] <-
  new_allmighty_table$neg_goodsgrna.bowtie[abs(new_allmighty_table$neg_lfc.bowtie) > new_allmighty_table$pos_lfc.bowtie]
new_allmighty_table$score.bowtie[abs(new_allmighty_table$neg_lfc.bowtie) > new_allmighty_table$pos_lfc.bowtie] <- 
  new_allmighty_table$neg_score.bowtie[abs(new_allmighty_table$neg_lfc.bowtie) > new_allmighty_table$pos_lfc.bowtie]

new_allmighty_table$goodsgrna.mageck <- new_allmighty_table$pos_goodsgrna.mageck
new_allmighty_table$score.mageck <- new_allmighty_table$pos_score.mageck
new_allmighty_table$goodsgrna.mageck[abs(new_allmighty_table$neg_lfc.mageck) > new_allmighty_table$pos_lfc.mageck] <-
  new_allmighty_table$neg_goodsgrna.mageck[abs(new_allmighty_table$neg_lfc.mageck) > new_allmighty_table$pos_lfc.mageck]
new_allmighty_table$score.mageck[abs(new_allmighty_table$neg_lfc.mageck) > new_allmighty_table$pos_lfc.mageck] <- 
  new_allmighty_table$neg_score.mageck[abs(new_allmighty_table$neg_lfc.mageck) > new_allmighty_table$pos_lfc.mageck]


###find goodsrna and score for each region
new_allmighty_table$goodsgrna <- new_allmighty_table$goodsgrna.bowtie
new_allmighty_table$score <- new_allmighty_table$score.bowtie
new_allmighty_table$num <- new_allmighty_table$num.bowtie

new_allmighty_table$goodsgrna[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)] <- 
    new_allmighty_table$goodsgrna.mageck[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)]
new_allmighty_table$score[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)] <-
    new_allmighty_table$score.mageck[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)]
new_allmighty_table$num[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)] <-
    new_allmighty_table$num.mageck[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)]


setDT(new_allmighty_table)

results_all_regions <- new_allmighty_table[,.(ID, Classification, Genes,
                                              Chromosome, Region_Start, Region_End,
                                              LFC, FDR, Score=score, Num=num, GoodsgRNA=goodsgrna,
                                              `Active Promoter`, `Flanking TSS downstream`,
                                              `Polycomb repressed regions`, `Quiescent/low`,
                                              `Strong Transcription`, `Weak enhancer`,
                                              `Repressed heterochromatin`, `Weak genic enhancer`,
                                              `Active Enhancer 1`, `Flanking TSS upstream`,
                                              `Active Enhancer 2`  
                                              )]


write.table(as.matrix(results_all_regions),
                sep = "\t", file = output_results_table_conclusion,
                quote = FALSE, row.names = FALSE, col.names = TRUE
    )


######designed guides
guides_designed <- matrix(scan(file = input_guides_design ,
          skip = 0,
          what = character(),
          sep = "\t"
        ),
        ncol = 3,
        byrow = TRUE
      )
guides_designed <- as.data.frame(guides_designed, stringsAsFactors= FALSE)
setDT(guides_designed)

guides_designed[like(V1, "-w"), .N]


###########sgRNAsummary
sgrna_bowtie_df <- read_new_file( input_sgrna_summary_bowtie )
sgrna_mageck_df <- read_new_file( input_sgrna_summary_mageck )
 
setDT(sgrna_bowtie_df)
colnames(sgrna_bowtie_df)
setnames(sgrna_bowtie_df, "Gene", "ID")
str(sgrna_bowtie_df)

setDT(sgrna_mageck_df)
colnames(sgrna_mageck_df)
setnames(sgrna_mageck_df, "Gene", "ID")
str(sgrna_mageck_df)

###find the corresponding sgrna for each region
new_allmighty_table$origin <- "bowtie"
new_allmighty_table$origin[abs(new_allmighty_table$LFC.mageck) > abs(new_allmighty_table$LFC.bowtie)] <- "mageck" 

key_table <- new_allmighty_table[,.(ID, origin),]

sgrna_bowtie_df <- merge(sgrna_bowtie_df, key_table, all = TRUE)
sgrna_mageck_df <- merge(sgrna_mageck_df, key_table, all = TRUE)

sgrna_bowtie_df[,.N,]
sgrna_mageck_df[,.N,]

sgrna_summary <- rbind(sgrna_bowtie_df[origin=="bowtie",,],
                       sgrna_mageck_df[origin=="mageck",,]
                       )

#setDT(sgrna_summary)
sgrna_summary[, origin:=NULL,]
sgrna_summary[order(ID),,]
setnames(sgrna_summary, "ID", "Gene")

write.table(as.matrix(sgrna_summary),
                sep = "\t", file = output_sgrna_summary_conclusion,
                quote = FALSE, row.names = FALSE, col.names = TRUE
    )



############overview

overview <- rbind(c("total regions", new_allmighty_table[,.N]),
      c("essential & FDR < 25%", new_allmighty_table[Classification == "essential" & FDR < 0.25, .N,]),
      c("essential & FDR < 10%", new_allmighty_table[Classification == "essential" & FDR < 0.10, .N,]),
      c("essential & FDR < 5%", new_allmighty_table[Classification == "essential" & FDR < 0.05, .N,]),
      c("non essential regions", new_allmighty_table[Classification == "non_essential", .N,]),
      c("essential regions with annotation", 
        new_allmighty_table[Classification == "essential" & FDR < 0.25 & !is.na(Genes), .N,]),
      c("essential regions without annotation", 
        new_allmighty_table[Classification == "essential" & FDR < 0.25 & is.na(Genes), .N,]),
      c("essential regions without annotation & with chromatin state Quiescent/low 100% ", 
        new_allmighty_table[Classification == "essential" & FDR < 0.25 & `Quiescent/low` == 100, .N,]),
      c("guides designed",
        guides_designed[like(V1, "-w"), .N]),
      c("guides analysed (in ctrl above treshold)",
        sgrna_summary[like(Gene, "-w"), .N])
      )

write.table(as.matrix(overview),
                sep = "\t", file = output_summary_table_overview,
                quote = FALSE, row.names = FALSE, col.names = FALSE
    )



####graphics

colors_venn <- c("#E41A1C" , "#7fbf7b")

mlist <- new_allmighty_table[Classification == "essential" & FDR.mageck < 0.25,(ID),]
blist <- new_allmighty_table[Classification == "essential" & FDR.bowtie < 0.25,(ID),]

venn_diagramm_data <- list(MAGeCK = mlist,
                           Bowtie2 = blist
                           )  


venn_diagram_title <- paste0("Regions lead to depletion (", new_allmighty_table[Classification == "essential", .N,], ")")

venn_category_names <- paste0(names(venn_diagramm_data)[1],
                          "(", 
                          new_allmighty_table[Classification == "essential" & FDR.mageck < 0.25,.N,] ,
                           ")"
                           )
for(i in 2:length(names(venn_diagramm_data))) {
  venn_category_names <- append(venn_category_names, 
                                paste0(names(venn_diagramm_data)[i],"(", 
                                new_allmighty_table[Classification == "essential" & FDR.bowtie < 0.25,.N,], ")"))
}

names(venn_diagramm_data) <- venn_category_names


png(file = output_compare_bowtie_mageck, width = 1200, height = 740)

venn(venn_diagramm_data
    , ilabels="counts"
    , zcolor = colors_venn
    , opacity = 0.6
    , ilcs = 3.5
    , sncs = 2.5
    , box = FALSE, ggplot = TRUE, ellipse = TRUE) +
labs(title = venn_diagram_title) +
theme(plot.title = element_text(size=34)) +
theme(legend.title = element_blank()) +
theme(plot.title = element_text(hjust = 0.5))

dev.off()

