
folder_args <- commandArgs(trailingOnly = TRUE)

input_file_essential <- folder_args[1]
outputfolder <- folder_args[2]

### set your files:
#setwd("C:/Users/MOA/Documents/")
#input_file_annotation_caterina <- "GitHub/lean_cho_go/data/lncrna.gtf"
#input_file_annotation_picrh <- "amolin/lean_cho_screen/resources/raw/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff"
#input_file_growth_enhancing <- "amolin/lean_cho_screen/resources/raw/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff"
#input_file_all <- "amolin/lean_cho_screen/results/gene_summary_all.txt"

#input_file_essential <- "negative_results_table.txt"
#input_file_all_genes <- "gene_summary_all.txt"


# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("orthogene")
# BiocManager::install("AnnotationHub")



library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(orthogene)
#library(purrr)
#library(tidyverse)
#library(clusterProfiler)
library(AnnotationHub)
#library(ggpubr)
#library(RColorBrewer)
#library(pathviewr)






plot_to_png <- function(plot_titel, outputfolder, run_name = project_run_name){
png(file = paste0(outputfolder, run_name, "_", deparse(substitute(plot_titel)),".png"),width = 1200, height = 740)
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


convert_to_mouse <- function(hamster_genes) {
  #' Convert hamster to mouse genes
  #'
  #' Converts the gene identifiers of the Chinese hamster to the orthologous
  #' gene of the house mouse.
  #'
  #' @param hamster_genes A character vector of hamster gene symbols
  #' @return A character vector of orthologous mouse genes.
  #'
  #' @examples
  #' convert_orthologs(c("Gapdh", "Fut8"))
  convert_orthologs(hamster_genes,
                    input_species = "Cricetulus griseus",
                    output_species = "Mus musculus",
                    method = "gprofiler",
                    non121_strategy = "drop_both_species"
  ) %>%
    rownames()
}



get_dotplot_for_go <- function(ont_xx, number_categories){
 ego <- enrichGO(gene = essential_genes_converted, 
                universe = all_genes_converted,
                # background genes. If missing, the all genes listed in
                # the database (eg TERM2GENE table) will be used as background.
                keyType = "SYMBOL", 
                OrgDb = org.Mm.eg.db,
                ont = ont_xx,      #"BP","MF","CC"
                pAdjustMethod = "none",
                pvalueCutoff  = 0.1,
                #minGSSize = 0,
                maxGSSize = 300000  # maximal size of genes annotated for testing
                #qvalueCutoff  = 0.05,
                #readable      = TRUE
)
 dotplot(ego, showCategory=number_categories) + 
    ggtitle( ont_xx)
}
 


essential_regions_df <- read_new_file(input_file_essential)
essential_genes_vec <- unique(unlist(strsplit(essential_regions_df$Genes, ", ")))
essential_genes_vec <- essential_genes_vec[!is.na(essential_genes_vec)]
essential_genes_converted <- convert_to_mouse(essential_genes_vec) 

all_regions_df <- read_new_file(input_file_all_genes)
all_genes_vec <- unique(unlist(strsplit(all_regions_df$Genes, ", ")))
all_genes_vec <- all_genes_vec[!is.na(all_genes_vec)]
all_genes_converted <- convert_to_mouse(all_genes_vec)


head(ego)
#dotplot(ego, showCategory=10) + 
#  ggtitle("all")

ego_simple <- clusterProfiler::simplify(ego)
head(ego_simple)
#dotplot(ego_simple, showCategory=10) + 
#  ggtitle("simple")


genes_mouse_all <- AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db), columns=c("SYMBOL", "GENENAME"))
# dim(genes)
# keys(org.Mm.eg.db)
# columns(org.Mm.eg.db)
# keytypes(org.Mm.eg.db)


#png(file = paste0(outputfolder, "plot_go_portrait_BP.png"),width = 440, height = 740)
#get_dotplot_for_go("BP", 8)
#dev.off()

#png(file = paste0(outputfolder, "plot_go_portrait_MF.png"),width = 440, height = 740)
#get_dotplot_for_go("MF", 8)
#dev.off()

#png(file = paste0(outputfolder, "plot_go_portrait_CC.png"),width = 440, height = 740)
#get_dotplot_for_go("CC", 8)
#dev.off()

png(file = paste0(outputfolder, "go_term_no_good.png"),width = 440, height = 740)
dev.off()











summary(as.factor(essential_genes_converted %in% genes_mouse_all$SYMBOL))
summary(as.factor(all_genes_converted %in% genes_mouse_all$SYMBOL))





annotation_snapshot <- "2023-10-21"
#annotation_id <- "AH114610" ## no record
annotation_id <- "AH117698" # org.Cricetulus_griseus.eg.sqlite


annotation_hub <- AnnotationHub()
snapshotDate(annotation_hub) <- annotation_snapshot


annot_db_Cg <- AnnotationHub::query(annotation_hub, "Cricetulus griseus")[[annotation_id]]

#converting essential genes to ENTREZID
universe_Crigri <- bitr(all_genes_vec, fromType = "SYMBOL", toType = "ENTREZID", annot_db_Cg)$ENTREZID
essential_Crigri <- bitr(essential_genes_vec, fromType = "SYMBOL", toType = "ENTREZID", annot_db_Cg)$ENTREZID


ego_annohub <- enrichGO(gene = essential_Crigri, 
#                universe = universe_Crigri,
                # background genes. If missing, the all genes listed in
                # the database (eg TERM2GENE table) will be used as background.
                keyType = "ENTREZID", 
                OrgDb = annot_db_Cg,
                ont = "all",      #"BP","MF","CC"
                pAdjustMethod = "none",
                pvalueCutoff  = 0.1,
                #minGSSize = 0,
                maxGSSize = 300000  # maximal size of genes annotated for testing
                #qvalueCutoff  = 0.05,
                #readable      = TRUE
)


head(ego_annohub)
#dotplot(ego_annohub, showCategory=10) + 
#  ggtitle("all_annohub")
#   if (ont_xx == "BP") {paste0("Biological Process")},
#   if (ont_xx == "MF") {paste0("Molecular Function")},
#   if (ont_xx == "CC") {paste0("Cellular Compartment")})

ego_annohub_simple <- clusterProfiler::simplify(ego_annohub)
head(ego_annohub_simple)
#dotplot(ego_annohub_simple, showCategory=10) + 
#  ggtitle("simple_annohub")


columns(annot_db_Cg)
genes_annotationhub_all <- AnnotationDbi::select(annot_db_Cg, keys = keys(annot_db_Cg), columns=c("SYMBOL", "GENENAME", "ENTREZID"))

summary(as.factor(essential_Crigri %in% genes_annotationhub_all$ENTREZID))
summary(as.factor(universe_Crigri %in% genes_annotationhub_all$ENTREZID))







