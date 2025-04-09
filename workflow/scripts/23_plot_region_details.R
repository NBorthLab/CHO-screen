folder_args <- commandArgs(trailingOnly = TRUE)

file_expressed_genes <- folder_args[1]
file_not_expressed_genes <- folder_args[2]
file_chromatin_states <- folder_args[3]
file_transcripts <- folder_args[4]
file_essentiality <- folder_args[5]
file_gtf <- folder_args[6]
input_selected_regions <- folder_args[7]
outputfile <- folder_args[8]
file_txdb <- folder_args[9]


library(karyoploteR)
library(GenomicFeatures)
library(dplyr)
library(stringr)

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

read_new_file_b <- function(path_to_file) {
  spaltenbeschriftung_newfile <-
    scan(
      file = path_to_file,
      what = character(),
      nlines = 1,
      sep = " "
    )
  newfile <-
    matrix(
      scan(
        file = path_to_file ,
        skip = 1,
        what = character(),
        sep = " "
      ),
      ncol = length(spaltenbeschriftung_newfile),
      byrow = TRUE
    )
  newfile <- as.data.frame(newfile, stringsAsFactors = FALSE)
  colnames(newfile) <- spaltenbeschriftung_newfile
  return(newfile)
}

change_values <- function(x, old_values, new_values) {
  old_chr_n <- setNames(new_values, old_values)
  return(unname(old_chr_n[x]))
}


colour_expressed <- "#009966"
colour_not_expressed <- "#990033"
colour_ignored <- "grey20"
colour_transcript <- "#000033"
filename_pdf_legende <- paste0(str_split_1(outputfile, ".pdf")[1], "_legend.pdf")


### read_regions_to plot:
selected_regions <- read_new_file(input_selected_regions)
selected_regions <- as.vector(selected_regions$id)

###table for expressed genes
data_expressed_genes <- scan(file_expressed_genes,
                             what = character(),
                             sep = "\t")

data_not_expressed_genes <- scan(file_not_expressed_genes,
                             what = character(),
                             sep = "\t")


data_essentiality <- read_new_file(file_essentiality)


data_essentiality <- subset(data_essentiality, select = c('ID', 'Chromosome', 'Region_Start', 'Region_End', 'classification', 'Genes', 'Quiescent/low'))

data_essentiality <- dplyr::select(data_essentiality, Chromosome, Region_Start, Region_End, everything())

region_essential <- subset(data_essentiality, classification == "essential" )
region_essential <- toGRanges(region_essential)


old_chr <- c( "chr1_0", "chr1_1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
              "chr8", "chr9", "chr10", "chrX", "unplaced_scaffold_1", "unplaced_scaffold_2")
new_chr <- c("NW_023276806.1", "NW_023276807.1", "NC_048595.1",
             "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1",
             "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1",
             "NC_048604.1", "NW_023276808.1", "NW_023276809.1")
plot_chr  <- c( "Chr 1_0", "Chr 1_1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr 7",
              "Chr 8", "Chr 9", "Chr 10", "Chr X", "unplaced scaffold 1", "unplaced scaffold 2")



data_transcript_15 <- read_new_file(file_transcripts)
data_transcripts_2 <- subset(data_transcript_15, select = c("GeneID",
                                    "Chr", "Start", "End", "Strand",
                                    "Length" ,"geo_mean", "any_transcripts"))


data_transcripts <- as.data.frame(data_transcripts_2, stringsAsFactors = FALSE)
data_transcripts <- subset(data_transcripts, select = c(GeneID, Chr, Start, End, geo_mean))
data_transcripts <- dplyr::select(data_transcripts, Chr, Start, End, everything())
data_transcripts$Chr <- factor(data_transcripts$Chr, levels = new_chr)

data_transcripts$Start <- as.numeric(data_transcripts$Start)
data_transcripts$End <- as.numeric(data_transcripts$End)

data_transcripts <- data_transcripts[order(data_transcripts$Chr, data_transcripts$Start), ]

data_transcripts$y <- as.numeric(data_transcripts$geo_mean)
data_transcripts$y_norm <- (data_transcripts$y - min(data_transcripts$y))/
  (max(data_transcripts$y) - min(data_transcripts$y))

data_transcripts$y0 <- min(data_transcripts$y_norm)
data_transcripts <- toGRanges(data_transcripts)


# Group by Chromosom
grouped_df <- group_by(data_essentiality, Chromosome)

##find min and max of start and end for each region
new_df <- as.data.frame(summarize(grouped_df,
                    Start = min(as.double(Region_Start)),
                    End = max(as.double(Region_End)))
)

genome.picrh <- toGRanges(new_df) ### not picrh more like lean_cho_genome
rm(new_df, grouped_df)

### data for annotation:
###check if txdb file exist and create if not:
if (!file.exists(file_txdb)){
  annotation_picrh_txdb <- makeTxDbFromGFF(file_gtf)
  saveDb(annotation_picrh_txdb, file_txdb)
}

##load txdb file for annotation:
annotation_picrh_txdb <- loadDb(file_txdb)


# chromatin_states as custom_cytoband:

###table for numbers of transcripts:
data_chromatin <- matrix(scan(file_chromatin_states,
                              what = character(), sep = "\t", skip = 1),
                         ncol = 9,byrow = TRUE)
data_chromatin <- as.data.frame(data_chromatin, stringsAsFactors=FALSE)
data_chromatin <- dplyr::select(data_chromatin, V1, V2, V3, V4, V9)
colnames(data_chromatin) <- c("Chr", "Start", "End", "name", "colour")
data_chromatin$gieStain <- data_chromatin$name
data_chromatin$Chr <- change_values(data_chromatin$Chr, old_chr, new_chr)

data_chromatin <- toGRanges(na.omit(data_chromatin))

### to set chromatin specific colours
chromatin_states_names <- unique(data_chromatin$name)
chromatin_states_colours <- unique(data_chromatin$colour)
colours_chromatin <- setNames(chromatin_states_colours, chromatin_states_names)

pdf(file = filename_pdf_legende, 
    paper = "a4",
    width = 11, 
    pointsize = 12)

### create legende:
legende_chromatin <- c("Repressed heterochromatin "
                       , "Quiescent/low"
                       , "Polycomb repressed regions"
                       , "Strong Transcription"
                       , "Weak genic enhancer"
                       , "Active Enhancer 1"
                       , "Active Enhancer 2"
                       , "Weak enhancer"
                       , "Active Promoter"
                       , "Flanking TSS upstream"
                       ,"Flanking TSS downstream"
)
clrs <- colours_chromatin[legende_chromatin]
ltype <- 1

legende_genes_names <- c("expressed"
                         , "not expressed"
                         , "indeterminate"
)

legende_genes_colours <- c(colour_expressed
                          , colour_not_expressed
                          , colour_ignored
                          )

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1,
     cex.main = 1.7,
     main = "Detail plots - legende")
legend("topleft",
       title="Chromatin states:",
       title.adj = c(0),
       text.font = 1,
       legend = legende_chromatin,
       col = clrs,
       fill = clrs,
       cex=1.3,
       bty='n' 
       )

legend("bottomleft",
       title="Gene status\n according to Riedl et al.:",
       title.adj = c(0),
       text.font = 1,
       legend = legende_genes_names,
       col = legende_genes_colours,
       fill = legende_genes_colours,
       cex=1.3,
       bty='n' 
       )
dev.off()


print("start plot function")

###Plot function###
plot_region <- function(zoom_region){
#prepare plot parameters
plot_params <- getDefaultPlotParams(plot.type = 2)
plot_params$data1height <- 80
plot_params$data1inmargin <- 40
plot_params$ideogramheight <- 80
plot_params$data2inmargin <- 40
plot_params$data2height <- 250
plot_params$leftmargin <- 0.15
plot_params$rightmargin <- 0.05

legende_y <- 1
legende_x <- 0.67
legende_cex <- 1  # text size for legend
legende_lwd <- 5 # Line width
label_margin = 0.04 # margin label to graph


#creates plot and chromatin-states:
kp <- plotKaryotype(genome = genome.picrh,
                    plot.params = plot_params,
                    cytobands = data_chromatin,
                     chromosomes = ifelse(exists("zoom_region"),
                                          as.character(zoom_region@seqnames@values[1]),
                                          "all"),
                    zoom = if(exists("zoom_region")) zoom_region else NULL,
                    plot.type = 2,
                    labels.plotter = if(exists("zoom_region")) NULL else kpAddChromosomeNames
                    )
if(exists("zoom_region")) kpAddLabels(kp, labels = paste0("chromatin\nstates"),
                                      data.panel = "ideogram",
                                      cex=legende_cex, label.margin = label_margin)
kpAddCytobands(kp, color.table = colours_chromatin)
kpAddBaseNumbers(kp, tick.dist = if(exists("zoom_region")) 50000 else 1000000)


title(main = if(exists("zoom_region"))
                   paste0(change_values(unlist(str_split(as.character(zoom_region), ":"))[1],
                                 new_chr,plot_chr)
                          , ": "
                          , unlist(str_split(as.character(zoom_region), ":"))[2])
                   else paste0("CHO-Genome"),
      adj = 0, line = 3, cex.main = legende_cex*1.5, font.main = 1)

###annotation:
#create genes data for annotations
genes_data <- makeGenesDataFromTxDb(txdb = annotation_picrh_txdb, karyoplot = kp)
#trim gene data to zoom region, to fix text-center-issue:
genes_data$genes <- restrict(genes_data$genes, start(zoom_region), end(zoom_region),
                             keep.all.ranges = TRUE)
genes_annotation_color <- ifelse(genes_data$genes$gene_id %in% data_expressed_genes, colour_expressed,
                               ifelse(genes_data$genes$gene_id %in% data_not_expressed_genes, 
                               colour_not_expressed, colour_ignored))

#tn <- 0
kpDataBackground(kp, data.panel = 1,
)
kpPlotGenes(kp, data=genes_data,
            data.panel = 1,
            plot.transcripts = FALSE,
            add.gene.names = TRUE,
            gene.col = genes_annotation_color,
          
            gene.name.cex = legende_cex * 0.65, gene.name.position = "top",
            clipping = FALSE, #print outside the plot, necessary for last gene name
            r0=0.1, r1=0.8,
            srt = 90  # turns text 90degree
            ) 
kpAddLabels(kp, labels = paste0("annotation"),
            data.panel = 1,
            srt=0,
            cex=legende_cex, label.margin = label_margin)


 tick_pos_text <- c(-4.5, seq(-3,7,1), 8.9)


tick_pos_scaliert <- 
  (tick_pos_text - min(tick_pos_text)) / (max(tick_pos_text) - min(tick_pos_text))

    
kpAxis( kp
        , r1 = 0.2
        , r0 = 0.8
        , col = "gray50"
        , cex = 0.5
        , tick.pos = tick_pos_scaliert
        , labels = tick_pos_text
        , data.panel = 2
)


    kpAbline(kp, h=c(0 - min(tick_pos_text)/(max(tick_pos_text)-min(tick_pos_text)))
         , col="gray60"
         , r1 = 0.2
         , r0 = 0.8
         , data.panel = 2
         , lty = 2)

kpBars(kp
       , data_transcripts
       , y0 = data_transcripts$y0
       , y1 = data_transcripts$y_norm
       , data.panel = 2
       , col = colour_transcript
       , border = colour_transcript
       , r1 = 0.2, r0 = 0.8
)

kpAddLabels(kp, labels = paste0("transcripts"),
            data.panel = 2,
            srt=0,
            cex=legende_cex, label.margin = label_margin)


}


plot_region_pdf <- function(vector_of_regions, filename_pdf) {
###create list of regions to plot:
df_regions_to_plot <- data_essentiality[data_essentiality$ID %in% vector_of_regions,]

###plot regions to file:
df_regions_to_plot$Chromosome <- gsub("\\s", "", df_regions_to_plot$Chromosome)
df_regions_to_plot$Region_Start <- gsub("\\s", "", df_regions_to_plot$Region_Start)
df_regions_to_plot$Region_End <- gsub("\\s", "", df_regions_to_plot$Region_End)

df_regions_to_plot <- rbind(df_regions_to_plot)

df_regions_to_plot$regions_to_plot <- with(df_regions_to_plot, paste0(Chromosome, ":", Region_Start, "-", Region_End))


###changed for debugg
pdf(file = paste0(filename_pdf), paper= "a4r", width = 11, pointsize = 12)

x <- 0
for (i in df_regions_to_plot$regions_to_plot) {
  x <- x+1
  print(x)
      zoom_region <- toGRanges(i)
      plot_region(toGRanges(i))
}

dev.off()

}



plot_region_pdf(selected_regions, outputfile)

