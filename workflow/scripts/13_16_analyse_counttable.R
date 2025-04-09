folder_args <- commandArgs(trailingOnly = TRUE)

path_to_input_file_counttable <- folder_args[1]
path_to_input_file_summary <- folder_args[2]
path_to_config_file <- folder_args[3]
path_to_output_folder <- folder_args[4]
set_for_mageck <- folder_args[5]
path_to_input_file_counttable_unfiltered <- folder_args[6]

  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("cluster")
  library("ggrepel")
  library("stringr")
  library("pheatmap")

#############definition of some colours & variables:

control_isch <- "burlywood3"
treated_isch <- "aquamarine3"
time0_isch <- "#FFD92F"

result_isch_light <- "lightcyan3" 
result_isch_middle <- "seashell"
result_isch_dark <- "cyan4"

black_isch <- "#999999"

name_for_treated <- "treated"
name_for_ctrl <- "ctrl"

colors=c( "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#762a83","#af8dc3",'#e7d4e8',
          '#f7f7f7','#d9f0d3','#7fbf7b','#1b7837')


# functions:
rownames_from_txt <- function(name_of_file) {
  return(scan(
    file = name_of_file,
    what = character(),
    flush = TRUE,
    sep = "\t"
  ))
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
  newfile <- as.data.frame(newfile, stringsAsFactors= FALSE)
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
  gene_names <- as.data.frame(count_table[,skip_columns], stringsAsFactors= FALSE)
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



#correlation plot:
cols = colorRampPalette(c(treated_isch, result_isch_middle))(31)
pal = colorRampPalette(cols)
cor_colors = data.frame(correlation = seq(0.7,1,0.01), 
                        correlation_color = pal(31)[1:31])  # assigns a color for each r correlation value
cor_colors$correlation_color = as.character(cor_colors$correlation_color)


panel.cor <- function(x, y, digits=2, prefix = "", cex.cor) {
par(usr = c(0, 1, 0, 1))
u <- par('usr') 
names(u) <- c("xleft", "xright", "ybottom", "ytop")
r <- abs(cor(x, y))
bgcolor = cor_colors[2+(-r+1)*100,2]    # converts correlation into a specific color
do.call(rect, c(col = bgcolor, as.list(u))) # colors the correlation box
txt <- format(c(r, 0.123456789), digits = digits)[1]
txt <- paste0(prefix, txt)
if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex.cor * r, col = ifelse(r > 0.5, 1, 2))
}

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = result_isch_light)
}


plot_samples_against_samples <- function(counttable_, samplenumbers){
  samples_to_analyse <- c("Gene", "sgRNA",samplenumbers)
  samplenames <- colnames(counttable_)
  samplenames_to_select <- ifelse(samplenames %in% samples_to_analyse, TRUE, FALSE)
  subset_for_plot <- subset(counttable_, select =  samplenames_to_select)
  
  if (ncol(subset_for_plot)>3){
    pairs(subset_for_plot[3:ncol(subset_for_plot)], 
          lower.panel = panel.smooth, 
          diag.panel = panel.hist, 
          upper.panel = panel.cor, 
          labels = gsub("_", " ", colnames(subset_for_plot[3:ncol(subset_for_plot)])))
  } else {
    pairs(data.frame(subset_for_plot[3],subset_for_plot[3]), 
          lower.panel = panel.smooth,
          diag.panel = panel.hist,
          upper.panel = panel.cor, 
          labels = gsub("_", " ", colnames(subset_for_plot[3:ncol(subset_for_plot)])))
  }
}


plot_to_png <- function(plot_titel, outputfolder, run_name = project_run_name){
png(file = paste0(outputfolder, run_name, "_", deparse(substitute(plot_titel)),".png"),width = 1200, height = 740)
}


## to set colored leaf labels to hierarchical clustering:
colLab <<- function(n) {
  if(is.leaf(n)) {
    a <- attributes(n)
    line=match(attributes(n)$label,h_clust_data[,1])
    treatment=h_clust_data[line,2];
    if(treatment==name_for_ctrl){col_treatment=control_isch};
    if(treatment==name_for_treated){col_treatment=treated_isch}
    attr(n, "nodePar") <- c(a$nodePar, list(cex=3, # dot size
                                            lab.cex=1,
                                            pch=20, # dot style
                                            col=col_treatment,lab.col=col_treatment,
                                            lab.font=1, # type style
                                            lab.cex=1))
  }
  return(n)
}


plot_mapping_ratio <- function (countSummary, Label = "Label", Reads = "Reads", Mapped = "Mapped", 
                                filename = NULL, width = 5, height = 4, ...) {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" is required. Please install it.", 
         call. = FALSE)
  }
  gg = data.frame(Label = rep(countSummary[, Label], 2),read = rep(
    countSummary[,Reads], 2), count = c(countSummary[, Mapped], countSummary[,Reads]- countSummary[, Mapped]),
    category = factor(rep(c("mapped", "unmapped"), c(nrow(countSummary), nrow(countSummary))),
                      levels = c("unmapped", "mapped")))
  gg$percent = paste0(round(gg$count * 100/gg$read, 1), "%")
  gg$pos = ceiling(gg$count/2)
  gg$pos[gg$category == "unmapped"] = ceiling(gg$pos
                                              [gg$category =="unmapped"] + gg$pos[gg$category == "mapped"] * 2)
  fill = c(result_isch_light, result_isch_dark)
  p <- ggplot(gg)
  p = p + geom_bar(aes_string(y = "count", x = "Label", fill = "category"), 
                   stat = "identity", width = 0.8, alpha = 0.9)
  p = p + geom_text(aes_string(x = "Label", y = "pos", label = "percent"), 
                    size = 7)                                                    # size of "%" within bars
  p = p + labs(x = NULL, y = "Reads", title = "Mapping ratio")
  p = p + scale_y_continuous(expand = c(0, 0))
  p = p + scale_fill_manual(values = fill)
  p = p + theme(legend.title = element_blank())
  p = p + theme_bw(base_size = 22)                                               # size of title, axistitle, aso.
  p = p + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(filename)) {
    ggsave(plot = p, filename = filename, units = "in", width = width, 
           height = height, ...)
  }
  return(p)
}


####data analyses####
count_table <- read_new_file(path_to_input_file_counttable)
count_table <- convert_df_to_numeric(count_table,c(1,2))

config_file <- scan(path_to_config_file, sep = ";", what = character())

project_run_name <- read_from_config_file("run_name: ")

if (set_for_mageck == TRUE) {
samplenumbers_control <- read_from_config_file("mageck_test_mageck_control: ")
samplenumbers_treated <- read_from_config_file("mageck_test_mageck_treated: ") 
samplenumbers_day_0_as_control <- read_from_config_file("mageck_test_mageck_day_0_as_control: ") 
} else {
samplenumbers_control <- read_from_config_file("mageck_test_d0_control: ")
samplenumbers_treated <- read_from_config_file("mageck_test_d0_treated: ")
samplenumbers_day_0_as_control <- read_from_config_file("mageck_test_day_0_as_control: ") 
}


###count table for analyses
samples_to_analyse <- c("Gene", "sgRNA",samplenumbers_day_0_as_control, samplenumbers_treated, samplenumbers_control)
samplenames <- colnames(count_table)
samplenames_to_select <- ifelse(samplenames %in% samples_to_analyse, TRUE, FALSE)
count_table <- subset(count_table, select =  samplenames_to_select)
rm(samples_to_analyse, samplenames, samplenames_to_select)


### normalisation to count per million:
cpm_count_table <- matrix(nrow = nrow(count_table), ncol = ncol(count_table))
for (i in (1:ncol(count_table))){
  if (is.character(count_table[2,i])){
    cpm_count_table[,i] <- count_table[,i]
  } else {
    cpm_count_table[,i] <- count_table[,i]/sum(count_table[,i])*10^6
  }
}
cpm_count_table <- as.data.frame(cpm_count_table)
colnames(cpm_count_table) <- colnames(count_table)

cpm_count_table <- convert_df_to_numeric(cpm_count_table, c(1,2))

cpm_count_table_long <- pivot_longer(cpm_count_table[,c(1,3:ncol(cpm_count_table))], !sgRNA, names_to =
                                       "samplename", values_to = "counts")
cpm_count_table_log <- cbind(cpm_count_table[,1:2], log2(as.matrix(cpm_count_table[,c(-1,-2)])+1))
cpm_count_table_long_log <- pivot_longer(cpm_count_table_log[,c(1,3:ncol(cpm_count_table_log))], !sgRNA, names_to =
                                           "samplename", values_to = "log2_counts")

cpm_count_table_long_log$treatment <- ifelse(cpm_count_table_long_log$samplename %in% samplenumbers_control, "ctrl",
                                             ifelse(cpm_count_table_long_log$samplename
                                                    %in% samplenumbers_treated, "treated", "time0"))
cpm_count_table_long$treatment <- ifelse(cpm_count_table_long$samplename %in% samplenumbers_control, "ctrl",
                                         ifelse(cpm_count_table_long$samplename 
                                                %in% samplenumbers_treated, "treated", "time0"))



###mapped_vs_unmapped_reads:
  summary_table <- read_new_file(path_to_input_file_summary)
  Mapped <- as.numeric(summary_table$Mapped)
  Reads <- as.numeric(summary_table$Reads)
  Label <- summary_table$Label

mapping_ratio_data <- data.frame(Label = Label,
                     Reads = Reads,
                     Mapped = Mapped)


png(file = paste0(path_to_output_folder, project_run_name, "_mapping_ratio.png"),width = 1200, height = 740)
plot_mapping_ratio(mapping_ratio_data)
dev.off()


#####BOXPLOTS:

if (ncol(cpm_count_table) > 2) {

plot_to_png(boxplot_cpm, path_to_output_folder)
cpm_count_table_long_log %>%
ggplot(aes(x=samplename, y=log2_counts, fill=treatment)) +
  geom_boxplot() +
  theme(legend.position="right", plot.title = element_text(size=14)) +
  labs(title = "Boxplot CPM-normalised") +
  theme(legend.title = element_blank()) +
  theme_bw(base_size = 22) +                                               # size of title, axistitle, etc
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c(time0_isch, treated_isch, control_isch ),  # manipulate the legend
                    name="Experimental\nCondition",
                    breaks=c( "time0", name_for_treated, name_for_ctrl),
                    labels=c("Plasmid", name_for_treated,name_for_ctrl ))
}
dev.off()


####violinplot:
if (ncol(cpm_count_table) > 2) {
plot_to_png(violinplot, path_to_output_folder )
ggplot(cpm_count_table_long_log,aes(x=samplename,y=log2_counts, 
                                    color=treatment)) +
    geom_violin(linewidth = 1.5, fill = "grey90" )+
    scale_color_manual(values = c(treated_isch, control_isch),
                      name="Experimental\nCondition",
                      breaks=c( name_for_treated, name_for_ctrl ),
                      labels=c( name_for_treated, name_for_ctrl ))+
    geom_boxplot(width=0.1) +
    theme(legend.position="right", plot.title = element_text(size=14)) +
    labs(title = "Violinplot CPM-normalised") +
    theme(legend.title = element_blank()) +
    theme_bw(base_size = 22) +                           # size of title, axistitle, etc
    theme(plot.title = element_text(hjust = 0.5)) 
}
dev.off()


###histogramm all
if (ncol(cpm_count_table) > 2) {
plot_to_png(histogramm, path_to_output_folder )
ggplot(cpm_count_table_long_log, aes(x = log2_counts)) +
  geom_histogram(binwidth=0.3, alpha = 0.8, position = "identity") +
  ylab("frequenzy")
}
dev.off()


#####distribution_of_counts unfiltered  
  
  count_table_unfiltered <- read_new_file(path_to_input_file_counttable_unfiltered)
  count_table_unfiltered <- convert_df_to_numeric(count_table_unfiltered,c(1,2))
  
  cpm_count_table_unfiltered <- matrix(nrow = nrow(count_table_unfiltered), ncol = ncol(count_table_unfiltered))
  for (i in (1:ncol(count_table_unfiltered))){
    if (is.character(count_table_unfiltered[2,i])){
      cpm_count_table_unfiltered[,i] <- count_table_unfiltered[,i]
    } else {
      cpm_count_table_unfiltered[,i] <- count_table_unfiltered[,i]/sum(count_table_unfiltered[,i])*10^6
    }
  }
  cpm_count_table_unfiltered <- as.data.frame(cpm_count_table_unfiltered)
  colnames(cpm_count_table_unfiltered) <- colnames(count_table_unfiltered)
  
  cpm_count_table_unfiltered <- convert_df_to_numeric(cpm_count_table_unfiltered, c(1,2))
  tabsmat=as.matrix(log2(cpm_count_table_unfiltered[,c(-1,-2)]+1))
  colnames(tabsmat)=colnames(cpm_count_table_unfiltered)[c(-1,-2)]

  samplecol= ifelse(dimnames(tabsmat)[[2]] %in% samplenumbers_control, control_isch,
                    ifelse(dimnames(tabsmat)[[2]]
                           %in% samplenumbers_treated, treated_isch, time0_isch))
  tgz=hist(tabsmat,breaks = 40, plot=F)  # setting breaks for histogram

  plot_to_png(distribution_of_counts_unfiltered, path_to_output_folder )
  if(ncol(tabsmat)>=1){
    histlist=lapply(1:ncol(tabsmat),function(X){ return (hist(tabsmat[,X],plot=F,breaks=tgz$breaks)) })
    xrange=range(unlist(lapply(histlist,function(X){X$mids})))
    yrange=range(unlist(lapply(histlist,function(X){X$counts})))
    hst1=histlist[[1]]
    plot(hst1$mids,hst1$counts,type='b',pch=20,xlim=c(0,xrange[2]*1.1),ylim=c(0,yrange[2]*1.2),
         xlab=bquote(~log[2] ~ "(Counts)"),ylab='Frequency',main='Distribution of read counts bevor filter',
         col = ifelse(dimnames(tabsmat)[[2]][1] %in% samplenumbers_control, control_isch,
                      ifelse(dimnames(tabsmat)[[2]][1] 
                             %in% samplenumbers_treated, treated_isch, time0_isch)), 
         cex.main = 2, cex.axis = 1.5, cex.lab = 1.7)
  }
  
  if(ncol(tabsmat)>=2){ 
    for(i in 2:ncol(tabsmat)){
      hstn=histlist[[i]]
      lines(hstn$mids,hstn$counts,type='b',pch=20,
            col=ifelse(dimnames(tabsmat)[[2]][i] %in% samplenumbers_control, control_isch,
                       ifelse(dimnames(tabsmat)[[2]][i] 
                              %in% samplenumbers_treated, treated_isch, time0_isch)))
    }
  }
  legend('topright',colnames(tabsmat),pch=20,lwd=1,col=samplecol, cex = 1.7)
dev.off()

  
#####distribution_of_counts filtered

if (nrow(cpm_count_table) > 2) {
tabsmat_2=as.matrix(log2(cpm_count_table[,c(-1,-2)]+1))
colnames(tabsmat_2)=colnames(cpm_count_table)[c(-1,-2)]

samplecol_2= ifelse(dimnames(tabsmat_2)[[2]] %in% samplenumbers_control, control_isch,
                  ifelse(dimnames(tabsmat_2)[[2]] 
                         %in% samplenumbers_treated, treated_isch, time0_isch))
  
tgz_2=hist(tabsmat_2,breaks = 40, plot=F)

plot_to_png(distribution_of_counts, path_to_output_folder )
if(ncol(tabsmat_2)>=1){
  histlist_2=lapply(1:ncol(tabsmat_2),function(X){ return (hist(tabsmat_2[,X],plot=F,breaks=tgz_2$breaks)) })
  xrange_2=range(unlist(lapply(histlist_2,function(X){X$mids})))
  yrange_2=range(unlist(lapply(histlist_2,function(X){X$counts})))
  hst1_2=histlist_2[[1]]
  plot(hst1_2$mids,hst1_2$counts,type='b',pch=20,xlim=c(0,xrange_2[2]*1.1),ylim=c(0,yrange_2[2]*1.2),
       xlab=bquote(~log[2] ~ "Counts"),ylab='Frequency',main='Distribution of read counts',
       col = ifelse(dimnames(tabsmat_2)[[2]][1] %in% samplenumbers_control, control_isch,
                    ifelse(dimnames(tabsmat_2)[[2]][1] 
                           %in% samplenumbers_treated, treated_isch, time0_isch)), 
       cex.main = 2, cex.axis = 1.5, cex.lab = 1.7)
}
 
if(ncol(tabsmat_2)>=2){ 
  for(i in 2:ncol(tabsmat_2)){
    hstn_2=histlist_2[[i]]
    lines(hstn_2$mids,hstn_2$counts,type='b',pch=20,
          col=ifelse(dimnames(tabsmat_2)[[2]][i] %in% samplenumbers_control, control_isch,
                     ifelse(dimnames(tabsmat_2)[[2]][i] 
                            %in% samplenumbers_treated, treated_isch, time0_isch)))
  }
}
legend('topright',colnames(tabsmat_2),pch=20,lwd=1,col=samplecol_2, cex = 1.7)
dev.off()
}



####PCA-plot:

if (nrow(cpm_count_table) > 2) {
 ctfit_tx<<-prcomp(t(cpm_count_table_log[,c(-1,-2)]),center=TRUE)
  treatment=ifelse(colnames(cpm_count_table_log)[3:length(colnames(cpm_count_table_log))] %in% samplenumbers_control, 
                   name_for_ctrl, ifelse(colnames(cpm_count_table_log)[3:length(colnames(cpm_count_table_log))] 
                          %in% samplenumbers_treated, name_for_treated, "time0"))

  if(length(ctfit_tx$sdev)==1){
    print('Skip PCA plot as there is only one sample.')
    return (0)
  }
   varpca=ctfit_tx$sdev^2
  varpca=varpca/sum(varpca)*100;
  if(length(varpca)>10){
    varpca=varpca[1:10];
  }
  
  pc1_titel <- paste0("PC1 (", round(varpca[1],0), "%)")
  pc2_titel <- paste0("PC2 (", round(varpca[2],0), "%)")
  pc3_titel <- paste0("PC3 (", round(varpca[3],0), "%)")

if(ncol(cpm_count_table_log[,c(-1,-2)])>1){
    pcareport=data.frame(PC1=ctfit_tx$x[,1],PC2=ctfit_tx$x[,2],sample=rownames(ctfit_tx$x))
    p <- ggplot(pcareport,aes(x=PC1,y=PC2,label=sample)) +
      geom_point(aes(colour=treatment), size = 4) +
      scale_colour_manual(values = c(control_isch, treated_isch))+
      xlim(min(pcareport$PC1)*1.2, max(pcareport$PC1)*1.2) + 
      ylim(min(pcareport$PC2)*1.2, max(pcareport$PC2)*1.2) +
      theme_bw(base_size = 22) +
      labs(title = "PCA-plot", x= pc1_titel, y= pc2_titel) +
      theme(plot.title = element_text(size=24)) +
      theme(legend.title = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(size = 6, segment.color = black_isch, segment.size = 0.3)
    
    plot_to_png(pca1vs2, path_to_output_folder )  
    print(p)
    dev.off()
}
  
if(ncol(cpm_count_table_log[,c(-1,-2)])>2){
        pcareport$PC3=ctfit_tx$x[,3]
      p<-ggplot(pcareport,aes(x=PC1,y=PC3,label=sample)) +
        geom_point(aes(colour=treatment), size = 4) +
        scale_colour_manual(values = c(control_isch, treated_isch))+
        xlim(min(pcareport$PC1)*1.2, max(pcareport$PC1)*1.2) + 
        ylim(min(pcareport$PC2)*1.2, max(pcareport$PC2)*1.2) +
        theme_bw(base_size = 22) +
        labs(title = "PCA-plot", x= pc1_titel, y= pc3_titel) +
        theme(plot.title = element_text(size=24)) +
        theme(legend.title = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_text_repel(size = 6, segment.color = black_isch, segment.size = 0.3)

      plot_to_png(pca1vs3, path_to_output_folder )  
      print(p)
      dev.off()
      
  p<-ggplot(pcareport,aes(x=PC2,y=PC3,label=sample)) +
        geom_point(aes(colour=treatment), size = 4) +
        scale_colour_manual(values = c(control_isch, treated_isch))+
        xlim(min(pcareport$PC1)*1.2, max(pcareport$PC1)*1.2) + 
        ylim(min(pcareport$PC2)*1.2, max(pcareport$PC2)*1.2) +
        theme_bw(base_size = 22) +
        labs(title = "PCA-plot", x= pc2_titel, y= pc3_titel) +
        theme(plot.title = element_text(size=24)) +
        theme(legend.title = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_text_repel(size = 6, segment.color = black_isch, segment.size = 0.3)

      plot_to_png(pca2vs3, path_to_output_folder )
      print(p)
      dev.off()
    }

  
####% variance of PCA plot:
  plot_to_png(pca_variance, path_to_output_folder )  
  plot(varpca,type='b',lwd=2,pch=20,xlab='PCs',ylab='% Variance explained');
  dev.off()
}



#plot correlation between samples:

if (nrow(cpm_count_table) > 2) {
samplenumbers_for_plot <- c(samplenumbers_day_0_as_control, samplenumbers_control)
plot_to_png(correlation_control, path_to_output_folder )
plot_samples_against_samples(cpm_count_table, samplenumbers_for_plot)
dev.off()

plot_to_png(correlation_treated, path_to_output_folder )
plot_samples_against_samples(cpm_count_table, samplenumbers_treated)
dev.off()

samplenumbers_for_plot <- (c(samplenumbers_day_0_as_control, samplenumbers_control, samplenumbers_treated))
plot_to_png(correlation_all, path_to_output_folder )
plot_samples_against_samples(cpm_count_table, samplenumbers_for_plot)
dev.off()
}

##########colored clustering
if (ncol(cpm_count_table) > 2) {

h_clust_data <- t(cpm_count_table)
samplenames <- unlist(colnames(cpm_count_table))
h_clust_data <- as.data.frame(h_clust_data)
colnames(h_clust_data) <- h_clust_data[1,]
h_clust_data <- h_clust_data[c(-1, -2),]
h_clust_data$samplenames <- samplenames[c(-1,-2)]
h_clust_data$treatment <- h_clust_data$samplenames

h_clust_data$treatment <- ifelse(h_clust_data$treatment %in% samplenumbers_control, "ctrl",
                                 ifelse(h_clust_data$treatment 
                                        %in% samplenumbers_treated, "treated", "time0"))
h_clust_data <- select(h_clust_data, samplenames, treatment, everything())

## Compute Euclidean distance between samples
dist=dist(h_clust_data[ , c(1:ncol(h_clust_data))] , diag=TRUE)

## Perform clustering with hclust
hc <- hclust(dist)
dhc <- as.dendrogram(hc)
dL <- dendrapply(dhc, colLab)

## plot
plot_to_png(hierarchical_clustering_fancy, path_to_output_folder )
plot(dL , main="Structure of the samples")
legend("topright", 
       legend = c(name_for_treated , name_for_ctrl ), 
       col = c(treated_isch, control_isch), 
       pch = c(20,20), bty = "n",  pt.cex = 3, cex = 0.8 , 
       text.col = black_isch, horiz = FALSE, inset = c(0, 0.1))

dev.off()
}

write.table(as.matrix(head(cpm_count_table)),
              sep = "\t", file = paste0(path_to_output_folder, project_run_name, "fakefile.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE
            )

