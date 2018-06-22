library("reshape2")
library("ggplot2")
library("plyr")
library(gridExtra)

#exon info
all_exon_bed <-  read.delim("/Users/david/Documents/data/exon_table_ucsc.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(all_exon_bed) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
all_exon_bed$ucsc_id <- unlist(lapply(all_exon_bed$info, function(x) strsplit(x, "_")[[1]][1]))

left_align_plots <- function(plots){
  grobs <- list()
  widths <- list()
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  do.call("grid.arrange", c(grobs, ncol = 1))
}

find_exon_lengths <- function(id){
  exons <- all_exon_bed[all_exon_bed$ucsc_id == id,]
  row.names(exons) <- 1:nrow(exons)
  exon_lengths <- (exons$stop - exons$start)
  return(list(exon_lengths, exons$strand[1]))
}
col7_exon_info <- find_exon_lengths("uc003ctz.3")
col7_exon_lengths <- col7_exon_info[[1]]
col7_strand <- col7_exon_info[[2]]
dmd_exon_lengths <- find_exon_lengths("uc004dda.2")
col17_exon_lengths <- find_exon_lengths("uc001kxr.4")
ttn_exon_lengths <- find_exon_lengths("uc031rqc.3")
gpi_exon_lengths <- find_exon_lengths("uc002nvg.3")

read_and_normalize <- function(fileDir, exon_lengths, strand){
  # read the files into a list of data.frames
  data.list <- lapply(list.files(path = fileDir, full.names = TRUE), function(x) read.csv(x, row.names = NULL))
  
  colnames(data.list[[7]]) <- colnames(data.list[[2]])
  colnames(data.list[[1]]) <- colnames(data.list[[2]])
  
  # concatenate into one big data.frame
  data.cat <- as.data.frame(do.call(rbind, data.list))
  all_samples <- subset(data.cat,!duplicated(data.cat$row.names))
  row.names(all_samples) <- all_samples$row.names
  all_samples$row.names <- c()
  colnms <- as.numeric(1:length(colnames(all_samples)))
  if(strand == "-"){
    colnms <- rev(colnms)
  }
  colnames(all_samples) <- colnms
  #correct for exon length:
  all_samples_n <- as.data.frame(t(t(all_samples) / exon_lengths))
  return(all_samples_n)
}
fileDirCol17 <- "/Users/david/Documents/data/cluster_output/col17"
fileDirCol7 <- "/Users/david/Documents/data/cluster_output/col7"
fileDirCol71ex <- "/Users/david/Documents/data/cluster_output/col7_1ex/"
fileDirDMD <- "/Users/david/Documents/data/cluster_output/dmd"
fileDirTTN <- "/Users/david/Documents/data/cluster_output/ttn"
fileDirGPI <- "/Users/david/Documents/data/cluster_output/gpi"
test <- "~/Documents/data/reverse_test/"
all_samples_COL17 <- read_and_normalize(fileDirCol17, col17_exon_lengths)
all_samples_COL7 <- read_and_normalize(fileDirCol7, col7_exon_lengths, col7_strand)
all_samples_COL7_1ex <- read_and_normalize(fileDirCol71ex, col7_exon_lengths, col7_strand)
all_samples_DMD <- read_and_normalize(fileDirDMD, dmd_exon_lengths)
all_samples_TTN <- read_and_normalize(fileDirTTN, ttn_exon_lengths)
all_samples_GPI <- read_and_normalize(fileDirGPI, gpi_exon_lengths)
all_test <- read_and_normalize(test, col7_exon_lengths)
#all_genes <- list(all_samples_COL7, all_samples_DMD)

make_bar_plot <- function(tbl, skippable_exons){
  multi_sample <- data.frame(colMeans(na.omit(tbl[row.names(tbl) %in% sample(row.names(tbl), round(length(row.names(tbl))/2)),])))
  for(i in 1:1000){
    multi_sample <- cbind(multi_sample ,colMeans(na.omit(tbl[row.names(tbl) %in% sample(row.names(tbl), round(length(row.names(tbl))/2)),])))
  }
  sampled_sds <- apply(multi_sample, 1, sd)
  barData <- data.frame(colMeans(na.omit(tbl)))
  colnames(barData) <- "dat"
  barData$exon <- as.numeric(row.names(barData))#factor(as.character(row.names(barData)), levels = unique(row.names(barData)))
  barData$values <- colSums(na.omit(tbl > 0))
  if(length(skippable_exons) > 1){ 
  barData$frame <- as.factor(ifelse(barData$exon %in% skippable_exons, "in frame", "out of frame"))
  }
  else{
    barData$frame <- "unknown"
  }
  bar <- ggplot(data=barData, aes(x=exon, y=dat, color=frame))
  bar <- bar + geom_bar(stat = "identity", position = "dodge")
  bar <- bar + theme(axis.text=element_text(size=8))
  bar <- bar + scale_color_manual(values=c("#5ab4ac", "#d8b365"))
  bar <- bar + geom_errorbar(aes(ymin=dat-sampled_sds, ymax=dat+sampled_sds), width=.2)
  bar <- bar + geom_text(aes(label=values), size=1.5, color="black", vjust = -0.5)
  
  return(bar)
}

make_stacked_bar_plot <-  function(tbl, skippable_exons){
  sums <- data.frame(colSums(na.omit(tbl)))
  colnames(sums) <- "stacked"
  sums$values <- colSums(na.omit(tbl > 0)) 
  sums$exon <- as.numeric(row.names(sums))
  if(length(skippable_exons) > 1){ 
    sums$frame <- as.factor(ifelse(sums$exon %in% skippable_exons, "in frame", "out of frame"))
  }
  else{
    sums$frame <- "unknown"
  }
  bar <- ggplot(data=sums, aes(x=exon, y=stacked, color=frame))
  bar <- bar + geom_bar(stat = "identity", position = "dodge")
  bar <- bar + theme(axis.text=element_text(size=8))
  bar <- bar + scale_color_manual(values=c("#5ab4ac", "#d8b365"))
  bar <- bar + geom_text(aes(label=values), size=1.5, color="black", vjust = -0.5)
  
  return(bar)
}

make_multi_bar_plot <- function(tbls, names){
  barplot_data <- data.frame()
  for(tbl_num in 1:length(tbls)){
    tbl_dat <- data.frame(colMeans(na.omit(tbls[[tbl_num]])))
    colnames(tbl_dat) <- "means"
    tbl_dat$exon <- rev(1:length(tbl_dat$means))#factor(as.character(1:length(tbl_dat$means)), levels = unique(1:length(tbl_dat$means)))
    tbl_dat$name <- names[tbl_num]
    barplot_data <- rbind(barplot_data, tbl_dat)
  }

  return(ggplot(barplot_data, aes(x=exon, y=means, fill=name)) + geom_bar(stat = "identity", position = "identity", alpha=.5)+scale_fill_brewer(palette="Dark2"))
}

sample_annotation <- read.delim("~/Documents/data/sample_annotation_24_03_2016_withFibroblasts.txt", stringsAsFactors = FALSE)

##---------------------------------Select tissues--------------------------------------------##
skin_studies <- sample_annotation[sample_annotation$annotation_organism_part == "skin",1]
brain_studies <- sample_annotation[sample_annotation$annotation_organism_part == "brain",1]
muscle_studies <- sample_annotation[sample_annotation$annotation_organism_part == "muscle",1]
blood_studies <- sample_annotation[sample_annotation$annotation_organism_part == "blood",1]
liver_studies <- sample_annotation[sample_annotation$annotation_organism_part == "liver",1]
heart_studies <- sample_annotation[sample_annotation$annotation_organism_part == "heart",1]
esophagus_studies <- sample_annotation[sample_annotation$annotation_organism_part == "esophagus",1]
intestine_studies <- sample_annotation[sample_annotation$annotation_organism_part == "intestine",1]

##----------------------------specify in/out of frame exons----------------------------------##
not_skip_exons_col7 <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
skippable_exons_col7 <- (1:118)[-not_skip_exons_col7] # exons that can be skipped (1 based)
not_skip_exons_dmd <- c(1,2,6,7,8,11,12,17,18,19,20,21,22,43,44,45,46,50,51,52,53,54,55,56,57,58,59,61,62,63,65,66,67,68,69,70,75,76,78,79)
skippable_exons_dmd <- (1:79)[-not_skip_exons_dmd]
skippable_exons_col17 <- c()
skippable_exons_ttn <- c()
skippable_exons_gpi <- c()

##---------------------------------specify gene to use---------------------------------------##

all_samples_n <- all_samples_COL7#all_samples_COL7_1ex#all_samples_GPI#all_samples_COL17#all_samples_TTN#all_samples_DMD#
skippable_exons <-   skippable_exons_col7#skippable_exons_gpi#skippable_exons_col17#skippable_exons_ttn#skippable_exons_dmd#
gene <- "COL7A1 (10 skips)"#"GPI"#"COL17A1"#"TTN"#"DMD"#  
##--------------------------------get data for tissues---------------------------------------##
skin_data <- all_samples_n[row.names(all_samples_n) %in% skin_studies,]
brain_data <- all_samples_n[row.names(all_samples_n) %in% brain_studies,]
muscle_data <- all_samples_n[row.names(all_samples_n) %in% muscle_studies,]
blood_data <- all_samples_n[row.names(all_samples_n) %in% blood_studies,]
liver_data <- all_samples_n[row.names(all_samples_n) %in% liver_studies,]
heart_data <- all_samples_n[row.names(all_samples_n) %in% heart_studies,]
esophagus_data <- all_samples_n[row.names(all_samples_n) %in% esophagus_studies,]
intestine_data <- all_samples_n[row.names(all_samples_n) %in% intestine_studies,]

##---------------------------------create ggplots--------------------------------------------##
b_all <- make_bar_plot(all_samples_n,skippable_exons) + labs(x="exon",y="mean read count", title=paste0(gene," skipped exons in all samples (n=", length(all_samples_COL7$`1`), ") (>0=", sum(na.omit(all_samples_COL7)>0), ")"))
b_skin <- make_bar_plot(skin_data,skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in skin samples (n=",length(skin_data$`1`), ") (>0=", sum(na.omit(skin_data)>0),")"))
b_brain <- make_bar_plot(brain_data, skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in brain samples (n=",length(brain_data$`1`), ") (>0=", sum(na.omit(brain_data)>0),")"))
b_muscle <- make_bar_plot(muscle_data, skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in muscle samples (n=",length(muscle_data$`1`), ") (>0=", sum(na.omit(muscle_data)>0),")"))
b_blood <- make_bar_plot(blood_data,skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in blood samples (n=",length(blood_data$`1`), ") (>0=", sum(na.omit(blood_data)>0),")"))
b_liver <- make_bar_plot(liver_data, skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in liver samples (n=",length(liver_data$`1`), ") (>0=", sum(na.omit(liver_data)>0),")"))
b_heart <- make_bar_plot(heart_data, skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in heart samples (n=",length(heart_data$`1`), ") (>0=", sum(na.omit(heart_data)>0),")"))
b_esophagus <- make_bar_plot(esophagus_data, skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in esophagus samples (n=",length(esophagus_data$`1`), ") (>0=", sum(na.omit(esophagus_data)>0),")"))
b_intestine <- make_bar_plot(intestine_data, skippable_exons) + labs(x="exon", y="mean read count", title=paste0(gene," skipped exons in intestine samples (n=",length(intestine_data$`1`), ") (>0=", sum(na.omit(intestine_data)>0),")"))
##---------------------------------display plots---------------------------------------------##
left_align_plots(list(b_skin, b_all, b_blood))
left_align_plots(list(b_brain, b_all, b_muscle))
left_align_plots(list(b_liver, b_all, b_heart))

make_multi_bar_plot(list(skin_data, blood_data, brain_data, muscle_data), c("skin", "blood", "brain", "muscle"))
make_multi_bar_plot(list(skin_data, esophagus_data, intestine_data), c("skin", "esophagus", "intestine"))
make_multi_bar_plot(list(skin_data, all_samples_COL17), c("skin", "all"))

##---------------------------------mutations scoring-----------------------------------------##

dominant_mutations <- read.csv("~/Documents/data/deb_central_dominant.csv", stringsAsFactors = FALSE)
recessive_mutations <- read.csv("~/Documents/data/deb_central_recessive.csv", stringsAsFactors = FALSE)
recessive_mutations$inheritance <- "recessive"
dominant_mutations$inheritance <- "dominant"
mutations <- rbind(recessive_mutations, dominant_mutations)

exonic_mutations <- mutations[mutations$Exon %in% unlist(lapply(1:118, function(x) paste0("exon " , x))),]
exonic_mutations$Exon <- unlist(lapply(exonic_mutations$Exon, function(x) as.numeric(strsplit(x, " ")[[1]][2])))
col7_means <- colMeans(na.omit(all_samples_COL7))
col7_means[exonic_mutations$Exon[6]]
hist(exonic_mutations$Exon,breaks=118, main = "mutations on deb-central", xlab = "exon")

hist(exonic_mutations[exonic_mutations$inheritance == "dominant",5], breaks = 120)
hist(exonic_mutations[exonic_mutations$inheritance == "recessive",5], breaks = 120)
ggplot(exonic_mutations[exonic_mutations$inheritance == "recessive",], aes(Exon)) + geom_histogram(bins = 118)
ggplot(exonic_mutations[exonic_mutations$inheritance == "dominant",], aes(Exon)) + geom_histogram(bins = 118)


