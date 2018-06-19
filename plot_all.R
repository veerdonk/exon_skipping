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
  return(exon_lengths)
}
col7_exon_lengths <- find_exon_lengths("uc003ctz.3")
dmd_exon_lengths <- find_exon_lengths("uc004dda.2")


read_and_normalize <- function(fileDir, exon_lengths){
  # read the files into a list of data.frames
  data.list <- lapply(list.files(path = fileDir, full.names = TRUE), function(x) read.csv(x, row.names = NULL))
  
  colnames(data.list[[7]]) <- colnames(data.list[[2]])
  colnames(data.list[[1]]) <- colnames(data.list[[2]])
  
  # concatenate into one big data.frame
  data.cat <- as.data.frame(do.call(rbind, data.list))
  all_samples <- subset(data.cat,!duplicated(data.cat$row.names))
  row.names(all_samples) <- all_samples$row.names
  all_samples$row.names <- c()
  
  colnames(all_samples) <- as.numeric(1:length(colnames(all_samples)))
  #correct for exon length:
  all_samples_n <- as.data.frame(t(t(all_samples) / exon_lengths))
  return(all_samples_n)
}
fileDirCol7 <- "/Users/david/Documents/data/cluster_output/col7"
fileDirDMD <- "/Users/david/Documents/data/cluster_output/dmd"
all_samples_COL7 <- read_and_normalize(fileDirCol7, col7_exon_lengths)
all_samples_DMD <- read_and_normalize(fileDirDMD, dmd_exon_lengths)
all_genes <- list(all_samples_COL7, all_samples_DMD)

# #row.names(data.cat) <- as.numeric(sapply(row.names(data.cat), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
# colnames(all_samples) <- as.numeric(sapply(colnames(all_samples), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
# skipped_exon_counts_mlt <- melt(t(all_samples))
# non_zero_mlt <- skipped_exon_counts_mlt[skipped_exon_counts_mlt$value > 0,]
# not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
# skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
# non_zero_mlt$skippable <- as.factor(ifelse(non_zero_mlt$Var1 %in% skippable_exons, "in frame", "out of frame"))
# non_zero_mlt <- non_zero_mlt[!is.na(non_zero_mlt$value),]
# non_zero_mlt <- non_zero_mlt[!is.infinite(non_zero_mlt$value),]
# all_sample_nn <- ggplot(non_zero_mlt, aes(x=Var1, y=value, group=Var1, color = skippable)) + labs(x="exon", y="read count", title="All samples* (not normalized)") + geom_boxplot() 
# 
# #normalized
# colnames(all_samples_n) <- as.numeric(sapply(colnames(all_samples_n), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
# skipped_exon_counts_mlt <- melt(t(all_samples_n))
# non_zero_mlt <- skipped_exon_counts_mlt[skipped_exon_counts_mlt$value > 0,]
# not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
# skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
# non_zero_mlt$skippable <- as.factor(ifelse(non_zero_mlt$Var1 %in% skippable_exons, "in frame", "out of frame"))
# non_zero_mlt <- non_zero_mlt[!is.na(non_zero_mlt$value),]
# non_zero_mlt <- non_zero_mlt[!is.infinite(non_zero_mlt$value),]
# all_sample_n <- ggplot(non_zero_mlt, aes(x=Var1, y=value, group=Var1, color = skippable)) + labs(x="exon", y="read count", title="All samples* corrected for exon length") + geom_boxplot() 
# 
# 
# barData <- data.frame(colMeans(na.omit(all_samples_n)))
# colnames(barData) <- "dat"
# barData$exon <- factor(as.character(row.names(barData)), levels = unique(row.names(barData)))
# barData$frame <- as.factor(ifelse(barData$exon %in% skippable_exons, "in frame", "out of frame"))
# bar <- ggplot(data=barData, aes(x=exon, y=dat, color=frame)) + geom_bar(stat = "identity") + theme(axis.text=element_text(size=4)) + labs(x="exon", y="mean read count", title="COL7 skipped exons") + scale_color_manual(values=c("#5ab4ac", "#d8b365"))
# 
# colsums <- data.frame(colSums(na.omit(all_samples_n)))
# colnames(colsums) <- "val"
# exon_nums <- 1:118
# colsums$exon <- factor(as.character(exon_nums), levels = unique(exon_nums))
# colsums$frame <- as.factor(ifelse(colsums$exon %in% skippable_exons, "in frame", "out of frame"))
# ggplot(data = colsums, aes(y=val, x=exon, color=frame)) + geom_bar(stat = "identity") + theme(axis.text=element_text(size=4)) + labs(x="exon", y="stacked counts", title="COL7 skipped exons") + scale_color_manual(values=c("#21f000", "#ff2d00"))
# 
# 
# total_jx_col7 <- read.csv("~/Documents/data/number_of_junctions.csv", header = FALSE)
# 
# binned <- count(na.omit(total_jx_col7$V2))
# 
# ggplot(binned, aes(x=x, y=freq)) + geom_bar(stat='identity', width=1) + labs(title="Number of junctions in COL7A1", x="Number of junctions", y="Number of samples")
# 
# 
# multi_sample <- data.frame(colMeans(na.omit(all_samples_n[row.names(all_samples_n) %in% sample(row.names(all_samples_n), round(length(row.names(all_samples_n))/2)),])))
# for(i in 1:1000){
#   multi_sample <- cbind(multi_sample ,colMeans(na.omit(all_samples_n[row.names(all_samples_n) %in% sample(row.names(all_samples_n), round(length(row.names(all_samples_n))/2)),])))
# }
# test <- data.frame(colMeans(na.omit(all_samples_n[row.names(all_samples_n) %in% sample(row.names(all_samples_n), round(length(row.names(all_samples_n))/50)),])))
# colnames(test) <- "val"
# exon_nums <- 1:118
# test$exon <- factor(as.character(exon_nums), levels = unique(exon_nums))
# test$frame <- as.factor(ifelse(test$exon %in% skippable_exons, "in frame", "out of frame"))
# ggplot(data = test, aes(y=val, x=exon, color=frame)) + geom_bar(stat = "identity") + theme(axis.text=element_text(size=4)) + labs(x="exon", y="stacked counts", title="COL7 skipped exons") + scale_color_manual(values=c("#21f000", "#ff2d00"))
# 
# sampled_sds <- apply(multi_sample, 1, sd)
# bar + geom_errorbar(aes(ymin=dat-sampled_sds, ymax=dat+sampled_sds), width=.2)

make_bar_plot <- function(tbl, skippable_exons){
  multi_sample <- data.frame(colMeans(na.omit(tbl[row.names(tbl) %in% sample(row.names(tbl), round(length(row.names(tbl))/2)),])))
  for(i in 1:1000){
    multi_sample <- cbind(multi_sample ,colMeans(na.omit(tbl[row.names(tbl) %in% sample(row.names(tbl), round(length(row.names(tbl))/2)),])))
  }
   
  
  sampled_sds <- apply(multi_sample, 1, sd)
  barData <- data.frame(colMeans(na.omit(tbl)))
  colnames(barData) <- "dat"
  barData$exon <- factor(as.character(row.names(barData)), levels = unique(row.names(barData)))
  barData$frame <- as.factor(ifelse(barData$exon %in% skippable_exons, "in frame", "out of frame"))
  bar <- ggplot(data=barData, aes(x=exon, y=dat, color=frame)) + geom_bar(stat = "identity", position = "dodge") + theme(axis.text=element_text(size=4)) + scale_color_manual(values=c("#5ab4ac", "#d8b365"))
  return(bar + geom_errorbar(aes(ymin=dat-sampled_sds, ymax=dat+sampled_sds), width=.2))
}

make_multi_bar_plot <- function(tbls, names){
  barplot_data <- data.frame()
  for(tbl_num in 1:length(tbls)){
    tbl_dat <- data.frame(colMeans(na.omit(tbls[[tbl_num]])))
    colnames(tbl_dat) <- "means"
    tbl_dat$exon <- 1:length(tbl_dat$means)
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


##---------------------------------specify gene to use---------------------------------------##

all_samples_n <- all_samples_COL7#all_samples_DMD#
skippable_exons <-   skippable_exons_col7#skippable_exons_dmd#
gene <- "COL7A1"#"DMD"#  
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

make_multi_bar_plot(list(skin_data, blood_data, brain_data), c("skin", "blood", "brain", "muscle"))
make_multi_bar_plot(list(skin_data, esophagus_data, intestine_data), c("skin", "esophagus", "intestine"))


##---------------------------------mutations scoring-----------------------------------------##

mutations <- read.csv("~/Documents/data/Mutations_2018-06-18_02_14_24.csv", stringsAsFactors = FALSE)

exonic_mutations <- mutations[mutations$Exon %in% unlist(lapply(1:118, function(x) paste0("exon " , x))),]
exonic_mutations$Exon <- unlist(lapply(exonic_mutations$Exon, function(x) as.numeric(strsplit(x, " ")[[1]][2])))
col7_means <- colMeans(na.omit(all_samples_COL7))
col7_means[exonic_mutations$Exon[6]]
hist(exonic_mutations$Exon,breaks=118, main = "mutations on deb-central", xlab = "exon")


