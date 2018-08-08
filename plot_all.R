##------------------------------------Load libraries----------------------------------------##
library("reshape2")
library("ggplot2")
library("plyr")
library(gridExtra)
##------------------------------------------------------------------------------------------##

##-----------------------------------read exon .bed ----------------------------------------##
all_exon_bed <-  read.delim("/Users/david/Documents/data/exon_table_ucsc.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(all_exon_bed) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
all_exon_bed$ucsc_id <- unlist(lapply(all_exon_bed$info, function(x) strsplit(x, "_")[[1]][1]))
##------------------------------------------------------------------------------------------##

##----------------------------Function for aligning plots-----------------------------------##
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
##------------------------------------------------------------------------------------------##

##---------------------------------Extracts exon lengths------------------------------------##
find_exon_lengths <- function(id){
  exons <- all_exon_bed[all_exon_bed$ucsc_id == id,]
  row.names(exons) <- 1:nrow(exons)
  exon_lengths <- (exons$stop - exons$start)
  return(list(exon_lengths, exons$strand[1]))
}
##------------------------------------------------------------------------------------------##

##--------------------------------Get lengths for genes-------------------------------------##
col7_exon_info <- find_exon_lengths("uc003ctz.3")
col7_exon_lengths <- col7_exon_info[[1]]
col7_strand <- col7_exon_info[[2]]
dmd_exon_lengths <- find_exon_lengths("uc004dda.2")
col17_exon_lengths <- find_exon_lengths("uc001kxr.4")
ttn_exon_lengths <- find_exon_lengths("uc031rqc.3")
gpi_exon_lengths <- find_exon_lengths("uc002nvg.3")
##------------------------------------------------------------------------------------------##

##------------------------function to read and normalize data-------------------------------##
read_and_normalize <- function(fileDir, exon_lengths, strand){
  # read the files into a list of data.frames
  data.list <- lapply(list.files(path = fileDir, full.names = TRUE), function(x) read.csv(x, row.names = NULL))
  
  #For the col7 data some colnames heve to be fixed
  colnames(data.list[[7]]) <- colnames(data.list[[2]])
  colnames(data.list[[1]]) <- colnames(data.list[[2]])
  
  # concatenate into one big data.frame
  data.cat <- as.data.frame(do.call(rbind, data.list))
  all_samples <- subset(data.cat,!duplicated(data.cat$row.names))
  row.names(all_samples) <- all_samples$row.names #fix row names
  all_samples$row.names <- c()
  colnms <- as.numeric(1:length(colnames(all_samples)))
  if(strand == "-"){ #if gene is on reverse strand reverse the colnames (exons)
    colnms <- rev(colnms)
  }
  colnames(all_samples) <- colnms
  
  # Currently not correcting for exon length
  # To correct for exon length switch the 2 lines below
  all_samples_n <- all_samples
  #all_samples_n <- as.data.frame(t(t(all_samples) / exon_lengths))
  
  return(all_samples_n)
}
##------------------------------------------------------------------------------------------##


##-----------------------------Define directories and run-----------------------------------##
fileDirCol7_no_scale <- "/Users/david/Documents/data/cluster_output/col7_no_scale"
fileDirCol7_new_scale <- "/Users/david/Documents/data/cluster_output/col7_new_scale"
fileDirCol7_false_positives <- "/Users/david/Documents/data/cluster_output/col7_false_positives"
fileDirCol7_false_positives_old <- "/Users/david/Documents/data/cluster_output/old_data/false_positives/old"
fileDirCol17 <- "/Users/david/Documents/data/cluster_output/col17"
fileDirCol7 <- "/Users/david/Documents/data/cluster_output/col7"
fileDirCol71ex <- "/Users/david/Documents/data/cluster_output/col7_1ex/"
fileDirDMD <- "/Users/david/Documents/data/cluster_output/dmd"
fileDirTTN <- "/Users/david/Documents/data/cluster_output/ttn"
fileDirGPI <- "/Users/david/Documents/data/cluster_output/gpi"
test <- "~/Documents/data/reverse_test/"
all_samples_COL7_new_scale <- read_and_normalize(fileDirCol7_new_scale, col7_exon_lengths, col7_strand)
all_samples_COL7_false_positives <- read_and_normalize(fileDirCol7_false_positives, col7_exon_lengths, col7_strand)
all_samples_COL7_false_positives_old <- read_and_normalize(fileDirCol7_false_positives_old, col7_exon_lengths, col7_strand)
all_samples_COL7_no_scale <- read_and_normalize(fileDirCol7_no_scale, col7_exon_lengths, col7_strand)
all_samples_COL17 <- read_and_normalize(fileDirCol17, col17_exon_lengths)
all_samples_COL7 <- read_and_normalize(fileDirCol7, col7_exon_lengths, col7_strand)
all_samples_COL7_1ex <- read_and_normalize(fileDirCol71ex, col7_exon_lengths, col7_strand)
all_samples_DMD <- read_and_normalize(fileDirDMD, dmd_exon_lengths)
all_samples_TTN <- read_and_normalize(fileDirTTN, ttn_exon_lengths)
all_samples_GPI <- read_and_normalize(fileDirGPI, gpi_exon_lengths)
all_test <- read_and_normalize(test, col7_exon_lengths)
##------------------------------------------------------------------------------------------##

##-----------------------------Creates ggplot barplots--------------------------------------##
make_bar_plot <- function(tbl, skippable_exons){
  # Randomly select half of the samples and calculate means
  multi_sample <- data.frame(colMeans(na.omit(tbl[row.names(tbl) %in% sample(row.names(tbl), round(length(row.names(tbl))/2)),])))
  # repeat 1000x
  for(i in 1:1000){
    multi_sample <- cbind(multi_sample ,colMeans(na.omit(tbl[row.names(tbl) %in% sample(row.names(tbl), round(length(row.names(tbl))/2)),])))
  }
  sampled_sds <- apply(multi_sample, 1, sd) #get sd
  # construct ggplot compliant data
  barData <- data.frame(colMeans(na.omit(tbl)))
  colnames(barData) <- "dat" 
  barData$exon <- rev(1:118)#factor(as.character(row.names(barData)), levels = unique(row.names(barData)))#as.numeric(row.names(barData))#factor(as.character(row.names(barData)), levels = unique(row.names(barData)))
  barData$values <- colSums(na.omit(tbl > 0))
  #determine whether exon is in- or out-of-frame
  if(length(skippable_exons) > 1){ 
  barData$frame <- as.factor(ifelse(barData$exon %in% skippable_exons, "in frame", "out of frame"))
  }
  else{
    barData$frame <- "unknown"
  }
  #build the plot object
  bar <- ggplot(data=barData, aes(x=exon, y=dat, color=frame)) #add data
  bar <- bar + geom_bar(stat = "identity", position = "dodge") #add type of plot
  bar <- bar + theme(axis.text=element_text(size=8)) #make text smaller
  bar <- bar + scale_color_manual(values=c("#5ab4ac", "#d8b365")) #choose bar colors
  bar <- bar + geom_errorbar(aes(ymin=dat-sampled_sds, ymax=dat+sampled_sds), width=.2) #add error bars
  bar <- bar + geom_text(aes(label=values), size=1.5, color="black", vjust = -0.5) #add number of samples above each bar
  
  return(bar)
}
##------------------------------------------------------------------------------------------##


##-------------------------------Make stacked bar plots-------------------------------------##
make_stacked_bar_plot <-  function(tbl, skippable_exons){
  #Works similarly to make_bar_plot but instead of means it uses stacked data
  sums <- data.frame(colSums(na.omit(tbl)))
  colnames(sums) <- "stacked"
  #create ggplot compliant data
  sums$values <- colSums(na.omit(tbl > 0)) 
  sums$exon <- as.numeric(row.names(sums))
  #determine in- or out of frame for each exon
  if(length(skippable_exons) > 1){ 
    sums$frame <- as.factor(ifelse(sums$exon %in% skippable_exons, "in frame", "out of frame"))
  }
  else{
    sums$frame <- "unknown"
  }
  
  #construct ggplot object
  bar <- ggplot(data=sums, aes(x=exon, y=stacked, color=frame)) #add data
  bar <- bar + geom_bar(stat = "identity", position = "dodge") #add plot type
  bar <- bar + theme(axis.text=element_text(size=8)) #set text size
  bar <- bar + scale_color_manual(values=c("#5ab4ac", "#d8b365")) #set colors
  bar <- bar + geom_text(aes(label=values), size=1.5, color="black", vjust = -0.5) #add number of samples above bars
  
  return(bar)
}
##------------------------------------------------------------------------------------------##


##-----------------------------Combine plots by alpha blending------------------------------##
make_multi_bar_plot <- function(tbls, names){
  barplot_data <- data.frame()
  # loops through all provided tables and adds them to a data frame
  for(tbl_num in 1:length(tbls)){
    tbl_dat <- data.frame(colMeans(na.omit(tbls[[tbl_num]])))
    colnames(tbl_dat) <- "means"
    tbl_dat$exon <- rev(1:length(tbl_dat$means))#factor(as.character(1:length(tbl_dat$means)), levels = unique(1:length(tbl_dat$means)))
    tbl_dat$name <- names[tbl_num]
    barplot_data <- rbind(barplot_data, tbl_dat)
  }

  return(ggplot(barplot_data, aes(x=exon, y=means, fill=name)) + geom_bar(stat = "identity", position = "identity", alpha=.5)+scale_fill_brewer(palette="Dark2"))
}
##------------------------------------------------------------------------------------------##

##-------------------------Read tissue annotation and extract tissues-----------------------##
sample_annotation <- read.delim("~/Documents/data/sample_annotation_24_03_2016_withFibroblasts.txt", stringsAsFactors = FALSE)

skin_studies <- sample_annotation[sample_annotation$annotation_organism_part == "skin",1]
brain_studies <- sample_annotation[sample_annotation$annotation_organism_part == "brain",1]
muscle_studies <- sample_annotation[sample_annotation$annotation_organism_part == "muscle",1]
blood_studies <- sample_annotation[sample_annotation$annotation_organism_part == "blood",1]
liver_studies <- sample_annotation[sample_annotation$annotation_organism_part == "liver",1]
heart_studies <- sample_annotation[sample_annotation$annotation_organism_part == "heart",1]
esophagus_studies <- sample_annotation[sample_annotation$annotation_organism_part == "esophagus",1]
intestine_studies <- sample_annotation[sample_annotation$annotation_organism_part == "intestine",1]

##----------------------------specify in/out of frame exons----------------------------------##
not_skip_exons_col7 <- c(1,2,3,4,6,7,24,25,27,113,118) # exons that cant be skipped (1 based)
skippable_exons_col7 <- (1:118)[-not_skip_exons_col7] # exons that can be skipped (1 based)
not_skip_exons_dmd <- c(1,2,6,7,8,11,12,17,18,19,20,21,22,43,44,45,46,50,51,52,53,54,55,56,57,58,59,61,62,63,65,66,67,68,69,70,75,76,78,79)
skippable_exons_dmd <- (1:79)[-not_skip_exons_dmd]
skippable_exons_col17 <- c()
skippable_exons_ttn <- c()
skippable_exons_gpi <- c()

##---------------------------specify gene to use for upcoming parts--------------------------##
all_samples_n <- all_samples_COL7#all_samples_COL7_new_scale#all_samples_COL7_false_positives#all_samples_COL7_no_scale#all_samples_GPI#all_samples_COL17#all_samples_TTN#all_samples_DMD#
skippable_exons <-   skippable_exons_col7#skippable_exons_gpi#skippable_exons_col17#skippable_exons_ttn#skippable_exons_dmd#
gene <- "COL7A1"#"GPI"#"COL17A1"#"TTN"#"DMD"#  
##-------------------------------------------------------------------------------------------##


##--------------------------------get data for tissues---------------------------------------##
skin_data <- all_samples_n[row.names(all_samples_n) %in% skin_studies,]
brain_data <- all_samples_n[row.names(all_samples_n) %in% brain_studies,]
muscle_data <- all_samples_n[row.names(all_samples_n) %in% muscle_studies,]
blood_data <- all_samples_n[row.names(all_samples_n) %in% blood_studies,]
liver_data <- all_samples_n[row.names(all_samples_n) %in% liver_studies,]
heart_data <- all_samples_n[row.names(all_samples_n) %in% heart_studies,]
esophagus_data <- all_samples_n[row.names(all_samples_n) %in% esophagus_studies,]
intestine_data <- all_samples_n[row.names(all_samples_n) %in% intestine_studies,]
##------------------------------------------------------------------------------------------##


##---------------------------------create ggplots-------------------------------------------##
b_all <- make_bar_plot(all_samples_n,skippable_exons) + labs(x="exon",y="mean normalized read count", title=paste0(gene," skipped exons in all samples"))# (n=", length(all_samples_n$`1`), ") (>0=", sum(na.omit(all_samples_n)>0), ")"))
b_skin <- make_bar_plot(skin_data,skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in skin samples")) #(n=",length(skin_data$`1`), ") (>0=", sum(na.omit(skin_data)>0),")"))
b_brain <- make_bar_plot(brain_data, skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in brain samples (n=",length(brain_data$`1`), ") (>0=", sum(na.omit(brain_data)>0),")"))
b_muscle <- make_bar_plot(muscle_data, skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in muscle samples (n=",length(muscle_data$`1`), ") (>0=", sum(na.omit(muscle_data)>0),")"))
b_blood <- make_bar_plot(blood_data,skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in blood samples (n=",length(blood_data$`1`), ") (>0=", sum(na.omit(blood_data)>0),")"))
b_liver <- make_bar_plot(liver_data, skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in liver samples (n=",length(liver_data$`1`), ") (>0=", sum(na.omit(liver_data)>0),")"))
b_heart <- make_bar_plot(heart_data, skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in heart samples (n=",length(heart_data$`1`), ") (>0=", sum(na.omit(heart_data)>0),")"))
b_esophagus <- make_bar_plot(esophagus_data, skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in esophagus samples")) #(n=",length(esophagus_data$`1`), ") (>0=", sum(na.omit(esophagus_data)>0),")"))
b_intestine <- make_bar_plot(intestine_data, skippable_exons) + labs(x="exon", y="mean normalized read count", title=paste0(gene," skipped exons in intestine samples (n=",length(intestine_data$`1`), ") (>0=", sum(na.omit(intestine_data)>0),")"))
##------------------------------------------------------------------------------------------##

##---------------------------------display plots--------------------------------------------##
left_align_plots(list(b_skin, b_all, b_blood))
left_align_plots(list(b_brain, b_all, b_muscle))
left_align_plots(list(b_liver, b_all, b_heart))

make_multi_bar_plot(list(skin_data, blood_data, brain_data, muscle_data), c("skin", "blood", "brain", "muscle"))
make_multi_bar_plot(list(skin_data, esophagus_data, intestine_data), c("skin", "esophagus", "intestine"))
make_multi_bar_plot(list(skin_data, all_samples_COL17), c("skin", "all"))
##------------------------------------------------------------------------------------------##


##-----------------------------SNPs and deb-cental mutations--------------------------------##
#install.packages("vcfR")
library("vcfR")

#read files for snps and deb-central mutations
col7_snps <- read.delim("/Users/david/Documents/data/all_col7a1snps.tsv")
col7_vcf <- read.vcfR("/Users/david/Documents/data/debcentral_10mei2017_variants_gavin.vcf")
#build 37 exon annotation because versions mismatch
col7_37 <- read.delim("/Users/david/Documents/data/col7_hg37_exons.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(col7_37) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
col7_37$ucsc_id <- unlist(lapply(col7_37$info, function(x) strsplit(x, "_")[[1]][1])) #extract ucsc ID

col7_37_exons <- col7_37[col7_37$ucsc_id == "uc003ctz.2",] #extract col7
col7_37_exons$exon <- rev(1:118)

col7_known_mutations <- data.frame(as.numeric(col7_vcf@fix[,"POS"])) #extract positions
colnames(col7_known_mutations) <- "position"
#extract cadd scores and exons for deb-central
col7_known_mutations$cadd <- unlist(lapply(col7_vcf@fix[,"INFO"], function(x) as.numeric(strsplit(strsplit(x, ";")[[1]][2], "=")[[1]][2])))
col7_known_mutations$exon <- unlist(lapply(col7_vcf@fix[,"INFO"], function(x) as.numeric(strsplit(strsplit(x, "|", fixed = TRUE)[[1]][9], "/", fixed = TRUE)[[1]][1])))

find_exon_for_pos <- function(pos){
  #use positions to find what exon mutation belongs to for SNPs
  for(start in col7_37_exons$start){
    if(start <= pos){
      for(stop in col7_37_exons$stop){
        if(stop >= pos){
          return(col7_37_exons[col7_37_exons$stop == stop,8])
        }
      }
    }
  }
  return(NA)
}

col7_snps$exon <- unlist(lapply(col7_snps$pos, find_exon_for_pos))#run function
mut_counts <- count(na.omit(col7_snps$exon)) #count number of occurences for each exon
# plot the counts in a bar plot
snp_plot <- ggplot(data=mut_counts, aes(y=freq, x=x)) + geom_bar(stat="identity") + labs(x="exon", y="number of SNPs", title="Distribution of SNPs in COL7A1")


deb_central_exons <- c()
# loop through the lines of the .vcf
for(i in 1:length(col7_vcf@fix[,1])){
  # extract exon numbers from vcf by splitting string
  deb_central_exons <- c(deb_central_exons, as.numeric(strsplit(strsplit(col7_vcf@fix[i,], "|", fixed = TRUE)$INFO[9], "/")[[1]][1]))
}
deb_counts <- count(deb_central_exons) #get exon frequency
# create plot
deb_plot <- ggplot(deb_counts, aes(x=x, y=freq)) + geom_bar(stat="identity") + labs(x="exon", y="number of mutations", title="Pathogenic mutations in deb-central")
# not all exons necessarily occur in the known mutations but to calculate correlation theyre needed
ns <- c(1:118) #118 exons
temp <- data.frame(ns[!(ns %in% deb_counts$x)]) #find which are missing
colnames(temp) <- "x" #match column name
temp$freq <- 0 #set frequency
all_deb_counts <- rbind(deb_counts, temp) #bind data frames
all_deb_counts <- all_deb_counts[order(all_deb_counts$x),] #order df
cor(all_deb_counts$x, colMeans(na.omit(all_samples_COL7_1ex)))#calculate correlation
##------------------------------------------------------------------------------------------##

##---------------------------------exon counts col7-----------------------------------------##
# retrieve the exon ranges used by recount - these are the "reduced" version
col7_ranges <- ranges(reduce(recount_exons$ENSG00000114270.16))

# reduced exon ranges combine some exons so labels for certain exons have to be combined
exon_labs <- c(1:118)
exon_labs[109:110] <- "109,110"
exon_labs[61:62] <- "61,62"
exon_labs[91:92] <- "91,92"
exon_labs[96:98] <- "96,97,98"
exon_labs <- unique(exon_labs)

# read the calculated exon counts
col7_exon_reads <- read.csv("/Users/david/Documents/data/cluster_output/col7_exons/col7_exon_reads.csv", header = FALSE)
col7_exon_reads$V1 <- NULL # remove sample IDs
col7_exon_reads <- as.data.frame(t(round(t(col7_exon_reads) / width(col7_ranges)))) #correct for exon length (longer exon = more reads)
colnames(col7_exon_reads) <- exon_labs #assign reduced exon labels
barData <- data.frame(rev(colSums(na.omit(col7_exon_reads)))) #construct ggplot compliant data
colnames(barData) <-  "dat"
barData$exon <-  factor(as.character(exon_labs), levels = unique(exon_labs)) #make sure each exon is numbered in plot
skippable_exons_ex <- c(skippable_exons, "109,110", "61,62", "91,92", "96,97,98") #add special labels to in/out-of-frame
barData$frame <- as.factor(ifelse(as.character(barData$exon) %in% skippable_exons_ex, "in frame", "out of frame"))#mark in-/out-of-frame

#create plot
exon_counts_plot <- ggplot(barData, aes(x=exon, y=dat, color=frame)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(size = 3.5, angle = 90, vjust = 0.5, hjust = 1)) + labs(x="exon", y="mean normalized read count", title="COL7A1 exon counts") + scale_color_manual(values=c("#5ab4ac", "#d8b365"))

# calculate means for the reduced exons to be able to calulate correlation with skip counts
col7_means <- rev(colMeans(na.omit(all_samples_n)))
col7_means[109:110] <- sum(col7_means[109:110])
col7_means[61:62] <- sum(col7_means[61:62])
col7_means[91:92] <- sum(col7_means[91:92])
col7_means[96:98] <- sum(col7_means[96:98])
col7_means_reduced <- col7_means[-c(110, 62, 92, 97,98)]

# calcualte correlation between exon counts and exon skip junction counts
cor(colMeans(na.omit(col7_exon_reads)), col7_means_reduced, method = "spearman")
##--------------------------------------------------------------------------------------------##

##----------------------------------Ordered bar chart-----------------------------------------##
# order exons by read count from high to low
ordered_exons <- rev(colMeans(na.omit(all_samples_COL7_1ex))[order(colMeans(na.omit(all_samples_COL7_1ex)))])

# construct ggplot complaint data
barData <- data.frame(unname(ordered_exons))
colnames(barData) <- "dat"
barData$exon <- factor(names(ordered_exons), levels = names(ordered_exons)) # make sure each exon is numbered as they wont be in order
# create plot
order_bar <- ggplot(data =barData, aes(x=exon, y=dat)) + geom_bar(stat = "identity", position = "dodge") #add data and type of plot
order_bar <- order_bar + theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1)) #Set x-axis text to be small enough to read exon numbers       
order_bar <- order_bar + labs(x="exon", y="mean normalized read count", title="COL7A1 exons ordered by exon skips") #add labels to plot
order_bar + geom_hline(yintercept = mean(colMeans(na.omit(all_samples_COL7_1ex))), linetype=2, color="red") #add line with mean skipping for COL7A1
##--------------------------------------------------------------------------------------------##

##---------------------------------Exon length vs read count----------------------------------##
#To see effect of exon length on number of junction reads
point_data <- data.frame(col7_exon_lengths) #collect lengths
colnames(point_data) <- "exon_lengths"
point_data$exon <- names(colSums(na.omit(all_samples_COL7_new_scale)))#collect exon numbers
point_data$read_count <- unname(colSums(na.omit(all_samples_COL7_new_scale)))#collect read counts
point_data$frame <- as.factor(ifelse(point_data$exon %in% skippable_exons, "in frame", "out of frame")) #add in-/out-of-frame data

#create plot to show correlation between exon length and junction count including pval of correlation
ggplot(point_data, aes(x=exon_lengths, y=read_count, color=frame)) + geom_point() + scale_color_manual(values=c("#5ab4ac", "#d8b365")) + labs(x="exon length", y="scaled read count", title="Number of reads for exon lengths") + annotate("text", x=300, y=max(point_data$read_count), label=paste0("pval: ", round(cor.test(point_data$read_count, point_data$exon_lengths)$p.val, 4)))
  
#calculate correlation between read count and exon length
cor(point_data$read_count, point_data$exon_lengths)
##------------------------------------------------------------------------------------------##

##------------------------------Average skipping of EB genes--------------------------------##
read_multi <- function(fileDir){
  ## function that reads all .csv files in a directory and concatenates them into 1 dataframe
  # read the files into a list of data.frames
  data.list <- lapply(list.files(path = fileDir, full.names = TRUE), function(x) read.csv(x, row.names = NULL))
  
  #set colnames to simple numbers to ensure tables can be merged
  for(i in 1:length(data.list)){
    colnames(data.list[[i]]) <- 1:length(colnames(data.list[[i]]))
  }
  
  # concatenate into one big data.frame
  data.cat <- as.data.frame(do.call(rbind, data.list))
  return(data.cat)
}
all_dir <- list.dirs("/Users/david/Documents/data/cluster_output/eb/") #all EB gene data saved in seperate folders here

eb_data <- list()
#loop through directories and exclude the starting directory
for(dir in all_dir){
  if(dir != "/Users/david/Documents/data/cluster_output/eb/"){
    eb_data[[dir]] <- read_multi(dir)#call concatenating function for each directory
  }
  
}

means <- c()
mean_per_gene <- c()
#loops through the list of EB genes
for(tbl in eb_data){
  colmns <- colMeans(na.omit(tbl[-1])) #calculate columns means. exclude first columns as it contains sample IDs
  # colmns[1] <- 0
  # colmns[length(colmns)] <- 0
  colmns[colmns == Inf] <- 0 #if any infinities are produced during calculations (due to NAs) set value to 0
  means <- c(means, colmns) #add column means to vector of all means
  mean_per_gene <- c(mean_per_gene, mean(colmns)) #add mean of all exons of entire mean to vector
}
means <- c(means, colMeans(na.omit(all_samples_COL7_1ex))) #add means of COL7a1
mean_per_gene <- c(mean_per_gene, mean(colMeans(na.omit(all_samples_COL7_1ex)))) #same

#create exon skipping plot including a horizontal line showing average skipping of EB genes
b_all + geom_hline(yintercept = median(means), linetype=2, color="red") + labs(title="Average exon skipping in EB genes")

#create table of all EB genes and their average skipping (for in report)
eb_table <- data.frame(mean_per_gene)
eb_table$symbol <- c(unlist(lapply(all_dir[-1], function(x) strsplit(x, "//", fixed = TRUE)[[1]][2])), "COL7A1")
eb_table$`mean skipped` <- eb_table$mean_per_gene
eb_table$mean_per_gene <- c()
##------------------------------------------------------------------------------------------##