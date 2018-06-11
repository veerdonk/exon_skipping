library("reshape2")
library("ggplot2")
fileDir <- "/Users/david/Documents/data/cluster_output/"

# read the files into a list of data.frames
data.list <- lapply(list.files(path = fileDir, full.names = TRUE), function(x) read.csv(x, row.names = NULL))

colnames(data.list[[7]]) <- colnames(data.list[[2]])
# concatenate into one big data.frame
data.cat <- as.data.frame(do.call(rbind, data.list))
all_samples <- subset(data.cat,!duplicated(data.cat$row.names))
row.names(all_samples) <- all_samples$row.names
all_samples$row.names <- c()

#correct for exon length:
all_samples_n <- as.data.frame(t(t(all_samples) / col7_exon_lengths))

#get skin samples:
skin_data <- all_samples_n[row.names(all_samples_n) %in% skin_samples,]

#row.names(data.cat) <- as.numeric(sapply(row.names(data.cat), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
colnames(all_samples) <- as.numeric(sapply(colnames(all_samples), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
skipped_exon_counts_mlt <- melt(t(all_samples))
non_zero_mlt <- skipped_exon_counts_mlt[skipped_exon_counts_mlt$value > 0,]
not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
non_zero_mlt$skippable <- as.factor(ifelse(non_zero_mlt$Var1 %in% skippable_exons, "in frame", "out of frame"))
non_zero_mlt <- non_zero_mlt[!is.na(non_zero_mlt$value),]
non_zero_mlt <- non_zero_mlt[!is.infinite(non_zero_mlt$value),]
all_sample_nn <- ggplot(non_zero_mlt, aes(x=Var1, y=value, group=Var1, color = skippable)) + labs(x="exon", y="read count", title="All samples* (not normalized)") + geom_boxplot() 

#normalized
colnames(all_samples_n) <- as.numeric(sapply(colnames(all_samples_n), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
skipped_exon_counts_mlt <- melt(t(all_samples_n))
non_zero_mlt <- skipped_exon_counts_mlt[skipped_exon_counts_mlt$value > 0,]
not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
non_zero_mlt$skippable <- as.factor(ifelse(non_zero_mlt$Var1 %in% skippable_exons, "in frame", "out of frame"))
non_zero_mlt <- non_zero_mlt[!is.na(non_zero_mlt$value),]
non_zero_mlt <- non_zero_mlt[!is.infinite(non_zero_mlt$value),]
all_sample_n <- ggplot(non_zero_mlt, aes(x=Var1, y=value, group=Var1, color = skippable)) + labs(x="exon", y="read count", title="All samples* corrected for exon length") + geom_boxplot() 


barData <- data.frame(colMeans(na.omit(all_samples_n)))
colnames(barData) <- "dat"
barData$exon <- factor(as.character(row.names(barData)), levels = unique(row.names(barData)))
barData$frame <- as.factor(ifelse(barData$exon %in% skippable_exons, "in frame", "out of frame"))
ggplot(data=barData, aes(x=exon, y=dat, color=frame)) + geom_bar(stat = "identity") + theme(axis.text=element_text(size=4)) + labs(x="exon", y="mean read count", title="COL7 skipped exons") + scale_color_manual(values=c("#21f000", "#ff2d00"))

colsums <- data.frame(colSums(na.omit(all_samples_n)))
colnames(colsums) <- "val"
exon_nums <- 1:118
colsums$exon <- factor(as.character(exon_nums), levels = unique(exon_nums))
colsums$frame <- as.factor(ifelse(colsums$exon %in% skippable_exons, "in frame", "out of frame"))
ggplot(data = colsums, aes(y=val, x=exon, color=frame)) + geom_bar(stat = "identity") + theme(axis.text=element_text(size=4)) + labs(x="exon", y="stacked counts", title="COL7 skipped exons") + scale_color_manual(values=c("#21f000", "#ff2d00"))


total_jx_col7 <- read.csv("~/Documents/data/number_of_junctions.csv", header = FALSE)

binned <- count(na.omit(total_jx_col7$V2))

ggplot(binned, aes(x=x, y=freq)) + geom_bar(stat='identity', width=1) + labs(title="Number of junctions in COL7A1", x="Number of junctions", y="Number of samples")
