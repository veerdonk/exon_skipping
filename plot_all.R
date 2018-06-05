library("reshape2")
library("ggplot2")
fileDir <- "/Users/david/Documents/data/cluster_output/"

# read the files into a list of data.frames
data.list <- lapply(list.files(path = fileDir, full.names = TRUE), function(x) read.csv(x, row.names = NULL))

colnames(data.list[[1]]) <- colnames(data.list[[2]])
# concatenate into one big data.frame
data.cat <- as.data.frame(do.call(rbind, data.list))
all_samples <- subset(data.cat,!duplicated(data.cat$row.names))
row.names(all_samples) <- all_samples$row.names
all_samples$row.names <- c()

#correct for exon length:
all_samples_n <- t(t(all_samples) / col7_exon_lengths)

#row.names(data.cat) <- as.numeric(sapply(row.names(data.cat), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
colnames(all_samples) <- as.numeric(sapply(colnames(all_samples), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
skipped_exon_counts_mlt <- melt(t(all_samples))
non_zero_mlt <- skipped_exon_counts_mlt[skipped_exon_counts_mlt$value > 0,]
not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
non_zero_mlt$skippable <- as.factor(ifelse(non_zero_mlt$Var1 %in% skippable_exons, "skippable", "not skippable"))
non_zero_mlt <- non_zero_mlt[!is.na(non_zero_mlt$value),]
non_zero_mlt <- non_zero_mlt[!is.infinite(non_zero_mlt$value),]
all_sample_nn <- ggplot(non_zero_mlt, aes(x=Var1, y=value, group=Var1, color = skippable)) + labs(x="exon", y="read count", title="All samples* (not normalized)") + geom_boxplot() 

#normalized
colnames(all_samples_n) <- as.numeric(sapply(colnames(all_samples_n), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
skipped_exon_counts_mlt <- melt(t(all_samples_n))
non_zero_mlt <- skipped_exon_counts_mlt[skipped_exon_counts_mlt$value > 0,]
not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
non_zero_mlt$skippable <- as.factor(ifelse(non_zero_mlt$Var1 %in% skippable_exons, "skippable", "not skippable"))
non_zero_mlt <- non_zero_mlt[!is.na(non_zero_mlt$value),]
non_zero_mlt <- non_zero_mlt[!is.infinite(non_zero_mlt$value),]
all_sample_n <- ggplot(non_zero_mlt, aes(x=Var1, y=value, group=Var1, color = skippable)) + labs(x="exon", y="read count", title="All samples* corrected for exon length") + geom_boxplot() 

bar_nn <- barplot(colMeans(na.omit(all_samples)), cex.names = 0.25)
bar_n <- barplot(colMeans(na.omit(all_samples_n)), cex.names = 0.25)

barData <- data.frame(colMeans(na.omit(all_samples)))
colnames(barData) <- "dat"
barData$exon <- as.numeric(row.names(barData))

means_nn <- ggplot(data=barData, aes(x=exon, y=dat)) + geom_bar(stat = "identity") + theme(axis.text=element_text(size=4)) + labs(x="exon", y="mean read count", title="means (not normalized)")

head(skipped_exon_counts_mlt)

skipped_exon_counts_mlt$Var2
