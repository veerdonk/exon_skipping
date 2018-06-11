##----------------------------------------libraries-----------------------------------------------------##
## Loading libraries
library("recount")
library("GenomicRanges")
library("GenomicFeatures")
library("ggplot2")
library("gridExtra")
library("plyr")
library(data.table)
library(DESeq2)
##------------------------------------------------------------------------------------------------------##
# "SRP050971" < skin study
study <-  "SRP067502"#"SRP065812"#"DRP000366"#
#col7_id <- "ENSG00000114270.16"
gene_id <- "ENSG00000185842.14"
ucsc_id <- "uc001how.3"
dnah14_id <- "ENSG00000185842.14"
#18 exons ENSG00000118900.14
#38 exons ENSG00000101605.12
#84 exons ENSG00000185842.14

all_studies <- abstract_search("",id_only = TRUE)
##-------------------------------------------COL7 exon data---------------------------------------------##
gene_exon_annotation <- function(){
  # loading ucsc exon annotation for COL7A1
  exon_bed <-  read.delim("/Users/david/Documents/data/dnah14_exons.bed", header = FALSE, stringsAsFactors = FALSE)
  # ucsc id
  #ucsc_col7_id <- "uc003ctz.3"
  colnames(exon_bed) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
  # split the 'info' column to get the ucsc ids and select the exons for COL7A1
  exon_bed$ucsc_id <- unlist(lapply(exon_bed$info, function(x) strsplit(x, "_")[[1]][1]))
  gene_exons <- exon_bed[exon_bed$ucsc_id == ucsc_id,]
  row.names(gene_exons) <- 1:nrow(gene_exons)
  return(gene_exons)
}
##------------------------------------------------------------------------------------------------------##

##-------------------------------------download/load study----------------------------------------------##
download_study_jx <- function(study){
  if(!file.exists(file.path(study, "rse_jx.Rdata"))){
    download_study(study, type = "rse-jx")
  }
  ## check and load data
  if(file.exists(file.path(study, "rse_jx.Rdata"))){
    load(file.path(study, "rse_jx.Rdata"))
  }
  #file.remove(file.path(study, "rse_jx.Rdata"))
  #file.remove(file.path(study))
  return(rse_jx)
}

download_jx_bed <- function(study){
  url <- paste0("http://duffel.rail.bio/recount/", study, "/", study, ".junction_id_with_transcripts.bed.gz")
  jx_bed_gz_loc <- paste0("/Users/david/Documents/data/studies/", study, ".junction_id_with_transcripts.bed.gz")
  if(!file.exists(file.path(jx_bed_gz_loc))){
    download.file(url, destfile = jx_bed_gz_loc)
  }
  
  jx_bed <- fread(input = paste0("zcat < " , jx_bed_gz_loc), header = FALSE, stringsAsFactors = FALSE, sep = "\t", colClasses = c("character", "numeric", "numeric","character", "numeric", "character"))
  colnames(jx_bed) <- c("chromosome", "start", "stop", "meta", "1000", "strand")
  jx_bed$id <- as.numeric(unlist(lapply(as.character(jx_bed$meta), function(x) strsplit(x, "|", fixed = TRUE)[[1]][1])))
  #file.remove(jx_bed_gz_loc)
  return(jx_bed)
}


##------------------------------------------------------------------------------------------------------##

prepare_recount_study <- function(rse_jx){
  
  rse <- scale_counts(rse_jx, by = 'mapped_reads', round = FALSE) #create a scaled rse
  #rse <- read_counts(rse_jx) #read counts for normalizing with TPM
  #rse <- rse_jx
  
  #normalize using DESeq2
  # dds <- DESeqDataSet(rse_jx, design = ~ 1)
  # dds <- estimateSizeFactors( dds )
  # sizeFactors(dds)
  # junction_counts <- log2( counts(dds, normalized=TRUE) + 1 ) # counts table for DESeq2 normalization
  # 
  junction_counts <- assays(rse)$counts #retrieve counts table
  single_id <- lapply(rowData(rse)$gene_id_proposed, function(x) x[1]) # select first annotated ID for each junction
  annotatedRows <- which(!is.na(unlist(single_id))) # select rows annotated with an ID
  annotated_junctions <- junction_counts[annotatedRows,] # select annotated rows from count table
  counts_with_ids <- cbind.data.frame(annotated_junctions, unlist(single_id[annotatedRows])) # add ID to table
  annotated_jx_ids <- rowData(rse)$junction_id[annotatedRows] # find the junction ids with annotated gene ids
  counts_with_ids <- cbind.data.frame(counts_with_ids, annotated_jx_ids) # add junction ids to counts table
  colnames(counts_with_ids) <- c(colnames(counts_with_ids)[1:(ncol(counts_with_ids)-2)], "gene", "junction_id") # assign column names
  
  # normalize to TPM
  # genes_in_recount <- as.data.frame(recount_genes)
  # rpks <- t(apply(counts_with_ids, 1, function(x) as.numeric(x[1:(length(x)-2)])/(genes_in_recount[genes_in_recount$gene_id == x["gene"], 7]/1000)))
  # tpm <- rpks/(colSums(rpks)/1e6)
  # 
  # normalized_with_id <- cbind.data.frame(tpm, counts_with_ids$gene, counts_with_ids$junction_id)
  # colnames(normalized_with_id) <- colnames(counts_with_ids)
  # 
  
  
  
  #change normalized_with_id to counts_with_ids to not use TPM normalized data.
  gene_jx <- as.data.frame(counts_with_ids[which(unlist(single_id[annotatedRows]) == gene_id),], stringsAsFactors = F)
  gene_jx$junction_id <- as.numeric(levels(gene_jx$junction_id))[gene_jx$junction_id]
  
  return(gene_jx)
}

find_skipped_exons <- function(study){
  ptm <- proc.time()
  print("Downloading study..")
  rse_jx <- download_study_jx(study)
  print("Done.")
  print("Downloading .bed files..")
  jx_bed <- download_jx_bed(study)
  print("Done.")
  print("Retrieving COL7A1 exons..")
  gene_exons <- gene_exon_annotation()
  print("Done.")
  print("Preparing and scaling..")
  gene_jx <- prepare_recount_study(rse_jx)
  print("Done.")
  print("Searching for skipped exons..")
  
  
  # retrieving all start and stop positions for col7a1
  start_stop <- as.data.frame(t(sapply(as.integer(gene_jx$junction_id), function(x) jx_bed[jx_bed$id == x,1:3])))
  start_stop <- as.data.frame(apply(start_stop, 2, unlist), stringsAsFactors = FALSE)
  
  if(length(start_stop) == 1){
    start_stop <- t(start_stop)
  }
  
  # binding chromosome + start/stop to col7 data frame
  gene_start_stop <- cbind(gene_jx, start_stop)
  
  gene_start_stop$start <- as.numeric(gene_start_stop$start)
  gene_start_stop$stop <- as.numeric(gene_start_stop$stop)
  print("pass")  
  alternative_length <- 1
  skipped_exons <- list()
  for(i in 1:length(gene_start_stop$start)){
    for(j in 1:nrow(gene_exons)){
      #if(col7_start_stop$start[i] >= col7_exons$start[1] && as.numeric(col7_start_stop$stop[i]) <= col7_exons$stop[118]){
      if((gene_exons$start[j] - alternative_length) > gene_start_stop$start[i] & (gene_exons$stop[j] + alternative_length) < gene_start_stop$stop[i]){
        if(exists(gene_exons$info[j], where = skipped_exons)){
          val <- skipped_exons[[gene_exons$info[j]]]
          skipped_exons[[gene_exons$info[j]]] <- c(val, gene_start_stop$junction_id[i])
        }
        else{
          skipped_exons[[gene_exons$info[j]]] <- gene_start_stop$junction_id[i]
        }
      }
      #}  
    }
  }
  
  # no_skips <- data.frame(count(unlist(unname(skipped_exons))))
  # no_skips_o <- no_skips[order(count(unlist(unname(skipped_exons)))$freq),]
  # no_skips_o$x <- as.factor(no_skips_o$x)
  # plot(sort(no_skips$freq), ylab="number of skips", xlab = "junction id (numbered)")
  
  good_jx_ids <- count(unlist(unname(skipped_exons)))[count(unlist(unname(skipped_exons)))$freq < 10,1]
  #jx_frequency <- count(unlist(unname(skipped_exons)))
  #freq_binned <- count(jx_frequency$freq)
  #p2 <- ggplot(freq_binned, aes(x=x, y=freq)) + geom_bar(stat='identity') + labs(x="number of skipped exons", y="number of junctions", title="Skips in COL7A1")
  #counts_per_skipped_exon <- sapply(names(skipped_exons), function(x) find_read_counts(skipped_exons[[x]], col7_jx))
  #count_table_skipped <- as.data.frame(t(sapply(names(skipped_exons), function(x) find_counts_table(skipped_exons[[x]], col7_jx, good_jx_ids))))
  
  # transformed cpount table to accomodate more samples
  
  count_table_skipped <- as.data.frame(sapply(gene_exons$info, function(x) find_counts_table(skipped_exons[[x]], gene_jx, good_jx_ids)))
  if(length(count_table_skipped) == 1){
    count_table_skipped <- as.data.frame(t(count_table_skipped))
    row.names(count_table_skipped) <- "sample1"
  }
  print(paste0("Done, time taken to process study: ", round((proc.time()-ptm)[3], 2), " S"))
  return(count_table_skipped)
  #return(counts_per_skipped_exon)
}

find_counts_table <- function(jx_ids, gene_jx, good_jx_ids){
  
  tempdf <- data.frame()
  for(id in jx_ids){
    if(id %in% good_jx_ids){
      tempdf <- rbind(tempdf, gene_jx[col7_jx$junction_id == id,1:(ncol(gene_jx)-2)])
    }
  }
  if(length(tempdf) == 0){
    
    tempdf <- as.data.frame(t(rep(0, (ncol(gene_jx)-2))))
    colnames(tempdf) <- colnames(gene_jx)[1:(ncol(gene_jx)-2)]
  }
  
  return(colSums(tempdf))
}

find_read_counts <- function(jx_ids, gene_jx){
  
  total_jx_rowcount <- 0
  for (id in jx_ids) {
    total_jx_rowcount <- total_jx_rowcount + rowSums(gene_jx[dnah14_jx$junction_id == id,1:(ncol(gene_jx)-2)])
  }
  
  return(total_jx_rowcount)
}

plot_preliminary <- function(skipped_exon_counts){
  ## TODO: make this compatible with a count table instead of rowSums
  
  
  # without sorting
  skipped_exon_counts_df <- as.data.frame(skipped_exon_counts)
  colnames(skipped_exon_counts_df) <- "junction_rowsum"
  skipped_exon_counts_df$exon <- as.numeric(sapply(row.names(skipped_exon_counts_df), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))
  not_skip_exons <- c(0,1,2,3,5,6,23,24,26,112,117) # exons that cant be skipped (0 based)
  skippable_exons <- (0:117)[-not_skip_exons] # exons that can be skipped (0 based)
  skipped_exon_counts_df$skippable <- as.factor(ifelse(skipped_exon_counts_df$exon %in% skippable_exons, "skippable", "not skippable"))
  p <-  ggplot(skipped_exon_counts_df, aes(x=exon, y=junction_rowsum, color = skippable))# + geom_point() + labs(x = "exon", y = "scaled read count", title="skipped exons")
  #q <- ggplot(skipped_exon_counts_df, aes)
  #
  
  #with sorting
  # skipped_ordered <- as.data.frame(skipped_exon_counts[order(skipped_exon_counts, decreasing = TRUE)])
  # colnames(skipped_ordered) <- "junction_rowsum"
  # skipped_ordered$exon <- as.numeric(sapply(row.names(skipped_ordered), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))
  # 
  # not_skip_exons <- c(0,1,2,3,5,6,23,24,26,112,117) # exons that cant be skipped (0 based)
  # skippable_exons <- (0:117)[-not_skip_exons] # exons that can be skipped (0 based)
  # skipped_ordered$skippable <- as.factor(ifelse(skipped_ordered$exon %in% skippable_exons, "skippable", "not skippable"))
  # 
  # skipped_ordered_no0 <- skipped_ordered[-1,]# remove to inlude exon 1
  # 
  # p <-  ggplot(skipped_ordered_no0, aes(x=exon, y=junction_rowsum, color = skippable)) + geom_point() + labs(x = "exon", y = "scaled read count", title="skipped exons")
  #
  return(p)
}

plot_count_table <- function(skipped_exon_count_table){
  row.names(skipped_exon_count_table) <- as.numeric(sapply(row.names(skipped_exon_count_table), function(x) strsplit(x, "_", fixed = TRUE)[[1]][3]))+1
  skipped_exon_counts_mlt <- melt(t(skipped_exon_count_table))
  not_skip_exons <- c(1,2,3,5,6,23,24,26,112,118) # exons that cant be skipped (1 based)
  skippable_exons <- (1:118)[-not_skip_exons] # exons that can be skipped (1 based)
  skipped_exon_counts_mlt$skippable <- as.factor(ifelse(skipped_exon_counts_mlt$Var2 %in% skippable_exons, "skippable", "not skippable"))
  
  p <- ggplot(skipped_exon_counts_mlt, aes(x=Var2, y=value, group=Var2, color=skippable)) + geom_boxplot()
  
  return(p)
}

skipped <- find_skipped_exons(study)
