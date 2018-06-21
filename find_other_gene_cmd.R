#!/usr/bin/env Rscript

##-----------------------------------------Args---------------------------------------------------------##
args <- commandArgs(trailingOnly = TRUE)
##------------------------------------------------------------------------------------------------------##

##----------------------------------------libraries-----------------------------------------------------##
## Loading libraries
library("recount")
library("plyr")
library(data.table)
##------------------------------------------------------------------------------------------------------##

##-----------------------------------------Constants----------------------------------------------------##
gene_id <- "ENSG00000114270.16" #COL7
ucsc_gene_id <- "uc003ctz.3" #COL7
all_studies <- abstract_search("",id_only = TRUE)
##------------------------------------------------------------------------------------------------------##

##-------------------------------------------gene exon data---------------------------------------------##
gene_exon_annotation <- function(exons_bed){
  # loading ucsc exon annotation for geneA1
  #"/Users/david/Documents/data/gene_exons.bed"
  exon_bed <-  read.delim(exons_bed, header = FALSE, stringsAsFactors = FALSE)
  # ucsc id
  
  colnames(exon_bed) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
  # split the 'info' column to get the ucsc ids and select the exons for gene
  exon_bed$ucsc_id <- unlist(lapply(exon_bed$info, function(x) strsplit(x, "_")[[1]][1]))
  gene_exons <- exon_bed[exon_bed$ucsc_id == ucsc_gene_id,]
  if(gene_exons$strand[1] == "-"){
  gene_exons$info <- rev(gene_exons$info)
  }
  rm(exon_bed)
  #row.names(gene_exons) <- 1:nrow(gene_exons)
  return(gene_exons)
}
##------------------------------------------------------------------------------------------------------##

##-------------------------------------download/load study----------------------------------------------##
download_study_jx <- function(study){
  if(!file.exists(paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study, "/rse_jx.Rdata"))){
    download_study(study, type = "rse-jx", outdir = paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study))
    print("rse downloaded..")
  }
  
  ## check and load data
  if(file.exists(paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study, "/rse_jx.Rdata"))){
    if(file.info(paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study, "/rse_jx.Rdata"))$size > 62914560){
      stop("file too large")
    }
    load(file.path(paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study, "/rse_jx.Rdata")))
    print("rse loaded..")
  }
  
  return(rse_jx)
}

download_jx_bed <- function(study){
  #construct url and download location
  url <- paste0("http://duffel.rail.bio/recount/", study, "/", study, ".junction_id_with_transcripts.bed.gz")
  jx_bed_gz_loc <- paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study, "/", study, ".junction_id_with_transcripts.bed.gz")
  if(!file.exists(jx_bed_gz_loc)){
    download.file(url, destfile = jx_bed_gz_loc)
  }
  #read the .bed file and assign column names
  jx_bed <- fread(input = paste0("zcat < " , jx_bed_gz_loc), header = FALSE, stringsAsFactors = FALSE, sep = "\t", colClasses = c("character", "numeric", "numeric","character", "numeric", "character"))
  colnames(jx_bed) <- c("chromosome", "start", "stop", "meta", "1000", "strand")
  #assign ids
  jx_bed$id <- as.numeric(unlist(lapply(as.character(jx_bed$meta), function(x) strsplit(x, "|", fixed = TRUE)[[1]][1])))
  
  return(jx_bed)
}
##------------------------------------------------------------------------------------------------------##


##----------------------------------------prepare/scale-------------------------------------------------##
prepare_recount_study <- function(rse_jx){
  # scale the rse by the number of mapped reads
  rse <- scale_counts(rse_jx, by = 'mapped_reads', round = FALSE) #create a scaled rse
  rm(rse_jx)
  junction_counts <- assays(rse)$counts #retrieve counts table
  single_id <- lapply(rowData(rse)$gene_id_proposed, function(x) x[1]) # select first annotated ID for each junction
  annotatedRows <- which(!is.na(unlist(single_id))) # select rows annotated with an ID
  annotated_junctions <- junction_counts[annotatedRows,] # select annotated rows from count table
  rm(junction_counts)
  counts_with_ids <- cbind.data.frame(annotated_junctions, unlist(single_id[annotatedRows])) # add ID to table
  annotated_jx_ids <- rowData(rse)$junction_id[annotatedRows] # find the junction ids with annotated gene ids
  rm(rse)
  gc()
  counts_with_ids <- cbind.data.frame(counts_with_ids, annotated_jx_ids) # add junction ids to counts table
  colnames(counts_with_ids) <- c(colnames(counts_with_ids)[1:(ncol(counts_with_ids)-2)], "gene", "junction_id") # assign column names
  
  #extract the junctions mapping to gene exons
  gene_jx <- as.data.frame(counts_with_ids[which(unlist(single_id[annotatedRows]) == gene_id),], stringsAsFactors = F)
  rm(counts_with_ids)
  gene_jx$junction_id <- as.numeric(levels(gene_jx$junction_id))[gene_jx$junction_id]
  rm(annotated_junctions)
  rm(single_id)
  
  return(gene_jx)
}
##------------------------------------------------------------------------------------------------------##

##---------------------------------controller/finding skips---------------------------------------------##
find_skipped_exons <- function(study, exons_bed){
  #run functions and print/time steps
  ptm <- proc.time()
  print("Downloading study..")
  rse_jx <- download_study_jx(study)
  gc()
  print("Done.")
  print("Downloading .bed files..")
  jx_bed <- download_jx_bed(study)
  gc()
  print("Done.")
  print("Retrieving gene exons..")
  gene_exons <- gene_exon_annotation(exons_bed)
  gc()
  print("Done.")
  print("Preparing and scaling..")
  gene_jx <- prepare_recount_study(rse_jx)
  gc()
  print("Done.")
  print("Searching for skipped exons..")
  
  # if no junction reads mapped to gene break operation
  if(length(gene_jx[,1]) == 0){
    return(FALSE)
  }
  
  #write.table(colSums(gene_jx[1:(length(gene_jx)-2)] > 0), file = "/groups/umcg-gcc//tmp03/umcg-dvdveerdonk/number_of_junctions.csv", append = TRUE, col.names = FALSE, sep = ",")
  
  # retrieving all start and stop positions for genea1
  start_stop <- as.data.frame(t(sapply(as.integer(gene_jx$junction_id), function(x) jx_bed[jx_bed$id == x,1:3])))
  start_stop <- as.data.frame(apply(start_stop, 2, unlist), stringsAsFactors = FALSE)
  
  # if only one junction is found transpose df so its usable
  if(length(start_stop) == 1){
    start_stop <- t(start_stop)
  }
  
  # binding chromosome + start/stop to gene data frame
  gene_start_stop <- cbind(gene_jx, start_stop)
  
  #converting factors to numeric
  gene_start_stop$start <- as.numeric(gene_start_stop$start)
  gene_start_stop$stop <- as.numeric(gene_start_stop$stop)
  
  # recount positions are off by 1
  alternative_length <- 1
  skipped_exons <- list() 
  # loop through all junctions and check whether an exon falls completely within its start/stop
  for(i in 1:length(gene_start_stop$start)){
    for(j in 1:nrow(gene_exons)){
      #if(gene_start_stop$start[i] >= gene_exons$start[1] && as.numeric(gene_start_stop$stop[i]) <= gene_exons$stop[118]){
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
  # check that junction does not 'skip' more than 10 exons (might be changed)
  good_jx_ids <- plyr::count(unlist(unname(skipped_exons)))[plyr::count(unlist(unname(skipped_exons)))$freq < 10,1]
  
  # transformed count table to accomodate more samples
  count_table_skipped <- as.data.frame(sapply(gene_exons$info, function(x) find_counts_table(skipped_exons[[x]], gene_jx, good_jx_ids)))
  
  # if only one sample is present the table needs to be transposed and given a rowname 
  if(length(count_table_skipped) == 1){
    count_table_skipped <- as.data.frame(t(count_table_skipped))
    row.names(count_table_skipped) <- study
  }
  # report the time taken to process the study
  print(paste0("Done, time taken to process study: ", round((proc.time()-ptm)[3], 2), " S"))
  print("")
  return(count_table_skipped)
  #return(counts_per_skipped_exon)
}
##------------------------------------------------------------------------------------------------------##

##-------------------------------------jx ids to count table--------------------------------------------##
find_counts_table <- function(jx_ids, gene_jx, good_jx_ids){
  #finds all counts for the junctions mapping to a single exon
  tempdf <- data.frame()
  for(id in jx_ids){
    if(id %in% good_jx_ids){
      tempdf <- rbind(tempdf, gene_jx[gene_jx$junction_id == id,1:(ncol(gene_jx)-2)])
    }
  }
  if(length(tempdf) == 0){
    
    tempdf <- as.data.frame(t(rep(0, (ncol(gene_jx)-2))))
    colnames(tempdf) <- colnames(gene_jx)[1:(ncol(gene_jx)-2)]
  }
  
  return(colSums(tempdf))
}
##------------------------------------------------------------------------------------------------------##

##------------------------------------loop through cmd args and run-------------------------------------##
# args <- "SRP065812"
# file ~/Documents/data/skipped/skipped_exons_multifile.csv
not_run_studies <- c()

time1 <- proc.time()
studies <- all_studies[args[3]:args[4]] # what studies have to be run
outfile <- args[1] # first argument should be desired output
exons_bed <- args[2] # second argument should be the location of the gene .bed file

io_error <- c()
no_gene <- c()
for(study in studies){
  if(!study %in% not_run_studies){
    print(study)
    skipped_exon_counts <- tryCatch(
      {find_skipped_exons(study, exons_bed)},
      error=function(e){
        io_error <- c(io_error, study)
        print(paste0(study, " produced an error.."))
      }
    )
    gc() # garbage collection to free up memory 
    if(length(skipped_exon_counts) > 1){
      if(!file.exists(outfile)){
        write.table(skipped_exon_counts, file = outfile, sep = ",", col.names = TRUE)
      }
      else{
        write.table(skipped_exon_counts, file = outfile, sep = ",", col.names = FALSE, append = TRUE)
      }
    }else{
      print("No data for geneA1 available, stopping..")
      no_gene <- c(no_gene, study)
    }
  }
}
print("Failed studies (IO):")
print(paste(io_error))
print("Failed studies (no gene):")
print(paste(no_gene))
print(proc.time() - time1)

##------------------------------------------------------------------------------------------------------##
