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
gene_id <- args[5]
ucsc_gene_id <- args[6]
all_studies <- abstract_search("",id_only = TRUE)
study_dir <- "/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/" #directory studies should be saved
##------------------------------------------------------------------------------------------------------##

##-------------------------------------------gene exon data---------------------------------------------##
gene_exon_annotation <- function(exons_bed){
  # loading ucsc exon annotation for gene
  exon_bed <-  read.delim(exons_bed, header = FALSE, stringsAsFactors = FALSE)
  colnames(exon_bed) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
  # split the 'info' column to get the ucsc ids and select the exons for gene
  exon_bed$ucsc_id <- unlist(lapply(exon_bed$info, function(x) strsplit(x, "_")[[1]][1])) #extract ids from string
  gene_exons <- exon_bed[exon_bed$ucsc_id == ucsc_gene_id,] #find the right lines
  if(gene_exons$strand[1] == "-"){
    gene_exons$info <- rev(gene_exons$info)
  }
  rm(exon_bed)
  return(gene_exons)
}
##------------------------------------------------------------------------------------------------------##

##-------------------------------------download/load study----------------------------------------------##
download_study_jx <- function(study){
  # download study if it is not already downloaded
  if(!file.exists(paste0(study_dir, study, "/rse_jx.Rdata"))){
    download_study(study, type = "rse-jx", outdir = paste0(study_dir, study))
    print("rse downloaded..")
  }
  
  ## check and load data
  if(file.exists(paste0(study_dir, study, "/rse_jx.Rdata"))){
    # If memory problems arise uncomment to limit allowed study size. NOTE: process excluded studies seperately
    #if(file.info(paste0(study_dir, study, "/rse_jx.Rdata"))$size > 62914560){
    #  stop("file too large")
    #}
    load(file.path(paste0(study_dir, study, "/rse_jx.Rdata")))
    print("rse loaded..")
  }
  
  return(rse_jx)
}

download_jx_bed <- function(study){
  #construct url and download location
  url <- paste0("http://duffel.rail.bio/recount/", study, "/", study, ".junction_id_with_transcripts.bed.gz") # construct download url
  jx_bed_gz_loc <- paste0(study_dir, study, "/", study, ".junction_id_with_transcripts.bed.gz") # construct drive location
  #download and save .bed if it is not found at location
  if(!file.exists(jx_bed_gz_loc)){
    download.file(url, destfile = jx_bed_gz_loc)
  }
  #read the .bed file and assign column names (bed can be large so fread is used)
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
  
  ## other way of scaling read count, this scaled to a total coverage of 40mil
  ## when using this uncomment block and comment scaling line above.
  #rse <- rse_jx
  #junction_counts <- assays(rse_jx)$counts
  #junction_counts <- t((t(junction_counts) / (colData(rse_jx)$mapped_read_count/1e6))*40)#scale all read counts to total coverage of 40mil reads.
  
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
  gene_jx$junction_id <- as.numeric(levels(gene_jx$junction_id))[gene_jx$junction_id] #assign the right junction IDs as numbers
  rm(annotated_junctions)
  rm(single_id)
  
  return(gene_jx)
}
##------------------------------------------------------------------------------------------------------##

##---------------------------------controller/finding skips---------------------------------------------##
find_skipped_exons <- function(study, exons_bed){
  #run functions and print/time steps for monitoring
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
    print("No junction reads map to chosen gene")
    return(FALSE)
  }
  
  ### code to get average read length to find false positives
  avg_read_length <- mean(colData(rse_jx)$avg_read_length)
  avg_read_sd <- sd(colData(rse_jx)$avg_read_length)
  possible_false_positives <- c()
  print(paste0("read length = ", avg_read_length))
  ###
  
  ## below line saves a different file in shich the total number of junctions for chosen gene are saved
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
  
  ## Complicated bit, this is the actual finding of the exons
  ## checks whether any junctions completely encompass an exon
  
  #loop through all junctions
  for(i in 1:length(gene_start_stop$start)){
    #for each junction loop through exons
    for(j in 1:nrow(gene_exons)){
      # check the start and stop positions of each exon against the start/stop of the junction
      if((gene_exons$start[j] - alternative_length) > gene_start_stop$start[i] & (gene_exons$stop[j] + alternative_length) < gene_start_stop$stop[i]){
        ## for finding possible false positives the start/stop positions and the average read length is used
        if(j >= 2){
          if(gene_start_stop$start[i]+((gene_exons$start[j+1]-gene_exons$stop[j])+(gene_exons$start[j]-gene_exons$stop[(j-1)])) + (avg_read_length+2*avg_read_sd) <= gene_start_stop$stop[i]){
            possible_false_positives <- c(possible_false_positives, gene_start_stop$junction_id[i])
        }}
        # if the exon is already in the list add the junction marking it as skipped to it
        if(exists(gene_exons$info[j], where = skipped_exons)){
          val <- skipped_exons[[gene_exons$info[j]]]
          skipped_exons[[gene_exons$info[j]]] <- c(val, gene_start_stop$junction_id[i])
        }
        # if exon isnt in the list make a new entry
        else{
          skipped_exons[[gene_exons$info[j]]] <- gene_start_stop$junction_id[i]
        }
      }  
    }
  }
  # some junctions mark many exons as skipped, this limits it to one.
  good_jx_ids <- count(unlist(unname(skipped_exons)))[count(unlist(unname(skipped_exons)))$freq == 1,1]
  
  ## uncomment to remove any reads that might be a false positive
  #good_jx_ids <- good_jx_ids[!good_jx_ids %in% possible_false_positives]

  # retrieve info from counts table for all skipped exons
  count_table_skipped <- as.data.frame(sapply(gene_exons$info, function(x) find_counts_table(skipped_exons[[x]], gene_jx, good_jx_ids)))
  
  # if only one sample is present the table needs to be transposed and given a rowname 
  if(length(count_table_skipped) == 1){
    count_table_skipped <- as.data.frame(t(count_table_skipped))
    row.names(count_table_skipped) <- study #study is used as name
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
  # loops through junction ids its given
  for(id in jx_ids){
    # checks whether it should be kept
    if(id %in% good_jx_ids){
      #retrieves junction count information and binds to temporary df
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
# exclude studies if wanted
not_run_studies <- c()
# start times
time1 <- proc.time()
studies <- all_studies[args[3]:args[4]] # what studies have to be run
outfile <- args[1] # first argument should be desired output
exons_bed <- args[2] # second argument should be the location of the gene .bed file

io_error <- c() #printed at the end to see what studies fail
no_gene <- c() # same but for studies not containing wanted gene
# loop trough studies given om command line
for(study in studies){
  if(!study %in% not_run_studies){
    print(study)
    # catch IO errors
    skipped_exon_counts <- tryCatch(
      {find_skipped_exons(study, exons_bed)}, # run skip finding functions
      error=function(e){
        io_error <- c(io_error, study)
        print(paste0(study, " produced an error.."))
      }
    )
    gc() # garbage collection to free up memory 
    if(length(skipped_exon_counts) > 1){ # if data is produced
      # write to file, create new if it doesnt exist, otherwise append
      if(!file.exists(outfile)){
        write.table(skipped_exon_counts, file = outfile, sep = ",", col.names = TRUE)
      }
      else{
        write.table(skipped_exon_counts, file = outfile, sep = ",", col.names = FALSE, append = TRUE)
      }
    }else{
      print("No data for gene available, stopping..")
      no_gene <- c(no_gene, study)
    }
  }
}
## print all the failed studies so they can be examined
print("Failed studies (IO):")
print(paste(io_error))
print("Failed studies (no gene):")
print(paste(no_gene))
print(proc.time() - time1)

##------------------------------------------------------------------------------------------------------##
