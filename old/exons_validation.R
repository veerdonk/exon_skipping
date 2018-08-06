library("recount")
study <- "SRP050971"
col7_id <- "ENSG00000114270.16"

url <- paste0("http://duffel.rail.bio/recount/", "v2/", study, "/rse_exon.Rdata")
exon_rse_loc <- paste0("/Users/david/Documents/data/studies/", study, "rse_exon.Rdata")

if(!file.exists(exon_rse_loc)){
  download.file(url, destfile = exon_rse_loc)
}
## check and load data
if(file.exists(exon_rse_loc)){
  load(exon_rse_loc)
}

rse_exon <- scale_counts(rse_exon, by = "mapped_reads")
exon_counts <- assays(rse_exon)$counts
#disjoint exons of col7
col7_recount_exons <- as.data.frame(exon_counts[rownames(exon_counts) == "ENSG00000114270.16",]) 
rownames(col7_recount_exons) <- c(1:length(rownames(col7_recount_exons)))

col7_disjoined_start_stop <- as.data.frame(ranges(disjoin(recount_exons$ENSG00000114270.16)))

exon_bed <-  read.delim("/Users/david/Documents/data/col7_exons.bed", header = FALSE, stringsAsFactors = FALSE)
# ucsc id
ucsc_col7_id <- "uc003ctz.3"
colnames(exon_bed) <- c("chromosome", "start", "stop", "info", "zero", "strand")#column names
# split the 'info' column to get the ucsc ids and select the exons for COL7A1
exon_bed$ucsc_id <- unlist(lapply(exon_bed$info, function(x) strsplit(x, "_")[[1]][1]))
col7_exons <- exon_bed[exon_bed$ucsc_id == ucsc_col7_id,]
row.names(col7_exons) <- 1:nrow(col7_exons)

get_exon_rows <- function(rownumber){
  if(length(rownumber) > 1){
    return(colSums(col7_recount_exons[rownumber,]))
  }
  else{
    return(col7_recount_exons[rownumber,])
  }
}

exons <- data.frame()
counter <- 0
col7_good_exons_counts <- data.frame()
goodrows <- c()
doubles <- c()
for(i in 1:length(col7_exons$start)){
  start <- col7_exons$start[i]
  stop <- col7_exons$stop[i]
  exon <- c()
  for(j in 1:length(col7_disjoined_start_stop$start)){
    if((col7_disjoined_start_stop$start[j]-1) >= start & col7_disjoined_start_stop$end[j] <= stop){
      exon <- c(exon, j)
      #print(exon)
    }
  }
  if(length(exon) > 1){
    doubles <- c(doubles, exon)
    row <- get_exon_rows(exon)
    
    col7_good_exons_counts <- rbind(col7_good_exons_counts, row)
    colnames(col7_good_exons_counts) <- colnames(col7_recount_exons)
    
    exons[i,1] <- col7_disjoined_start_stop$start[exon[1]]
    exons[i,2] <- col7_disjoined_start_stop$end[tail(exon, n=1)]
    exons[i,3] <- paste(exon, collapse = ",")
  }
  else{
    print(exon)
    row <-  get_exon_rows(exon)
    col7_good_exons_counts <- rbind(col7_good_exons_counts, row)
    exons[i,1] <- col7_disjoined_start_stop$start[exon]
    exons[i,2] <- col7_disjoined_start_stop$end[exon]
    exons[i,3] <- paste(exon, collapse = ",")
    counter <- counter + 1
  }
}
row.names(col7_good_exons_counts) <- c(1:118)

plot(col7_good_exons_counts$SRR1698134)

#-----------------------------------------------------------------------------------------------------#
#---------------------------------correcting for exon length------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

col7_exon_lengths <- (col7_exons$stop - col7_exons$start)

divided <- col7_good_exons_counts/col7_exon_lengths
#niet goed...
#norm <- col7_good_exons_counts$SRR1698134 - mean(col7_good_exons_counts$SRR1698134)/(max(col7_good_exons_counts$SRR1698134 - min(col7_good_exons_counts$SRR1698134)))

melted_exons <- melt(t(divided))

melted_exons$skippable <- as.factor(ifelse(melted_exons$Var2 %in% (skippable_exons0+1), "skippable", "not skippable"))

# p_exons = boxplot van exon counts
p_exons <- ggplot(melted_exons, mapping = aes(x=Var2, y=value, group=Var2, color=skippable)) + geom_boxplot()

# niet op lengte gecorrigeerd:
melted_exons_nn <- melt(t(col7_good_exons_counts))

melted_exons_nn$skippable <- as.factor(ifelse(melted_exons$Var2 %in% (skippable_exons0+1), "skippable", "not skippable"))

p_exon_nn <- ggplot(melted_exons_nn, mapping = aes(x=Var2, y=value, group=Var2, color=skippable)) + geom_boxplot()

#-----------------------------------------------------------------------------------------------------#
#---------------------------------correcting for G/C percentage---------------------------------------#
#-----------------------------------------------------------------------------------------------------#

library(stringr)
exonic_seq <- read.csv("~/Downloads/ExonsSpreadsheet-Homo_sapiens_Transcript_Exons_ENST00000328333.csv", stringsAsFactors = FALSE)
# fromJSON(exonic_seq$Sequence[2])

# grep("[ACTG]*", unlist(strsplit(exonic_seq$Sequence[2], " ")))
# regexpr("([ACTG])*", exonic_seq$Sequence[2])

exonic_seq$seq <- apply(exonic_seq, 1, function(x) str_extract(x[8], "[ACTG]{5,}"))
only_exonic <- exonic_seq[!is.na(exonic_seq$No.),]
only_exonic$gc <- apply(only_exonic, 1, function(x) round((str_count(x[9], "[GC]")/nchar(x[9]))*100, 1))
only_exonic$skippable <- as.factor(ifelse(only_exonic$No. %in% (skippable_exons0+1), "skippable", "not skippable"))

# p_gc = barplot van gc percentages
p_gc <- ggplot(only_exonic, mapping = aes(x=No., y=gc, color=skippable)) + geom_col() + geom_hline(yintercept = 35) + geom_hline(yintercept = 65)

#-----------------------------------------------------------------------------------------------------#
#--------------------------------------left aligns plots----------------------------------------------#
#-----------------------------------------------------------------------------------------------------#
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
#-----------------------------------------------------------------------------------------------------#

gc_labels <- sapply(only_exonic$gc, function(x) ifelse(x < 60 && x > 50, x, ""))
p_exons + stat_summary(geom = 'text', label = only_exonic$gc, fun.y = max, vjust = -1)
