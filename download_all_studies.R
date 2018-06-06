library("recount")
#find all study accession numbers
all_studies <- abstract_search("",id_only = TRUE)

counter = 0
for(study in all_studies){
  url <- paste0("http://duffel.rail.bio/recount/", study, "/", study, ".junction_id_with_transcripts.bed.gz")
  jx_bed_gz_loc <- paste0(getwd(), "/", study, "/", study, ".junction_id_with_transcripts.bed.gz")
  tryCatch(
  if(!file.exists(file.path(jx_bed_gz_loc))){
    download.file(url, destfile = jx_bed_gz_loc)
  }  
  )
  
  tryCatch(
  if(!file.exists(paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study, "/"))){
    download_study(study, type = "rse-jx", outdir = paste0("/groups/umcg-gcc/tmp03/umcg-dvdveerdonk/studies/", study))
  }
  )
  if(counter %% 10){
    gc()
  }
  counter = counter + 1
}

