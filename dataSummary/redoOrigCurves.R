library(tidyverse)
library(synapser)
synLogin()
##now get the manifest from synapse
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
  dplyr::rename(common_name='Sample')


##PDX contain list of files
pdx<-manifest|>
  dplyr::select(common_name,PDX_Drug_Data)|>
  # left_join(pdx_samps)|>
  distinct()|>
  subset(!is.na(PDX_Drug_Data))


##MTS contain lists of directories
mts<-manifest|>
  dplyr::select(common_name,MicroTissueDrugFolder)|>
  #left_join(org_samps)|>
  distinct()|>
  subset(!is.na(MicroTissueDrugFolder))

# Modify the extract_date_hour function to return a named vector
extract_date_hour <- function(experiment_id) {
  pattern <- "(\\d{6})_?(\\d{2,3})?"
  matches <- str_match(experiment_id, pattern)
  date <- matches[, 2]
  hour <- matches[, 3]
  date[is.na(date)] <- NA  # Replace with NA instead of blank
  hour[is.na(hour)] <- 48  # Replace with 48 instead of blank (default)
  return(list(date = date, hour = hour))
}
##first function to get children from parentId
getDrugDataByParent<-function(parid,sampleId){
  qtab<-synTableQuery(paste('select id,name,experimentalCondition,parentId from syn21993642 where parentId=\'',parid,'\''))$asDataFrame()|>
    as.data.frame()|>
    subset(!is.na(experimentalCondition))|>
    dplyr::select(id,name,experimentalCondition)|>
    subset(name!='synapse_storage_manifest.csv')
  ##now we need to parse the metadatda table get the info
  #  print()
  res<-do.call(rbind,lapply(qtab$id,function(x){
    # print(x)
    sname <- subset(qtab,id==x)
    #print(sname)
    time <- sname$experimentalTimePoint
    sname <-extract_date_hour(sname$name)
    #print(x)
    #print(sname)
    annotes <- synGetAnnotations(x)
    data <- data.table::fread(synGet(x)$path)|>
      dplyr::filter(response_type=='percent viability')|>
      mutate(improve_sample_id=sampleId,
             DOSE=(10^dosage)*1000000, ##dosage is log(M), need to move to micromolar
             GROWTH=response, #/100,
             source = "NF Data Portal",
             #CELL = improve_sample_id,
             Drug = compound_name,
             study = paste0('MT ',sname$date,' exp'),
             time = annotes$experimentalTimepoint,
             time_unit = 'hours') %>%
      dplyr::select(improve_sample_id,DOSE,GROWTH,source,chem_name,study,time)

    return(data)
  }))
  return(res)
}

##now loop through manifest to get all the files
mts_fold <- mts |>
  tidyr::separate_rows(MicroTissueDrugFolder,sep=',')

#mts_fold <- mts_fold[which(!mts_fold$V1%in%c("NA",NA)),]

alldrugs<-do.call(rbind,lapply(mts_fold$MicroTissueDrugFolder,function(x){
  samp<-subset(mts_fold,MicroTissueDrugFolder==x)
  print(samp$common_name)
  res<-getDrugDataByParent(x,samp$common_name)
  #print(res)
  return(res)
}))

readr::write_tsv(alldrugs,file='origDrugResponse.tsv')
