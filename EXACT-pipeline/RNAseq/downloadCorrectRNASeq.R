##functions to analyze RNASeq data


if (file.exists('map.csv')) {
  map <- read.csv('map.csv')
} else {
  library(biomaRt)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  map<-getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters ='ensembl_gene_id', values =unique(rnaseqdata$ensembl_gene_id), mart = ensembl)
  write.csv(map, 'map.csv', row.names = FALSE)
}


rnaSeqTab<-getRNASeq<-function(){
  ##metadata query first queries file view for all files - in this case i do counts.sf files, but you can also query by data type
  ##Yang needs to update this!!
  manifest <- synapser::synTableQuery("select * from syn52369043 where name like '%quant.genes.sf'")$filepath
  manifest <- readr::read_csv(manifest)|>dplyr::select(name,id,specimenID)

  batch1metadata<-readr::read_csv(synGet('syn64608368')$path)|>
    dplyr::select(specimenID,individualID,Drug,Timepoint,Date,Proceeded_id)|>
    distinct()
  batch2metadata<-readr::read_csv(synGet('syn66302373')$path)|>
    dplyr::select(specimenID,individualID,Drug,Timepoint,Date,Proceeded_id)|>
    distinct()
  batch3metadata <- read.csv(synGet("syn66050299")$path)|> # was syn6559564 before 20250415
    dplyr::select(specimenID,individualID,Drug,Timepoint,Date,Proceeded_id)|>
    distinct()
  batch3metadata$specimenID <- sub("-", "_", batch3metadata$specimenID)
  batch4metadata<-readxl::read_excel(synGet('syn68080392')$path)|>
    dplyr::select(`Label on Cap`,Drug,Timepoint,`Individual ID`)|>
    distinct()
  batch4metadata$individualID <- batch4metadata$`Individual ID`
  batch4metadata$Timepoint <- sub("HR","h",batch4metadata$Timepoint)
  batch4metadata$Timepoint <- sub(" h","h",batch4metadata$Timepoint)
  colnames(batch4metadata)[1] <- "specimenID"
  batch4metadata$name <- paste0(batch4metadata$specimenID,"_quant.genes.sf")
  batch4metadata$Date <- "B4"
#  batch4metadata$Full_name <- batch4metadata$Drug

  ##METDATA - replace with annotations
  mn2_b1<-manifest|>
    dplyr::right_join(batch1metadata)|>
    dplyr::select(id,specimenID,individualID,Drug,Timepoint,Date)

  mn2_b2<-manifest|>
    tidyr::separate(name,into=c('Proceeded_id','sf'),sep='_quant.genes.')|>
    dplyr::select(-specimenID)|>
    dplyr::right_join(batch2metadata)|>
    dplyr::select(id,specimenID,individualID,Drug,Timepoint,Date)

  jh2002_b3 <- manifest|>
    dplyr::filter(stringr::str_starts(name, "B3_"))
  jh2002_b3$specimenID <- sub(".quant.genes.sf","",jh2002_b3$name) # filenames used to end in _quant.genes.sf before April 2, 2025
  jh2002_b3 <- merge(jh2002_b3, batch3metadata, by="specimenID", all=TRUE) |>
    dplyr::select(id,specimenID,individualID,Drug,Timepoint,Date)

  wu225_b4 <- manifest|>
    dplyr::filter(stringr::str_starts(name, "B4-"))
  wu225_b4$individualID <- "WU-225"
  wu225_b4$specimenID <- sub("_quant.genes.sf","",wu225_b4$name)
  wu225_b4 <- merge(wu225_b4, as.data.frame(batch4metadata), by=c("name","specimenID","individualID"), all=TRUE) |>
    dplyr::select(id,specimenID,individualID,Drug,Timepoint,Date)


  allmeta<-rbind(mn2_b1,mn2_b2,jh2002_b3,wu225_b4)
  if (any(is.na(allmeta$id))) {
    warning("missing synIDs for ", paste(allmeta[is.na(allmeta$id),]$specimenID, collapse=" "))
    allmeta <- allmeta[!is.na(allmeta$id),]
  }

  mtab<-do.call(rbind,lapply(allmeta$id,function(x){
    readr::read_tsv(synapser::synGet(x)$path)|>
      dplyr::mutate(id=x)
  }))|>left_join(allmeta)
  mtab$Drug <- sub(" ","",mtab$Drug)
#  mtab$Drug<-unlist(sapply(mtab$Drug,function(x) return(ifelse(x%in%names(drugnames),drugnames[[x]],x))))

  return(mtab)
}


addGeneNames <- function(rnaseqdata){

  rnaseqdata <- rnaseqdata |>
    tidyr::separate(Name,into=c('ensembl_gene_id','version'),remove=FALSE)

  # change "MN2" to "MN-2", "untreated" to "DMSO", "VIR" to "Vorinostat",
  # "IRI" to "Irinotecan", "RMC-4630" to "RMC4630", "Capmetinib" to "Capmatinib"
  #rnaseqdata[rnaseqdata$individualID=="MN2",]$individualID <- "MN-2"
  #rnaseqdata[rnaseqdata$Drug=="untreated",]$Drug <- "DMSO"
  #rnaseqdata[rnaseqdata$Drug=="VIR",]$Drug <- "Vorinostat"
  #rnaseqdata[rnaseqdata$Drug=="IRI",]$Drug <- "Irinotecan"
  #rnaseqdata[rnaseqdata$Drug=="RMC-4630",]$Drug <- "RMC4630"
  #rnaseqdata[rnaseqdata$Drug=="Capmetinib",]$Drug <- "Capmatinib"

  print(paste("Got RNASeq from",length(unique(rnaseqdata$id)),'samples across',length(unique(rnaseqdata$Drug)),'drugs and',length(unique(rnaseqdata$Timepoint)),'timepoints'))

  if (file.exists('map.csv')) {
    map <- read.csv('map.csv')
  } else {
    library(biomaRt)
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    map<-getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters ='ensembl_gene_id', values =unique(rnaseqdata$ensembl_gene_id), mart = ensembl)
    write.csv(map, 'map.csv', row.names = FALSE)
  }

  rnaseqgenes <- rnaseqdata |> left_join(map) |>
    subset(hgnc_symbol!="")|>
    dplyr::group_by(hgnc_symbol,Date,id,Drug,Timepoint,individualID)|>
    summarize(TPM=sum(TPM),NumReads=sum(NumReads))

  return(rnaseqgenes)
}

batchCorrect <- function(reduced.data, sampmeta){

tpm<-reduced.data|>
  dplyr::select(Name,id,TPM)|>
  tidyr::pivot_wider(names_from=id,values_from=TPM)|>
  tibble::column_to_rownames('Name')|>
  as.matrix()

numReads<-reduced.data|>
  dplyr::select(Name,id,NumReads)|>
  tidyr::pivot_wider(names_from=id,values_from=NumReads)|>
  tibble::column_to_rownames('Name')|>
  as.matrix()


counts<-reduced.data|>
  dplyr::select(Name,id,NumReads)|>
  tidyr::pivot_wider(names_from=id,values_from=NumReads)|>
  tibble::column_to_rownames('Name')|>
  as.matrix()


dds<-DESeqDataSetFromMatrix(countData = round(counts),
                              colData = sampmeta,
                              design = ~ Drug+Date+Timepoint )




dds<-DESeq(dds)

#lets look at both rlog and vsd as transformations
vsd <- vst(dds, blind=FALSE)

mat <- assay(vsd)
mm <- model.matrix(~Timepoint+Drug, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Date, design=mm)
assay(vsd) <- mat



##cortab
md=colData(vsd)[,c('Drug','Timepoint','Date')]
corvsd=cor(assay(vsd),use='pairwise.complete.obs')
#pheatmap::pheatmap(corvsd,annotation_row=md,annotation_col=md)

##now add the VSD counts to the table
vtab<-assay(vsd)|>as.data.frame()|>
  tibble::rownames_to_column('Name')|>
  tidyr:::pivot_longer(starts_with('syn'),names_to='id',values_to='VSD_batchCorrected')

rnaseqdata<-rnaseqdata|>
  left_join(vtab)

readr::write_csv(rnaseqdata,file='rnaseqcounts.csv')

synFile <- synapser::File('rnaseqcounts.csv', parent = 'syn64004220')
synapser::synStore(synFile)

}

doPCA<-function(rnaseqdata){

  sampmeta<-rnaseqdata|>
    dplyr::select(id,specimenID,Timepoint,Drug,Date,individualID)|>
    distinct()
  sampmeta$indDate <- paste0(sampmeta$individualID, ", ",sampmeta$Date)

  pcs<-prcomp(log2(0.0001+tpm),scale=TRUE,center=TRUE)

  pcvals<-pcs$rotation|>
    as.data.frame()|>
    dplyr::select(PC1,PC2)|>
    tibble::rownames_to_column('id')|>
    left_join(sampmeta)
  pcvals$Timepoint <- factor(pcvals$Timepoint, c("0h", "8h", "24h"))

  ggplot(pcvals,aes(x=PC1,y=PC2,col=Drug,shape=Timepoint))+geom_point(size=4)+
    facet_grid(~indDate)+scale_color_discrete()
  ggplot(pcvals,aes(x=PC1,y=PC2,col=Drug,shape=Timepoint))+geom_point(aes(size=4))+facet_grid(~Date)+scale_color_manual(values=cols)

  ggsave('batched_pca.png',width=12,height=6)

  ggplot(pcvals,aes(x=PC1,y=PC2,col=Drug,shape=Timepoint))+geom_point(size=4)+
    facet_grid(~individualID)+scale_color_discrete()
  ggsave('batched_pca_noDates.png',width=12,height=6)

  ggplot(pcvals,aes(x=PC1,y=PC2,col=Drug,alpha=Timepoint, shape=indDate))+geom_point(size=4)+
    scale_alpha_discrete(range=c(0.1,1))+scale_color_discrete()
  ggsave('batched_pca_noFacet.png',width=8,height=8)

  drug.info <- list("Palbociclib" = "CDK inhibitor", "Ribociclib" = "CDK inhibitor",
                    "Trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                    "Decitabine" = "DNMT inhibitor", "Vorinostat" = "HDAC inhibitor",
                    "Mirdametinib" = "MEK inhibitor", "Selumetinib" = "MEK inhibitor",
                    "Trametinib" = "MEK inhibitor", "Capmatinib" = "MET inhibitor",
                    "Olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                    "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "TOP inhibitor",
                    "Irinotecan" = "TOP inhibitor", "Verteporfin" = "YAP inhibitor")
  drugCol <- colorRampPalette(RColorBrewer::brewer.pal(length(drug.info),"Set3"))(length(drug.info)) # need color ramp for >12 colors
  drugCol <- c("black", drugCol)
  names(drugCol) <- c("DMSO", names(drug.info))
  ggplot(pcvals,aes(x=PC1,y=PC2,col=Drug,alpha=Timepoint, shape=individualID))+geom_point(size=4)+
    scale_alpha_discrete(range=c(0.2,1))+scale_color_discrete() + labs(shape="MPNST") + theme_minimal() +
    scale_color_manual(values=drugCol, breaks=names(drugCol))
  ggsave('batched_pca_noFacet_noDates.png',width=8,height=8)
  ggplot(pcvals,aes(x=PC1,y=PC2,col=Drug,size=Timepoint, shape=individualID))+geom_point()+
    scale_color_discrete() + labs(shape="MPNST") + theme_minimal() +
    scale_color_manual(values=drugCol, breaks=names(drugCol))
  ggsave('batched_pca_noFacet_noDates_sizeTime.png',width=8,height=8)

  ggplot(pcvals[pcvals$Date %in% c("B", "C"),],aes(x=PC1,y=PC2,col=Drug,alpha=Timepoint, shape=indDate))+geom_point(size=4)+
    scale_alpha_discrete(range=c(0.1,1))+scale_color_discrete()
  ggsave('batched_pca_noFacet_datesBC.png',width=8,height=6)
}
