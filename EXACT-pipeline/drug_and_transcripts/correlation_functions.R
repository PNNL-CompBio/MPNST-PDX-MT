###correlation functions

mem.maxVSize(vsize=32762*2)
################################### data

##get sample data
samples <- readr::read_csv(synGet('syn72249842')$path)

###ask jermey to fix this
#samples$other_id <- stringr::str_replace_all(samples$other_id,'MN_2','MN-2')
#samples$other_id <- stringr::str_replace_all(samples$other_id,'JH_2_002','JH-2-002')
#samples$other_id <- stringr::str_replace_all(samples$other_id,'WU_225','WU-225')

library(synapser)

synLogin()
#get drug data
edrugs <- c('cudc-907','carboplatin','carbozantinib', 'rapamycin','eribulin','etoposide',
            'shp-099','dacarbazine','sapanisertib','everolimus', 'bortezomib', 'cerdulatinib',
            'brigatinib','erlotinib','defactinib', 'digoxin','epirubicine','pexidartinib','verteporphin',
            'irinotecan', 'ripretinib', 'sn-38','nazartinib')

ddrugs <- readr::read_tsv(synGet('syn72249846')$path) |>
  subset(chem_name %in% c('dmso','iag933','cudc-9087','sn-38','pexidartinib','verteporphin','decitabine',
                          tolower(names(drug.class))))


##get experimental data

target_study <- "EXACT"
dr_metric = "auc"

sing_exp <- readr::read_tsv(synGet('syn72331824')$path) |>
  subset(dose_response_metric == dr_metric & study == target_study) |>
  mutate(treat_time = as.factor(time)) |>
  left_join(samples) |>
  left_join(ddrugs) |>
  subset(!is.na(chem_name)) |>
  dplyr::select(common_name, improve_sample_id, other_id,
                study,
                drug = chem_name, treat_time,
                auc = dose_response_value) |>
  distinct()



###combination data
exp <- readr::read_tsv(synGet('syn72249845')$path) |>
  mutate(treat_time = as.factor(time)) |>
  left_join(samples) |>
  tidyr::separate(other_id, into=c('samp','initial_drug','time','treat','type'),sep='_')|>
  left_join(ddrugs)|>
  dplyr::select(common_name, improve_sample_id, initial_drug,
                drug=chem_name,expr_time=time, treat_time,
                delta_auc=dose_response_value)|>
  mutate(initial_drug=tolower(initial_drug))|>
  distinct()

#exp$initial_drug = stringr::str_replace_all(exp$initial_drug,'capmetinib','capmatinib')

exp <- exp |>
  subset(initial_drug != 'irinotecan')


#get gene expression

##get gene expression data
get_harmonized_gex <- function(){
  message("Pulling gene expression")
  genes <- readr::read_csv(synGet('syn72249894')$path)

  full_trans <- readr::read_csv(synGet('syn72249843')$path) |>
    left_join(genes) |>
    subset(!is.na(gene_symbol)) |>
    dplyr::select(transcriptomics,improve_sample_id, gene_symbol) |>
    unique()


  ##now let's filter for transcripts expressed across all
  tcounts <- full_trans |>
    group_by(gene_symbol) |>
    summarize(nsamps = n_distinct(improve_sample_id)) |>
    subset(nsamps == length(unique(full_trans$improve_sample_id)))

  full_trans <- full_trans |>
    subset(gene_symbol %in% tcounts$gene_symbol)

  trans <- full_trans |>
    #  filter(improve_sample_id %in% sing_exp$improve_sample_id)|>
    distinct() |>
    group_by(improve_sample_id, gene_symbol) |>
    summarize(exp = mean(transcriptomics, na.rm = T))


  tmat <- trans |>
    tidyr::pivot_wider(names_from=improve_sample_id,values_from=exp)|>
    tibble::column_to_rownames('gene_symbol') |>
    drop_na()

  ##filter out hi nas
  navals <- apply(tmat,1,function(x) length(which(is.na(x)))/length(x))
  tmat <- tmat[which(navals < 0.25),]

  ##now we can median center and log2 norom

  return(tmat)
}

################################### correlation data


##define function to compute correlation
compute_cor<-function(exp, tmat, drug,etime,ttime){
 # print(drug)
  dexp <- exp[exp$drug == drug &
                   exp$expr_time == etime &
                   exp$treat_time == ttime,]|>
    drop_na()|>
    tibble::column_to_rownames('improve_sample_id')
  if(nrow(dexp) > 0){
    samps <- intersect(rownames(dexp), colnames(tmat))
    av <- dexp[samps,'auc']
    print(paste('computing correlation across',length(samps),'samples'))
    cval <- suppressMessages(cor(av,
                                 t(tmat[,samps]),
                                 use = 'p',
                                 method = 'spearman'))

    return(data.frame(corval=t(cval), gene=colnames(cval),
                      drug=drug, ttime=ttime, etime=etime))
  }
}

#compute significance separately to avoid overloading memory
compute_sig <- function(exp, tmat, drug, etime, ttime){
  dexp <- subset(exp,drug == drug &
                   expr_time == etime &
                   treat_time == ttime)
  if(nrow(dexp) > 0) {
    cp <- apply(tmat[,as.character(dexp$improve_sample_id)],1,
                function(x){
                  y <- try(pv <- cor.test(as.numeric(x),
                                          dexp$auc,use ='p', method='spearman'), silent=TRUE)
                  if(!inherits(y, "try-error"))
                    return(pv$p.value)
                  else(return(1.0))
                })
    return(data.frame(sig = t(cp), gene = colnames(cval),
                      drug=drug, ttime=ttime, etime=etime,
                      improve_sample_id = unique(dexp$improve_sample_id)))
  }
}
#update memory size just in case
##for each expr time point and second drug we calculate correlation

get_ssgsea<-function(tmat){
  library(GSVA)

  gsetslst <- list(INNATE_RESPONSE=c("AIM2", "ALPK1", "AP3B1"),
                   ADAPTIVE_RESPONSE=c("CD27", "CD70", "EBAG9"))
  gsetslst
  gsetsgsc <- geneIdsToGeneSetCollection(gsetslst)

  ##now run ssGSEA on entire matrix

}

