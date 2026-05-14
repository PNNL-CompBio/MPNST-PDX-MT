###correlation functions

mem.maxVSize(vsize = 32762*2)
################################### data

library(synapser)
library(GSVA)
library(msigdbr)
library(patchwork)
library(data.table)
library(dplyr)
library(tidyr)
synLogin()

##get sample data
samples <- readr::read_csv(synGet('syn72249842')$path)

#get drug data
edrugs <- c('cudc-907','carboplatin','carbozantinib', 'rapamycin','eribulin','etoposide',
            'shp-099','dacarbazine','sapanisertib','everolimus', 'bortezomib', 'cerdulatinib',
            'brigatinib','erlotinib','defactinib', 'digoxin','epirubicine','pexidartinib','verteporfin',
            'irinotecan', 'ripretinib', 'sn-38','nazartinib')

print('loading drugs')
ddrugs <- readr::read_tsv(synGet('syn72249846')$path) |>
  subset(chem_name %in% c('dmso','iag933','cudc-9087','sn-38','pexidartinib','verteporphin','decitabine',
                          tolower(names(drug.class))))


##get experimental data
ddrugs$chem_name[ddrugs$chem_name=='verteporphin']<-'verteporfin'

target_study <- "EXACT"
dr_metric = "auc"

print('loading single agent dose response')
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


print('loading combination data')
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

md <- exp|>group_by(drug,treat_time) |> summarize(count =n())
#get gene expression


get_gex_for_samples <- function(){

}


##deduplication of RNA takes a long time
#so here i store it locally
all_rnaseq <- function(){

  if(file.exists('gex.rda')){
    message('Reading gene expression from file')
    load('gex.rda')
  }else{
  message("Loading gene expression")
  genes <- readr::read_csv(synGet('syn72249894')$path)

  all_trans <- readr::read_csv(synGet('syn72249843')$path) |>
    left_join(genes) |>
    subset(!is.na(gene_symbol)) |>
    dplyr::select(transcriptomics,improve_sample_id, gene_symbol) |>
    unique()
  save(all_trans,file='gex.rda')

  }
  return(all_trans)
}

##get nz gex across all samples
get_harmonized_gex <- function(){

  all_trans <- all_rnaseq()

    ##now let's filter for transcripts expressed across all
  #tcounts <- all_trans |>
  #  group_by(gene_symbol) |>
  #  summarize(nsamps = n_distinct(improve_sample_id))

#  all_samps <- all_trans |>
#      subset(nsamps == length(unique(full_trans$improve_sample_id)))

  #full_trans <- all_trans #all_samps |>
  #  subset(gene_symbol %in% tcounts$gene_symbol)

  trans <- all_trans |>
    #  filter(improve_sample_id %in% sing_exp$improve_sample_id)|>
    #distinct() |>
    group_by(improve_sample_id, gene_symbol) |>
    summarize(exp = mean(transcriptomics, na.rm = T))

  tmat <- trans |>
    tidyr::pivot_wider(names_from=improve_sample_id,values_from=exp)|>
    tibble::column_to_rownames('gene_symbol') #|>
    #drop_na()

  ##filter out hi nas
  navals <- apply(tmat,1,function(x) length(which(is.na(x)))/length(x))
 # tmat <- tmat[which(navals < 0.25),]

  ##log2
  tmat <- 0.0001+log2(tmat)
  ##now median center by sample so we can better compare
  tmat <- apply(tmat,2,function(x) x-median(x,na.rm=T))

  #mvals <- apply(tmat,1,mean,na.rm = T)
  #for(i in 1:nrow(tmat))
  #  tmat[i,which(is.na(tmat[i,]))] <- mvals[i]
  ##now we can median center and log2 norom

  return(tmat)
}

get_gsea_scores <- function(tmat){


  m_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  gene_sets <- split(x = m_df$gene_symbol, f = m_df$gs_name)

  # Define the ssGSEA parameters
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = as.matrix(tmat),
    geneSets = gene_sets
  )

  # Run GSVA using the ssGSEA method
  ssgsea_scores <- GSVA::gsva(param = ssgsea_param)

  return(ssgsea_scores)
}

################################### correlation data


##define function to compute correlation
compute_cor<-function(exp, tmat, drug,etime,ttime){
 # print(drug)
  dexp <- exp[exp$drug == drug &
                   exp$expr_time == etime &
                   exp$treat_time == ttime,] |>
    drop_na() |>
    tibble::column_to_rownames('improve_sample_id')
  if(nrow(dexp) > 0){
    samps <- intersect(rownames(dexp), colnames(tmat))
    av <- dexp[samps,'auc']
  #  print(paste('computing correlation across',length(samps),'samples'))
    cval <- suppressMessages(cor(av,
                                 t(tmat[,samps]),
                                 use = 'p',
                                 method = 'p'))

    return(data.frame(corval=t(cval), gene=colnames(cval),
                      drug=drug, ttime=ttime, etime=etime))
  }
}

#compute significance separately to avoid overloading memory
compute_sig <- function(exp, tmat, drug, etime, ttime){
#  print(drug)
  dexp <- exp[exp$drug == drug &
                exp$expr_time == etime &
                exp$treat_time == ttime,]|>
    drop_na()|>
    tibble::column_to_rownames('improve_sample_id')
 # print(dim(dexp))
  if(nrow(dexp) > 0) {
    samps <- intersect(rownames(dexp), colnames(tmat))
    av <- dexp[samps,'auc']
    cp <- apply(tmat[,samps],1,
                function(x){
                  y <- try(pv <- cor.test(as.numeric(x),
                                          av,use ='p', method='p'), silent=TRUE)
                  if(!inherits(y, "try-error"))
                    return(pv$p.value)
                  else(return(1.0))
                })
    return(data.frame(sig = cp, gene = names(cp), qval = p.adjust(cp),
                      drug=drug, ttime=ttime, etime=etime))
  }
}

#quick function to recreate bubble plots for ssgsea
bubble_plot <- function(merged_cor, sig_thresh = 0.05){

  #setOrder <- mean.results[order(mean.results$NES),]$Gene_set
  merged_cor$Pathway = gsub('HALLMARK_','',merged_cor$gene)
  drug.order <- intersect(names(classCol),merged_cor$drug)
  merged_cor$InvCor =  -1 * merged_cor$corval
  #filter out pathways that are significant in at least one samp;e
  sig_paths <- merged_cor|>
      group_by(Pathway) |>
      summarize(sigval = min(sig) < sig_thresh) |>
      subset(sigval)

  dot.df <- merged_cor[merged_cor$Pathway  %in% sig_paths$Pathway,]#[gsea.df$Gene_set %in% topSets,]
  dot.df$Time <- factor(dot.df$ttime, levels=c(48,120))
  if (any(dot.df$sig == 0)) {
    dot.df[dot.df$sig == 0,]$sig <- 0.0001
  }
  dot.df$minusLogP <- -log(dot.df[,"sig"], base = 10)
  dot.df$minusLogQ <- -log(dot.df[,"qval"], base = 10)
  maxAbsNES <- ceiling(max(abs(dot.df$corval)))
  dot.plot <- ggplot2::ggplot(
    na.omit(dot.df),
    ggplot2::aes(
      x = drug, y = Pathway, color = corval,
      size = minusLogP
    )
  ) +
    ggplot2::geom_point() +
   # ggplot2::scale_y_discrete(limits = setOrder) +
    ggplot2::scale_x_discrete(limits=drug.order) +
    scale_color_gradient2(low="blue",high="red", mid="grey",
                          limits=c(-maxAbsNES, maxAbsNES)) +
    #viridis::scale_color_viridis() +
    theme_minimal(base_size=12) +
    ggplot2::labs(color = "Correlation", size = "-log(P)") +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, vjust=1, hjust=1,
                                   color = classCol[drug.order])) +
    geom_point(data = subset(dot.df, sig < sig_thresh), col = "black", stroke = 1.5, shape = 21) +
    facet_grid(etime ~ ttime)
  dot.plot
  # most are only in prot and not in RNA
  #ggplot2::ggsave("enrichedPathways_dotPlot.pdf", dot.plot, width=12, height=49)
}
