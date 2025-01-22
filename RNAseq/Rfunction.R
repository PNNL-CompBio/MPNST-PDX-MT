

# Function needed
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

PCA_Plot_3<- function (data,Annotation,VAR,Color) {
  # logcountdata row:genes,column: samples
  pca <- prcomp(data) 
  pca_out<-as.data.frame(pca$x)
  df_out<- pca_out %>%tibble::rownames_to_column(var=VAR)  %>% left_join(., Annotation) 
  #df_out<- merge (pca_out,Annotation,by.x=0,by.y=0) 
  
  # label_color<- factor(df_out[,group])
  ggplot(df_out,aes_string(x="PC1",y="PC2")) +geom_point(aes_string(colour = Color))
}

PCA_Plot_3_shape<- function (data,Annotation,VAR,Color,Shape) {
  # logcountdata row:genes,column: samples
  pca <- prcomp(data) 
  pca_out<-as.data.frame(pca$x)
  df_out<- pca_out %>%tibble::rownames_to_column(var=VAR)  %>% left_join(., Annotation) 
  #df_out<- merge (pca_out,Annotation,by.x=0,by.y=0) 
  
  # label_color<- factor(df_out[,group])
  ggplot(df_out,aes_string(x="PC1",y="PC2")) +geom_point(aes_string(colour = Color,shape=Shape),size=5)
}



Drug_count_plot=function(logcount,sample_anno){
  
  logcount_2=logcount%>%t() %>%as.data.frame() %>%mutate(Sample_name=rownames(.))
  Level=sample_anno %>%arrange(Drug)%>%.$Sample_name
  data_mod = reshape2::melt(logcount_2, id.vars="Sample_name", 
                            measure.vars=rownames(logcount)) %>%left_join(sample_anno)%>%arrange(Drug)
  
  data_mod$Sample_name=factor(data_mod$Sample_name,levels =Level )
  #return(logcount_2)
  ggplot(data_mod)+geom_boxplot(aes(x=Sample_name, y=value, color=Drug,fill=Drug))+xlab("Sample")+ylab("log2(Count+1)")+theme(axis.text.x=element_blank())
  
  
}

selectGenes_BY_group <- function(Counts, min.count=10, N=0.90,Sample_Anno){
  Group_list=list()
  for (i in unique(Sample_Anno$Group)) {
    sample_name=Sample_Anno%>%dplyr::filter(Group==i)%>%.$Sample_name
    counts=Counts[, sample_name]
    lib.size <- base::colSums(counts)
    MedianLibSize <- median(lib.size)
    CPM.Cutoff <- min.count / MedianLibSize*1e6
    CPM <- edgeR::cpm(counts,lib.size=lib.size)
    
    min.samples <- round(N * ncol(counts))
    
    f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
    flist <- genefilter::filterfun(f1)
    Group_list[[i]] <- genefilter::genefilter(CPM, flist)%>%as.data.frame()%>%setNames(paste0(i  ,"_","TF"))
  }
  ## the same as:
  #keep <- apply(CPM, 1, function(x, n = min.samples){
  #  t = sum(x >= CPM.Cutoff) >= n
  #  t
  #})
  
  Output=do.call("cbind",Group_list)
  keep.exprs_2=Output%>%mutate(Sum=apply(.,1,sum))%>%dplyr::filter(Sum==length(unique(Sample_Anno$Group)))
  myFilt <-Counts[rownames(keep.exprs_2),]
  
}

selectGenes <- function(counts, min.count=10, N=0.90){
  
  lib.size <- base::colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- edgeR::cpm(counts,lib.size=lib.size)
  
  min.samples <- round(N * ncol(counts))
  
  f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  flist <- genefilter::filterfun(f1)
  keep <- genefilter::genefilter(CPM, flist)
  
  ## the same as:
  #keep <- apply(CPM, 1, function(x, n = min.samples){
  #  t = sum(x >= CPM.Cutoff) >= n
  #  t
  #})
  
  return(keep)
}


selectGenes_BY_group <- function(Counts, min.count=10, N=0.90,Sample_Anno){
  Group_list=list()
  for (i in unique(Sample_Anno$Group)) {
    sample_name=Sample_Anno%>%dplyr::filter(Group==i)%>%.$Sample_name
    counts=Counts[, sample_name]
    lib.size <- base::colSums(counts)
    MedianLibSize <- median(lib.size)
    CPM.Cutoff <- min.count / MedianLibSize*1e6
    CPM <- edgeR::cpm(counts,lib.size=lib.size)
    
    min.samples <- round(N * ncol(counts))
    
    f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
    flist <- genefilter::filterfun(f1)
    Group_list[[i]] <- genefilter::genefilter(CPM, flist)%>%as.data.frame()%>%setNames(paste0(i  ,"_","TF"))
  }
  ## the same as:
  #keep <- apply(CPM, 1, function(x, n = min.samples){
  #  t = sum(x >= CPM.Cutoff) >= n
  #  t
  #})
  
  Output=do.call("cbind",Group_list)
  keep.exprs_2=Output%>%mutate(Sum=apply(.,1,sum))%>%dplyr::filter(Sum==length(unique(Sample_Anno$Group)))
  myFilt <-Counts[rownames(keep.exprs_2),]
  
}


Volcano_plot_rnaseq_simple=function(All_res_Data,Pathway_list,pathway_name){
  
  cols <- c("DOWN"="#26b3ff", "NO"= "grey","UP"= "#ffad73","up" = "red", "down" = "blue") 
  sizes <- c( "DOWN"=0.5,"UP"=0.5 ,"up" = 3.5, "down" = 3.5, "NO" = 0.5) 
  alphas <- c("UP" = 1, "DOWN" = 1, "NO" = 0.5)
  
  pathway_genes=Pathway_list [[pathway_name]]%>%as.character()
  
  library(ggrepel)
  # Data=All_res_Data%>%dplyr::filter(!gene_symbol%in%pathway_genes)%>%mutate(diffexpressed=ifelse(padj>0.05|is.na(padj),"NO",ifelse(log2FoldChange>0,"UP","DOWN")))%>%mutate(delabel=NA)
  
  pathway_Data=All_res_Data%>%.[.$gene_symbol%in%pathway_genes,]%>%mutate(diffexpressed=ifelse(padj>0.05|is.na(padj),"NO",ifelse(log2FoldChange>0,"up","down")))%>%mutate(delabel=ifelse(diffexpressed!="NO"&(log2FoldChange>0.27|log2FoldChange<as.numeric( "-0.27")),gene_symbol,NA))
  # add a column of NAs
  #de=rbind(Data,pathway_Data)
  de=pathway_Data
  
  MAX_y=max(-log10(de$padj))
  # de$delabel <- NA
  # 
  # de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    geom_point(show.legend = F) + ggtitle(pathway_name)+
    theme_minimal() +
    geom_text_repel(size=3)+scale_color_manual(values=cols)+scale_size_manual(values=sizes)+scale_alpha_manual(values=alphas)+
    geom_vline(xintercept=c(-0.27, 0.27), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+theme_bw()+  theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.grid.major = element_blank()) 
  
}

Volcano_plot_rnaseq_simple_2 <- function(All_res_Data, Pathway_list = NULL, pathway_name = NULL, top_n = 10) {
  library(ggplot2)
  library(ggrepel)
  
  # Define colors, sizes, and alphas for plotting
  cols <- c("DOWN" = "#26b3ff", "NO" = "grey", "UP" = "#ffad73", "up" = "red", "down" = "blue")
  sizes <- c("DOWN" = 0.5, "UP" = 0.5, "up" = 3.5, "down" = 3.5, "NO" = 0.5)
  alphas <- c("UP" = 1, "DOWN" = 1, "NO" = 0.5)
  
  # Check if a pathway list is provided
  if (!is.null(Pathway_list)) {
    pathway_genes <- Pathway_list[[pathway_name]] %>% as.character()
    
    # Filter data to include only pathway genes and assign labels
    pathway_Data <- All_res_Data %>%
      filter(gene_symbol %in% pathway_genes) %>%
      mutate(diffexpressed = ifelse(padj > 0.05 | is.na(padj), "NO", ifelse(log2FoldChange > 0, "up", "down"))) %>%
      mutate(delabel = ifelse(diffexpressed != "NO" & abs(log2FoldChange) > 0.27, gene_symbol, NA))
    
    de <- pathway_Data
    
  } else {
    
    # Process data without pathway filtering
    de <- All_res_Data %>%
      mutate(diffexpressed = ifelse(padj > 0.05 | is.na(padj), "NO", ifelse(log2FoldChange > 0, "up", "down"))) %>%
      mutate(delabel = ifelse(diffexpressed != "NO" & abs(log2FoldChange) > 0.27, gene_symbol, NA))
  }
  
  # Select top significant genes for both up and down
  significant_genes_up <- de %>%
    dplyr::filter(diffexpressed == "up") %>%
    arrange(padj) %>%
    head(top_n)
  
  significant_genes_down <- de %>%
    dplyr::filter(diffexpressed == "down") %>%
    arrange(padj) %>%
    head(top_n)
  
  significant_genes <- bind_rows(significant_genes_up, significant_genes_down)
  
  # Generate the volcano plot
  ggplot(data = de, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
    geom_point(show.legend = F) + 
    geom_point(data = significant_genes, aes(x = log2FoldChange, y = -log10(padj)), color = "black", size = 3) +
    geom_text_repel(data = significant_genes, aes(label = delabel), size = 3) + 
    ggtitle(ifelse(is.null(pathway_name), "Volcano Plot", pathway_name)) +
    theme_minimal() +
    scale_alpha_manual(values = alphas) +
    scale_color_manual(values = cols) +
    scale_size_manual(values = sizes) +
    geom_vline(xintercept = c(-0.27, 0.27), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
}




Volcano_plot_rnaseq_add_label_update=function(All_res_Data,Pathway_list,pathway_name){
  
  cols <- c("DOWN"="#26b3ff", "NO"= "grey","UP"= "#ffad73","up" = "red", "down" = "blue") 
  sizes <- c( "DOWN"=0.5,"UP"=0.5 ,"up" = 3.5, "down" = 3.5, "NO" = 0.5) 
  alphas <- c("UP" = 1, "DOWN" = 1, "NO" = 0.5)
  
  pathway_genes=Pathway_list [[pathway_name]]%>%as.character()
  
  library(ggrepel)
  # Data=All_res_Data%>%dplyr::filter(!gene_symbol%in%pathway_genes)%>%mutate(diffexpressed=ifelse(padj>0.05|is.na(padj),"NO",ifelse(log2FoldChange>0,"UP","DOWN")))%>%mutate(delabel=NA)
  
  pathway_Data=All_res_Data%>%.[.$gene_symbol%in%pathway_genes,]%>%mutate(diffexpressed=ifelse(padj>0.05|is.na(padj),"NO",ifelse(log2FoldChange>0,"up","down")))%>%mutate(delabel=ifelse(diffexpressed!="NO"&(log2FoldChange>0.27|log2FoldChange<as.numeric( "-0.27")),gene_symbol,NA))
  # add a column of NAs
  #de=rbind(Data,pathway_Data)
  de=pathway_Data
  
  MAX_y=max(-log10(de$padj))
  # de$delabel <- NA
  # 
  # de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
    geom_point(show.legend = F) + ggtitle(pathway_name)+
    theme_minimal() +
    geom_text_repel(size=3,force_pull= 0,max.overlaps  = Inf,
                    data = subset(de, log2FoldChange >2),
                    nudge_x= 6 - data$log2FoldChange,
                    segment.size= 0.2,
                    segment.color= "grey50",
                    direction = "y",
                    hjust         = 0.5 ) +geom_text_repel(size=3,force_pull= 0,max.overlaps  = Inf,
                                                           data = subset(de, log2FoldChange >0.27&log2FoldChange<2),
                                                           nudge_x= 4.5 - subset(de, log2FoldChange >1&log2FoldChange<2)$log2FoldChange,
                                                           segment.size= 0.2,
                                                           segment.color= "grey50",
                                                           direction = "y",
                                                           hjust         = 0.5 ) +geom_text_repel(size=3,force_pull= 0,max.overlaps  = Inf,
                                                                                                  data = subset(de, log2FoldChange < as.numeric(-0.27)),
                                                                                                  nudge_x= -5 - subset(de, log2FoldChange <as.numeric(-1))$log2FoldChange,
                                                                                                  segment.size= 0.2,
                                                                                                  segment.color= "grey50",
                                                                                                  direction = "y",
                                                                                                  hjust         = 0.5 ) +scale_color_manual(values=cols)+scale_size_manual(values=sizes)+scale_alpha_manual(values=alphas)+
    geom_vline(xintercept=c(-0.27, 0.27), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")+theme_bw()+  theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.grid.major = element_blank()) 
  
  
  
  
}

Gene_pathway_sig_enrich=function(gene) {
  
  library(msigdbr)
  m_df <- msigdbr(species = "Homo sapiens")
  Cat=c("C2","C4","C5","C6","H")
  emlist=list()
  for (i in Cat) {
    
    m_t2g <- msigdbr(species = "Homo sapiens", category = i)%>% dplyr::select(gs_name, gene_symbol)
    #MSigDb over-presentaton analysis
    em= enricher(gene, TERM2GENE=m_t2g)
    
    emlist[[i]]=em@result%>%as.data.frame()%>%dplyr::filter(pvalue<0.05)
  }
  
  out=do.call("rbind",emlist)
  
}

Bar_plot_RES <- function(RES_data, FC_cutoff = 0.67,Ylim=2000) {
  library(ggplot2)
  
  # Prepare the data for plotting with the specified FC cutoff
  plot_data <- RES_data %>%
    mutate(diffexpressed = ifelse(padj > 0.05 | is.na(padj), "NO", 
                                  ifelse(log2FoldChange > FC_cutoff, "up", 
                                         ifelse(log2FoldChange < -FC_cutoff, "down", "NO")))) %>%
    dplyr::filter(diffexpressed %in% c("up", "down")) %>%  # Focus only on "up" and "down"
    group_by(diffexpressed) %>%
    summarise(count = n()) %>%  # Count the number of up and down regulated genes
    ungroup()
  
  # Create the bar plot
  ggplot(plot_data, aes(x = diffexpressed, y = count, fill = diffexpressed)) +
    geom_bar(stat = "identity", width = 0.7) +ylim(c(0,Ylim))+  # Create a bar plot
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +  # Set colors for up and down
    labs(x = "Regulation", y = "Number of Genes", 
         title = paste("Number of Up and Down Regulated Genes (FC cutoff =", FC_cutoff, ")")) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
}

PCA_Plot_3_shape<- function (data,Annotation,VAR,Color,Shape) {
  # logcountdata row:genes,column: samples
  pca <- prcomp(data) 
  pca_out<-as.data.frame(pca$x)
  df_out<- pca_out %>%tibble::rownames_to_column(var=VAR)  %>% left_join(., Annotation)%>%dplyr::filter(!is.na(Drug))
  #df_out<- merge (pca_out,Annotation,by.x=0,by.y=0) 
  
  # label_color<- factor(df_out[,group])
  ggplot(df_out,aes_string(x="PC1",y="PC2")) +geom_point(aes_string(colour = Color,shape=Shape),size=5)
}


TPM_heatmap=function(TPM_Data,Gene_names,Annotation_col,Annotation_color,gap_set,Title) {
  Input_data=TPM_Data%>%.[Gene_names,]
  pheatmap:: pheatmap(Input_data,cluster_cols = F,cluster_rows = F,
                      scale = "row",  # Scale rows (genes)
                      clustering_method = "complete",  # Hierarchical clustering method
                      color = colorRampPalette(c("blue", "white", "red"))(50),  # Color palette
                      main = Title,  # Title
                      gaps_col = gap_set,annotation_col = Annotation_col,
                      annotation_colors =Annotation_color,
                      fontsize = 8,  # Font size for labels
                      cellwidth = 15,  # Width of each cell
                      cellheight = 15,  # Height of each cell
                      show_rownames = T,  # Show row names (gene names)
                      show_colnames = TRUE  # Show column names (sample names)
  )
}


Enrich_dotplot_Single=function(FGSEA_out){
  library(ggbreak)
  FGSEA_out$NES=round(FGSEA_out$NES,2)
  #mycolor=RColorBrewer::brewer.pal(3,"YlOrRd")
  data=FGSEA_out%>%mutate(Norm_padj=-log10(padj))%>%mutate(State=ifelse(NES>0&padj<0.05,"Activated",ifelse(NES<0&padj<0.05,"Inhibited", ifelse(NES>0&padj<0.07,"P_Activated",ifelse(NES<0&padj<0.07,"P_Inhibited","no_change")))))%>%arrange(NES)
  data$pathway=factor(data$pathway,levels =unique(data$pathway ))
  MAX_NES=max(data$NES)%>%ceiling()
  min_NES=min(data$NES)%>%floor()
  # ggplot(data) +geom_point(aes(x = NES, y = pathway, color = -log10(padj),size=gene_counts))+scale_size("gene_counts") +
  color_set= setNames(c('red','green',"wheat","turquoise" ,'grey'),c("Activated", "Inhibited","P_Activated","P_Inhibited","no_change"))
  
  P= ggplot(data) +geom_point(aes(x = NES, y = pathway, color =State ,size=Norm_padj)) +
    theme_bw() +labs(x="Normalized enrichment score",y="")+
    theme(axis.text.x = element_text(size=rel(1.15)),
          axis.title = element_text(size=rel(1.15)),plot.title = element_text(hjust=0.5, 
                                                                              face = "bold"),legend.title = element_text(size=rel(1.15),hjust=0.5, face="bold"),axis.text.y = element_text(size =15))+
    scale_size_continuous (name="-log10(padj)",range=c(2,12))+scale_color_manual(name="State",values=color_set)
  # if (min_NES>0|MAX_NES<0) {
  #   return(P)
  #   } else {
  #   Breaks=c(min_NES,-1,1,MAX_NES)
  #  return(P+scale_x_break(c(-1, 1), scales = 1.5))
  # xlab("Normalized enrichment score") +ylab("")+scale_x_break(c(-1, 1), scales = 1.5)
  #   }
}

GSEA_heatmap=function(GSEA_Data,Annotation_col,Annotation_color,Title) {
  pheatmap:: pheatmap(GSEA_Data,cluster_cols = F,
                      scale = "row",  # Scale rows (genes)
                      clustering_method = "complete",  # Hierarchical clustering method
                      color = colorRampPalette(c("blue", "white", "red"))(50),  # Color palette
                      main = Title,  # Title
                      gaps_col = c(5,5),annotation_col = Annotation_col,
                      annotation_colors =Annotation_color,
                      fontsize = 8,  # Font size for labels
                      cellwidth = 15,  # Width of each cell
                      cellheight = 15,  # Height of each cell
                      show_rownames = T,  # Show row names (gene names)
                      show_colnames = TRUE  # Show column names (sample names)
  )
}


singcore_GSEA=function (data,genelist, is.na=T) {
  nm=names(genelist)
  plist=list()
  Data_prepare=singscore::rankGenes (data)
  
  for (i in seq_along(genelist)) { 
    Genes=genelist[[i]]%>%as.character()
    scoredf<- singscore::simpleScore (Data_prepare, upSet=unique(Genes))%>%.[1]
    #colnames(scoredf)=paste0(nm[i],"_singscore")
    colnames(scoredf)=paste0(nm[i])
    plist[[i]]=scoredf
  }
  output=do.call("cbind",plist)
}


Singscore_boxplot_timepoint_drug <- function(TPM_filter, Pathway_list, Drug_DMSO_data,Drug_colors) {
  # Step 1: Calculate Singscore
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  singscore_df <- singcore_GSEA(TPM_filter, Pathway_list) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Pathway")
  
  # Step 2: Transform to long format
  singscore_long <- singscore_df %>%
    tidyr::pivot_longer(
      cols = -Pathway, 
      names_to = "Sample_ID", 
      values_to = "Score"
    )
  
  # Step 3: Merge with annotation data
  merged_data <- singscore_long %>%
    left_join(Drug_DMSO_data, by = "Sample_ID") 
  
  # Step 4: Generate boxplot
  p <- ggplot(merged_data, aes(x = Timepoint, y = Score, fill = Drug)) +
    geom_boxplot(
      aes(group = interaction(Timepoint, Drug)),
      position = position_dodge(width = 0.5), # Controls box spacing
      outlier.shape = NA, alpha = 0.7
    ) +
    geom_jitter(
      aes(color = Drug),
      position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), # Closer jitter points
      size = 1
    ) +
    facet_wrap(~Pathway, scales = "free_y", ncol = 3)+scale_fill_manual(values = Drug_colors)+
    theme_minimal() +
    labs(
      x = "Timepoint",
      y = "Singscore",
      fill = "Drug",
      color = "Drug"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate x-axis labels
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.spacing = unit(1, "lines"), # Add space between facets
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add black border
      strip.background = element_rect(fill = "white", color = "black"), # Frame facet labels
      strip.text = element_text(face = "bold", size = 10) # Bold facet labels
    )
  
  # Step 5: Return the plot
  return(p)
}

Singscore_lineplot_timepoint_drug <- function(TPM_filter, Pathway_list, Drug_DMSO_data,Drug_colors,Title=NULL) {
  if(is.null(Title)) {
    Title="Singscore line plot"
    
  }
  # Load necessary libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # Step 1: Calculate Singscore
  singscore_df <- singcore_GSEA(TPM_filter, Pathway_list) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Pathway")
  
  # Step 2: Transform to long format
  singscore_long <- singscore_df %>%
    tidyr::pivot_longer(
      cols = -Pathway, 
      names_to = "Sample_ID", 
      values_to = "Score"
    )
  
  # Step 3: Merge with annotation data
  merged_data <- singscore_long %>%
    left_join(Drug_DMSO_data, by = "Sample_ID") 
 
  pval_data <- merged_data %>%
    group_split(Pathway, Timepoint) %>% # Split data by Pathway and Timepoint
    purrr::map_dfr(function(df) {
      cat("Processing Pathway:", unique(df$Pathway), 
          " Timepoint:", unique(df$Timepoint), "\n")
      
      control_group <- "DMSO" # Specify control group
      other_drugs <- unique(df$Drug[df$Drug != control_group])%>%as.character()
      
      if (length(other_drugs) == 0) {
        cat("No drugs found for comparison.\n")
        return(NULL)
      }
      
      # Perform t-tests for each drug vs. control group
      test_results <- purrr::map(other_drugs, function(drug) {
        subset_df <- df %>% filter(Drug %in% c(control_group, drug))
        if (nrow(subset_df) == 0) {
          cat("Subset is empty for Drug:", drug, "\n")
        }
        
        # Ensure exactly 2 levels of Drug
        if (n_distinct(subset_df$Drug) == 2) {
          test <- t.test(Score ~ Drug, data = subset_df)
          data.frame(
            Drug = drug,
            Pathway = unique(df$Pathway),
            Timepoint = unique(df$Timepoint),
            p_value = test$p.value
          )
        } else {
          NULL
        }
      })
      
      # Combine results for all drugs
      bind_rows(test_results)%>%
        mutate(Significance = ifelse(p_value < 0.05, "*", ifelse(p_value < 0.01, "**",""))) # Add significance markers
    })
  
   
  # Step 5: Aggregate data (calculate mean and variance)
  aggregated_data <- merged_data %>%
    group_by(Pathway, Drug, Timepoint) %>%
    summarise(
      Mean_Score = mean(Score, na.rm = TRUE),
      Variance = var(Score, na.rm = TRUE),
      .groups = "drop"
    )
  pval_data2=pval_data%>%left_join(aggregated_data)%>%mutate(y_pos=Mean_Score+Variance)
  
  # Step 6: Generate line plot with vertical error bars
  p <- ggplot(aggregated_data, aes(x = Timepoint, y = Mean_Score, group = Drug, color = Drug)) +
    geom_line(size = 1) + # Line for mean scores
    geom_point(size = 2) + # Points for mean scores
    geom_errorbar(aes(ymin = Mean_Score - Variance, ymax = Mean_Score + Variance), 
                  width = 0.2, size = 0.4) + # Vertical error bars
    facet_wrap(~Pathway, scales = "free_y", ncol = 3) +
    theme_minimal() +
    labs(
      title = Title,
      x = "Timepoint",
      y = "Mean Singscore (± Variance)",
      color = "Drug"
    ) +scale_color_manual(values = Drug_colors)+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate x-axis labels
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.spacing = unit(1, "lines"), # Add space between facets
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add black border
      strip.background = element_rect(fill = "white", color = "black"), # Frame facet labels
      strip.text = element_text(face = "bold", size = 10) # Bold facet labels
    )
  
  # Step 7: Add significance annotations (* for p-value < 0.05)
  p <- p + geom_text(
    data = pval_data2 %>% filter(Significance != ""),
    aes(x = Timepoint, y = y_pos, label = Significance, group = Drug),
    inherit.aes = FALSE,
    color = "black",
    size = 8
  )
  return(p)
}


Self_define_term=function(Pathway_names,Pathway_Data,Gene_input){
  necessary_cols <- c("Term_Description", "lowest_p", "Up_regulated", 
                      "Down_regulated", "no_change","group")
  out=list()
  for(i in Pathway_names){
    Gene_in_pathway=Pathway_Data[[i]]%>%as.character()
    Gene_list=Gene_input%>%dplyr::filter(Gene.symbol%in%Gene_in_pathway)%>%mutate(Change=ifelse(logFC>0,"Up_regulated","Down_regulated"))
    
    Up_Genes=Gene_list%>%dplyr::filter(Change=="Up_regulated")
    DN_Genes=Gene_list%>%dplyr::filter(Change!="Up_regulated")
    Gene_anno=data.frame(Term_Description=i,lowest_p=NA,Up_regulated=toString  (Up_Genes$Gene.symbol, collapse = ", "),
                         Down_regulated=toString  (DN_Genes$Gene.symbol, collapse = ", "))
    out[[i]]=Gene_anno
    
  }
  out_all=do.call("rbind",out)
}

DEGs_heatmap <- function(pathway_data) {
  # Extract gene names
  genes <- unique(unlist(strsplit(pathway_data$Up_regulated, ", ")))
  genes <- union(genes, unlist(strsplit(pathway_data$Down_regulated, ", ")))
  
  # Initialize list to store expression matrices for each sample
  expression_matrix <- list()
  
  # Fill expression matrix with expression values for each sample
  for (i in 1:nrow(pathway_data)) {
    sample <- pathway_data[i, ]
    upregulated <- unlist(strsplit(sample$Up_regulated, ", "))
    downregulated <- unlist(strsplit(sample$Down_regulated, ", "))
    no_change <- unlist(strsplit(sample$no_change, ", "))
    
    for (I in sample$group ) {
      DATA=sample%>%dplyr::filter(group==I)
      # Create matrix for current sample
      sample_matrix <- matrix(0, nrow = length(genes), ncol = 1)
      rownames(sample_matrix) <- genes
      colnames(sample_matrix) <-I  # Assign group names as column names
      upregulated_genes=unlist(strsplit(DATA$Up_regulated, ", "))
      downregulated_genes <- unlist(strsplit(DATA$Down_regulated, ", "))
      # Fill sample matrix with expression values
      for (gene_index in 1:length(genes)) {
        gene <- genes[gene_index]
        if (gene %in% upregulated_genes) {
          sample_matrix[gene, I] <- 1
        } else if (gene %in% downregulated_genes) {
          sample_matrix[gene, I] <- -1
        }
      }
      
    } 
    
    # Store sample matrix in the list
    expression_matrix[[i]] <- sample_matrix
  }
  
  # Merge expression matrices using cbind
  merged_matrix <- do.call(cbind, expression_matrix)
  
  # Create heatmap
  OUT=pheatmap::pheatmap(merged_matrix,
                         cluster_rows = TRUE,  # Cluster rows (genes)
                         cluster_cols = FALSE,
                         color = c("blue", "grey", "red"),
                         breaks = c(-0.9, -0.1, 0.9),
                         scale = "none",cellheight =10,   treeheight_row = 5,legend = F,
                         main = paste("Heatmap for", unique(pathway_data$Term_Description)))
}

pheatmap_pathway_3<- function (expression_data, Gene_select,pathway_name,Annotation) {
  # Gene_select<-Pathway_list[,pathway_name]%>%unique()
  Gene_select_anno= expression_data[Gene_select,] 
  pheatmap::pheatmap(Gene_select_anno, cellwigermline= 2, cellheight =10, scale = "row",
                     treeheight_row = 5,main=pathway_name,
                     show_rownames = T,show_colnames = F,
                     annotation_col= Annotation,  # Title
                     gaps_col = c(3,3,8),
                     # annotation_row=Annotation,
                     #annotation_legend=Label_def,
                     cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean")
  
}


pheatmap_pathway_4<- function (expression_data, Gene_select,pathway_name,Annotation,cutcol) {
  # Gene_select<-Pathway_list[,pathway_name]%>%unique()
  Gene_select_anno= expression_data[Gene_select,] 
  pheatmap::pheatmap(Gene_select_anno, cellwigermline= 2, cellheight =10, scale = "row",
                     treeheight_row = 5,main=pathway_name,
                     show_rownames = T,show_colnames = F,
                     annotation_col= Annotation,  # Title
                     gaps_col = cutcol,
                     # annotation_row=Annotation,
                     #annotation_legend=Label_def,
                     cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean")
  
}

merge_and_rename_genes <- function(pathway_out_1, pathway_out_2,group1_name="Group1",group2_name="Group2") {
  Out_2 <- data.frame()
  
  # Add missing genes column and rename for pathway_out_1
  for (pathway_name in pathway_out_1$Term_Description) {
    # Extract genes for the current pathway
    Data1 <- pathway_out_1 %>% 
      dplyr::filter(Term_Description == pathway_name) %>% 
      mutate(group = group1_name)
    
    Data2 <- pathway_out_2 %>% 
      dplyr::filter(Term_Description == pathway_name) %>% 
      mutate(group = group2_name)
    
    # Split Up_regulated and Down_regulated by ","
    Data1_Up_regulated <- stringr::str_split(Data1$Up_regulated, ", ")%>%unlist() %>%as.character()
    Data1_Down_regulated <- stringr::str_split(Data1$Down_regulated, ", ")%>%unlist() %>%as.character()
    Data2_Up_regulated <- stringr::str_split(Data2$Up_regulated, ", ",)%>%unlist() %>%as.character()
    Data2_Down_regulated <- stringr::str_split(Data2$Down_regulated, ", ", )%>%unlist() %>%as.character()
    
    # Merge genes for pathway_out_1 and pathway_out_2 for the current pathway
    merged_genes_1 <- c(Data1_Up_regulated, Data1_Down_regulated)
    merged_genes_2 <- c(Data2_Up_regulated, Data2_Down_regulated)
    all_genes=unique(c(merged_genes_1,merged_genes_2))
    
    missing_1 <- setdiff( all_genes, merged_genes_1)
    Data1$no_change <- toString  (missing_1, collapse = ", ")
    
    missing_genes_2 <- setdiff( all_genes, merged_genes_2)
    Data2$no_change <- toString(missing_genes_2,collapse = ", ")
    
    out <- rbind(Data1, Data2)
    Out_2 <- rbind(Out_2, out)
  }
  Out_2$Term_Description=toupper(Out_2$Term_Description)
  return(Out_2)
}


Term_gene_graph_3_update=function (result_df, num_terms = 10, layout = "stress", use_description = FALSE, 
                                   node_size = "num_genes", node_colors = c("#E5D7BF", "green", "red", "grey"),facet_variable = NULL) {
  # For paired plot
  
  if (!is.numeric(num_terms) & !is.null(num_terms)) {
    stop("`num_terms` must either be numeric or NULL!")
  }
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }
  ID_column <- ifelse(use_description, "Term_Description", 
                      "ID")
  val_node_size <- c("num_genes", "p_val")
  if (!node_size %in% val_node_size) {
    stop("`node_size` should be one of ", paste(dQuote(val_node_size), 
                                                collapse = ", "))
  }
  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", 
                      "Down_regulated", "no_change","group")  # Include "no_change"
  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "), 
                 "must be present in `results_df`!"), collapse = " "))
  }
  if (length(node_colors) != 4) {  # Corrected length check
    stop("`node_colors` must contain exactly 4 colors")  # Updated error message
  }
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) {
      num_terms <- NULL
    }
  }
  Group=unique(result_df$group)
  
  ####################
  glist=list()
  input=result_df
  for (G in Group)    {
    result_df<- input[order(input$Term_Description, decreasing = FALSE),  ]%>%dplyr::filter(group==G)
    if (!is.null(num_terms)) {
      result_df <- result_df[1:num_terms, ]
    }
    graph_df <- data.frame()
    for (i in base::seq_len(nrow(result_df))) {
      up_genes <- unlist(strsplit(result_df$Up_regulated[i], 
                                  ", "))
      down_genes <- unlist(strsplit(result_df$Down_regulated[i], 
                                    ", "))
      
      no_change= unlist(strsplit(result_df$no_change[i], 
                                 ", "))
      
      for (gene in c(up_genes, down_genes,no_change)) {
        graph_df <- rbind(graph_df, data.frame(Term = result_df[i, 
                                                                ID_column], Gene = gene))%>%arrange(Gene)
      }
    }
    up_genes <- unlist(lapply(result_df$Up_regulated, function(x) strsplit(x, ", ")))
    DN_genes=  unlist(lapply(result_df$Down_regulated, function(x) strsplit(x, ", ")))
    
    g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
    cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
    cond_up_gene <- names(igraph::V(g)) %in% up_genes
    cond_down_gene=  names(igraph::V(g)) %in% DN_genes
    node_type <- ifelse(cond_term, "term", 
                        ifelse(cond_up_gene, "up", 
                               ifelse(cond_down_gene, "down", "no_change")))
    node_type <- factor(node_type, levels = c("term", "up", "down", "no_change"))
    node_type <- factor(node_type, levels = c("term", "up", "down","no_change"))
    igraph::V(g)$type <- node_type
    type_descriptions <- c(term = "enriched term", up = "up-regulated gene", 
                           down = "down-regulated gene", no_change = "no_change gene")  # Include "no_change" description
    type_descriptions <- type_descriptions[levels(node_type)]
    names(node_colors) <- c("term", "up", "down", "no_change")  # Added "no_change"
    node_colors <- node_colors[levels(node_type)]
    if (node_size == "num_genes") {
      sizes <- igraph::degree(g)
      sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
      size_label <- "# genes"
    }
    else {
      idx <- match(names(igraph::V(g)), result_df[, ID_column])
      sizes <- -log10(result_df$lowest_p[idx])
      sizes[is.na(sizes)] <- 2
      size_label <- "-log10(p)"
    }
    
    igraph::V(g)$size <- sizes
    igraph::V(g)$label.cex <- 0.5
    igraph::V(g)$frame.color <- "gray"
    glist[[G]]=g
  }    
  
  ##################################################
  
  plist=list()
  
  for (s in Group) {
    g=glist[[s]]
    #xy <- graphlayouts::layout_with_stress(glist[[1]])
    p <- ggraph::ggraph(g, layout = layout)
    
    p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
    p <- p + ggraph::geom_node_point(ggplot2::aes(color = .data$type, 
                                                  size = .data$size))
    p <- p + ggplot2::scale_size(range = c(5, 10), breaks = round(seq(round(min(igraph::V(g)$size)), 
                                                                      round(max(igraph::V(g)$size)), length.out = 4)), name = size_label)
    p <- p + ggplot2::theme_void()
    p <- p + suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name), 
                                                     nudge_y = 0.2, repel = TRUE, max.overlaps = 20))
    p <- p + ggplot2::scale_color_manual(values = node_colors, 
                                         name = NULL, labels = type_descriptions)
    if (is.null(num_terms)) {
      p <- p + ggplot2::ggtitle("Term-Gene Graph")
    }
    else {
      p <- p + ggplot2::ggtitle("Term-Gene Graph", subtitle = paste(c("Top", 
                                                                      num_terms, "terms"), collapse = " "))
    }
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                            plot.subtitle = ggplot2::element_text(hjust = 0.5))+ggtitle(s)
    plist[[s]]=p
  }
  
  return(plist)
}



Enrich_dotplot=function(FGSEA_out){
  library(ggbreak)
  FGSEA_out$NES=round(FGSEA_out$NES,2)
  #mycolor=RColorBrewer::brewer.pal(3,"YlOrRd")
  data=FGSEA_out%>%mutate(Norm_padj=-log10(padj))%>%mutate(State=ifelse(NES>0&padj<0.05,"Activated",ifelse(NES<0&padj<0.05,"Inhibited", ifelse(NES>0&padj<0.07,"P_Activated",ifelse(NES<0&padj<0.07,"P_Inhibited","no_change")))))%>%arrange(NES)
  data$pathway=factor(data$pathway,levels =unique(data$pathway ))
  # Up_data=data%>%dplyr::filter(NES>0)%>%arrange(NES)
  #   DN_data=data%>%dplyr::filter(NES<0)
  #   data=rbind(DN_data,Up_data)
  # data$pathway=factor(data$pathway,levels =data$pathway )
  MAX_NES=max(data$NES)%>%ceiling()
  min_NES=min(data$NES)%>%floor()
  # ggplot(data) +geom_point(aes(x = NES, y = pathway, color = -log10(padj),size=gene_counts))+scale_size("gene_counts") +
  color_set= setNames(c('red','green',"wheat","turquoise" ,'grey'),c("Activated", "Inhibited","P_Activated","P_Inhibited","no_change"))
  
  P= ggplot(data) +geom_point(aes(x = NES, y = pathway, color =State ,size=Norm_padj)) +
    theme_bw() +labs(x="Normalized enrichment score",y="")+
    theme(axis.text.x = element_text(size=rel(1.15)),
          axis.title = element_text(size=rel(1.15)),plot.title = element_text(hjust=0.5, 
                                                                              face = "bold"),legend.title = element_text(size=rel(1.15),hjust=0.5, face="bold"),axis.text.y = element_text(size =15))+
    scale_size_continuous (name="-log10(padj)",range=c(2,12))+scale_color_manual(name="State",values=color_set)+facet_wrap(~group,nrow=1)
  # if (min_NES>0|MAX_NES<0) {
  #   return(P)
  #   } else {
  #   Breaks=c(min_NES,-1,1,MAX_NES)
  #  return(P+scale_x_break(c(-1, 1), scales = 1.5))
  # xlab("Normalized enrichment score") +ylab("")+scale_x_break(c(-1, 1), scales = 1.5)
  #   }
}


Enrich_dotplot_update=function(FGSEA_out){
  library(ggbreak)
  FGSEA_out$NES=round(FGSEA_out$NES,2)
  #mycolor=RColorBrewer::brewer.pal(3,"YlOrRd")
  data <- FGSEA_out %>%
    mutate(
      Norm_padj = -log10(padj),  # Transform p-adjusted values
      State = case_when(
        NES > 0 & padj < 0.05 ~ "Activated",
        NES < 0 & padj < 0.05 ~ "Inhibited",
        TRUE ~ "no_change" # Default
      )
    ) %>%
    arrange(NES)
 
  data$pathway=factor(data$pathway,levels =unique(data$pathway ))
  MAX_NES=max(data$NES)%>%ceiling()
  min_NES=min(data$NES)%>%floor()
  # ggplot(data) +geom_point(aes(x = NES, y = pathway, color = -log10(padj),size=gene_counts))+scale_size("gene_counts") +
  color_set= setNames(c('red','green','grey'),c("Activated", "Inhibited","no_change"))
    P= ggplot(data) +geom_point(aes(x = NES, y = pathway, color =State ,size=Norm_padj)) +
    theme_bw() +labs(x="Normalized enrichment score",y="")+
    theme(axis.text.x = element_text(size=rel(1.15)),
          axis.title = element_text(size=rel(1.15)),plot.title = element_text(hjust=0.5, 
                                                                              face = "bold"),legend.title = element_text(size=rel(1.15),hjust=0.5, face="bold"),axis.text.y = element_text(size =15))+
    scale_size_continuous (name="-log10(padj)",range=c(2,12))+scale_color_manual(name="State",values=color_set)+facet_wrap(~group,nrow=1)

}


Common_sig_pathway_NES <- function(data, Pathway_key = "MAPK") {
  library(dplyr)
  library(tidyr)
  
  # Use Pathway_key in grepl to filter pathways dynamically
  data_pathway <- lapply(data, function(x) {
    x %>% dplyr::filter(grepl(Pathway_key, pathway))
  })
  
  # Extract unique significant pathway names (padj < 0.05) from all datasets
  all_sig_pathway_name <- lapply(data_pathway, function(x) {
    x %>%
      as.data.frame() %>%
      dplyr::filter(padj < 0.05) %>%
      .$pathway
  }) %>%
    unlist() %>%
    unique()
  
  # Initialize a list to store the NES data for significant pathways
  NES_sig_pathway_pickup <- list()
  
  # Loop over each dataset in data
  for (i in names(data)) {
    Nes_data <- data[[i]]  # Get the NES data for the current group
    
    # Join with the significant pathway names and handle missing values
    NES_sig_pathway_pickup[[i]] <- data.frame(pathway = all_sig_pathway_name) %>%
      left_join(Nes_data, by = "pathway") %>%
      mutate_at("NES", ~ tidyr::replace_na(., 0)) %>%  # Replace missing NES values with 0
      mutate_at("padj", ~ tidyr::replace_na(., 1)) %>%  # Replace missing padj values with 1
      mutate(Group = i)  # Add group identifier
  }
  NES_sig_pathway_pickup2=do.call("rbind", NES_sig_pathway_pickup)
  
}
run_pathview_analysis <- function(test_data, pathway_name, pathway_id, output_dir,outname) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Map gene symbols to Entrez IDs
  test_data$entrez_id <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, 
    keys = test_data$gene_symbol, 
    column = "ENTREZID", 
    keytype = "SYMBOL", 
    multiVals = "first"
  )
  
  # Remove rows with missing Entrez IDs
  test_data <- test_data %>% dplyr::filter(!is.na(entrez_id))
  
  # Prepare gene data for pathview
  gene_data <- test_data$log2FoldChange
  names(gene_data) <- test_data$entrez_id
  
  # Set output filename based on pathway name
  
  # Run pathview for the given pathway
  pathview_result <- pathview::pathview(
    gene.data = gene_data, 
    pathway.id = pathway_id,  # Provide pathway ID (e.g., "04110")
    species = "hsa", 
    out.dir = output_dir, 
    bins = list(gene = 5, cpd = 5),  # Adjust bins if needed
    kegg.native = TRUE
  )
  
  outfie_name=file.path(".",paste0(pathway_id,".pathview.png"))
  newname=file.path(output_dir,outname)
  file.rename(outfie_name,newname)
  png=file.path(".",paste0(pathway_id,".png"))
  xml=file.path(".",paste0(pathway_id,".xml"))
  file.remove(png)
  file.remove(xml)
  # Return the pathview result object
  
}


run_pathview_gene_symbol <- function(test_data, pathway_name, pathway_id, output_dir,outname) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Remove rows with missing Entrez IDs
  test_data <-  test_data[,c("gene_symbol","log2FoldChange")]
  
  # Prepare gene data for pathview
  gene_data <- test_data$log2FoldChange
  names(gene_data) <- test_data$gene_symbol
  
  # Set output filename based on pathway name
  
  # Run pathview for the given pathway
  pathview_result <- pathview::pathview(
    gene.data = gene_data, 
    pathway.id = pathway_id,  # Provide pathway ID (e.g., "04110")
    species = "hsa", 
    out.dir = output_dir, 
    bins = list(gene = 5, cpd = 5),  # Adjust bins if needed
    kegg.native = TRUE,gene.idtype = "SYMBOL"
  )
  
  outfie_name=file.path(".",paste0(pathway_id,".pathview.png"))
  newname=file.path(output_dir,outname)
  file.rename(outfie_name,newname)
  png=file.path(".",paste0(pathway_id,".png"))
  xml=file.path(".",paste0(pathway_id,".xml"))
  file.remove(png)
  file.remove(xml)
  # Return the pathview result object
  
}

run_pathview_group_analysis <- function(test_data1, test_data2, pathway_name, pathway_id, output_dir, outname) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Helper function to map gene symbols to Entrez IDs
  map_entrez_ids <- function(data) {
    data$entrez_id <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, 
      keys = data$gene_symbol, 
      column = "ENTREZID", 
      keytype = "SYMBOL", 
      multiVals = "first"
    )
    data <- data %>% filter(!is.na(entrez_id))
    return(data)
  }
  
  # Process both datasets
  data1 <- map_entrez_ids(test_data1)
  data2 <- map_entrez_ids(test_data2)
  
  # Prepare merged gene data (columns are groups)
  merged_data <- merge(data1[, c("entrez_id", "log2FoldChange")], 
                       data2[, c("entrez_id", "log2FoldChange")], 
                       by = "entrez_id", suffixes = c("_group1", "_group2"),all=T)
  
  merged_data[is.na(merged_data)] <- 0
 
  gene_data <- as.matrix(merged_data[, -1])  # Convert to matrix
  rownames(gene_data) <- merged_data$entrez_id
  
  # Run pathview to visualize both groups in the same pathway
  pathview_result <- pathview::pathview(
    gene.data = gene_data,
    pathway.id = pathway_id,
    species = "hsa",
    out.dir = output_dir,
    bins = list(gene = 5, cpd = 5),
    kegg.native = TRUE,
    multi.state = TRUE,   # Handle multiple states (group1 and group2)
    same.layer = FALSE    # Display different groups on different layers
  )
  
  # Rename and clean up output files
  png_name <- file.path(".", paste0(pathway_id, ".pathview.multi.png"))
  new_name <- file.path(output_dir, outname)
  file.rename(png_name, new_name)
  
  # Remove extra files generated during visualization
  file.remove(file.path(".", paste0(pathway_id, ".png")))
  file.remove(png_name)
  file.remove(file.path(".", paste0(pathway_id, ".xml")))
  
  # Return the pathview result

}



run_pathview_group_gene<- function(test_data1, test_data2, pathway_name, pathway_id, output_dir, outname) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Helper function to map gene symbols to Entrez IDs
  map_entrez_ids <- function(data) {
    data$entrez_id <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, 
      keys = data$gene_symbol, 
      column = "ENTREZID", 
      keytype = "SYMBOL", 
      multiVals = "first"
    )
    data <- data %>% filter(!is.na(entrez_id))
    return(data)
  }
  
  # Process both datasets
  data1 <- map_entrez_ids(test_data1)
  data2 <- map_entrez_ids(test_data2)
  
  # Prepare merged gene data (columns are groups)
  merged_data <- merge(data1[, c("entrez_id", "log2FoldChange")], 
                       data2[, c("entrez_id", "log2FoldChange")], 
                       by = "entrez_id", suffixes = c("_group1", "_group2"),all=T)
  
  merged_data[is.na(merged_data)] <- 0
  
  gene_data <- as.matrix(merged_data[, -1])  # Convert to matrix
  rownames(gene_data) <- merged_data$entrez_id
  
  # Run pathview to visualize both groups in the same pathway
  pathview_result <- pathview::pathview(
    gene.data = gene_data,
    pathway.id = pathway_id,
    species = "hsa",
    out.dir = output_dir,
    bins = list(gene = 5, cpd = 5),
    kegg.native = TRUE,
    multi.state = TRUE,   # Handle multiple states (group1 and group2)
    same.layer = FALSE    # Display different groups on different layers
  )
  
  # Rename and clean up output files
  png_name <- file.path(".", paste0(pathway_id, ".pathview.multi.png"))
  new_name <- file.path(output_dir, outname)
  file.rename(png_name, new_name)
  
  # Remove extra files generated during visualization
  file.remove(file.path(".", paste0(pathway_id, ".png")))
  file.remove(png_name)
  file.remove(file.path(".", paste0(pathway_id, ".xml")))
  
  # Return the pathview result
  
}

run_pathview_group_gene_symbol <- function(test_data1, test_data2, pathway_name, pathway_id, output_dir, outname) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Merge both datasets by gene symbol, filling missing values with 0
  merged_data <- merge(
    test_data1[, c("gene_symbol", "log2FoldChange")], 
    test_data2[, c("gene_symbol", "log2FoldChange")], 
    by = "gene_symbol", 
    suffixes = c("_group1", "_group2"), 
    all = TRUE
  )
  
  merged_data[is.na(merged_data)] <- 0  # Set missing values to 0
  
  # Prepare the gene data matrix
  gene_data <- as.matrix(merged_data[, -1])  # Convert to matrix
  rownames(gene_data) <- merged_data$gene_symbol  # Use gene symbols as rownames
  
  # Run pathview visualization
  pathview_result <- pathview::pathview(
    gene.data = gene_data,
    pathway.id = pathway_id,
    species = "hsa",
    out.dir = output_dir,
    bins = list(gene = 5, cpd = 5),
    kegg.native = TRUE,
    multi.state = TRUE,   # Display both groups in the same pathway
    same.layer = FALSE,   # Separate layers for each group
    gene.idtype = "SYMBOL"  # Use gene symbols directly
  )
  
  # Rename and clean up the output files
  png_name <- file.path(".", paste0(pathway_id, ".pathview.multi.png"))
  new_name <- file.path(output_dir, outname)
  file.rename(png_name, new_name)
  
  # Remove extra files generated during visualization
  file.remove(file.path(".", paste0(pathway_id, ".png")))
  file.remove(png_name)
  file.remove(file.path(".", paste0(pathway_id, ".xml")))
  
  # Return the pathview result object
 
}


# Define the function
plot_pathway_heatmap <- function(data, pathway_col, padj_col, NES_col, group_col, drug_col) {
  # Ensure the input data is a data frame
  if (!is.data.frame(data)) stop("Input data must be a data frame")
  
  # Check if the required columns exist
  required_cols <- c(pathway_col, padj_col, NES_col, group_col, drug_col)
  if (!all(required_cols %in% colnames(data))) {
    stop("Specified columns do not exist in the data frame")
  }
  
  # Generate all possible combinations of drugs, including single drugs
  unique_drugs <- unique(data[[drug_col]])
  
  # Compute metrics for each drug combination
  metrics <- data %>%
    group_by(.data[[pathway_col]]) %>%
    summarise(
      # Get unique drugs associated with the pathway
      drugs = paste(sort(unique(.data[[drug_col]])), collapse = "+"),  # Concatenate drugs
      sig_count = sum(.data[[padj_col]] < 0.05),                       # Count significant padj values
      NES_sum = n(),NES_value=sum(.data[[NES_col]]),Var=prod(.data[[NES_col]])                           # Sum of absolute NES values
    ) %>%
    ungroup()
  
  
 
  # Sort pathways based on significance and NES
  sorted_pathways <- metrics %>%
    arrange(desc(NES_sum), drugs,desc(Var),desc(sig_count),desc(NES_value)) %>% # Sort by significance count and NES sum
    mutate(!!sym(pathway_col) := factor(.data[[pathway_col]], levels = unique(.data[[pathway_col]])))
 
   # Merge sorted pathways with the input data
  plot_data <- data %>%
    mutate(
      color_val = ifelse(.data[[padj_col]] > 0.05, NA, .data[[NES_col]]) # Color grey if padj > 0.05
    )
  

  
  plot_data$pathway=factor( plot_data$pathway,levels =rev(sorted_pathways$pathway ))
  heatmap_plot <- ggplot(plot_data, aes(
    x = .data[[group_col]], 
    y = .data[[pathway_col]], 
    fill = color_val, 
    size = -log10(.data[[padj_col]])
  )) +
    geom_point(shape = 21, color = "black") +  # Circle with border
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-6, 6),
      na.value = "grey"  # Set NA values to grey
    ) +
    scale_size(range = c(2, 10)) +  # Point size based on -log10(padj)
    labs(
      title = "Pathway Heatmap",
      x = "Group",
      y = "Pathway",
      fill = "NES",
      size = "-log10(padj)"
    ) +
    # Enhanced theme
    theme_minimal(base_size = 12) +  # Clean theme with base size for text
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Title in center and bold
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),      # Rotate x-axis labels
      axis.text.y = element_text(size = 10),                            # Y-axis text size
      panel.grid.major = element_blank(),                               # Remove major gridlines
      panel.grid.minor = element_blank(),                               # Remove minor gridlines
      panel.background = element_rect(fill = "white"),                  # White background
      legend.position = "right",                                        # Legend to the right
      legend.title = element_text(size = 12),                           # Legend title size
      legend.text = element_text(size = 10)                             # Legend text size
    )
  
  # Return the plot
  return(heatmap_plot)
}


# Define the function
get_synapse_annotations <- function(syn_id) {
  # Login to Synapse (ensure credentials are set up)
 
  
  # Get the folder entity
  folder <- synGet(syn_id)
  
  # Check if it's a folder
  if (folder$entityType != "org.sagebionetworks.repo.model.Folder") {
    stop("The provided Synapse ID is not a folder.")
  }
  
  # Get the list of files in the folder
  files <- synGetChildren(syn_id)
  
  # Initialize a list to store annotations
  annotations_list <- list()
  
  # Loop through each file in the folder
  for (file in files) {
    # Get file details
    entity <- synGet(file$id)
    
    # Retrieve annotations
    annotations <- synGetAnnotations(entity)
    
    # Combine file ID, name, and annotations into a data frame
    annotations_list[[file$name]] <- data.frame(
      FileName = file$name,
      SynID = file$id,
      Annotations = sapply(annotations, toString) # Flatten annotations
    )
  }
  
  # Combine all annotations into one data frame
  annotations_df <- bind_rows(annotations_list)
  
  # Return the result
  return(annotations_df)
}

