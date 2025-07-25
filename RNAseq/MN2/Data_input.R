

# Define your own config

getwd()
config <- config::get(file = 'config.yml')
# 
Working_path=config$Path
base_output_dir <- paste0(Working_path,"/","Report")
email=config$Email
Authtoken=config$Token_file
knitr::opts_chunk$set(echo = F,warnings =F, message=F,include=F)
options(knitr.table.format = 'markdown',encoding = 'UTF-8',warnings =F, message=F ) 
# Check and create output directory if it doesn't exist
output_dir <- file.path(Working_path, "Data_preparation")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
Reference_path=config$Reference_path

dir(Reference_path)

getwd()
# load function
source(file.path("/Volumes/lyu.yang/MAC_1/Github/MPNST-PDX-MT/RNAseq","Rfunction.R" ))
 

if (! file.exists( file.path(Reference_path,"Human_GSEA","Protein_coding_genelist.rdata"))){
  syID="syn63828384"
  synapser::synGet(syID,downloadLocation = Reference_path,ifcollision="keep.local")
  
}
load(file.path(Reference_path,"Human_GSEA","Protein_coding_genelist.rdata"))



# Download raw counts, TPM and annoation data from synapse

# Batch1 data dowloading from synapse
Rawdata_path=file.path(Working_path,"Rawdata_update")
Processed_path=file.path(Working_path,"Processed_data")
dir.create(Processed_path)
Rawdata_path_b1=file.path(Rawdata_path,"b1")
# login synapse
Token=read.table(Authtoken)%>%as.character()
synapser::synLogin(email=email,authToken=Token)
# download files from synapse
Rawcounts_syID_b1="syn64608365"
RawTPM_syID_b1="syn64608366"
Annotation_syID_b1="syn64608368"

for (i in c(Rawcounts_syID_b1,RawTPM_syID_b1,Annotation_syID_b1)){
  synapser::synGet(i,downloadLocation = Rawdata_path_b1,ifcollision="keep.local")
  
  
}

Annotation_b1=read.table(file.path(Rawdata_path_b1,"01172025_Batch1_Samplemap.csv"),sep=",",header = T)
Annotation_b1$Drug=factor(Annotation_b1$Drug,levels = c("DMSO","MIR"))
Annotation_b1$Timepoint=factor(Annotation_b1$Timepoint,levels=c("0h","8h","24h"))
# Download Batch2 data from synapse

Rawdata_path=file.path(Working_path,"Rawdata_update")
Rawdata_path_b2=file.path(Rawdata_path,"b2")
# login synapse
Token=read.table(Authtoken)%>%as.character()
synapser::synLogin(email=email,authToken=Token)
# download files from synapse
Rawcounts_syID_b2="syn64608369"
RawTPM_syID_b2="syn64608370"
#Annotation_syID_b2="syn64608372" # DMSO_batch_version1
Annotation_syID_b2="syn66302373"

for (i in c(Rawcounts_syID_b2,RawTPM_syID_b2,Annotation_syID_b2)){
  synapser::synGet(i,downloadLocation = Rawdata_path_b2,ifcollision="keep.local")
  
  
}



Annotation_b2=read.table(file.path(Rawdata_path_b2,"01172025_Batch2_Samplemap_swap_DMSO.csv"),sep=",",header = T)%>%dplyr::mutate(Drug=ifelse(Full_name=="Trametinib","TRAM",Drug))
dir(Rawdata_path_b2)
Annotation_b2$Drug=factor(Annotation_b2$Drug)
Annotation_b2$Timepoint=factor(Annotation_b2$Timepoint,levels=c("0h","8h","24h"))
Drug_table_syID="syn64608501"

synapser::synGet(Drug_table_syID,downloadLocation = Rawdata_path,ifcollision="keep.local")

All_drug_name=unique(c(Annotation_b1$Drug,Annotation_b2$Drug))

# human GMT files download
GMT_folder_id="	syn63828408"
GSEA_path=file.path(Reference_path,"Human_GSEA")
download_synapse_folder(folder_id=GMT_folder_id, destination_path=GSEA_path)

Drug_anno=readxl::read_excel(file.path(Rawdata_path,"Drug_annotation.xlsx"))%>%arrange(Drug_class)
Drug_name = Drug_anno$Drug %>% unique()%>%na.omit()
Drug_name=factor(Drug_name,levels = c("DMSO" ,"PAL",  "RIB" , "DEC" , "VOR" , "TRAM"  ,"MIR"  ,"SEL"  ,"OLA",  "TNO",  "RMC" , "CAP"  ,"TRAB", "DOX"  ,"IFO"  ,"IRI"))
setdiff(All_drug_name,Drug_name)


length(Drug_name)
# Ensure you have 14 unique Drug names (1 for black, 13 other colors)
if (length(Drug_name) != 16) {
  stop("The number of unique drugs should be 16, including the one colored black.")
}

# Use RColorBrewer to get 13 colors from the Paired palette
#paired_colors =RColorBrewer:: brewer.pal(n = 12, name = "Paired")
paired_colors <- c(
  "#A6CEE3",  # light blue
  "#1F78B4",  # dark blue
  "#B2DF8A",  # light green
  "#33A02C",  # medium green
  "#FB9A99",  # light red/pink
  "#E31A1C",  # bright red
  "#FDBF6F",  # orange-yellow
  "#FF7F00",  # orange
  "#CAB2D6",  # lavender
  "#6A3D9A",  # dark purple
  "#1B9E77",  # teal
  "#B15928"   # brown
)
#extra_color =RColorBrewer:: brewer.pal(n = 8, name = "Set3")[1:3]  # Take one color from another palette (e.g., Set3)

extra_color =c(
  "#8DD3C7",  # light cyan
  "#7CAE00",  # lime green
  "#BEBADA"   # light purple
)
# Combine black and the 13 colors
drug_colors = c("black", paired_colors, extra_color)

# Ensure the color vector matches the order of Drug_name
names(drug_colors) = sort(Drug_name)

Timepoint_colors= RColorBrewer::brewer.pal(5,"Set1")[2:4]
names(Timepoint_colors)=c("0h","8h","24h")

## load counts and TPM
set.seed(12345)
dir(Rawdata_path_b1)
Raw_counts_b1=read.table(file.path(Rawdata_path_b1,"salmon.merged.gene_counts_corrected_id.csv"),sep=",",header = T,row.names = 1)
TPM_counts_b1=read.table(file.path(Rawdata_path_b1,"salmon.merged.gene_tpm_corrected_id.csv"),sep=",",header = T,row.names = 1)
Annotation_b1=read.table(file.path(Rawdata_path_b1,"01172025_Batch1_Samplemap.csv"),sep=",",header = T)
Annotation_b1$Drug=factor(Annotation_b1$Drug,levels = c("DMSO","MIR"))
Annotation_b1$Timepoint=factor(Annotation_b1$Timepoint,levels=c("0h","8h","24h"))
TPM_b1_2=TPM_counts_b1%>%.[,c(2:ncol(.))]%>%dplyr::filter(gene_name%in%Protein_coding_genelist)%>%group_by(gene_name)%>%summarise_all(mean)%>%as.data.frame()%>%matrix.please()%>%as.data.frame()
Counts_b1_2=Raw_counts_b1%>%.[,c(2:ncol(.))]%>%dplyr::filter(gene_name%in%Protein_coding_genelist)%>%group_by(gene_name)%>%summarise_all(mean)%>%as.data.frame()%>%matrix.please()%>%as.data.frame()%>%mutate_all(.,as.integer)

Raw_counts_b2=read.table(file.path(Rawdata_path_b2,"salmon.merged.gene_counts_corrected_id.csv"),sep=",",header = T)
TPM_counts_b2=read.table(file.path(Rawdata_path_b2,"salmon.merged.gene_tpm.tsv_corrected_ID.csv"),sep=",",header = T)


Annotation_b2$Drug=factor(Annotation_b2$Drug)
Annotation_b2$Timepoint=factor(Annotation_b2$Timepoint,levels=c("0h","8h","24h"))
TPM_b2_2=TPM_counts_b2%>%.[,c(2:ncol(.))]%>%dplyr::filter(gene_name%in%Protein_coding_genelist)%>%group_by(gene_name)%>%summarise_all(mean)%>%as.data.frame()%>%matrix.please()%>%as.data.frame()
Counts_b2_2=Raw_counts_b2%>%.[,c(2:ncol(.))]%>%dplyr::filter(gene_name%in%Protein_coding_genelist)%>%group_by(gene_name)%>%summarise_all(mean)%>%as.data.frame()%>%matrix.please()%>%as.data.frame()%>%mutate_all(.,as.integer)

