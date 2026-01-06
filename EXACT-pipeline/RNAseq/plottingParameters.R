##setting plotting options to ensure consistency across figures


###get colors for drug names, drug classes, and MOA
library(RColorBrewer)

##this script sets 3 color palletes: one for drug names, one for classes,
#and one for mechamisms of action as defined by the CCLE

###drug MOA
drug.info <- list(DMSO = "DMSO", "Palbociclib" = "CDK inhibitor", "Ribociclib" = "CDK inhibitor",
                  "Trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "Decitabine" = "DNMT inhibitor", "Vorinostat" = "HDAC inhibitor",
                  "Mirdametinib" = "MEK inhibitor", "Selumetinib" = "MEK inhibitor",
                  "Trametinib" = "MEK inhibitor", "Capmatinib" = "MET inhibitor",
                  "Olaparib" = "PARP inhibitor", "RMC-4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "TOP inhibitor")#,
#                  "Irinotecan" = "TOP inhibitor", "Verteporfin" = "YAP inhibitor")


moaCol <- colorRampPalette(RColorBrewer::brewer.pal(length(unique(drug.class)),"Paired"))(length(unique(drug.class))) # need color ramp for >12 colors
names(moaCol) <- unique(drug.class)
moaCol <-  moaCol[unlist(drug.class)] ##spread it out
names(moaCol) <- names(drug.class)


drug.class <- list(Trametinib = "MEK1/2i", Mirdametinib="MEK1/2i",
                   Selumetinib = "MEK1/2i", TNO155="SHP2i", `RMC-4630`="SHP2i",
                   Palbociclib= "CDK4/6i",Ribociclib= "CDK4/6i",Decitabine= "DNMTi",
                   Olaparib= "PARPi",Vorinostat= "HDACi",Capmatinib= "c-Meti",
                   #Verteporfin	= "TEADi",IAG933= "TEADi",Pexidartinib= "CSFRi",
                   Doxorubicin= "chemotherapy",Ifosfamide= "chemotherapy",`SN-38`	= "chemotherapy",
                   Trabectedin	= "chemotherapy")


classCol <- colorRampPalette(RColorBrewer::brewer.pal(length(unique(drug.class)),"Paired"))(length(unique(drug.class))) # need color ramp for >12 colors
names(classCol) <- unique(drug.class)
classCol <-  classCol[unlist(drug.class)] ##spread it out
names(classCol) <- names(drug.class)

drugnames <- intersect(names(drug.info),names(drug.class))

drugCol <- colorRampPalette(RColorBrewer::brewer.pal(length(unique(drugnames)),"Set3"))(length(unique(drugnames))) # need color ramp for >12 colors
names(drugCol) <- unique(drugnames)


##very cofusing: change drug.info to data frame!
drug.df <- data.frame(Drug=drugnames, MOA = unlist(drug.info[drugnames]),
                      Class=unlist(drug.class[drugnames]))

