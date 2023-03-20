# check if used pathways match
# MelBrainSys_PhDProject_TKraft/DataSets/Annotations/geneset_with_cancerGenes_merged.txt & Co
# compared to annotation/pathway-definitions.ods

library(xlsx)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

# further pathways compiled by Nathan/Theresa
dat = read.csv(file = "/data/bcu_projects/MelBrainSys_PhDProject_TKraft/DataSets/Annotations/geneset_with_cancerGenes_merged.txt", header = T,sep = "\t",stringsAsFactors = F)
furtherPWs = NULL

dat = read.csv(file = "/data/bcu_projects/MelBrainSys_PhDProject_TKraft/DataSets/Annotations/geneset_with_metabolics_merged.txt", header = T,sep = "\t",stringsAsFactors = F)
head(dat,2)
for(i in 8:ncol(dat)) {
    n = colnames(dat)[i]
    furtherPWs[n] = list(dat$Gene[ which(dat[,i]==1)])
}

dat = read.csv(file = "/data/bcu_projects/MelBrainSys_PhDProject_TKraft/DataSets/Annotations/geneset_with_pathways_merged.txt", header = T,sep = "\t",stringsAsFactors = F)
head(dat,2)
for(i in 8:ncol(dat)) {
    n = colnames(dat)[i]
    furtherPWs[n] = list(dat$Gene[ which(dat[,i]==1)])
}

sapply(furtherPWs,length)
#furtherPWs$TF.Cofactor

# load pathways of interest:
pathways = read.csv(file = "../../../annotation/Rebekka-pathways-2021-06-14/Rebekkas-pathways-gene-annot.tsv",header = T,sep = "\t",stringsAsFactors = F)
head(pathways,3)
pathwayList = sapply(sort(unique(pathways$pathwayName)), function(pw) sort(unique(pathways[ pathways$pathwayName==pw,"gene_symbol"])))
pathwayList

# remove K_biosyn, only 2 genes, made problems before
pathwayList = c(pathwayList,furtherPWs)
goodPathwayList = pathwayList
wh = which(names(goodPathwayList) == "K_biosyn")
goodPathwayList = goodPathwayList[-wh]

metabolicPW = c("glycolysis", "citrate_cycle", "pyruvate_meta", "pentP", "oxidative_phospho", "fatty_acids", 
                "inositolP", "pyrimidine", "purine", "ADE_meta", "W_meta", "RP_meta", "F_meta", 
                "QD_meta", "FYW_biosyn", "CM_meta", "H_meta", "GST_meta", "VIL_biosyn", "Y_meta")
    # = all of Theresa's geneset_with_metabolics_merged.txt, omit K_biosyn 
signallingPW = c("Melanoma", "ppar", "mapk", "erbb", "cytokine", "cell_cycle", "p53", "mtor", 
                 "pi3k_akt", "apoptosis", "wnt", "tgfb", "vegf", "focal_adh", "ECM", "adh_jun", "jak_stat", 
                 "notch", "hedgehog", "replication", "BER", "NER", "HR", "NHEJ", "mismatch", "telomere",
                "ECM-receptor interaction")
    # = all of Theresa's geneset_with_pathways_merged.txt   &  geneset_with_cancerGenes_merged.txt   & Rebekka's Melanome
immunePW = c("PD-L1 expression and PD-1 checkpoint pathway in cancer", "Toll-like receptor signaling pathway", 
             "RIG-I-like receptor signaling pathway","NOD-like receptor signaling pathway", 
             "Natural killer cell mediated cytotoxicity", "T cell receptor signaling pathway",
             "Antigen processing and presentation", "B cell receptor signaling pathway", 
             "Th1 and Th2 cell differentiation", "Th17 cell differentiation", "IL-17 signaling pathway", 
             "Leukocyte transendothelial migration", "Cytokine-cytokine receptor interaction")
    # all of Rebekka's immune except "Melanome"
metabolicPW; signallingPW; immunePW
# now get gene lists for each name:
metabolicPW = goodPathwayList[metabolicPW]
signallingPW = goodPathwayList[signallingPW]
immunePW = goodPathwayList[immunePW]

pwCategories = c(signallingPW=list(signallingPW), metabolicPW=list(metabolicPW), immunePW=list(immunePW))
pwCategories

cancerSignPWs = read.csv(file = "annotation/pathway-definitions-cancer-signaling.csv", header = T,sep = "\t",stringsAsFactors = F )
colnames(cancerSignPWs) = gsub("\\."," ",colnames(cancerSignPWs))
head(cancerSignPWs,3)

immunePWs = read.csv(file = "annotation/pathway-definitions-immune.csv", header = T,sep = "\t",stringsAsFactors = F )
colnames(immunePWs) = gsub("\\."," ",colnames(immunePWs))
head(immunePWs,3)

metabolicPWs = read.csv(file = "annotation/pathway-definitions-metabolic.csv", header = T,sep = "\t",stringsAsFactors = F )
colnames(metabolicPWs) = gsub("\\."," ",colnames(metabolicPWs))
head(metabolicPWs,3)

tmp = read.csv(file = "annotation/Nathan-pathway-name-mapping.csv", header = T,sep = "\t",stringsAsFactors = F)
head(tmp)
pwNameMapping = tmp$Konrad
names(pwNameMapping) = tmp$Nathan
head(pwNameMapping)

names(pwNameMapping)
names(pwCategories$metabolicPW)
all(names(pwCategories$metabolicPW) %in% names(pwNameMapping))
all(names(pwCategories$immunePW) %in% names(pwNameMapping))
all(names(pwCategories$signallingPW) %in% names(pwNameMapping))

categ = names(pwCategories)[1]
categ
for (pwName in names(pwCategories[[categ]])) {
    pwName = gsub("-"," ",pwName)
    pwNamesNice = pwNameMapping[pwName]
    if( ! pwNamesNice %in% colnames(cancerSignPWs)) {
        cat(pwNamesNice," not inside\n")
    } else {
        cat(length(sort(cancerSignPWs[ cancerSignPWs[,pwNamesNice]==1,"Gene"]))," ",
            length(sort(pwCategories[[categ]][[pwName]]))," - ",pwName,"\n")
    }
}

categ = names(pwCategories)[2]
categ
for (pwName in names(pwCategories[[categ]])) {
    pwName = gsub("-"," ",pwName)
    pwNamesNice = pwNameMapping[pwName]
    if( ! pwNamesNice %in% colnames(metabolicPWs)) {
        cat(pwNamesNice," not inside\n")
    } else {
        cat(length(sort(metabolicPWs[ metabolicPWs[,pwNamesNice]==1,"Gene"]))," ",
            length(sort(pwCategories[[categ]][[pwName]]))," - ",pwName,"\n")
    }
}

categ = names(pwCategories)[3]
categ
for (pwName in names(pwCategories[[categ]])) {
    pwNameHyphen = pwName
    pwName = gsub("-"," ",pwName)
    cat(length(sort(immunePWs[ immunePWs[,pwName]==1,"Gene"]))," ",
        length(sort(pwCategories[[categ]][[pwNameHyphen]]))," - ",pwName,"\n")
}
#  ->  0..5 genes missing!

# TODO: 
remove H_metab or renoma to His metab
add Melanoma into Supp Tab 4 sheet cancer signalling
