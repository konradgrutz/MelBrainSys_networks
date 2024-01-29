library(xlsx)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

immunePWs = read.xlsx(file = "annotation/pathway-definitions.xls", header = T, sheetName = "ImmunePathways")
colnames(immunePWs) = gsub("\\."," ",colnames(immunePWs))
head(immunePWs,3)

metabolicPWs = read.xlsx(file = "annotation/pathway-definitions.xls", header = T, sheetName = "MetabolomePathways")
colnames(metabolicPWs) = gsub("\\."," ",colnames(metabolicPWs))
head(metabolicPWs,3)

cancerSignPWs = read.xlsx(file = "annotation/pathway-definitions.xls", header = T, sheetName = "CancerSignalingPathways")
colnames(cancerSignPWs) = gsub("\\."," ",colnames(cancerSignPWs))
head(cancerSignPWs,3)

# metab: Lys metab is added compared to the done analyses
# cancer sign: DNA replication added, ECM missing, Melanoma missing, replication missing
# immune: ECM receptor interaction added

pwCategories = list(signalingPWs = NULL, metabolicPWs = NULL, immunePWs = NULL)

for (pw in colnames(metabolicPWs[,-c(1:7)])) {
    pwCategories$metabolicPWs = 
        append(pwCategories$metabolicPWs,
               list(metabolicPWs$Gene[ metabolicPWs[,pw]==1]))
}
names(pwCategories$metabolicPWs) = colnames(metabolicPWs[,-c(1:7)])

for (pw in colnames(cancerSignPWs[,-c(1:7)])) {
    pwCategories$signalingPWs = 
        append(pwCategories$signalingPWs,
               list(cancerSignPWs$Gene[ cancerSignPWs[,pw]==1]))
}
names(pwCategories$signalingPWs) = colnames(cancerSignPWs[,-c(1:7)])

pw = colnames(immunePWs[,-c(1:7)])[1]
pw
for (pw in colnames(immunePWs[,-c(1:7)])) {
    pwCategories$immunePWs = 
        append(pwCategories$immunePWs,
               list(immunePWs$Gene[ immunePWs[,pw]==1]))
}
names(pwCategories$immunePWs) = colnames(immunePWs[,-c(1:7)])

pwCategories

sort(sapply(unlist(pwCategories, recursive = F),length))

# remove Lys metbaolism, only 2 genes
pwCategories$metabolicPWs$`Lys biosyn` 
pwCategories$metabolicPWs = pwCategories$metabolicPWs[ names(pwCategories$metabolicPWs) != "Lys biosyn"]

names(pwCategories$metabolicPWs)

saveRDS(object = pwCategories, file = "pathway-categories.rds")
