# defines and prepares altered genes, i.e. differential promoter methylation + expression tendency
# saved in altered-genes-per-patient.rds
# overlap of these between metastases pairs

library(pheatmap)
library(ggplot2)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

# define met. pairs, and annotation colors 
tmp = readRDS(file = "annotation/samplePairs-annotation-colors-clusters.rds")
patient_colors = tmp$patient_colors
tissue_colors = tmp$tissue_colors
neededSamples = tmp$neededSamples
sampleMapping = tmp$sampleMapping
samplePairs = tmp$samplePairs
samplePairPerSubgroup = tmp$samplePairPerSubgroup
subgroupPerSamplePair = tmp$subgroupPerSamplePair
subgroupNames = tmp$subgroupNames
colAnnot = tmp$colAnnot
annotColors = tmp$annotColors

expr = read.csv(file = "data/MelBrainSys_RNA_HMM_output.txt",
                header = T,sep = "\t",stringsAsFactors = F)
colnames(expr)[colnames(expr)=="Genes"] = "EnsemblID"
head(expr,3)
length(which(duplicated(expr$Gene))) # no duplicated genes

methyl = read.csv(file = "data/MelBrainSys_diffMethylatedGenes_expressionTrend.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(methyl) = gsub("\\.","_",colnames(methyl))
colnames(methyl)[colnames(methyl)=="P16_Blun"] = "P16_BLun"
colnames(methyl)[colnames(methyl)=="P03_Blun"] = "P03_BLun"
head(methyl,3)

# remove duplicated genes:
length(which(duplicated(methyl$Gene)))
g = methyl$Gene[which(duplicated(methyl$Gene))]
wh= which(methyl$Gene %in% g)
g
length(wh)
methyl = methyl[-wh,]
length(which(duplicated(methyl$Gene)))

# gene symbol / EnsgID mapping
dat = read.csv(file = "annotation/EnsgEntrez-hg19.txt.gz",header = T,sep = "\t",stringsAsFactors = F)
head(dat,3)
wh = which(duplicated(dat$Gene.stable.ID))
length(wh)
dat = dat[-wh,]
wh = which(duplicated(dat$Gene.name))
length(wh)
dat = dat[-wh,]
ensg2gene = dat$Gene.name
names(ensg2gene) = dat$Gene.stable.ID
head(ensg2gene)

# only common genes between expr and methyl table:
# first, find gene symbol for EnsgIDs of expression table:
wh = which(expr$EnsemblID %in% names(ensg2gene))
length(wh)
round(length(which(expr$EnsemblID %in% names(ensg2gene)))/nrow(expr)*100) # 98% of genes known
expr = expr[wh,]
expr$Gene = ensg2gene[ expr$EnsemblID ]
head(expr,2)

setdiff(colnames(expr), colnames(methyl))
setdiff(colnames(methyl), colnames(expr))
# 1 sample is not in methyl data: P04_BSki_2
commonPatSamplePair = sort(setdiff(intersect(colnames(expr),colnames(methyl)),c("Gene","EnsemblID")))
commonPatSamplePair

# for each patient, a list of diff expr + methyl genes
pat = commonPatSamplePair[1]
regGenesPerPat = NULL  # regulated genes (= altered genes)
upRegGenesPerPat = NULL # upregulated genes
downRegGenesPerPat = NULL # downregulated
for(pat in commonPatSamplePair) {
    upExprGenes = expr$Gene[which(expr[,pat]=="+")]
    downExprGenes = expr$Gene[which(expr[,pat]=="-")]
    upMethylGenes = methyl$Gene[which(methyl[,pat]=="1")]
    downMethylGenes = methyl$Gene[which(methyl[,pat]=="-1")]
    upRegGenes = intersect(upExprGenes, downMethylGenes)
    downRegGenes = intersect(downExprGenes, upMethylGenes)
    regGenesPerPat[[length(regGenesPerPat)+1]] = c(upRegGenes, downRegGenes)
    upRegGenesPerPat[[length(upRegGenesPerPat)+1]] = upRegGenes
    downRegGenesPerPat[[length(downRegGenesPerPat)+1]] = downRegGenes
}
names(regGenesPerPat) = commonPatSamplePair
names(upRegGenesPerPat) = commonPatSamplePair
names(downRegGenesPerPat) = commonPatSamplePair
regGenesPerPat

saveRDS(object = list(regGenesPerPat = regGenesPerPat,
                      upRegGenesPerPat = upRegGenesPerPat,
                      downRegGenesPerPat = downRegGenesPerPat),
                      file = "altered-genes-per-patient.rds")

# matrix with pair-wise overlapps:
ovlpMat = matrix(data = 0,nrow = length(commonPatSamplePair), ncol = length(commonPatSamplePair),
                 dimnames = c(list(commonPatSamplePair),list(commonPatSamplePair)))
for(n1 in commonPatSamplePair) {
    for(n2 in commonPatSamplePair) {
        ovlpMat[n1,n2] = length(intersect(regGenesPerPat[[n1]],regGenesPerPat[[n2]]))
        if(n1==n2){ovlpMat[n1,n2]=NA}
    }
}

ovlpMat

# diff. expr/methyl. genes are all in our data?
ourExpr = read.csv(file = "regNet/Data/MelBrainSys-expression.csv",header = T,sep = "\t",stringsAsFactors = F)
head(ourExpr,2)

allGenes = intersect(sort(unique(unname(unlist(regGenesPerPat)))),ourExpr$geneSymbol)
length(allGenes); head(allGenes)

# matrix with pair-wise overlapps:
genesPatients = matrix(data = 0,nrow = length(allGenes), ncol = length(samplePairs),
                 dimnames = c(list(allGenes),list(samplePairs)))

for(n1 in samplePairs) {
    genesPatients[ intersect(allGenes,regGenesPerPat[[n1]]),n1] = 1
}

hmGenesPat = pheatmap(mat = genesPatients,annotation_col=colAnnot, legend = F,color = c("white","deepskyblue3"),
                      show_rownames = F, angle_col = 45,
                   annotation_colors = annotColors)

svg(filename = "FiguresTables/SupplFigure-2-genes-vs-samplePairs-part1.svg")
hmGenesPat
dev.off()

rs = rowSums(genesPatients)
wh = which(rs > 3)
length(wh)
hmGenesPatSub = pheatmap(mat = genesPatients[wh,],annotation_col=colAnnot, legend = F,
                         color = c("white","deepskyblue3"),show_rownames = T, angle_col = 45,
                   annotation_colors = annotColors)

svg(filename = "FiguresTables/SupplFigure-2-genes-vs-samplePairs-part2.svg")
hmGenesPatSub
dev.off()

rs = as.numeric(table(rowSums(genesPatients)))
data = data.frame(num_met_pairs = rowSums(genesPatients), gene =rownames(genesPatients))

barplotGenesShared = ggplot(data, aes(x=num_met_pairs)) + geom_histogram(binwidth=1, fill = c("deepskyblue3"), col="black") + 
    xlab("number of metastases pairs") + ylab("number of genes")  + 
    annotate("text", x=1:5,y=rs, label=rs, vjust=-1)
barplotGenesShared
svg(filename = "FiguresTables/SupplFigure-2-genes-vs-samplePairs-part3.svg")
barplotGenesShared
dev.off()

histoGenesShared = ggplot(data[ data$num_met_pairs>=3,], aes(x = reorder(gene,-num_met_pairs), y=num_met_pairs)) + 
    geom_bar(stat="identity", fill = c("deepskyblue3"), col="darkgrey", width = 0.7) + 
    xlab("gene") + ylab("number of metastases pairs") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
histoGenesShared

svg(filename = "FiguresTables/SupplFigure-2-genes-vs-samplePairs-part4.svg")
histoGenesShared
dev.off()
