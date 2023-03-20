# calculate auto-correlation of the expression of neighboring genes 
#  based on TCGA expression values, Figure S1

library(reshape2) # for melt()
library(ggplot2)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

outDirectory = paste0(basePath,"FiguresTables/")

geneExprTable = read.csv(file = "regNet/Data/TCGA-expression.csv", 
                         header = T, sep="\t",stringsAsFactors = F)
head(geneExprTable,3); dim(geneExprTable)

# make sure genes are sorted by chr/position
class(geneExprTable$chr); class(geneExprTable$pos)
o = order(geneExprTable$chr,geneExprTable$pos)
head(o)
table(o[2:length(o)]-o[1:(length(o)-1)]) # already in desired order

# double check that there are no NAs inside:
table(apply(geneExprTable[,-c(1:3)],2,function(c) length(which(is.na(c)))))
# no NA in any column

# speed up by having a double matrix of gene expr values:
geneExprTableM = as.matrix(geneExprTable[,-c(1:3)])
class(geneExprTableM); typeof(geneExprTableM); dim(geneExprTableM)

# for each chr, get indexes of respective genes in expression array:
chrIdx = sapply(unique(geneExprTable$chr), function(chr) which(geneExprTable$chr == chr))

length(chrIdx)
sapply(names(chrIdx), function(chr) head(chrIdx[[chr]]))
# OK

# calculate auto-correlation using acf() of stat package
# iterate over chromosomes separately
# iterate over all 270 TCGA samples
# allAcfs is a list of chromosomes, each containing a matrix of dimension maximum lag x samples, 
#  matrix entries are the correlations of gene expressions of genes with "lag" distance in the sample on that chromosome
allAcfs = lapply (names(chrIdx), function(chr) {
    sapply(1:ncol(geneExprTableM), function(sampleIdx) {
        acfRes = acf(geneExprTableM[chrIdx[[chr]],sampleIdx], lag.max = 200,plot = F)
        acfRes$acf[,1,1]
    })
})
names(allAcfs) = names(chrIdx)

dim(allAcfs$chr21)
allAcfs$chr21[1:4,1:4]

neededChr = setdiff(names(chrIdx),"chrY")
o=order(as.numeric(gsub("chr(\\d+.*)","\\1",neededChr)))
neededChr = neededChr[o]
neededChr

neededQuantiles = c(0 ,0.1 ,0.25 ,0.5 ,0.75, 0.9, 1 )

# quantiles over all samples, separately for all combinations of all chr and lags
allAcfQuantiles = NULL
for(chr in neededChr) {
    q = t(sapply(1:nrow(allAcfs[[chr]]),function(lag) quantile(allAcfs[[chr]][lag,],neededQuantiles)))
    colnames(q) = gsub("\\%","",paste0("quant_",colnames(q)))
    allAcfQuantiles = rbind(allAcfQuantiles, data.frame(chr=chr,lag = 1:nrow(allAcfs[[chr]]), q))
}

# condense the data of all chromosomes: 
allAcfQuantilesChrCondensed = NULL
for(lag in sort(unique(allAcfQuantiles$lag))) {
    allAcfQuantilesChrCondensed = 
        rbind(allAcfQuantilesChrCondensed,
              data.frame(lag=lag,t(colMeans(allAcfQuantiles[ allAcfQuantiles$lag==lag,-c(1,2)]))))
}

head(allAcfQuantilesChrCondensed,3)

fillCols = c("0%-100%" = "#f3e43f", "10%-90%" = "#e99d63",  "25%-75%" = "#e96363")
lineType = c("10%-90%" = "dashed", "0%-100%" = "dotted")
p = ggplot(data = allAcfQuantilesChrCondensed, mapping=aes(x=lag)) +
    geom_ribbon(aes(ymin = quant_0, ymax = quant_100, fill = "0%-100%"), alpha = .5) + 
    geom_ribbon(aes(ymin = quant_10, ymax = quant_90, fill = "10%-90%"), alpha = .5) + 
    geom_ribbon(aes(ymin = quant_25, ymax = quant_75, fill = "25%-75%"), alpha = .5) + 
    geom_line(mapping = aes(y = quant_50), col="red") + 
    geom_vline(xintercept = seq(0,200,50), col="grey" ) + 
    geom_vline(xintercept = seq(0,200,10), col="grey",linetype="dotted" ) + 
        xlim(c(2,200)) + ylim(c(-0.25,0.35))  + 
     labs(y = "correlation") + 
    scale_fill_manual(values = fillCols, name = "correlation quantiles") #+ 
p

png(filename = paste0(outDirectory, "SupplFigure-S1-expr-auto-correlation.png"),
    width = 2000,height=1200, res = 340)
p
dev.off()
