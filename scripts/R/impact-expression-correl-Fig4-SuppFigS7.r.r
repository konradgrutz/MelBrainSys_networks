# calculates correlations between average log expression ratio and 
#  average log impact ratio of the altered genes on all pathways of each metastases pair
# creates SupplFigure 5 with scatter plots (SupplFigure-S5-impact-expression-correl.png)
# and Figure 4 with barplots of the correlations (Figure-4-plot-impact-expression-correl-barplot.png)

library(ggplot2)
library(grid) # for textGrob for common axis titles in grid.arrange
library(gridExtra)
library(pheatmap)
library(scales) # for wrap_format, to wrap labels
library(stringr) # edit names in pheatmap
library(pvclust) # for clustering stability analyses
library(xlsx)

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

outDirectory = paste0(basePath,"FiguresTables/")

# pathway definitions
pwCategories = readRDS(file = "pathway-categories.rds")
head(pwCategories$signalingPWs$PPAR)

genesPerPathway = unlist(unname(pwCategories),recursive = F)
lapply(head(genesPerPathway,3),head)

# pre-calculated impact log-ratios on pathways from each altered gene
metPairImpacts = readRDS(file = "metPairs-impactRatios-onPathways.rds")
head(metPairImpacts,3)
wh = which(!is.na(metPairImpacts$log2_mean_impact_ratio))
length(wh); nrow(metPairImpacts)
metPairImpactsNoNA = metPairImpacts[wh,]

# load log expression values
expr = read.csv(file = "regNet/Data/MelBrainSys-expression.csv", header = T,sep = "\t",stringsAsFactors = F)
head(expr,3)

impactExprTable = NULL
for (samplePair in names(sampleMapping)) {
    sampleBr = sampleMapping[[samplePair]][[1]]
    sampleExC = sampleMapping[[samplePair]][[2]]
    for (pw in names(genesPerPathway)) {
        wh = which(expr$geneSymbol %in% genesPerPathway[[pw]])
        exprRatio = mean(expr[ wh,sampleBr] - expr[wh,sampleExC])
        impRatio = mean(metPairImpactsNoNA$log2_median_impact_ratio[ 
            metPairImpactsNoNA$pathway == pw & 
            metPairImpactsNoNA$samplePair == samplePair])
        impactExprTable = 
            rbind(impactExprTable,
                  data.frame(stringsAsFactors = F,
                      sample_pair = samplePair, pathway = pw, 
                      "impact_ratio" = impRatio,
                      "expression_ratio" = exprRatio,
                      subgroup = subgroupNames[subgroupPerSamplePair[samplePair]],
                      "subgroup ID" = subgroupPerSamplePair[samplePair]))
    }
}
impactExprTable = impactExprTable[ order(impactExprTable$sample_pair,
                                         impactExprTable$pathway),]

head(impactExprTable,3)

samplePairPerSubgroup

cors = NULL
for(sp in names(sampleMapping)) {
    tmp = cor.test(x = impactExprTable[ impactExprTable$sample_pair==sp,"impact_ratio"], 
            y = impactExprTable[ impactExprTable$sample_pair==sp,"expression_ratio"], method = "pearson")
    cors = rbind(cors, 
                 data.frame(sample_pair = sp, r = round(tmp$estimate,2), p = tmp$p.value))
}
cors$FDR = p.adjust(cors$p)
cors

corLabels = data.frame(stringsAsFactors = F,sample_pair =  cors$sample_pair,
           label=paste0("r = ",round(cors$r,2),
                        ", p = ",format(cors$p,digits=2),
                        ", FDR = ",format(cors$p,digits=2)
                       ))
corLabels

samplePairPerSubgroup$SG1

# split by subgroup (higher, balanced, lower impact in brain):
p1 = ggplot(data=impactExprTable[impactExprTable$sample_pair %in% samplePairPerSubgroup$SG1,], 
           mapping = aes(x=impact_ratio,y=expression_ratio)) + geom_point() + 
             facet_wrap(~sample_pair,nrow = 1) + 
             scale_y_continuous(limits=c(-1.5,2.5),labels = wrap_format(10), name = "") + 
             scale_x_continuous(name= "") + ggtitle("SG1, higher impact in brain") + 
            geom_smooth(method=lm, fullrange=F) + # regression line + confidence
             geom_text(data = corLabels[ corLabels$sample_pair %in% samplePairPerSubgroup$SG1,],  
                       mapping = aes(x = -Inf, y = -Inf, label = label), size=4, 
                       hjust   = -0.05,  vjust   = -1, col="#222222") +
            geom_hline(yintercept = 0, col=alpha("red",0.3)) + geom_vline(xintercept = 0, col=alpha("red",0.3)) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1

p2 = ggplot(data=impactExprTable[impactExprTable$sample_pair %in% samplePairPerSubgroup$SG2,], 
           mapping = aes(x=impact_ratio,y=expression_ratio)) + geom_point() + 
             facet_wrap(~sample_pair,nrow = 1) + 
             scale_y_continuous(limits=c(-1.5,2.5),labels = wrap_format(10), name = "") + 
             scale_x_continuous(name= "") + ggtitle("SG2, lower impact in brain") + 
            geom_smooth(method=lm, fullrange=F) + # regression line + confidence
            geom_text(data = corLabels[ corLabels$sample_pair %in% samplePairPerSubgroup$SG2,],  
                  mapping = aes(x = -Inf, y = -Inf, label = label), size=4, 
                   hjust   = -0.05,  vjust   = -1, col="#222222") +
            geom_hline(yintercept = 0, col=alpha("red",0.3)) + geom_vline(xintercept = 0, col=alpha("red",0.3)) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2

p3 = ggplot(data=impactExprTable[impactExprTable$sample_pair %in% samplePairPerSubgroup$SG3,], 
           mapping = aes(x=impact_ratio,y=expression_ratio)) + geom_point() + 
             facet_wrap(~sample_pair,nrow = 1) + 
             scale_y_continuous(limits=c(-1.5,2.5),labels = wrap_format(10), name = "") + 
             scale_x_continuous(name= "") + ggtitle("SG3,slightly lower impact in brain") + 
            geom_smooth(method=lm, fullrange=F) + # regression line + confidence
            geom_text(data = corLabels[ corLabels$sample_pair %in% samplePairPerSubgroup$SG3,],    
                       mapping = aes(x = -Inf, y = -Inf, label = label), size=4, 
                       hjust   = -0.05,  vjust   = -1, col="#222222") +
            geom_hline(yintercept = 0, col=alpha("red",0.3)) + geom_vline(xintercept = 0, col=alpha("red",0.3)) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3

png(filename = paste0(outDirectory,"SupplFigure-S7-impact-expression-correl.png"), 
    width = 3000, height = 3000, res=220)
grid.arrange(grobs=list(p1,p3,p2), nrow=3, # widths=c(1,1,1.1), 
            left = textGrob(expression(expression~log[2]~ratio), rot = 90, vjust = 1),
                         bottom = textGrob(expression(impact~log[2]~ratio), vjust = 0.4))
dev.off()

#impExprVoomCorrelSamplePair
correls = data.frame(stringsAsFactors = F,
                     correlation = cors$r, p = cors$p, p_adjusted = cors$FDR, 
                sample_pair = cors$sample_pair)
correls$subgroupID = subgroupPerSamplePair[ correls$sample_pair]
correls$subgroup = as.factor(subgroupNames[ correls$subgroupID])
correls$subgroupNum = as.numeric(gsub("SG","",correls$subgroupID))
head(correls,3)

t(correls[,c(1,4)])

# signif star labels:
o = order(correls$subgroupNum)
whSign = which(correls$p_adjusted[o] < 0.05)
whY = correls$correlation[o][whSign]-0.025

p = ggplot(data = correls,mapping = aes(x=reorder(sample_pair,subgroupNum),
                                        y=correlation, fill=subgroupID)) + 
    geom_bar(stat="identity") + 
    xlab("") + ylab("Correlation between impact and expression") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    scale_x_discrete(position = "top")  + 
    scale_fill_manual(values=c(annotColors$subgroup["SG1"], annotColors$subgroup["SG2"], 
                               annotColors$subgroup["SG3"]), 
                     breaks=c("SG1","SG2","SG3")) + 
    #theme(legend.title = element_blank()) + 
    annotate('text', x = whSign, y =whY, label = '*', size = 4)
p

png(filename = paste0(outDirectory,"Figure-4-plot-impact-expression-correl-barplot.png"),
    width = 1400, height = 1600,res = 240)
p
dev.off()
