# produces
# Figure 5: circos plots of source and target genes
# SupplTable S10 target genes mean impact and expr. ratios per group
# SupplFigure S8 target gene candidates scatterplot ompact vs. expression

library(reshape2)
library(ggplot2)
library(circlize)
library(scales) # for alpha semi-transparency
library(parallel)
library(pheatmap)
#library(tidytext) # reorder values within each facet of ggplot
library(gridExtra)
library(xlsx)
library(readxl)
library(ggrepel)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

outDirectory = paste0(basePath,"FiguresTables/")

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

chrSizes = read.csv("annotation/hg19.chrom.sizes", header = F,sep = "\t", stringsAsFactors = F, col.names = c("chr","size"))
chrSizes = chrSizes[ chrSizes$chr %in% c(paste0("chr",1:22),"chrX"),]
o = order(as.numeric(gsub("chr(.+)","\\1",chrSizes$chr)))
chrSizes = chrSizes[o,]

load("annotation/Homo_sapiens.GRCh37.75-chrRename-noHaplo.RData")

head(geneAnnot,2)
dim(geneAnnot)
geneAnnot = geneAnnot[ geneAnnot$V3=="gene",]
geneAnnot$gene = gsub(".*gene_name ([^;]+).*","\\1",geneAnnot$V9)
wh = which(duplicated(geneAnnot$gene))
length(wh)
dupGeneNames = sort(unique(geneAnnot$gene[wh]))
geneAnnot = geneAnnot[ -which(geneAnnot$gene %in% dupGeneNames),]
dim(geneAnnot)
rownames(geneAnnot) = geneAnnot$gene

geneAnnot$pos = round((geneAnnot$V5+geneAnnot$V4)/2) # gene pos = middle of start and end
geneAnnot$pos = as.integer(geneAnnot$pos)
colnames(geneAnnot)[1]="chr"
geneAnnot$chr = as.character(geneAnnot$chr)

head(geneAnnot,3)

tmp = readRDS(file = "altered-genes-per-patient.rds")
regGenesPerPat = tmp$regGenesPerPat
upRegGenesPerPat = tmp$upRegGenesPerPat
downRegGenesPerPat= tmp$downRegGenesPerPat

sort(regGenesPerPat$P04_BSki_1)

sapply(regGenesPerPat,head)

altGenes = NULL
for(sampleName in names(regGenesPerPat)) {
    upGenes = upRegGenesPerPat[[sampleName]]; downGenes = downRegGenesPerPat[[sampleName]]
    altGenes = rbind(altGenes,
                 data.frame(samplePair = sampleName, geneAnnot[upGenes,],direction=1,y=1),
                 data.frame(samplePair = sampleName, geneAnnot[downGenes,],direction=2,y=1))
}
altGenes$pos = as.integer(altGenes$pos)
altGenes$y = as.integer(altGenes$y)    
head(altGenes,2)

dat = read_excel(path = "annotation/pathway-definitions.xls", sheet = 5)

head(dat,2)

TFs = dat$Gene[ dat[,"TF/Cofactor"]==1]
oncogene = dat$Gene[ dat$Oncogene==1]
tumorsuppr = dat$Gene[ dat$TumorSuppressorGene==1]
cancerCensus = dat$Gene[ dat$CancerCensusGene==1]

# pathway definitions
pwCategories = readRDS(file = "pathway-categories.rds")
genesPerPathway = unlist(unname(pwCategories),recursive = F)

signalingGenes = unique(sort(unname(unlist(pwCategories$signalingPWs))))
immuneGenes = unique(sort(unname(unlist(pwCategories$immunePWs))))
metabolicGenes = unique(sort(unname(unlist(pwCategories$metabolicPWs))))
geneList = list("immuneGenes"=immuneGenes, "metabolicGenes"=metabolicGenes, 
                "signalingGenes"=signalingGenes, "cancerCensus"=cancerCensus,
                "oncogene"=oncogene,"tumorsuppr"=tumorsuppr,"TFs"= TFs)

cancerGenes = sort(unique(c(cancerCensus, tumorsuppr, oncogene)))
metabRestGenes = setdiff(metabolicGenes,cancerGenes)
signalRestGenes = setdiff(signalingGenes, cancerGenes)
TFrestGenes = setdiff(TFs, cancerGenes)

getGeneCateg = function(genes) {
    categs = rep("gene",length(genes))
    wh = which(genes %in% cancerGenes); categs[wh] = "cancer"
    wh = which(genes %in% metabRestGenes); categs[wh] = "metabolism"
    wh = which(genes %in% signalRestGenes); categs[wh] = "signaling"
    wh = which(genes %in% TFrestGenes); categs[wh] = "TF"
    categs
}

categColors = c("TF" = "#4444FF","metabolism"="orange",
            "cancer"="red","signaling"="green","gene"="grey")
geneColors = c("up"="red","down"="#008800","target"="#6666FF", "mixed"="#884400")

samplePairGroupNiceNames = c("SG1"="higher","SG3" = "slightly lower", "SG2"="lower")

tmp = readRDS(file="metPairs-impactRatios-onAllTargetGenes.rds")
logMedianImpRatios = tmp$logMedianImpRatios
logMedianImpRatios[[1]][1:3,1:4]

samplePairPerSubgroup[["SG1"]]

perGroupMeanLogImp

# individual and mean group log impact ratios of the top genes:
perGroupLogImpMat = NULL
for( sg in names(samplePairPerSubgroup)) {
    logImpMat = NULL
    for(sp in samplePairPerSubgroup[[sg]]) {
        logImpMat = rbind(logImpMat, colMeans(logMedianImpRatios[[sp]], na.rm = T))
    }
    #  remove genes which have any NA entry:
    wh = which(apply(logImpMat,2,function(c) all(!is.na(c)) ))
    logImpMat = logImpMat[,wh]
    perGroupLogImpMat[ length(perGroupLogImpMat)+1] = list(logImpMat)
}
names(perGroupLogImpMat) = names(samplePairPerSubgroup)

perGroupLogImpMat[[1]][1:3,1:4]

perGroupMeanLogImp = NULL
for (sg in names(samplePairPerSubgroup)) {
    # for each sample pair, mean impact from DEGs on each target:
    samplePairGroupMeanLogImp = sapply(samplePairPerSubgroup[[sg]],
                           function(sp) colMeans(logMedianImpRatios[[sp]],na.rm = T))
    # only genes present in all sample pairs:
    wh = which(apply(samplePairGroupMeanLogImp,1,function(r) length(which(is.na(r))))==0)
    samplePairGroupMeanLogImp = samplePairGroupMeanLogImp[wh,]
    # now mean over all sample pairs for each target gene:
    meanGroupMeanLogImp = apply(samplePairGroupMeanLogImp,1,mean,na.omit=T)
    perGroupMeanLogImp[length(perGroupMeanLogImp)+1] = list(meanGroupMeanLogImp)
}
names(perGroupMeanLogImp) = names(samplePairPerSubgroup)

topInteractions = NULL
for(sampleName in names(regGenesPerPat)) {
    logRatioLong = melt(logMedianImpRatios[[sampleName]], 
                        varnames = c("source","target"), value.name = "log_impact_ratio")
    logRatioLong$source = as.character(logRatioLong$source)
    logRatioLong$target = as.character(logRatioLong$target)
    # remove NaN, NA values:
    wh = which(is.nan(logRatioLong$log_impact_ratio) | is.na(logRatioLong$log_impact_ratio))
    if(length(wh)>0) {
        logRatioLong = logRatioLong[-wh,]
    }
    topUp = logRatioLong[ order(logRatioLong$log_impact_ratio, decreasing = T),]
    topUp = topUp[ which(topUp$log_impact_ratio>0),]
    topDwn = logRatioLong[ order(logRatioLong$log_impact_ratio, decreasing = T),]
    topDwn = topDwn[ which(topDwn$log_impact_ratio<0),]
    topInteractions = 
        rbind(topInteractions,
              data.frame(sample_pair = sampleName, rbind(topUp,topDwn)))
}

head(topInteractions)

topInteractCoord = cbind(sample_pair = topInteractions$sample_pair,
     source = topInteractions$source, geneAnnot[ topInteractions$source,c("chr","pos")],
      target = topInteractions$target, geneAnnot[ topInteractions$target,c("chr","pos")],
      log_impact_ratio = topInteractions$log_impact_ratio)
colnames(topInteractCoord) = c("sample_pair","source_gene","source_chr","source_pos",
                               "target_gene","target_chr","target_pos","log_impact_ratio")
topInteractCoord$direction = ifelse(topInteractCoord$log_impact_ratio>0,1,2)

head(topInteractCoord)

topConsistent = NULL
for(gr in names(samplePairPerSubgroup)) {
    m = perGroupLogImpMat[[gr]] # log impact ratios of this group
    m = m[, which(apply(m,2,function(c) all(c>0))) ]# only cols with positive log imp
    if(ncol(m)<1) { next } # n.b. lower and slightly lower have no results
    meanLogImp = apply(m,2,mean)
    m = meanLogImp[ order(meanLogImp,decreasing = T)][1:10]
    m = m[which(m>0)]
    cat(gr,"  ", length(m),"  genes\n")
    if(length(m)==0){next}
    tmp = data.frame(group = gr, gene = names(m), 
                     mean_impact = unname(m), stringsAsFactors = F)
    topConsistent = rbind(topConsistent,tmp)
}
topConsistent = topConsistent[order(topConsistent$group,topConsistent$mean_impact, decreasing = T),]
topConsistent
# no genes for lower and slightly lower group

bottomConsistent = NULL
for(gr in names(samplePairPerSubgroup)) {
    m = perGroupLogImpMat[[gr]] # log impact ratios of this group
    m = m[, which(apply(m,2,function(c) all(c<0))) ] # only cols with negative log imp
    if(ncol(m)<1) { next } # 
    meanLogImp = apply(m,2,mean)
    m = meanLogImp[ order(meanLogImp,decreasing = F)][1:10]
    m = m[which(m<0)]
    cat(gr,"  ", length(m),"  genes\n")
    if(length(m)==0){next}
    tmp = data.frame(group = gr, gene = names(m), 
                     mean_impact = unname(m), stringsAsFactors = F)
    bottomConsistent = rbind(bottomConsistent,tmp)
}
bottomConsistent = bottomConsistent[order(bottomConsistent$group,bottomConsistent$mean_impact),]
bottomConsistent
# 

plotGroupCircosConsistent = function(groupName = "higher") {
    
    # top and bottom 10 target genes with highest/lowest mean log impact in that group:
    topGenes = topConsistent$gene[ topConsistent$group == groupName]
    bottomGenes = bottomConsistent$gene[ bottomConsistent$group == groupName]
    length(topGenes)
    length(bottomGenes)
    
    # all interactions in that group to top and bottom target genes:
    # n.b. topInteractCoord contains all interactions for all sample pairs, sorted by log impact ratio
    groupInteractions = topInteractCoord[ topInteractCoord$sample_pair %in% samplePairPerSubgroup[[groupName]],]
    # reduce to top and bottom target genes:
    groupInteractions = groupInteractions[ groupInteractions$target_gene %in% c(topGenes, bottomGenes), ]
    head(groupInteractions)
    
    # for each source gene, add num of sample pairs in which it is DEG:
    altGenesSub = altGenes[ altGenes$samplePair %in% samplePairPerSubgroup[[groupName]] ,]
    tmp = sort(table(altGenesSub$gene), decreasing = T)
    frqPerGene = as.numeric(tmp); names(frqPerGene) = names(tmp)
    groupInteractions$frqDEG = frqPerGene[groupInteractions$source_gene]
    
    # mean log imp ratio between all source and target genes in this group:
    meanGroupInteractions = aggregate(groupInteractions$log_impact_ratio, 
                                      by = list(source_gene = groupInteractions$source_gene, 
                                                target_gene=groupInteractions$target_gene), FUN = mean)
    colnames(meanGroupInteractions)[3]="mean_log_impact_ratio"
    meanGroupInteractions$frqDEG = frqPerGene[meanGroupInteractions$source_gene]
    meanGroupInteractions = cbind(meanGroupInteractions,geneAnnot[ meanGroupInteractions$source_gene,c("chr","pos")])
    colnames(meanGroupInteractions)[c(5,6)]=c("source_chr","source_pos")
    meanGroupInteractions = cbind(meanGroupInteractions,geneAnnot[ meanGroupInteractions$target_gene,c("chr","pos")])
    colnames(meanGroupInteractions)[c(7,8)]=c("target_chr","target_pos")
    
    # sort by frqDEG (numer of sample pairs in which it is DEG) and abs(mean_log_impact_ratio):
    meanGroupInteractions = meanGroupInteractions[ order(-meanGroupInteractions$frqDEG, -abs(meanGroupInteractions$mean_log_impact_ratio)),]
    meanGroupInteractions$direction = ifelse(meanGroupInteractions$mean_log_impact_ratio>0,1,2)
    # remove chrY interactions:
    wh = which(meanGroupInteractions$source_chr=="chrY" | meanGroupInteractions$target_chr=="chrY")
    if(length(wh)>0) {
        meanGroupInteractions = meanGroupInteractions[-wh,]
    }
    # take the abs top 100 interactions:
    topMeanGroupInteractions = head(meanGroupInteractions,100)
    
    altGenesSub = altGenes[ altGenes$samplePair %in% samplePairPerSubgroup[[groupName]] & altGenes$gene %in% topMeanGroupInteractions$source_gene,]
    altGenesSub
    # aggregate expr. direction of group's altGenes:
    #     direction: 1 higher, 2 lower, 1,2  = mixed between sample pairs
    altGenesSubAggr = aggregate(altGenesSub$direction, by = list(gene = altGenesSub$gene), 
        FUN = function(x) paste0(sort(unique(x)),collapse=","))
    altGenesSubAggr = cbind(altGenesSubAggr, geneAnnot[ altGenesSubAggr$gene,c("chr","pos")], y=1)
    colnames(altGenesSubAggr)[2] = "direction"
                            
    # combine source and target genes for plotting:
    targetGannot = geneAnnot[geneAnnot$gene %in% topMeanGroupInteractions$target_gene,]
    targetSourceGene = rbind(altGenesSubAggr[,c("chr","pos","gene","direction")],
                            cbind(targetGannot[,c("chr","pos","gene")], direction="0"))
    targetSourceGene$direction = ifelse(targetSourceGene$direction=="0","target",
                                        ifelse(targetSourceGene$direction==1,"up",
                                               ifelse(targetSourceGene$direction==2,"down","mixed")))
    targetSourceGene$categ = "gene"; targetSourceGene$categ = getGeneCateg(genes = targetSourceGene$gene)

    # for source genes, add in how many SPs it is a DEG
    targetSourceGene$frqDEG[targetSourceGene$direction != "target"] = 
        sapply(targetSourceGene$gene[ targetSourceGene$direction!="target"], 
            function(gene) unique(meanGroupInteractions$frqDEG[ meanGroupInteractions$source_gene == gene]))
    targetSourceGene$geneLabel = targetSourceGene$gene
    wh = which(!is.na(targetSourceGene$frqDEG))
    targetSourceGene$geneLabel[wh] = paste0(targetSourceGene$geneLabel[wh]," (",targetSourceGene$frqDEG[wh],")")

    # n.b., no sense to add in how many sample pairs the target gene is target gene
    #  because its always in almost all sample pairs, but not necessarily from the top source genes
               
    # save log impact for later:
    logImp = round(perGroupMeanLogImp[[groupName]][ targetSourceGene$gene[ targetSourceGene$direction=="target"] ],3)
    targetSourceGene$meanLogImpact = logImp[ targetSourceGene$gene ]
    
    # plot:
    circos.par(points.overflow.warning = FALSE)
    circos.par("track.height" = 0.05)
    # initialize with max extension for each chromosome (0..chr size):
    circos.initializeWithIdeogram(species = "hg19", chromosome.index = chrSizes$chr, plotType="")

                            
    # labels and ticks for target and source genes:
    circos.labels(sectors = targetSourceGene$chr, x = targetSourceGene$pos, labels = targetSourceGene$geneLabel, 
                                  side = "outside", connection_height=mm_h(5), cex=0.7, 
                          line_col=alpha(geneColors[targetSourceGene$direction],0.3),
                         col=geneColors[targetSourceGene$direction])

    circos.track(chrSizes$chr, ylim = c(0,1),
        panel.fun = function(x, y) {
            circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(3), gsub("chr","",CELL_META$sector.index), cex=0.9)
            circos.axis(labels=F, major.tick=F)
    })
                            
    circos.par(track.margin=c(0,0)) 

    # add colored dot for target source & genes in catagory
    circos.trackPoints(sectors = targetSourceGene$chr, x = targetSourceGene$pos, y = rep(.5,nrow(targetSourceGene)), 
                       col = categColors[targetSourceGene$categ], pch = 16, cex = 0.5)

    # links between source and target genes:
    for(i in 1:nrow(topMeanGroupInteractions)) {
        l = topMeanGroupInteractions[i,]
        circos.link(l$source_chr, l$source_pos, l$target_chr, l$target_pos, h = 0.7, directional = 1, 
            col = c(alpha("darkorange",0.5),alpha("cornflowerblue",0.5))[l$direction],
            #col = c(rgb(1,0.7,0.5,0.5),rgb(0.76,0.86,1,0.5))[l$direction],# lighter colors
            arr.length=0.15)  
    }
    title(main = samplePairGroupNiceNames[groupName], cex.main=0.9)
                                                                             
    # return target genes of the group for literature research
    list(targetGenes = targetSourceGene[ targetSourceGene$direction=="target",],
         sourceGenes = targetSourceGene[ targetSourceGene$direction!="target",])
}

l = plotGroupCircosConsistent(groupName = "SG1")
targetGenesHigher = l$targetGenes
sourceGenesHigher = l$sourceGenes

l = plotGroupCircosConsistent(groupName = "SG3")
targetGenesSlightlyLower = l$targetGenes
sourceGenesSlightlyLower = l$sourceGenes

l = plotGroupCircosConsistent(groupName = "SG2")
targetGenesLower = l$targetGenes
sourceGenesLower = l$sourceGenes

png(filename = paste0(outDirectory,"Figure-5-circos-topImpactedGenes.png"),
    width = 4000,height = 4000, pointsize = 15, res=400)
par(mfrow = c(2,2))
plotGroupCircosConsistent(groupName = "SG1")
plotGroupCircosConsistent(groupName = "SG3")
plotGroupCircosConsistent(groupName = "SG2")
dev.off()

head(geneAnnot)

exprImpForSave = NULL
for(g in c(targetGenesHigher$gene, targetGenesSlightlyLower$gene, targetGenesLower$gene)) {
    exprImpForSave = rbind(exprImpForSave, data.frame(stringsAsFactors = F, gene = g))
}
exprImpForSave$cluster = c(rep("higher",nrow(targetGenesHigher)), 
                rep("slightly lower",nrow(targetGenesSlightlyLower)),
                rep("lower",nrow(targetGenesLower)))
exprImpForSave$mean_log_impact = c(targetGenesHigher$meanLogImpact, targetGenesSlightlyLower$meanLogImpact, targetGenesLower$meanLogImpact)
exprImpForSave$chr = c(targetGenesHigher$chr, targetGenesSlightlyLower$chr, targetGenesLower$chr)
exprImpForSave = cbind(exprImpForSave,
           geneAnnot[ c(targetGenesHigher$gene, targetGenesSlightlyLower$gene, targetGenesLower$gene),
                         c("V4","V5","V7")])
exprImpForSave = exprImpForSave[,c(1,4:7,2:3)]
colnames(exprImpForSave)[3:5] = c("start","end","strand")
head(exprImpForSave,3)

exprImpForSavePerGroup = list()
for ( sg in names(samplePairPerSubgroup)) {
    sg_nice = samplePairGroupNiceNames[sg]
    # for each sample pair, mean impact from each altered gene on each target gene:
    samplePairGroupMeanLogImp = sapply(samplePairPerSubgroup[[sg]],
                       function(sp) colMeans(logMedianImpRatios[[sp]],na.rm = T))
    impPerSP = 
        samplePairGroupMeanLogImp[ exprImpForSave[ exprImpForSave$cluster == sg_nice,"gene"],]
    exprImpForSavePerGroup[[sg]] = 
        data.frame(exprImpForSave[ exprImpForSave$cluster == sg_nice,],
                   sd_log_impact = apply(impPerSP, 1, sd), 
                  impPerSP)
}

head(exprImpForSavePerGroup$SG1,3)

# calc expression ratios as in network-propagation-diff-methyl-expr-2021-11-10.r.ipynb 
voomExpr = read.csv(file = "regNet/Data/MelBrainSys-expression.csv", header = T, sep = "\t", stringsAsFactors = F)
#head(voomExpr[,-c(1:3)],4)
quantile(as.matrix(voomExpr[,-c(1:3)])) # 
# calc log-ratios for each sample pair:
all(unlist(sampleMapping) %in% colnames(voomExpr))
voomExprRatios = sapply(sampleMapping, function(pair) voomExpr[,pair[1]]-voomExpr[,pair[2]])
rownames(voomExprRatios) = voomExpr$geneSymbol
head(voomExprRatios,3)

for ( sg in names(samplePairPerSubgroup)) {
    genes = exprImpForSavePerGroup[[sg]]$gene
    groupVoomExpr =  voomExprRatios[genes,samplePairPerSubgroup[[sg]]]
    groupVoomExprMean = rowMeans(voomExprRatios[genes,samplePairPerSubgroup[[sg]]])
    groupVoomExprSD = apply(voomExprRatios[genes,samplePairPerSubgroup[[sg]]],1,sd)
    exprImpForSavePerGroup[[sg]] = 
        data.frame(exprImpForSavePerGroup[[sg]], 
                   mean_log_expression = groupVoomExprMean, sd_log_expression = groupVoomExprSD,
                  groupVoomExpr)
    # remove .1 .2 .. from sample names
    colnames(exprImpForSavePerGroup[[sg]]) = gsub("\\.\\d+$","",colnames(exprImpForSavePerGroup[[sg]]))
}

exprImpForSavePerGroup$SG1

write.xlsx(exprImpForSavePerGroup$SG1, 
           file = paste0(outDirectory,"SupplTable-S10-target-gene-candidates.xls"),
           sheetName = "higher cluster", append = F, row.names = F)
write.xlsx(exprImpForSavePerGroup$SG3, 
           file = paste0(outDirectory,"SupplTable-S10-target-gene-candidates.xls"),
           sheetName = "sligthly lower cluster", append = T, row.names = F)
write.xlsx(exprImpForSavePerGroup$SG2, 
           file = paste0(outDirectory,"SupplTable-S10-target-gene-candidates.xls"),
           sheetName = "lower cluster", append = T, row.names = F)

neededCols = c("gene","chr","start","end","strand","cluster","mean_log_impact","mean_log_expression")
exprImpForPlot = rbind(exprImpForSavePerGroup$SG1[,neededCols],
                      exprImpForSavePerGroup$SG3[,neededCols],
                      exprImpForSavePerGroup$SG2[,neededCols])
colnames(exprImpForPlot)[6] = "subgroup"
exprImpForPlot$subgroup = ifelse(exprImpForPlot$subgroup=="higher","SG1",
                                 ifelse(exprImpForPlot$subgroup=="lower","SG2","SG3"))

head(exprImpForPlot,3)
tail(exprImpForPlot,3)

exprImpForPlot

# get middle of cluster group impact ratios to draw intersection lines
o = order(exprImpForPlot$mean_log_impact)
cbind(exprImpForPlot$mean_log_impact[o],exprImpForPlot$subgroup[o])
intersect1 = mean(exprImpForPlot$mean_log_impact[o][c(10:11)])
intersect2 = mean(exprImpForPlot$mean_log_impact[o][c(21:22)])
intersect1
intersect2

p = ggplot(data = exprImpForPlot, mapping = aes(x = mean_log_impact, y = mean_log_expression, 
                                                col=subgroup)) + #, label=gene)) + 
    geom_vline(xintercept = intersect1, col="grey") + geom_vline(xintercept = intersect2, col="grey") + 
    geom_smooth(method='lm', col="darkgrey",mapping = aes(x = mean_log_impact, y = mean_log_expression, group=1) ) +
    geom_point(size=1) + 
    geom_text_repel(hjust=-0.1, size=2, col="#222222", segment.colour="darkgrey", fontface="italic", 
                    mapping = aes(label=gene)) + 
    scale_color_manual(values = annotColors$subgroup) + 
        xlab(expression(impact~log["2"]~ratio)) + ylab(expression(expression~log["2"]~ratio))
   # annotate(geom = "text",x = -1.5, y = 2, label=paste0("p = ",round(testRes$p.value,3),", R = ",round(testRes$estimate,2)))
p

png(filename = paste0(outDirectory,"SupplFigure-S8-target-gene-candidates-scatterplot.png"),
    width = 1600,height = 1200, res=300)
p
dev.off()
