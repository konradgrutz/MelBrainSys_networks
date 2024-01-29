# violin plots of impact log2-ratios, Suppl. Figure S4
# cluster heatmaps of log-ratio impacts on pathways, Figure 3
# Suppl. Table S7: log-impact and log-expression per pathway
# Suppl. Table S9: expression impact quadrants shared pathways

library(ggplot2)
library(gridExtra)
library(pheatmap)
library(scales) # for wrap_format, to wrap labels
library(stringr) # edit names in pheatmap
library(pvclust) # for clustering stability analyses
#library(readxl)
library(xlsx)

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

# pathway definitions
pwCategories = readRDS(file = "pathway-categories.rds")
head(pwCategories$signalingPWs$PPAR)

# pre-calculated impact log-ratios on pathways from each altered gene
metPairImpacts = readRDS(file = "metPairs-impactRatios-onPathways.rds")
head(metPairImpacts,3)
wh = which(!is.na(metPairImpacts$log2_mean_impact_ratio))
length(wh); nrow(metPairImpacts)
metPairImpactsNoNA = metPairImpacts[wh,]

# load log expression values

expr = read.csv(file = "regNet/Data/MelBrainSys-expression.csv", header = T,sep = "\t",stringsAsFactors = F)
head(expr,3)

genesPerPathway = unlist(unname(pwCategories),recursive = F)
lapply(head(genesPerPathway,3),head)

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
                      "impact log2 ratio" = impRatio,
                      "expression log2 ratio" = exprRatio,
                      subgroup = subgroupNames[subgroupPerSamplePair[samplePair]],
                      "subgroup ID" = subgroupPerSamplePair[samplePair]))
    }
}
impactExprTable = impactExprTable[ order(impactExprTable$sample_pair,
                                         impactExprTable$pathway),]

head(impactExprTable)

write.xlsx(x = impactExprTable, row.names = F,
           file = paste0(outDirectory,"SupplTableS7-log-impact-log-expression-per-pathway.xls"))

samplePairsSorted = c(samplePairPerSubgroup[[1]],samplePairPerSubgroup[[2]],samplePairPerSubgroup[[3]])
allPWs = names(genesPerPathway)

emptyPresenceTable = matrix(dimnames = c(list(allPWs), list(samplePairsSorted)),
    data = 0,nrow = length(allPWs), ncol = length(samplePairsSorted))
emptyPresenceTable[1:2,1:10]

presenceTables = list("lower impact higher expression" = emptyPresenceTable,
                     "higher impact lower expression" = emptyPresenceTable,
                     "higher impact higher expression" = emptyPresenceTable,
                     "lower impact lower expression" = emptyPresenceTable)

lapply(presenceTables,function(t) t[1:2,1:4])

impDir = sign(impactExprTable$impact.log2.ratio)
exprDir = sign(impactExprTable$expression.log2.ratio)
for (i in 1:nrow(impactExprTable)) {
    type = paste0(ifelse(impDir[i] > 0,"higher","lower"), " impact ",
                  ifelse(exprDir[i] > 0,"higher","lower")," expression", collapse="")
    #message(paste0(round(impactExprTable[i,"impact.log2.ratio"],2), " ",
    #       round(impactExprTable[i,"expression.log2.ratio"],2)," -> ",type))
    SP = impactExprTable$sample_pair[i]
    pw  = impactExprTable$pathway[i]
    #message(type, " ", SP," ", pw, paste0(impactExprTable[i,], collapse=" "))
    presenceTables[[type]][pw,SP]=1
}

impactExprTable[ sample(1:nrow(impactExprTable),10),]

presenceTables$`higher impact lower expression`["inositol phosphate","P08_BSof_1"]

lapply(presenceTables,function(t) t[1:2,1:4])

write.xlsx(x = presenceTables[[1]], sheetName = names(presenceTables)[1], append = F, 
           file = paste0(outDirectory,"SupplTable-S9-expression-impact-shared-pathways.xls"))
write.xlsx(x = presenceTables[[2]], sheetName = names(presenceTables)[2], append = T, 
           file = paste0(outDirectory,"SupplTable-S9-expression-impact-shared-pathways.xls"))
write.xlsx(x = presenceTables[[3]], sheetName = names(presenceTables)[3], append = T, 
           file = paste0(outDirectory,"SupplTable-S9-expression-impact-shared-pathways.xls"))
write.xlsx(x = presenceTables[[4]], sheetName = names(presenceTables)[4], append = T, 
           file = paste0(outDirectory,"SupplTable-S9-expression-impact-shared-pathways.xls"))

# beatify the tables in excel
# and manually make 5th sheet where sheet 1-4 are merged

p1 = ggplot(data = metPairImpactsNoNA, mapping = aes(x=samplePair, y=log2_median_impact_ratio)) + 
    geom_violin(alpha=1,fill="lightblue",draw_quantiles = c(0.25, 0.5, 0.75))+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    geom_hline(yintercept = 0,col="red",alpha=0.5) + 
    xlab("metastases pair") + ylab(expression(median~impact~log[2]~ratio)) 

p1

png(filename = paste0(outDirectory,"SupplFigure-S6-log-impact-ratio-violins.png"),
    width = 2000,height = 2000, res= 400)
p1
dev.off()

# one heatmap scale for all 3 plots
# therefore load all data before and get their max range
rg = c(0,0)
for(i in 1:length(pwCategories)) {
   # pwCategNice = pwCatNamesNice[i]
    pwCateg = names(pwCategories)[i]
    meanLogRatios = sapply(unique(metPairImpactsNoNA$samplePair), function(sp) 
        sapply(names(pwCategories[[pwCateg]]), function(pw) 
            mean(metPairImpactsNoNA$log2_median_impact_ratio[ 
                metPairImpactsNoNA$pathway == pw & metPairImpactsNoNA$samplePair == sp])))
    rg[1] = min(rg[1],min(meanLogRatios))
   rg[2] = max(rg[2],max(meanLogRatios))
}
rg

annotColors

# heatmap of mean impact log-ratio per pathway/sample
heatmaps=NULL
for(i in 1:length(pwCategories)) {
    pwCateg = names(pwCategories)[i]
    meanLogRatios = sapply(unique(metPairImpactsNoNA$samplePair), function(sp) 
        sapply(names(pwCategories[[pwCateg]]), function(pw) 
            mean(metPairImpactsNoNA$log2_median_impact_ratio[ 
                metPairImpactsNoNA$pathway == pw & metPairImpactsNoNA$samplePair == sp])))
    drows = dist(meanLogRatios, method = "euclidean")
    dcols = dist(t(meanLogRatios), method = "euclidean")
    showAnnotLegend = ifelse(i==4,T,F)
       ph = pheatmap(mat = meanLogRatios,clustering_distance_rows = drows, clustering_distance_cols = dcols,                     
                     angle_col = 45,annotation_col = colAnnot, annotation_colors = annotColors, 
                  color = colorRampPalette(c("blue","white","red"))(100),
                  breaks = c(seq(rg[1],0,length.out = 51),seq(0,rg[2],length.out = 50)[-1]),
             main = pwCateg, annotation_legend = showAnnotLegend, legend = showAnnotLegend,
                  labels_row=str_wrap(rownames(meanLogRatios),width = 30))
   heatmaps[[length(heatmaps)+1]] = ph[[4]] # for later plotting side by side with grid.arrange
    cat(unname(ph$tree_row$labels[ph$tree_row$order]), sep = ", ")
}

# plot last PW heatmap double with the legend
ph = pheatmap(mat = meanLogRatios,clustering_distance_rows = drows, clustering_distance_cols = dcols,                     
                 angle_col = 45,annotation_col = colAnnot, annotation_colors = annotColors, 
              color = colorRampPalette(c("blue","white","red"))(100),
              breaks = c(seq(rg[1],0,length.out = 51),seq(0,rg[2],length.out = 50)[-1]),
         main = paste0(pwCateg), annotation_legend = T, legend = T,
              labels_row=str_wrap(rownames(meanLogRatios),width = 30))
heatmaps[[length(heatmaps)+1]] = ph[[4]] # for later plotting side by side with grid.arrange


length(heatmaps)

svg(filename = paste0(outDirectory,"Figure-3-medianLogImpRatios-cluster-heatmap.svg"), width = 15, height = 15)
grid.arrange(grobs=heatmaps, ncol=2)#, widths=c(1.05,1,1.25))
dev.off()

# one heatmap scale for all 3 plots
# therefore load all data before and get their max range
rg = c(0,0)
for(i in 1:length(pwCategories)) {
   # pwCategNice = pwCatNamesNice[i]
    pwCateg = names(pwCategories)[i]
    meanLogRatios = sapply(unique(metPairImpactsNoNA$samplePair), function(sp) 
        sapply(names(pwCategories[[pwCateg]]), function(pw) 
            mean(metPairImpactsNoNA$log2_median_impact_ratio[ 
                metPairImpactsNoNA$pathway == pw & metPairImpactsNoNA$samplePair == sp])))
    rg[1] = min(rg[1],min(meanLogRatios))
   rg[2] = max(rg[2],max(meanLogRatios))
}
rg

clusterStapPlots = NULL
for(i in 1:length(pwCategories)) {
    pwCateg = names(pwCategories)[i]
    meanLogRatios = sapply(unique(metPairImpactsNoNA$samplePair), function(sp) 
        sapply(names(pwCategories[[pwCateg]]), function(pw) 
            mean(metPairImpactsNoNA$log2_median_impact_ratio[ 
                metPairImpactsNoNA$pathway == pw & metPairImpactsNoNA$samplePair == sp])))
    clusterStab = pvclust(data = meanLogRatios, method.dist="euclidean", 
                          method.hclust="complete", nboot=10000, parallel=as.integer(30))
    plot(clusterStab, main = pwCateg )      
    clusterStapPlots = append(clusterStapPlots,list(clusterStab))
}
names(clusterStapPlots) = names(pwCategories)

length(clusterStapPlots)

par(mfrow = c(2,2))
plot(clusterStapPlots[[1]], main = names(clusterStapPlots)[1] )      
plot(clusterStapPlots[[2]], main = names(clusterStapPlots)[2] )      
plot(clusterStapPlots[[3]], main = names(clusterStapPlots)[3] )      


svg(filename = paste0(outDirectory,"Figure-3-medianLogImpRatios-cluster-heatmap-stability.svg"), 
    width = 15, height = 15)
par(mfrow = c(2,2))
plot(clusterStapPlots[[1]], main = names(clusterStapPlots)[1] )      
plot(clusterStapPlots[[2]], main = names(clusterStapPlots)[2] )      
plot(clusterStapPlots[[3]], main = names(clusterStapPlots)[3] )      

dev.off()
