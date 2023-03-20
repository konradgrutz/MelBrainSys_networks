library(parallel)
library(pheatmap)
library(xlsx)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/validation-cohort/")
setwd(myPath)

localRlibs = paste0(basePath,"conda/lib/R/library/")
library(regNet)

outDirectory = paste0(basePath,"FiguresTables/")

# define met. pairs, and annotation colors 
tmp = readRDS(file = "../../annotation/samplePairs-annotation-colors-clusters.rds")
patient_colors = tmp$patient_colors
tissue_colors = tmp$tissue_colors
discSampleMapping = tmp$sampleMapping
discSubgroupPerSamplePair = tmp$subgroupPerSamplePair
colAnnot = tmp$colAnnot
annotColors = tmp$annotColors

valiSampleMapping = list(P78_BIn = c('P78_Br_1a','P78_In_3d'), P77_BLy = c('P77_Br_e','P77_Ly_1d'),
                  P74_BLy = c('P74_Br_1a','P74_Ly_1a'), P13_BLy = c('P13_Br_2','P13_Ly_e'),
                  P111_BLy = c('P111_Br_b','P111_Ly_1b'), P108_BLy = c('P108_Br_a','P108_Ly_2a'),
                  P101_BLi = c('P101_Br_a','P101_Li_b'),P106_BLy = c('P106_Br_ns1','P106_Ly_1b'),
                   P107_BLu = c('P107_Br_c','P107_Lu_a'))
valiSampleMapping

nbCPUs = 25
localGeneCutoff = 30
pValCutoff = 0.001
colSumsThreshold = 0.001
numOfNWs = 2
networkName = "TcgaMelanomaExprMeth"
dataSetName = "MelBrainSys"

tmp = readRDS(file = "../../altered-genes-per-patient.rds")
regGenesPerPat = tmp$regGenesPerPat
formerRegGenes = unique(sort(unname(unlist(regGenesPerPat))))
head(formerRegGenes)
length(formerRegGenes)
sapply(regGenesPerPat,length)

myFullDataSet = loadGeneExpressionAndCopyNumberDataSet(
    geneExpressionFile = "MelBrainSys_ExpressionData_2022_allNeededGenes_regNet.txt", 
    geneCopyNumberFile = "MelBrainSys_MethylationData_2022_allNeededGenes_regNet.txt", 
    path = "Data/" )
allGenesNWs = myFullDataSet$genes
numGenes = length(allGenesNWs)

calc1sampleImpacts = function(paramSets=paramSets,DGEsDMPs=DGEsDMPs,nbCPUs=nbCPUs) {
    allImpacts1 = mclapply(X = paramSets, mc.cores = nbCPUs, FUN = function(params) {
        sample = params["sample"]
        nw = params["nw"]
        path = paste0(myPath,"TrainNetwork-",nw,"/")
        impacts = getImpacts_PatientSpecificAbsoluteImpacts(
            patient = sample, sourceGenes = DGEsDMPs, targetGenes = allGenesNWs, 
            dataSetName = dataSetName,networkName = networkName, 
            pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, 
            colSumsThreshold = colSumsThreshold, path = path, 
            outputFile =  paste0("tmp-",nw,"-",sample), output=F )
        f = file(description = logFile,open = "a")
        writeLines(con = f,text = paste0("netw ",nw,", sample ",sample))
        close(f)
        rownames(impacts) = DGEsDMPs
        as.matrix(impacts[,-1])
    })
}

# extract for each sample pair, each sample, each one of the 25 NWs, DGEsDMPs impacts on all other genes:
logFile="log-getImpacts.txt";
if(file.exists(logFile)){file.remove(logFile)}
startt = Sys.time()
allImpacts = NULL
for( samplePair in names(valiSampleMapping)) {
    message(samplePair)
    samples = valiSampleMapping[[samplePair]]
    message("  ",samples[1])
    paramSets = lapply(samples[1], function(sample) lapply(1:numOfNWs, function(nw) c(sample=sample, nw=nw)))
    paramSets = unlist(paramSets,recursive = F)
    DGEsDMPs = formerRegGenes
    DGEsDMPs = intersect(DGEsDMPs, allGenesNWs) # only DGEs that are in our data
    message("    ",length(DGEsDMPs), " DGEsDMPs")
    brainImpacts = calc1sampleImpacts(paramSets = paramSets, DGEsDMPs = DGEsDMPs,nbCPUs = nbCPUs)

    message("  ",samples[2])
    paramSets = lapply(samples[2], function(sample) lapply(1:numOfNWs, function(nw) c(sample=sample, nw=nw)))
    paramSets = unlist(paramSets,recursive = F)
    nonBrainImpacts = calc1sampleImpacts(paramSets = paramSets, DGEsDMPs = DGEsDMPs,nbCPUs = nbCPUs)
    tmp = list(brainImpacts, nonBrainImpacts)
    names(tmp) = samples
    allImpacts[[length(allImpacts)+1]] = tmp
}
endt = Sys.time()
round(difftime(endt, startt,units="min"))
names(allImpacts) = names(valiSampleMapping)
# ca 5 min for 2 networks and 10 met. pairs

# calc. median impacts for each sample over all networks
#  and from these calc. log median impact ratios for each sample pair
nbCPUs=8
if(file.exists(logFile)){file.remove(logFile)}
startt=Sys.time()
logMedianImpRatios = NULL

for(sp in names(allImpacts)) {
    message(sp)
    medianImpacts = NULL
    for( sampleName in names(allImpacts[[sp]])) {
        message("   ",sampleName)
        impacts = allImpacts[[sp]][[sampleName]]
        message("   ",length(impacts))
        numSourcGenes = nrow(impacts[[1]])
        #numSourcGenes
        arr = array(data = NA,dim = c(numOfNWs,numSourcGenes,numGenes))
        for(j in 1:length(impacts)) {
            arr[j,,]=impacts[[j]]
        }
        medianSampleImpacts = mclapply(1:numGenes,mc.cores = nbCPUs,function(c) {
            sapply(1:numSourcGenes, function(r)
                median(x=arr[,r,c], na.rm=T)
            )
        })
        medianSampleImpacts = do.call(cbind,medianSampleImpacts)
        rownames(medianSampleImpacts)=rownames(impacts[[1]])
        colnames(medianSampleImpacts)=allGenesNWs
        medianImpacts[[length(medianImpacts)+1]] = medianSampleImpacts
        f = file(description = logFile,open = "a")
        writeLines(con = f,text = paste0("sample ",sampleName))
        close(f)
    }
    names(medianImpacts) = names(allImpacts[[sp]])
    # calc. log median impact ratios 
    logMedianImpRatios[[length(logMedianImpRatios)+1]] = 
            log2(medianImpacts[[1]]/medianImpacts[[2]])
}
# 1-3 min
names(logMedianImpRatios) = names(allImpacts)
endt=Sys.time()
round(difftime(endt, startt,units="min"))

logMedianImpRatios[[1]][1:2,1:5]
# for each metastases pair source genes x target genes matrix with impact log ratios on these target genes

# and now mean of log median imp ratio over source genes
meanLogMedianImpRatios=NULL
for(sp in names(logMedianImpRatios)) {
    meanLogMedianImpRatios[[length(meanLogMedianImpRatios)+1]] = 
        apply(logMedianImpRatios[[sp]], 2,function(m) mean(m,na.rm=T))
}
names(meanLogMedianImpRatios) = names(logMedianImpRatios)
names(meanLogMedianImpRatios)[1:3]
head(meanLogMedianImpRatios[[1]],4)
meanLogMedianImpRatiosMat = sapply(meanLogMedianImpRatios,function(v) v)
head(meanLogMedianImpRatiosMat,3)
dim(meanLogMedianImpRatiosMat)

meanLogMedianImpRatiosMatValiCoh = meanLogMedianImpRatiosMat

tmp = readRDS(file="../../metPairs-impactRatios-onAllTargetGenes.rds")
meanLogMedianImpRatiosMatDiscCoh = tmp$meanLogMedianImpRatiosMat
head(meanLogMedianImpRatiosMatDiscCoh,2)

commonGenes = intersect(rownames(meanLogMedianImpRatiosMatDiscCoh), 
                        rownames(meanLogMedianImpRatiosMatValiCoh))
length(commonGenes)

impRatiosDiscValid = cbind(meanLogMedianImpRatiosMatValiCoh[ commonGenes,],
                                  meanLogMedianImpRatiosMatDiscCoh[ commonGenes,])
head(impRatiosDiscValid,3)

valiSamples = names(valiSampleMapping)
discSamples = names(discSampleMapping)

# assign subgroups to met. pairs for valid. cohort.
# do this manually depending an where they cluster in the heatmap below
#  I didn't do this here, because it will be different anyway when using 25 instead of 2 networks
valiSubgroupPerSamplePair = c("SG2","SG2","SG3","SG3","SG3","SG3","SG3","SG3","SG3")
names(valiSubgroupPerSamplePair) = c('P13_BLy','P108_BLy','P78_BIn','P77_BLy','P74_BLy','P111_BLy',
                         'P101_BLi','P106_BLy','P107_BLu')
valiSubgroupPerSamplePair

colAnnotNew = data.frame( 
    cohort = c(rep("validation cohort",length(valiSamples)), 
               rep("discovery cohort",length(discSamples)) ), 
        row.names =  c(valiSamples, discSamples ),
    subgroup = c(discSubgroupPerSamplePair[discSamples], valiSubgroupPerSamplePair[valiSamples]) )
# use colors of the circos plot, better group colors: green orange
colAnnotColorNew = list(
    cohort = c("discovery cohort"="#eebb77", "validation cohort"="#22CC77"),
    subgroup = c("SG1"="tomato1","SG2"="cornflowerblue","SG3"="lightskyblue1"))

drows = dist(impRatiosDiscValid, method = "euclidean", )
dcols = dist(t(impRatiosDiscValid), method = "euclidean")
rg=range(impRatiosDiscValid, na.rm=T)
ph = pheatmap(mat = impRatiosDiscValid,clustering_distance_rows = drows, 
              cellwidth = 23, #cellheight=0.01,
    annotation_colors = colAnnotColorNew,
    annotation_col = colAnnotNew,
    cluster_rows = F, show_rownames = F,
    clustering_distance_cols = dcols, display_numbers = F, angle_col = 45, 
    color = colorRampPalette(c("blue","white","red"))(100),
    breaks = c(seq(rg[1],0,length.out = 51),seq(0,rg[2],length.out = 50)[-1]))

png(file = paste0(outDirectory, "SupplFigure-S7-MeanLogImp-clustering-discovery-validation-cohort.png"),
    res = 250,width = 2400,height = 2000)
ph
dev.off()

valiCohGroups = list(SG2 = c(names(valiSubgroupPerSamplePair[ which(valiSubgroupPerSamplePair == "SG2")])), 
                    SG3 = c(names(valiSubgroupPerSamplePair[ which(valiSubgroupPerSamplePair == "SG3")])))
valiCohGroups

SPs = valiCohGroups[['SG2']]
m = meanLogMedianImpRatiosMatValiCoh[,SPs]
m = m[which(apply(m,1,function(r) all(r<0))),]
dim(m)
o=order(rowMeans(m))
valiCohSG2bottomGenes = rowMeans(m[o,])
head(round(valiCohSG2bottomGenes,3))

SPs = valiCohGroups[['SG3']]
m = meanLogMedianImpRatiosMatValiCoh[,SPs]
m = m[which(apply(m,1,function(r) length(which(r<0))>=6)),] 
#   this is the relaxed condition for SG3: 
#   in at least 6 of the 7 met. pairs the impact ratio must be consistent
dim(m)
o=order(rowMeans(m))
valiCohSG3bottomGenesRelaxed = rowMeans(m[o,])
head(round(valiCohSG3bottomGenesRelaxed,3))

# discovery cohort bottom target genes
discCohSG2bottomGenes = read.xlsx(file = "../../FiguresTables/SupplTable-S10-target-gene-candidates.xls",
                          sheetName = "lower cluster")
discCohSG3bottomGenes = read.xlsx(file = "../../FiguresTables/SupplTable-S10-target-gene-candidates.xls",
                          sheetName = "sligthly lower cluster")

# in the original analyses all candidate genes of SG2 and SG3 had negative mean impact log-ratio
#  so also restrict here to only negative ones
discCohSG2bottomGenes = discCohSG2bottomGenes[ which(discCohSG2bottomGenes$mean_log_impact<0),]
discCohSG3bottomGenes = discCohSG3bottomGenes[ which(discCohSG3bottomGenes$mean_log_impact<0),]
head(discCohSG2bottomGenes,3)
head(discCohSG3bottomGenes,3)

# how many genes in discovery cohort are found in validation cohort
length(intersect(discCohSG2bottomGenes$gene, names(valiCohSG2bottomGenes))) # 7 of 10 found again
nrow(discCohSG2bottomGenes)
length(intersect(discCohSG3bottomGenes$gene, names(valiCohSG3bottomGenesRelaxed))) # 5 of 10 found again
nrow(discCohSG3bottomGenes)

# where do target genes of validation cohort rank in discovery cohort
pos = unlist(sapply(discCohSG2bottomGenes$gene, function(g) which(names(valiCohSG2bottomGenes) == g)))
round(pos/length(valiCohSG2bottomGenes)*100) # at the top n-percentile of the list

pos = unlist(sapply(discCohSG3bottomGenes$gene, function(g) which(names(valiCohSG3bottomGenesRelaxed) == g)))
round(pos/length(valiCohSG3bottomGenesRelaxed)*100) # at the top n-percentile of the list



SG2candidates = intersect(discCohSG2bottomGenes$gene, names(valiCohSG2bottomGenes))
SG3candidates = intersect(discCohSG3bottomGenes$gene, names(valiCohSG3bottomGenesRelaxed))
SG2candidates
SG3candidates

SPs = valiCohGroups[["SG2"]]
m = meanLogMedianImpRatiosMatValiCoh[,SPs]
m[ SG2candidates ,]
meanLogIR = rowMeans(m[SG2candidates,])
meanLogIR
pos = unlist(sapply(SG2candidates, function(g) which(names(valiCohSG2bottomGenes) == g)))
percentiles = round(pos/length(valiCohSG2bottomGenes)*100,1) # at the top n-percentile of the list
percentiles
outSG2 = data.frame(gene = SG2candidates, meanLogIR,percentiles,
                    m[ SG2candidates ,])
outSG2

SPs = valiCohGroups[["SG3"]]
m = meanLogMedianImpRatiosMatValiCoh[,SPs]
m[ SG3candidates ,]
meanLogIR = rowMeans(m[SG3candidates,])
meanLogIR
pos = unlist(sapply(SG3candidates, function(g) which(names(valiCohSG3bottomGenesRelaxed) == g)))
percentiles = round(pos/length(valiCohSG3bottomGenesRelaxed)*100,1) # at the top n-percentile of the list
percentiles
outSG3 = data.frame(gene = SG3candidates, meanLogIR,percentiles,
                    m[ SG3candidates ,])
outSG3

write.xlsx(x = outSG2,file = paste0(outDirectory,"Suppl-Table-S12-target-gene-ranking.xls"),
           sheetName = "SG2", row.names = F)
write.xlsx(x = outSG3,file = paste0(outDirectory,"Suppl-Table-S12-target-gene-ranking.xls"),
           sheetName = "SG3", row.names = F, append = T)
