# calculates for each metastases pair 
#   average impacts log-ratios (brain/extra-cranial) on pathways
#   median impact log-ratios of all altered genes on all genes
#   average impacts log-ratios on each gene
#   
# 
# from the altered genes of all met. pairs
# as basis for other scripts
# saved for later: metPairs-impactRatios-onPathways.rds,
#   logMedianImpRatios, metPairs-impactRatios-onAllTargetGenes.rds

library(parallel)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/")

setwd(basePath)

localRlibs = paste0(basePath,"conda/lib/R/library/")
library(regNet)

localGeneCutoff = 30
pValCutoff = 0.001
numOfNWs = 2
networkName = "TcgaMelanomaExprMeth"
dataSetName = "MelBrainSys"
nbCPUs = 10 
colSumsThreshold = 0.001

pwCategories = readRDS(file = "pathway-categories.rds")
pwCategories$signalingPWs$PPAR

genePerPathwayList = unlist(unname(pwCategories),recursive = F)

lapply(head(genePerPathwayList,3), head)

tmp = readRDS(file = "altered-genes-per-patient.rds")
regGenesPerPat = tmp$regGenesPerPat

# define met. pairs, and annotation colors 
tmp = readRDS(file = "annotation/samplePairs-annotation-colors-clusters.rds")
sampleMapping = tmp$sampleMapping

getRealImpacts = function(
    samplePair, sample, sourceGenes, targetGenes, pValCutoff, localGeneCutoff, colSumsThreshold, numOfNWs) {
    nwSubdirPrefix = 'TrainNetwork-'
    impacts = matrix(data = 0,nrow = length(sourceGenes),ncol = 0, dimnames = list(c(sourceGenes)))
    for(networkNum in 1:numOfNWs){
        path = paste0(myPath,nwSubdirPrefix,networkNum,"/")
        impactsTmp = getAverageImpacts_PatientSpecificAbsoluteImpacts(
            patient = sample,sourceGenes = sourceGenes,targetGenes = targetGenes, dataSetName = dataSetName, 
            networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, 
            colSumsThreshold = colSumsThreshold, path = path, outputFile = "tmp")
        impacts = cbind(impacts,NA)
        impacts[impactsTmp$SourceGene,ncol(impacts)] = impactsTmp$AverageImpactOnTargetGenes
    }
    impacts
}

allSamplImp = NULL
nbTotalIt = length(sampleMapping)*2
currIt=0
startt = Sys.time()
for(samplePair in names(sampleMapping)) {
    message(samplePair)
    samplePairList = NULL
    for (sampleIdx in 1:2) {
        currIt = currIt + 1
        sample = sampleMapping[[samplePair]][sampleIdx]
        message("  sampleIdx ",sampleIdx, " ",sample)
        impPerPW = mclapply(X = genePerPathwayList,  mc.cores = nbCPUs,FUN = function(pwGenes) {
            getRealImpacts(samplePair = samplePair, sample = sample,
                                 sourceGenes = regGenesPerPat[[samplePair]], targetGenes = pwGenes, 
                                 pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, 
                                 colSumsThreshold = colSumsThreshold, numOfNWs = numOfNWs)
        })
        samplePairList = append( samplePairList, list(impPerPW) )
        names(samplePairList)[sampleIdx] = sample

        endt = Sys.time()
        totalTime = difftime(endt,startt, units = "min")*nbTotalIt/currIt
        restTime = round(totalTime * (nbTotalIt-currIt)/nbTotalIt)
        message("    ",currIt," of ",nbTotalIt,", ",restTime," min still needed")
    }
    allSamplImp = append( allSamplImp, list(samplePairList) )
}
names(allSamplImp) = names(sampleMapping)
endt = Sys.time()
message(difftime(endt,startt,units = "min")," min in total needed")
# 15 min in total for 2 networks 10 metastases pairs with 10 CPUs 

# sample pair - sample mate - pathway - matrix of genes x networks: mean impact of gene on pathway in sample mate
head(allSamplImp[[1]][[1]][[2]])
class(allSamplImp[[1]][[1]][[37]])
dim(allSamplImp[[1]][[1]][[37]])
head(allSamplImp[[1]][[1]][[37]])

# transform values to 1 data.frame for later
# for each met pair and pathway
# 1 column for impacts in brain / non-brain, log-ratio

metPairImpacts = NULL
oldWarningOpt <- getOption("warn")
options(warn = -1)
for (samplePair in names(sampleMapping)) {
    message(samplePair)
    dataAllPw = NULL
    for (pw in names(genePerPathwayList)) {
        sample1 = sampleMapping[[samplePair]][[1]]
        sample2 = sampleMapping[[samplePair]][[2]]
        realImp1 = allSamplImp[[samplePair]][[sample1]][[pw]]
        realImp2 = allSamplImp[[samplePair]][[sample2]][[pw]]
        dataImp = data.frame(
            median_brain_impact = apply(realImp1,1,function(r) median(na.omit(r))),
            median_nonBrain_impact = apply(realImp2,1,function(r) median(na.omit(r))),
                        stringsAsFactors=F )
        dataImp$log2_median_impact_ratio = 
                        log2(dataImp$median_brain_impact/dataImp$median_nonBrain_impact)
        dataImp$log2_mean_impact_ratio = 
            log2(apply(realImp1,1,function(r) 
                mean(na.omit(r)))/apply(realImp2,1,function(r) mean(na.omit(r))))
        data1pw = data.frame(
            samplePair=samplePair, gene = rownames(allSamplImp[[samplePair]][[1]][[pw]]), pathway = pw, 
            dataImp, stringsAsFactors = F)
        metPairImpacts = rbind(metPairImpacts, data1pw)
    }
}
options(warn = oldWarningOpt)
#names(metPairImpacts) = names(sampleMapping)[1:length(metPairImpacts)]

dim(metPairImpacts)
head(metPairImpacts)

saveRDS(object = metPairImpacts, file = "metPairs-impactRatios-onPathways.rds")

myFullDataSet = loadGeneExpressionAndCopyNumberDataSet(
    geneExpressionFile = "MelBrainSys-expression.csv", 
    geneCopyNumberFile = "MelBrainSys-methylation.csv", path = paste0(myPath,"Data/") )

allGenesNWs = myFullDataSet$genes
numGenes = length(allGenesNWs)
numGenes; head(allGenesNWs)

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
for( samplePair in names(sampleMapping)) {
    message(samplePair)
    samples = sampleMapping[[samplePair]]
    message("  ",samples[1])
    paramSets = lapply(samples[1], function(sample) lapply(1:numOfNWs, function(nw) c(sample=sample, nw=nw)))
    paramSets = unlist(paramSets,recursive = F)
    DGEsDMPs = regGenesPerPat[[samplePair]]
    DGEsDMPs = intersect(DGEsDMPs, allGenesNWs) # only DGEs that are in our data
    message("    ",length(DGEsDMPs), " DGEsDMPs")
    brainImpacts = calc1sampleImpacts(paramSets = paramSets, DGEsDMPs = DGEsDMPs,nbCPUs = nbCPUs)

    message("  ",samples[2])
    paramSets = lapply(samples[2], function(sample) lapply(1:numOfNWs, function(nw) c(sample=sample, nw=nw)))
    paramSets = unlist(paramSets,recursive = F)
    DGEsDMPs = regGenesPerPat[[samplePair]]
    DGEsDMPs = intersect(DGEsDMPs, allGenesNWs) # only DGEs that are in our data
    message("    ",length(DGEsDMPs), " DGEsDMPs")
    nonBrainImpacts = calc1sampleImpacts(paramSets = paramSets, DGEsDMPs = DGEsDMPs,nbCPUs = nbCPUs)
    tmp = list(brainImpacts, nonBrainImpacts)
    names(tmp) = samples
    allImpacts[[length(allImpacts)+1]] = tmp
}
endt = Sys.time()
round(difftime(endt, startt,units="min"))
names(allImpacts) = names(sampleMapping)
# ca 5 min for 2 networks and 10 met. pairs

# calc. median impacts for each sample over all 25 networks
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

logMedianImpRatios$P03_BLun[1:4,1:3]
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
head(meanLogMedianImpRatiosMat)
dim(meanLogMedianImpRatiosMat)

saveRDS(object = list(logMedianImpRatios = logMedianImpRatios,
                      meanLogMedianImpRatiosMat = meanLogMedianImpRatiosMat), 
        file="metPairs-impactRatios-onAllTargetGenes.rds")
