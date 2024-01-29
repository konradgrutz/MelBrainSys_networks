# predicts expression of the test data with the trained and random networks
# predicts expression of MelBrainSys data
# creates correlation plots of these predicted and original expression:
#   Figure-2-gene-expr-prediction-correl.png

library(ggplot2)
library(gridExtra)
library(parallel)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/")

localRlibs = paste0(basePath,"conda/lib/R/library/")
library(regNet)

setwd(myPath)

outDirectory = paste0(basePath,"FiguresTables/")

output=T
networkName = "TcgaMelanomaExprMeth"

localGeneCutoff = 30
pValCutoff = 0.001
nbCPUs = 10 
numOfNWs = 2 # for test purposes only 2 networks. Hence note, the results can deviate from the paper.
numRandNWinstances = 2
dataSetsPath = "Data/"

# random network instance is saved
# under NetworkModel/WholeNetwork/ using a pre-defined file-naming convention with random
# network name prefix paste0( "RandomNetwork_", randomNetworkInstance, 
#  "_PVal-ueCutoff_", pValCutoff, "_BasedOn_", networkName )

create1randomNW = function(parameterSet) {
    system(paste("echo 'calculating set:",parameterSet$runningNumber,"'")) 
     #  -> progress output in console where you started jupyter notebook
    determineRandomNetworkWithFilteringForSignificantPredictors(output = T,
        networkName = parameterSet$networkName, pValCutoff =  parameterSet$pValCutoff, 
        randomNetworkInstance = parameterSet$nwInstanceNb, 
        path = parameterSet$path )
}

# make parameter sets for creation of random nw:
randNWinstancesRange = 1:numRandNWinstances
nwRange = 1:numOfNWs
parameterSets = NULL
i=0
for(nwInstanceNb in randNWinstancesRange ) {
    for(nwSubdir in paste0("TrainNetwork-",nwRange)) {
        i = i + 1
        parameterSets[[i]] = list(networkName=networkName,pValCutoff=pValCutoff,
                                path = nwSubdir, runningNumber=i, nwInstanceNb = nwInstanceNb)
    }
}

startt = Sys.time()
tmp = mclapply(X = parameterSets, create1randomNW, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt
# 1 min with 10 cpus for 2 nw instances, 2 trained networks

# n.b. already parallelized, don't over-parallelize, and uses 1.5% of RAM
predict1nw = function(parameterSet) {
    system(paste("echo 'calculating set:",parameterSet$runningNumber,"'"))
    
    geneExFile = parameterSet$dataSetExpr
    geneMthFile = parameterSet$dataSetMethyl
    system(paste0("echo 'using: ",dataSetsPath,geneExFile," and ",dataSetsPath,geneMthFile,"'"))

    if(! file.exists(paste0(dataSetsPath,geneExFile))) {
        system(paste0("echo 'does not exist: ",dataSetsPath,geneExFile,"'"))
    }
    if(! file.exists(paste0(dataSetsPath,geneMthFile))) {
        system(paste("echo 'does not exist:",dataSetsPath,geneMthFile,"'"))
    }
    testDataSet = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = geneExFile,
        geneCopyNumberFile = geneMthFile, path = dataSetsPath )
    predictGeneExpression(data = testDataSet, dataSetName = parameterSet$dataSetName, 
                          networkName = parameterSet$networkName, pValCutoff = parameterSet$pValCutoff, 
                          localGeneCutoff = parameterSet$localGeneCutoff, path = parameterSet$trainNWpath )
}

# create parameter sets for prediction with trained networks
createNormalPredictionSets = function(networkName,dataSetName,
                                pValCutoffs,localGeneCutoffs,nwRange) {
    parameterSets = NULL
    i=0
    for (pValCutoff in pValCutoffs) {
        for (localGeneCutoff in localGeneCutoffs) {
            for (networkNum in nwRange) {
                i = i + 1
                geneExFile = paste0(dataSetName, "_ExpressionData_regNet_Run_",networkNum,".txt")
                geneMthFile = paste0(dataSetName , "_MethylationData_regNet_Run_",networkNum,".txt")
                parameterSets[[i]] = list(dataSetExpr = geneExFile, dataSetMethyl = geneMthFile,
                                          dataSetName = dataSetName, networkName=networkName,
                                          pValCutoff=pValCutoff,localGeneCutoff = localGeneCutoff, 
                                         trainNWpath = paste0("TrainNetwork-",networkNum),runningNumber=i)
            }
        }

    }
    return(parameterSets)
}

# create parameter sets for prediction with random networks
## for each of the trained networks there will be several random networks
createRandomPredictionSets = function(networkName,randomNameBase,dataSetName,
                                pValCutoffs,localGeneCutoffs,numNWinstanceRange,nwRange) {
    parameterSets = NULL
    i=0
    for (nwInstanceNb in numNWinstanceRange ) {
        for (pValCutoff in pValCutoffs) {
            for (localGeneCutoff in localGeneCutoffs) {
                for (networkNum in nwRange) {
                    i = i + 1
                    randomNetworkName = paste0(randomNameBase,nwInstanceNb,"_PValueCutoff_",pValCutoff,
                                               "_BasedOn_",networkName)
                    geneExFile = paste0(dataSetName, "_ExpressionData_regNet_Run_",networkNum,".txt")
                    geneMthFile = paste0(dataSetName , "_MethylationData_regNet_Run_",networkNum,".txt")
                    parameterSets[[i]] = list(dataSetExpr = geneExFile, dataSetMethyl = geneMthFile,
                                              dataSetName = dataSetName, networkName=randomNetworkName,
                                              pValCutoff=pValCutoff,localGeneCutoff = localGeneCutoff, 
                                             trainNWpath = paste0("TrainNetwork-",networkNum),runningNumber=i)
                }
            }
           
        }
    }
    return(parameterSets)
}

# parameter sets for prediction of test data with normal networks:
parameterSetsNormal = 
    createNormalPredictionSets(networkName = networkName,dataSetName = "TestSet", 
                               pValCutoffs = 0.001, localGeneCutoffs = 30, nwRange = 1:numOfNWs)

# parameter sets for prediction of random nteworks:
networkNameBase = "RandomNetwork_" 
parameterSetsRandom = 
    createRandomPredictionSets(networkName = networkName, randomNameBase = "RandomNetwork_",
                         dataSetName = "TestSet", pValCutoffs =  0.001,
                         localGeneCutoffs = 30, numNWinstanceRange = 1:numRandNWinstances, nwRange = 1:numOfNWs)

# parameter sets for prediction of MelBrain Data with TCGA trained:
parameterSetsMelBrainSys = 
    createNormalPredictionSets(networkName = networkName,
                         dataSetName = "MelBrainSys", pValCutoffs = 0.001, 
                         localGeneCutoffs = 30, nwRange = 1:numOfNWs)

# rename melBrainSys data sets, they don't follow the TrainNetwok number:
parameterSetsMelBrainSys[[1]]$dataSetExpr
parameterSetsMelBrainSys[[1]]$dataSetMethyl

for(i in 1:length(parameterSetsMelBrainSys)) {
    parameterSetsMelBrainSys[[i]]$dataSetExpr = "MelBrainSys-expression.csv"
    parameterSetsMelBrainSys[[i]]$dataSetMethyl = "MelBrainSys-methylation.csv"
}

# predict test data with trained network:
startt = Sys.time()
tmp = mclapply(X = parameterSetsNormal, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt
# 35 sec with 19 cpus for 2 data sets (2 NWs each 1 NW instances)

# predict test data with random networks:
startt = Sys.time()
tmp = mclapply(X = parameterSetsRandom, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt
# 39 sec with 19 cpus for 2 data sets (2 NWs each 2 NW instances)

# predict MelBrainSys data with trained networks:
startt = Sys.time()
tmp = mclapply(X = parameterSetsMelBrainSys, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt
# 24 sec with 19 cpus for 2 data sets (2 NWs each 1 NW instances)

# parameter sets for prediction of MelBrain BRAIN Data with TCGA trained:
parameterSetsMelBrainSysBrain = 
    createNormalPredictionSets(networkName = networkName,
                         dataSetName = "MelBrainSysBrain", pValCutoffs = 0.001, 
                         localGeneCutoffs = 30, nwRange = 1:numOfNWs)
# rename data sets:
for(i in 1:length(parameterSetsMelBrainSysBrain)) {
    parameterSetsMelBrainSysBrain[[i]]$dataSetExpr = "MelBrainSysBrainSamples-expression.csv"
    parameterSetsMelBrainSysBrain[[i]]$dataSetMethyl = "MelBrainSysBrainSamples-methylation.csv"
}

# parameter sets for prediction of MelBrain EXTRACRANIAL Data with TCGA trained:
parameterSetsMelBrainSysECM = 
    createNormalPredictionSets(networkName = networkName,
                         dataSetName = "MelBrainSysECM", pValCutoffs = 0.001, 
                         localGeneCutoffs = 30, nwRange = 1:numOfNWs)
# rename data sets:
for(i in 1:length(parameterSetsMelBrainSysECM)) {
    parameterSetsMelBrainSysECM[[i]]$dataSetExpr = "MelBrainSysNonBrainSamples-expression.csv"
    parameterSetsMelBrainSysECM[[i]]$dataSetMethyl = "MelBrainSysNonBrainSamples-methylation.csv"
}

# predict MelBrainSys Brain sample data with trained networks:
startt = Sys.time()
tmp = mclapply(X = parameterSetsMelBrainSysBrain, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt

# predict MelBrainSys extracranial sample data with trained networks:
startt = Sys.time()
tmp = mclapply(X = parameterSetsMelBrainSysECM, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt

# parameter sets for prediction of MelBrain Validation:
parameterSetsMelBrainSysValidation = 
    createNormalPredictionSets(networkName = networkName, 
                         dataSetName = "MelBrainSysValidation", pValCutoffs = 0.001, 
                         localGeneCutoffs = 30, nwRange = 1:numOfNWs)
# rename validation data sets:
for(i in 1:length(parameterSetsMelBrainSysValidation)) {
    parameterSetsMelBrainSysValidation[[i]]$dataSetExpr = "MelBrainSys-validation-expression.csv"
    parameterSetsMelBrainSysValidation[[i]]$dataSetMethyl = "MelBrainSys-validation-methylation.csv"
}

# predict MelBrainSys validation sample data:
startt = Sys.time()
tmp = mclapply(X = parameterSetsMelBrainSysValidation, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt

# parameter sets for prediction of MelBrain Validation BRAIN Data with TCGA trained:
parameterSetsMelBrainSysValidationBrain = 
    createNormalPredictionSets(networkName = networkName, 
                         dataSetName = "MelBrainSysValidationBrain", pValCutoffs = 0.001, 
                         localGeneCutoffs = 30, nwRange = 1:numOfNWs)
# rename validation data sets:
for(i in 1:length(parameterSetsMelBrainSysValidationBrain)) {
    parameterSetsMelBrainSysValidationBrain[[i]]$dataSetExpr = "MelBrainSysBrain-validation-expression.csv"
    parameterSetsMelBrainSysValidationBrain[[i]]$dataSetMethyl = "MelBrainSysBrain-validation-methylation.csv"
}

# predict MelBrainSys BRAIN validation sample data with trained networks:
startt = Sys.time()
tmp = mclapply(X = parameterSetsMelBrainSysValidationBrain, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt

# parameter sets for prediction of MelBrain Validation EXTRACRANIAL Data with TCGA trained:
parameterSetsMelBrainSysValidationECM = 
    createNormalPredictionSets(networkName = networkName, 
                         dataSetName = "MelBrainSysValidationECM", pValCutoffs = 0.001, 
                         localGeneCutoffs = 30, nwRange = 1:numOfNWs)
# rename validation data sets:
for(i in 1:length(parameterSetsMelBrainSysValidationECM)) {
    parameterSetsMelBrainSysValidationECM[[i]]$dataSetExpr = "MelBrainSysECM-validation-expression.csv"
    parameterSetsMelBrainSysValidationECM[[i]]$dataSetMethyl = "MelBrainSysECM-validation-methylation.csv"
}

# predict MelBrainSys EXTRACRANIAL validation sample data with trained networks:
startt = Sys.time()
tmp = mclapply(X = parameterSetsMelBrainSysValidationECM, predict1nw, mc.cores = nbCPUs)
endt = Sys.time()
endt - startt

readInNormalNWPredictions = function(dataSetName,networkName, pValCutoffs,
                             localGeneCutoffs,nwRange) {
    # first build list of files and check which files are there:
    dataSets = NULL
    for(localGeneCutoff in localGeneCutoffs) {
        message("localGeneCutoff ",localGeneCutoff,"\n")
        for (pValCutoff in pValCutoffs) {
            message("reading files of pValCutoff ",pValCutoff,"...\n")
            for (networkNum in nwRange) {
                resultFile = 
                    paste0("TrainNetwork-",networkNum,"/NetworkPredictions/",dataSetName,
                           "_PredictionOfGeneExpressionBasedOn_",
                           networkName,"_PValueCutoff_",
                           pValCutoff,"_LocalGeneCutoff_",localGeneCutoff,".txt")
                dataSets = rbind(dataSets,
                                 data.frame(file=resultFile,networkNum=networkNum,pValCutoff=pValCutoff,
                                           localGeneCutoff=localGeneCutoff, 
                                         stringsAsFactors = F))
            }
        }
    }
    message("class of file ",class(dataSets$file))
    wh = which(! file.exists(dataSets$file))
    if (length(wh)>0) {
        message( length(wh), " files do not exist: \n", paste0(dataSets$file[wh],collapse = "\n"))
        return(NULL)
    }
    predResults=NULL
    for(i in 1:nrow(dataSets)) {
        message("reading in: ",dataSets$file[i],"\n")
        dat = read.csv(file = dataSets$file[i],header = T,sep = "\t",stringsAsFactors = F)
        predResults = 
            rbind(predResults,data.frame(localGeneCutoff=dataSets$localGeneCutoff[i],
                                         pValCutoff=dataSets$pValCutoff[i],
                                         networkNum=dataSets$networkNum[i],dat))
    }
    return(predResults)
}

readInRandomNWPredictions = function(dataSetName,networkNameBase,networkName, pValCutoffs,
                             localGeneCutoffs,NWinstances,nwRange) {
    # first build list of files and check which files are there:
    dataSets = NULL
    for(localGeneCutoff in localGeneCutoffs) {
        message("localGeneCutoff ",localGeneCutoff,"\n")
        for (pValCutoff in pValCutoffs) {
            message("  pValCutoff ",pValCutoff,"...\n")
            for(nwInstanceNb in NWinstances) {
                for (networkNum in nwRange) {
                    resultFile = 
                        paste0("TrainNetwork-",networkNum,"/NetworkPredictions/",dataSetName,
                               "_PredictionOfGeneExpressionBasedOn_",
                               networkNameBase,nwInstanceNb,"_PValueCutoff_",pValCutoff,
                               "_BasedOn_",networkName,"_PValueCutoff_",
                               pValCutoff,"_LocalGeneCutoff_",localGeneCutoff,".txt")
                    dataSets = 
                        rbind(dataSets,
                              data.frame(file=resultFile,networkNum=networkNum,nwInstanceNb=nwInstanceNb,
                                         pValCutoff=pValCutoff,localGeneCutoff=localGeneCutoff, stringsAsFactors = F))
                }
            }
        }
    }
    message("class of file ",class(dataSets))

    wh = which(! file.exists(dataSets$file))
    if (length(wh)>0) {
        message( length(wh), " files do not exist: \n", paste0(dataSets$file[wh],collapse = "\n"))
        return(NULL)
    }
    predResults=NULL
    for(i in 1:nrow(dataSets)) {
        message(class(dataSets$file[i]))
        message("reading in: ",dataSets$file[i],"\n")
        dat = read.csv(file = dataSets$file[i],header = T,sep = "\t",stringsAsFactors = F)
        predResults = 
            rbind(predResults,
                  data.frame(localGeneCutoff = dataSets$localGeneCutoff[i],pValCutoff = dataSets$pValCutoff[i],
                             networkNum = dataSets$networkNum[i],nwInstanceNb = dataSets$nwInstanceNb[i],dat))
    }
    return(predResults)
}

normalPredTestResults = 
    readInNormalNWPredictions(dataSetName = "TestSet",networkName = networkName,
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

melBrainPredResults = 
    readInNormalNWPredictions(dataSetName = "MelBrainSys",networkName = networkName, 
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

randomPredTestResults = 
    readInRandomNWPredictions(dataSetName = "TestSet",networkNameBase = "RandomNetwork_",
                      networkName = networkName, pValCutoffs = "0.001",
                      localGeneCutoffs = 30, nwRange=1:numOfNWs , NWinstances = 1:numRandNWinstances)

melBrainMBMPredResults = 
    readInNormalNWPredictions(dataSetName = "MelBrainSysBrain",networkName = networkName, 
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

melBrainECMPredResults = 
    readInNormalNWPredictions(dataSetName = "MelBrainSysECM",networkName = networkName, 
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

melBrainValidationPredResults = 
    readInNormalNWPredictions(dataSetName = "MelBrainSysValidation",networkName = networkName, 
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

melBrainValidationMBMPredResults = 
    readInNormalNWPredictions(dataSetName = "MelBrainSysValidationBrain",networkName = networkName, 
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

melBrainValidationECMPredResults = 
    readInNormalNWPredictions(dataSetName = "MelBrainSysValidationECM",networkName = networkName, 
                              pValCutoffs = "0.001", localGeneCutoffs = 30, nwRange = 1:numOfNWs)

head(melBrainPredResults,3)
head(normalPredTestResults,3)
head(randomPredTestResults,3)
head(melBrainMBMPredResults,3)
head(melBrainECMPredResults,3)
head(melBrainValidationPredResults,3)
head(melBrainValidationMBMPredResults,3)
head(melBrainValidationECMPredResults,3)

# now get median for each gene in each random network (nwInstanceNb) of networkNum networks
allGenes = unique(sort(randomPredTestResults$Gene))
length(allGenes)
medianrandomPredTestResults = NULL
for (networkNum in 1:numOfNWs) {
    dat = randomPredTestResults[randomPredTestResults$localGeneCutoff == localGeneCutoff & 
                      randomPredTestResults$pValCutoff == pValCutoff & 
                                randomPredTestResults$networkNum == networkNum,]
    median_correlation = sapply(allGenes, function(g) median(dat[ dat$Gene == g,"Correlation" ]))
    median_pvalue = sapply(allGenes, function(g) median(dat[ dat$Gene == g,"P.Value" ]))
    medianrandomPredTestResults = 
        rbind(medianrandomPredTestResults,
              data.frame(gene=allGenes,localGeneCutoff=localGeneCutoff,pValCutoff=pValCutoff,
                         networkNum=networkNum,median_correlation=median_correlation, 
                         median_pvalue=median_pvalue))
}
# takes ca. 3 minutes, maybe we could programm it better or parallelize

head(medianrandomPredTestResults,3)

# put all results into 1 data.frame for ggplot
mergedPredResultsOwnNW = 
    data.frame(stringsAsFactors = F, type="Random networks",
               medianrandomPredTestResults[,c("localGeneCutoff","pValCutoff","networkNum","gene","median_correlation","median_pvalue")])

colnames(mergedPredResultsOwnNW) = c("type","localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")
tail(mergedPredResultsOwnNW,3)
mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
          data.frame(stringsAsFactors = F, type="Discovery cohort",
                     melBrainPredResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
          data.frame(stringsAsFactors = F, type="TCGA test data",
                     normalPredTestResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
      data.frame(stringsAsFactors = F, type="Validation cohort",
         melBrainValidationPredResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
      data.frame(stringsAsFactors = F, type="Discovery cohort extracranial",
         melBrainECMPredResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
      data.frame(stringsAsFactors = F, type="Discovery cohort brain",
         melBrainMBMPredResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
      data.frame(stringsAsFactors = F, type="Validation cohort brain",
         melBrainValidationMBMPredResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

mergedPredResultsOwnNW = 
    rbind(mergedPredResultsOwnNW,
      data.frame(stringsAsFactors = F, type="Validation cohort extracranial",
         melBrainValidationECMPredResults[,c("localGeneCutoff","pValCutoff","networkNum","Gene","Correlation", "P.Value")]))
tail(mergedPredResultsOwnNW,3)

# only results without NA values
mergedPredResultsOwnNWNoNA = mergedPredResultsOwnNW[which(!is.na(mergedPredResultsOwnNW$Correlation)),]
nrow(mergedPredResultsOwnNWNoNA); nrow(mergedPredResultsOwnNW)
table(mergedPredResultsOwnNWNoNA$type)

# median correlation and p-value for each:
sapply(unique(mergedPredResultsOwnNWNoNA$type), function(type) 
    c(cor = median(mergedPredResultsOwnNWNoNA[ mergedPredResultsOwnNWNoNA$type == type,"Correlation"]),
    pval = median(mergedPredResultsOwnNWNoNA[ mergedPredResultsOwnNWNoNA$type == type,"P.Value"]))
       )

# now calculate mean correlations over all trained networks:
mergedPredResultsOwnNWNoNAmean = NULL
for(type in unique(mergedPredResultsOwnNWNoNA$type)) {
    message(type)
    tmp = mergedPredResultsOwnNWNoNA[ mergedPredResultsOwnNWNoNA$pValCutoff == pValCutoff & 
                                     mergedPredResultsOwnNWNoNA$localGeneCutoff == localGeneCutoff &
                      mergedPredResultsOwnNWNoNA$type == type,c("Gene","Correlation","networkNum")]
    tmpWide = reshape(tmp, idvar = "Gene", timevar = "networkNum", direction = "wide")
    meanCor = rowMeans(tmpWide[,-1], na.rm = T)
    mergedPredResultsOwnNWNoNAmean = 
        rbind(mergedPredResultsOwnNWNoNAmean,
                  data.frame(type=type,localGeneCutoff=localGeneCutoff,
                             pValCutoff=pValCutoff,Gene=names(meanCor),Correlation=meanCor))
}

correlRange = c(-1,1)

table(mergedPredResultsOwnNWNoNAmean$type)

### random NW + TCGA test + MelBrain discovery + validation

wh = which(mergedPredResultsOwnNWNoNAmean$type %in% 
           c("Random networks","TCGA test data","Discovery cohort",
                                "Validation cohort"))
p1 = ggplot(data = mergedPredResultsOwnNWNoNAmean[wh,],
           mapping = aes(x=Correlation, fill=type)) + 
    geom_histogram(position="identity",alpha=0.4, bins = 100) + 
    scale_fill_manual(breaks =c("Random networks","TCGA test data","Discovery cohort",
                                "Validation cohort"),
                      values = c("#666666","aquamarine4","blue",
                                "darkorchid1")) +
    scale_y_sqrt(breaks=c(10,100,250,500,1000,2000), name="Number of genes", limits=c(0,3000))+
    scale_x_continuous( name="Correlation: predicted vs. original expression", limits=c(-1,1)) + 
    #theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) + # axis labels perpendicular
    theme(legend.title=element_blank(), legend.position=c(.8,.8),legend.text=element_text(size=7),
          legend.spacing.x = unit(0.25, 'cm')) # more space between legend label and colored box
p1

# now plot mean correlation between predicted and real gene expression

# discovery brain + ecm + validation brain + ecm

wh = which(mergedPredResultsOwnNWNoNAmean$type %in% 
           c("Discovery cohort brain", "Discovery cohort extracranial",
                               "Validation cohort brain", "Validation cohort extracranial"))

p2 = ggplot(data = mergedPredResultsOwnNWNoNAmean[wh,],
           mapping = aes(x=Correlation, fill=type)) + 
    geom_histogram(position="identity",alpha=0.4, bins = 100) + 
    scale_fill_manual(breaks = c("Discovery cohort brain", "Discovery cohort extracranial",
                                "Validation cohort brain", "Validation cohort extracranial"),
                      values = c("brown2", "cornflowerblue","gold2" , "aquamarine4")) +
    scale_y_sqrt(breaks=c(10,100,250,500,1000,2000), name="Number of genes", limits=c(0,3000))+
    scale_x_continuous( name="Correlation: predicted vs. original expression", limits=c(-1,1)) + 
    #labs(tag = "A") + theme(plot.tag.position = c(-0.5, 2000))  + 
    #theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) + # axis labels perpendicular
    theme(legend.title=element_blank(), legend.position=c(.75,.8),legend.text=element_text(size=7),
          legend.spacing.x = unit(0.25, 'cm')) # more space between legend label and colored box
p2

allplots = grid.arrange(p1,p2,nrow=1)

png(filename = paste0(outDirectory,"Figure-2-gene-expr-prediction-correl.png"),
    width = 4000,height = 1600, res=370)
plot(allplots)
dev.off()
