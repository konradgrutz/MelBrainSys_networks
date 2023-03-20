# computes network flow/impact matrix of each gene on each other gene for each sample
# outputs in TrainNetwork-*/NetworkPropagation/NetworkFlow/

library(parallel)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/")

localRlibs = paste0(basePath,"conda/lib/R/library/")
library(regNet)

setwd(myPath)

localGeneCutoff = 30
pValCutoff = 0.001
numOfNWs = 2
networkName = "TcgaMelanomaExprMeth"

# remaining samples (those that are not in a tissue pair):
neededSamples = c("P03_Br","P03_Lu", "P04_Br", "P04_Sk_GA", "P08_Br", "P08_St_GA", "P08_St_BA", "P08_St_YA", 
                     "P16_Br", "P16_Lu", "P18_Br", "P18_Lu_GA", "P18_Lu_YA", "P39_Br", "P39_Lu", 
                     "P42_Br_GA", "P42_Ly_GA", "P42_Ly_YA")
length(neededSamples)

myFullDataSet = loadGeneExpressionAndCopyNumberDataSet(
    geneExpressionFile = "MelBrainSys-expression.csv", 
    geneCopyNumberFile = "MelBrainSys-methylation.csv", path = "Data/" )

# calc. all patient/sample specific flow matrixes for all networks
# 90 min for 10 samples and 2 networks, 160 CPUs in total

localGeneCutoff = 30
colSumsThreshold = 0.001
nbTotalIt = numOfNWs*length(neededSamples)*1
currIt = 0
totalStartt = Sys.time()
for (patient in neededSamples) {
    message(patient)
    for (netwNum in 1:numOfNWs) {
        nwSubdir = paste0("TrainNetwork-",netwNum)
        message("  ",nwSubdir)
        path = paste0(myPath,nwSubdir)
        computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts(
            patient = patient, data = myFullDataSet, dataSetName = "MelBrainSys", 
            networkName = networkName, pValCutoff = pValCutoff, 
            localGeneCutoff = localGeneCutoff, path = path)
        computeNetworkFlow_PatientSpecificAbsoluteImpacts(
            patient = patient, dataSetName = "MelBrainSys", 
            networkName = networkName, pValCutoff = pValCutoff, 
            localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, 
            path = path, output=T)
        endt = Sys.time()
        currIt = currIt + 1
        message("  ",endt-totalStartt," so far needed")
        restTime = round((nbTotalIt-currIt)*difftime(endt,totalStartt, units = "min")/currIt)
        restTimeH = round((nbTotalIt-currIt)*difftime(endt,totalStartt, units = "h")/currIt,1)
        message(restTime," min (= ", restTimeH," hours) still needed")
    }
}
difftime(time1 = endt, time2 = totalStartt) 
