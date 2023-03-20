library(parallel)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/validation-cohort/")

localRlibs = paste0(basePath,"conda/lib/R/library/")
library(regNet)

setwd(myPath)

localGeneCutoff = 30
pValCutoff = 0.001
colSumsThreshold = 0.001
numOfNWs = 2
networkName = "TcgaMelanomaExprMeth"

# create standard folder structure for regNet
projectNames = paste0("TrainNetwork-",1:numOfNWs,"/")
projectNames

# Create basic folder structure
path = "/data/bcu-projects/gruetzmko/MelBrainSys/regNet/new-sample-pairs/"
for(projectName in projectNames) {
    createBasicFolderStructure( projectName = projectName, path = myPath,output = F)    
}

samplePairs = list(P78_BIn = c('P78_Br_1a','P78_In_3d'), P77_BLy = c('P77_Br_e','P77_Ly_1d'),
                  P74_BLy = c('P74_Br_1a','P74_Ly_1a'), P13_BLy = c('P13_Br_2','P13_Ly_e'),
                  P111_BLy = c('P111_Br_b','P111_Ly_1b'), P108_B_Ly = c('P108_Br_a','P108_Ly_2a'),
                  P101_BLi = c('P101_Br_a','P101_Li_b'),P106_BLy = c('P106_Br_ns1','P106_Ly_1b'),
                   P107_BLu = c('P107_Br_c','P107_Lu_a'))
usedSamples = unname(unlist(samplePairs))
usedSamples

myFullDataSet = loadGeneExpressionAndCopyNumberDataSet(
    geneExpressionFile = "MelBrainSys_ExpressionData_2022_allNeededGenes_regNet.txt", 
    geneCopyNumberFile = "MelBrainSys_MethylationData_2022_allNeededGenes_regNet.txt", path = "Data/" )

nwFlow1sample = function(sample, logFile) {
    starttotal = Sys.time()
    for (netwNum in 1:numOfNWs) {
        startt = Sys.time()
        nwSubdir = paste0("TrainNetwork-",netwNum)
        f = file(logFile,open = "a")
        writeLines(text = paste0(sample," ",netwNum," started"), con = f)
        close(f)
        
        computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts(
            patient = sample, data = myFullDataSet, dataSetName = "MelBrainSys", 
            networkName = networkName, pValCutoff = pValCutoff, 
            localGeneCutoff = localGeneCutoff, path = nwSubdir)
        computeNetworkFlow_PatientSpecificAbsoluteImpacts(
            patient = sample, dataSetName = "MelBrainSys", 
            networkName = networkName, pValCutoff = pValCutoff, 
            localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, 
            path = nwSubdir, output=T)
        
        endt = Sys.time()
        f = file(logFile,open = "a")
        writeLines(text = paste0(
            "  ",as.character(round(difftime(endt, startt, units = "min"),1)),
            " min needed for NW ",netwNum," of sample ",sample),
                   con = f)
        close(f)
        
    }
    endt = Sys.time()
    f = file(logFile,open = "a")
    writeLines(text = paste0(sample," DONE  ",
                            as.character(round(difftime(endt, starttotal, units = "min"))),
                            "needed in total"), con = f)
    close(f)
}

nbCPUs = 3 # already parallelized internally, don't overdo
logFile = "nwFlowParallel.log"
if(file.exists(logFile)) { file.remove(logFile)}

startt = Sys.time()
log_mclapply = mclapply(X = usedSamples, mc.cores = nbCPUs, FUN = function(sample) nwFlow1sample(sample,logFile))
endt = Sys.time()

paste0(as.character(round(difftime(endt, startt, units = "min")))," min needed")
# about 90 min for 10 samples, 2 networks, 160 CPUs in total


