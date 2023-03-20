# from redmine (end of ticket): https://imbredmine.medizin.tu-dresden.de/redmine/issues/187

# for each of the 1-n networks: 
#  this will create the basic regNet folder structure
#  and train a network with the corresponding training data (expression and methylation files)
#  the results are split up into single parts for each network, which then are combined

myPath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"

library(parallel)

localRlibs = paste0(myPath,"conda/lib/R/library/")
library(regNet)

# we train 25 networks in total
# here you set the range of which networks to train in the current run
startRun = 1
endRun = 2

dataSubPath = "Data/"
regNetPath = paste0(myPath, "regNet/")
output=T
setwd(regNetPath)
loadPath = paste0( regNetPath, dataSubPath )

# get data file directions:
geneExprFiles = paste0("TrainSet_ExpressionData_regNet_Run_",startRun:endRun,".txt")
geneExprFiles
if(all(file.exists(paste0(dataSubPath,geneExprFiles)))) {
    cat("good: all gene expression training set files exist\n")
} else {
    cat("bad: not all gene expression training set files exist, please check:\n",
        paste0(dataSubPath,geneExprFiles,collapse = ", "),"\n")
    q("no")
}

geneMethFiles = paste0("TrainSet_MethylationData_regNet_Run_",startRun:endRun,".txt")
geneMethFiles
if(all(file.exists(paste0(dataSubPath,geneMethFiles)))) {
    cat("good: all gene methylation training set files exist\n")
} else {
    cat("bad: not all gene methylation training set files exist, please check:\n",
        paste0(dataSubPath,geneMethFiles,collapse = ", "),"\n")
    q("no")
}



networkName = "TcgaMelanomaExprMeth"
# set low totalNumberOfJobs 
totalNumberOfJobs = 20 # number of regNetJobs
# count number of jobs/genes in input data:
nbParallelJobs = 10 #  how many jobs in parallel with mclapply
# when totalNumberOfJobs == nbParallelJobs, all regNet jobs are calculated at the same time
# if you have less CPUs and RAM, lower both totalNumberOfJobs and nbParallelJobs

data$loc[1:4]

# start this with nohup and RScript, not in jupyter

starttTotal = Sys.time()

for( idx in 1:length(geneExprFiles)) {
    projectName = paste0("TrainNetwork-",(startRun:endRun)[idx])
    cat("project",projectName,"\n")
    cat("  ",geneExprFiles[idx],"\n")
    cat("  ",geneMethFiles[idx],"\n")
    
    startt = Sys.time()
    projectPath = createBasicFolderStructure( projectName = projectName, path = regNetPath, output = output )
    
    data = loadGeneExpressionAndCopyNumberDataSet( 
        geneExpressionFile = geneExprFiles[idx], 
        geneCopyNumberFile = geneMethFiles[idx], 
        path = loadPath )

    # start parallel calculation of this training set:
    tmp = mclapply(1:totalNumberOfJobs,mc.cores = nbParallelJobs, function(j) 
        learnNetwork_ParallelComputation( data = data, networkName = networkName, cores = totalNumberOfJobs, 
            job = j, path = projectPath, nfolds = 10, cvReplicates = 10, output = output ))

    #save.image(file = paste0("TrainNetwork-",(startRun:endRun)[j],"-image.RData"))
    endt = Sys.time()
    timeDiff = difftime(endt,startt,units = "min")
    message(timeDiff," needed")
}

endtTotal = Sys.time()
cat("all network trainings done. Needed.\n")
endtTotal - starttTotal

# to see the progress, check the subdirectories "TrainNetwork-*" that appear in the regNet directory

myPath = paste0(myPath,"regNet/")
setwd(myPath)

# combine single network parts from network-training.r to whole network, for each of the 1-n networks

startt = Sys.time()
for(nwSubdir in paste0("TrainNetwork-",startRun:endRun)) {
    path = paste0(myPath, nwSubdir)
    message(path,"\n")
    combineSingleJobs(networkName = networkName, cores = totalNumberOfJobs, path = path, output = output)
}
endt = Sys.time()
endt-startt# 1.5 min for 2 Networks
# creates TcgaMelanomaExprMeth_NetworkCreator_CVStatistics.txt and TcgaMelanomaExprMeth_NetworkCreator.Rout
#  in the TrainNetwork-*/NetworkModel/WholeNetwork/  subdirectories
