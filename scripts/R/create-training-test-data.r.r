setwd("/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/")

dataPath = "regNet/Data/"

trainingSetRatio = 3/4 # how many samples used for training, the remainder are used for test data

# full data set gene expression
fullExprDataFile = paste0( dataPath, "TCGA-expression.csv" )
# Gene, Chromsome, Location, Samples
fullExprData = read.delim( file = fullExprDataFile, header = TRUE, sep = "\t", stringsAsFactors = F)
colnames(fullExprData) = gsub("\\.","-",colnames(fullExprData))
head(fullExprData,3)

# full methylation data set
#
fullMethDataFile = paste0( dataPath, "TCGA-methylation.csv" )
#Gene, Chromsome, Location, Samples
fullMethData = read.delim( file = fullMethDataFile, header = TRUE, sep = "\t" )
colnames(fullMethData) = gsub("\\.","-",colnames(fullMethData))
head(fullMethData,3)

totalSamples = 4:ncol( fullExprData )
nbSamplesTotal = length( totalSamples )
nbTrainSamples = round( nbSamplesTotal * trainingSetRatio, dig = 0 )
cat( nbSamplesTotal, "samples in total,", nbTrainSamples ,"for training,",
    nbSamplesTotal-nbTrainSamples,"for testing \n")

createTraining_and_TestDataSets = function( dataSetNb ) {
    cat("data set number",dataSetNb,"\n")
    trainSamples = sample( x = totalSamples, size = nbTrainSamples  )
    testSamples  = setdiff( totalSamples, trainSamples )

    cat("training samples\n", trainSamples ,"\n")
    cat("test samples\n", testSamples  ,"\n")

    # Expression Data:
    # save training data set
    trainDataFile = paste0( dataPath, "TrainSet_ExpressionData_regNet_Run_", dataSetNb, ".txt" )
    write.table( fullExprData[ , c( 1:3, trainSamples ) ], file = trainDataFile, row.names = FALSE, col.names = TRUE, quote = FALSE, dec = ".", sep = "\t" )
    cat("saved to ",trainDataFile,"\n")
    
    # save test data set
    testDataFile = paste0( dataPath, "TestSet_ExpressionData_regNet_Run_", dataSetNb, ".txt" )
    write.table( fullExprData[ , c( 1:3, testSamples ) ], file = testDataFile, row.names = FALSE, col.names = TRUE, quote = FALSE, dec = ".", sep = "\t" )
    cat("saved to ",testDataFile,"\n")
    
    # Methylation Data:
    # save training data set
    trainDataFile = paste0( dataPath, "TrainSet_MethylationData_regNet_Run_", dataSetNb, ".txt" )
    write.table( fullMethData[ , c( 1:3, trainSamples ) ], file = trainDataFile, row.names = FALSE, col.names = TRUE, quote = FALSE, dec = ".", sep = "\t" )
    cat("saved to ",trainDataFile,"\n")
    
    # save test data set
    testDataFile = paste0( dataPath, "TestSet_MethylationData_regNet_Run_", dataSetNb, ".txt" )
    write.table( fullMethData[ , c( 1:3, testSamples ) ], file = testDataFile, row.names = FALSE, col.names = TRUE, quote = FALSE, dec = ".", sep = "\t" )
    cat("saved to ",testDataFile,"\n")
}

# create data set 1..10
nbDataSets = 10
for( i in 1:nbDataSets ) {
    createTraining_and_TestDataSets( dataSetNb = i)
}










