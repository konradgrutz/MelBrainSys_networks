myPath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"

setwd(myPath)

install.packages("bin/covTest/", repos = NULL, type="source")

library(covTest)

library(devtools) # for function install()

.libPaths()

regNetParentDir = "bin/regNet-master/" 
setwd(regNetParentDir)
localRlibs = paste0(myPath,"conda/lib/R/library/")
install( pkg = "regNet", args = c( paste0("--library=",localRlibs ) ) )

getwd()
