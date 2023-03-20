basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

neededSamples = c("P03_Br","P03_Lu", "P04_Br", "P04_Sk_GA", "P08_Br", "P08_St_GA", "P08_St_BA", "P08_St_YA", 
                     "P16_Br", "P16_Lu", "P18_Br", "P18_Lu_GA", "P18_Lu_YA", "P39_Br", "P39_Lu", 
                     "P42_Br_GA", "P42_Ly_GA", "P42_Ly_YA")

# mapping of patient samples to metastases pairs 
sampleMapping = list("P03_BLun" = c("P03_Br","P03_Lu"), 
                  "P04_BSki_1" = c("P04_Br", "P04_Sk_GA"), 
                  "P08_BSof_1" = c("P08_Br", "P08_St_GA"), 
                  "P08_BSof_2" = c("P08_Br", "P08_St_BA"), 
                  "P08_BSof_3" = c("P08_Br", "P08_St_YA"), 
                  "P16_BLun" = c("P16_Br", "P16_Lu"), 
                  "P18_BLun_1" = c("P18_Br", "P18_Lu_GA"), 
                  "P18_BLun_2" = c("P18_Br", "P18_Lu_YA"), 
                  "P39_BLun" = c( "P39_Br", "P39_Lu"), 
                  "P42_BLym_1" = c("P42_Br_GA", "P42_Ly_GA"), 
                  "P42_BLym_2" = c("P42_Br_GA", "P42_Ly_YA"))
samplePairs = names(sampleMapping)

patient_colors = c(P03 = '#0dbfbf', P04= '#ff6db6ff', P08 = '#006ddbff',
                   P16 = '#c2e0fcff', P18 = '#db6d00ff', P39 = '#25f923ff',
                   P42 = '#fffc6dff')

tissue_colors = c(BLun = '#00bb00ff', BLym = '#bb2cbbff', 
                  BSki = '#999999ff', BSof = '#e5e5e5ff')

# sample pair cluster colors
colClusterHi = "tomato1"
colClusterLow = "cornflowerblue"
colClusterSlightLow = "lightskyblue1"
samplePairPerSubgroup = list("SG1"=c('P04_BSki_1','P08_BSof_2','P16_BLun','P42_BLym_1'), 
                     "SG3"=c('P03_BLun','P08_BSof_3','P18_BLun_1','P42_BLym_2'),
                     "SG2" = c('P08_BSof_1','P18_BLun_2','P39_BLun'))
subgroupPerSamplePair = c('SG1','SG1','SG1','SG1','SG3','SG3','SG3','SG3','SG2','SG2','SG2')
names(subgroupPerSamplePair) = c('P04_BSki_1','P08_BSof_2','P16_BLun','P42_BLym_1','P03_BLun','P08_BSof_3','P18_BLun_1','P42_BLym_2','P08_BSof_1','P18_BLun_2','P39_BLun')
subgroupNames = c("SG1"="higher","SG2"="lower","SG3"="slightly lower")

subgroup_colors = c("SG1"=colClusterHi, "SG2"=colClusterLow, "SG3"=colClusterSlightLow)

# plot colors for patients tissues
colAnnot <- data.frame(patient = gsub("^(P\\d+).+","\\1",samplePairs), 
                       tissue=gsub("^P\\d+_([A-Za-z]+).*","\\1",samplePairs))
row.names(colAnnot)=samplePairs
colAnnot$subgroup = subgroupPerSamplePair[ rownames(colAnnot)]
colAnnot
annotColors = list(patient = patient_colors, tissue = tissue_colors, 
                   subgroup = subgroup_colors)
annotColors

saveRDS(object = list(
    neededSamples = neededSamples,
    sampleMapping = sampleMapping,
    samplePairs = samplePairs,
    samplePairPerSubgroup = samplePairPerSubgroup,
    subgroupPerSamplePair = subgroupPerSamplePair,
    subgroupNames = subgroupNames,
    colAnnot = colAnnot,
    annotColors = annotColors),
        file = "annotation/samplePairs-annotation-colors-clusters.rds")
