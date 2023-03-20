# make overrepresentation analyses
# create
#  SupplTable-S8-overrepresentation-analysis.xls

library(org.Hs.eg.db) # gene annotation, for mapping symbols to entrez/ensembl
library(clusterProfiler) # GO enrichment analysis
library(xlsx)
library(parallel)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/")
setwd(basePath)

outDirectory = paste0(basePath,"FiguresTables/")

# define met. pairs, and annotation colors 
tmp = readRDS(file = "annotation/samplePairs-annotation-colors-clusters.rds")
samplePairPerSubgroup = tmp$samplePairPerSubgroup
subgroupNames = tmp$subgroupNames

tmp = readRDS(file="metPairs-impactRatios-onAllTargetGenes.rds")
impRatios = tmp$meanLogMedianImpRatiosMat
head(impRatios)

allGenesNWs = rownames(impRatios)

# mapping all genes to Entrez
# "ENTREZID"
allGenesNWsMapped <- bitr(allGenesNWs, fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Hs.eg.db)
# 5.2% failed mapping (when only mapping on ENSEMBL/ENTREZ it is 5.16 and 5.18%)
wh = which(duplicated(allGenesNWsMapped$SYMBOL))
length(wh) # 324 duplicated, remove the duplicates
allGenesNWsMapped = allGenesNWsMapped[-wh,]
rownames(allGenesNWsMapped) = allGenesNWsMapped$SYMBOL
head(allGenesNWsMapped,2)

# GO ID / descr mapping:
goIDdescrMap = function(res) {
    tmp = res[,c("ID","Description")]
    wh = which(duplicated(tmp$ID))
    if(length(wh)>0) {tmp = tmp[-wh,]}
    goIDdescr = tmp$Description; names(goIDdescr) = tmp$ID
    return(goIDdescr)
}

calcOverrepr1set = function( geneSet, allGenesSet, maxQval) {
    geneSetEnsg = na.omit(allGenesSet[geneSet, "ENSEMBL"])
    geneSetEntr = na.omit(allGenesSet[geneSet, "ENTREZID"]) 
    res = NULL
    for (goCateg in c("MF","CC","BP")) {
        go = enrichGO(gene = geneSetEnsg, universe = allGenesSet$ENSEMBL,
                    OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = goCateg, 
                    pAdjustMethod = "BH", pvalueCutoff  = 1, qvalueCutoff  = maxQval, readable = T)
        if(class(go)!="NULL" & nrow(go)>0) {
            res = rbind(res,data.frame(stringsAsFactors=F, categ = goCateg,num_genes=length(geneSetEnsg),go))
        }
    }
    kegg = enrichKEGG(gene = geneSetEntr, universe = allGenesSet$ENTREZID, qvalueCutoff = maxQval,
             keyType = 'ncbi-geneid', organism = 'hsa',pvalueCutoff = 1)
    if(class(kegg)!="NULL" & nrow(kegg)>0) {
        res = rbind(res,data.frame(stringsAsFactors=F, categ = "KEGG" ,num_genes=length(geneSetEntr) ,kegg))
    }
    res
}

calcOverReprWcutoff = function(impRatios, decreasing, nbCPUs, percentile) {
    GOres = mclapply(mc.cores = nbCPUs,X = colnames(impRatios), FUN = function(samplePair) {
            topGenes = head(sort(impRatios[,samplePair], decreasing = decreasing),
                            round(percentile*nrow(allGenesNWsMapped)/100))
            if(decreasing) { 
                topGenes = topGenes[ which(topGenes > 0)] 
            } else {
                topGenes = topGenes[ which(topGenes < 0)] 
            }
            message(samplePair," ",length(topGenes)," topGenes")
            if (length(topGenes)>0) {
                topGenes = names(topGenes)
                tmp = calcOverrepr1set(geneSet = topGenes, allGenesSet = allGenesNWsMapped, maxQval = 0.1)
                if(class(tmp)!="NULL") {
                    return(data.frame(stringsAsFactors = F,samplePair = samplePair,tmp))
                } else { return(NULL) }
            } else { return(NULL) }
        })
    do.call(rbind,GOres)
}

set.seed(seed = 42) # reproducible results

nbCPUs = 1 # parallelization doesn't work, enrichment functions break with error "database disk image is malformed"
allRevMeanGOres = calcOverReprWcutoff(impRatios = impRatios, decreasing = F,
                                      nbCPUs = nbCPUs,percentile = 5 )
# takes 15-30min

allFwdMeanGOres = calcOverReprWcutoff(impRatios = impRatios, decreasing = T,
                                      nbCPUs = nbCPUs,percentile = 5 )
# takes 15-30min

head(allFwdMeanGOres,2)
tail(allFwdMeanGOres,2)

# save for paper
write.xlsx(allFwdMeanGOres, sheetName="upper 5%",append=F, row.names = F,
           file = paste0(outDirectory,"SupplTable-S8-overrepresentation-analysis.xls"))
# save for paper
write.xlsx(allRevMeanGOres, sheetName="lower 5%", append=T,row.names = F,
           file = paste0(outDirectory,"SupplTable-S8-overrepresentation-analysis.xls"))

goIDdescr = goIDdescrMap(res = rbind(allFwdMeanGOres, allRevMeanGOres))

FwdAndRevResults = list("upper 5%" = allFwdMeanGOres, 
                     "lower 5%"=allRevMeanGOres)

sharedGOs = NULL # contains all results that are in >= 2 samples
anyGOs = 
    list("upper 5%" = NULL, "lower 5%"=NULL) 
    # anyGOs: for each list (upper/lower), and each SG, which GO terms are overrepresented
    # needed later to get exclusive GO terms
for(geneList in names(FwdAndRevResults)) {
    message(geneList)
    goRes = FwdAndRevResults[[geneList]]
    for(SG in names(subgroupNames)) {
        SPs = samplePairPerSubgroup[[SG]]
        anyGOs[[geneList]][[SG]] = 
            unique(sort(goRes$ID[ goRes$samplePair %in% SPs]))
        frqGO = sort(table(goRes$ID[ goRes$samplePair %in% SPs]), decreasing = T)
        frqGO = frqGO[ frqGO > 1]
        message("subgroup ",paste0(SG,": ",paste0(SPs, collapse=" ")),
                ", ",length(frqGO)," shared categories")
        if(length(frqGO)>0) {
            tmp = data.frame(stringsAsFactors = F, gene_list= geneList, 
                        subgroup = SG, num_sample_pairs = length(SPs),
                        num_sample_pairs_sharing = as.numeric(frqGO),
                             GO_ID = names(frqGO),description = goIDdescr[names(frqGO)])
            out = ""
            for(i in 1:min(10,nrow(tmp))) { out = paste0(out,paste0(paste0(" | ",tmp[i,],collapse=""),collapse=" | ")," |\n") }
            message(out)
            sharedGOs = rbind(sharedGOs,tmp)
        }
    }
}

head(sharedGOs,3)

exclSharedGOs = NULL
for (geneList in unique(sharedGOs$gene_list)) {
    message(geneList)
    for(SG in unique(sharedGOs$subgroup)) {
        otherSGs = setdiff(names(subgroupNames),SG)
        unwantedGOs = unique(unlist(anyGOs[[geneList]][otherSGs]))
        tmp = sharedGOs[ sharedGOs$gene_list==geneList & sharedGOs$subgroup==SG,]
        a = length(unique(tmp$GO_ID))
        tmp = tmp[ !tmp$GO_ID %in% unwantedGOs,]
        b = length(unique(tmp$GO_ID))
        message("  ",SG," -> other SGs: ", paste0(otherSGs, collapse=" "), ", ",
                length(unwantedGOs)," unwanted GOs, ", a, " before ",b," after removal" )
        exclSharedGOs = rbind(exclSharedGOs,tmp)
    }
}

head(exclSharedGOs,3)

table(sharedGOs$subgroup, sharedGOs$gene_list)

table(exclSharedGOs$subgroup, exclSharedGOs$gene_list)

exclSharedGOs$subgroup_ID = exclSharedGOs$subgroup
exclSharedGOs$subgroup = paste(subgroupNames[ exclSharedGOs$subgroup],"in brain")

exclSharedGOs = exclSharedGOs[,c('gene_list','subgroup','subgroup_ID',
                                 'num_sample_pairs','num_sample_pairs_sharing','GO_ID','description')]
head(exclSharedGOs)
table(exclSharedGOs$subgroup, exclSharedGOs$subgroup_ID)

write.xlsx(exclSharedGOs, sheetName="exclusive categories", append=T, row.names = F,
           file = paste0(outDirectory,"SupplTable-S8-overrepresentation-analysis.xls"))
