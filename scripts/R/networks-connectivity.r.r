# creates 
#  - connectivity table for paper - SupplTable-S5-connectivity-stats.xls
#  - DotPlot of connections - SupplFigure-S5-source-target-interactions.png
#  - known melanoma driver gene ranking in networks - SupplTable-S6-key-driver-ranks.xls

library(ggplot2)
library(gridExtra)
library(parallel)
library(xlsx)
library(pheatmap)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
myPath = paste0(basePath,"regNet/")

localRlibs = paste0(basePath,"conda/lib/R/library/")
library(regNet)

setwd(myPath)

outDirectory = paste0(basePath,"FiguresTables/")

output=T
projectName = "MyNetwork-1"
networkName = "TcgaMelanomaExprMeth"
numOfNWs = 2

myPath

localGeneCutoff = 30
pValCutoff = 0.001
nbCPUs = 10 

startt = Sys.time()
for(nwSubdir in paste0("TrainNetwork-",1:numOfNWs)) {
    path = paste0(myPath, nwSubdir)
    message("calculating for ",nwSubdir)
    getNetworkConnections( networkName = networkName, pValCutoff = nwSubdir,
                                 localGeneCutoff = localGeneCutoff, path = path )
}
endt = Sys.time()
endt-startt # 5 min for 2 networks

connectFiles = paste0("TrainNetwork-",1:numOfNWs,
                      "/NetworkConnectivity/",networkName,"_GeneSpecificNetworkConnections_PValueCutoff_",pValCutoff,
                      "_LocalGeneCutoff_30.txt")
all(file.exists(connectFiles))

connectStats = data.frame()
for( f in connectFiles) {
    conn = read.csv(file = f, header = T,sep = "\t",stringsAsFactors = F)
    connectStats = 
        rbind(connectStats, 
              data.frame(inRep = sum(conn$IncomingRepressorLinks), 
                         outRep = sum(conn$OutgoingRepressorLinks),
                         inAct = sum(conn$IncomingActivatorLinks), 
                         ourAct = sum(conn$OutgoingActivatorLinks),
                         inLinks = sum(conn$IncomingLinks), 
                         outLinks = sum(conn$OutgoingLinks)))
}

head(connectStats,3)
round(apply(connectStats,2,mean))
sum(round(apply(connectStats,2,mean))[c(1,3)]) # sum incoming
sum(round(apply(connectStats,2,mean))[c(2,4)]) # sum outgoing
round(mean(apply(connectStats[c(2,4)],1,sum))) # mean of sum of outgoing
round(sd (apply(connectStats[,c(2,4)],1,sum))) # s.d. of sum of outgoing
round(mean(connectStats[,2]/connectStats[,4]*100)) # mean ratio of outgoing repr. vs. act
round(sd(connectStats[,2]/connectStats[,4]*100)) # s.d. ratio of outgoing repr. vs. act

all(connectStats$inLinks > connectStats$outLinks)
# there are more inbound than outbound links

# for inLinks the promoter methylation is counted extra, while for outLinks not
conn = read.csv(file = connectFiles[1], header = T,sep = "\t",stringsAsFactors = F)
head(conn,2)
length(grep("GeneSpecificCopyNumber",conn$IncomingRepressorLinkGenes))
length(grep("GeneSpecificCopyNumber",conn$IncomingActivatorLinkGenes))
length(grep("GeneSpecificCopyNumber",conn$OutgoingRepressorLinkGenes))
length(grep("GeneSpecificCopyNumber",conn$OutgoingActivatorLinkGenes))
# only incomming links marked with methylation ("copynumber" in the original regNet paper)

promMethPerNW = NULL
for (f in connectFiles) {
    conn = read.csv(file = f, header = T,sep = "\t",stringsAsFactors = F)
    promMethPerNW = c(promMethPerNW,
                      length(grep("GeneSpecificCopyNumber",conn$IncomingRepressorLinkGenes)) + 
                      length(grep("GeneSpecificCopyNumber",conn$IncomingActivatorLinkGenes)))
}

promMethPerNW
round(mean(promMethPerNW)); round(sd(promMethPerNW))
length(promMethPerNW)
connectStats
all(connectStats$inLinks == connectStats$outLinks + promMethPerNW)
# incoming links are as many more as there are promoter methyl. links

# for each gene and each of the 25 NWs, how many outgoing links
outgoingLinksPerNW = NULL
for( f in connectFiles) {
    tmp = read.csv(file = f, header = T,sep = "\t",stringsAsFactors = F)
    if(is.null(outgoingLinksPerNW)) { 
        outgoingLinksPerNW = matrix(data = tmp$OutgoingLinks, ncol = 1)
        rownames(outgoingLinksPerNW) = tmp$Gene
    } else {
        rownames(tmp) = tmp$Gene
        outgoingLinksPerNW = cbind(outgoingLinksPerNW,tmp[ rownames(outgoingLinksPerNW),"OutgoingLinks"])
    }
}

head(outgoingLinksPerNW)
outgoingLinksMedian = apply(outgoingLinksPerNW,1,median)
outgoingLinksMean = apply(outgoingLinksPerNW,1,mean)
outgoingLinksSD = apply(outgoingLinksPerNW,1,sd)
head(outgoingLinksMedian)
head(outgoingLinksMean)
head(outgoingLinksSD)
quantile(outgoingLinksMean, 0:20/20)

genesAndAliasesForPaper = list(
    "KIT" = c("SCFR","CD117","C-Kit","P145","MASTC"),
    "BRAF" = c("BRAF1", "BRAF-1","RAFB1","B-RAF1","B-Raf","NS7","P94"), 
    "NRAS" = c("N-Ras","ALPS4","NRAS1","HRAS1","CMNS","NCMS","NS6"), 
    #"KRAS" = c("KRAS","C-Ki-Ras2","C-K-RAS","K-RAS2A","K-RAS2B","K-RAS4A","K-RAS4B","C-K-Ras","KI-RAS","K-Ras","CFC2","RALD","NS3","OES","NS"), 
    "NF1" = c("NFNS","VRNF","WSS"), 
    "TP53" = c("P53","LFS1","BMFS5","TRP53","BCC7"), 
    "DNMT3A" = c("DNMT3A2","HESJAS","TBRS"), 
    "DNMT3B" = c("FSHD4","ICF1","ICF"), 
    "TET1" = c("BA119F7.1","KIAA1676","CXXC6"), 
    "IDH1" = c("PICD","IDH","IDCD","IDPC","IDP"), 
    "ARID2" = c("BAF200","KIAA1557","SMARCF3","FLJ30619","CSS6","P200"), 
    "CTCF" = c("CFAP108","FAP108","MRD21"), 
    "CDKN2A" = c("P14ARF","CDK4I","MTS1","ARF","MTS1","CDKN2","CMM2","INK4","P16","P19","P16INK4a","MLM","P19Arf","INK4a","TP16"), 
    "PTEN" = c("MMAC1","TEP1","PTEN1","MHAM","BZS","PTENepsilon","PTENbeta","CWS1","GLM2","DEC"), 
    "ATF3" = c("GC01P211461","GC01P208625"),
    "MITF"=c("BHLHe32","WS2A","BHLHE32","WS2","COMMAD","MITF-A","CMM8"),
    "CREB" = c("CREB1","CREB-1"),
    "MYC" = c("BHLHe39","C-Myc","MYCC","BHLHE39","MRTL"),
    "FOS" = c("cFOS","AP-1","C-Fos","P55","G0S7"),
    "SOX10" = c("WS2E","DOM","WS4","PCWH","WS4C","SOX-10","SRY-Box 10","SRY"),
    "BRN2" = c("POU3F2","OCT7","POUF3","OTF7","OTF-7","Brain-2","Brn-2","Oct-7","N-Oct3"),
    "TET2" = c("KIAA1546","MDS","IMD75","FLJ20032"),
    "EHMT2" = c("KMT1C","G9A","C6orf30","BAT8","NG36","G9a","NG36/G9a","GAT8"),
    "TERT" = c("HEST2","TCS1","EST2","TRT","TP2","PFBMFT1","DKCA2","DKCB4","CMM9","HTRT"),
    "SDHD" = c("CybS","CII-4","PGL1","QPs3","SDH4","PGL","MC2DN3","CBT1","CWS3"))

g = names(genesAndAliasesForPaper)[1]
g
aliases = c(g,genesAndAliasesForPaper[[g]])
aliases
wh = which(toupper(aliases) %in% toupper(names(outgoingLinksMean)))
wh

outDF = NULL
for(g in names(genesAndAliasesForPaper)) {
    aliases = c(g,genesAndAliasesForPaper[[g]])
    wh = which(toupper(aliases) %in% toupper(names(outgoingLinksMean)))
    if(length(wh)>0) {
        
        # find average percentile of the gene's number of outgoing links via ecdf() function
        mean_percentile = mean(
            sapply(1:ncol(outgoingLinksPerNW), function(i) {
                n = outgoingLinksPerNW[aliases[wh[1]],i]
                100*ecdf(x = outgoingLinksPerNW[,i])(n)
            })
        )

        outDF = rbind(outDF,
                      data.frame(stringsAsFactors = F,gene=g, aliases=paste0(aliases[-1],collapse=","), 
                        mean_outgoing_links = outgoingLinksMean[aliases[wh[1]]], 
                                 sd_outgoing_links = outgoingLinksSD[aliases[wh[1]]], 
                                mean_percentile_outgoing_links = mean_percentile))
    } else {
        outDF = rbind(outDF,data.frame(stringsAsFactors = F,gene=g, 
               aliases=paste0(aliases[-1],collapse=","), mean_outgoing_links=NA, sd_outgoing_links=NA, 
                                mean_percentile_outgoing_links = NA))
    }
}

rownames(outDF) = outDF$gene
outDF = outDF[ order(-outDF$mean_percentile_outgoing_links,outDF$gene),]
outDF

write.xlsx(outDF, file = paste0(outDirectory,"SupplTable-S6-key-driver-ranks.xls"))

conn = read.csv(file = connectFiles[1], header = T,sep = "\t",stringsAsFactors = F)
activatorConn = matrix(data=0,nrow=nrow(conn),ncol=nrow(conn),
                       dimnames = c(list(conn$Gene),list(conn$Gene)))
repressorConn = activatorConn
fileCount=0
for(f in connectFiles) {
    fileCount = fileCount + 1
    message(fileCount," of ",length(connectFiles))
    conn = read.csv(file = f, header = T,sep = "\t",stringsAsFactors = F)
    for(i in 1:nrow(conn)) {
        #if(i %% 1000 == 0) { message("   line ",i," of ",nrow(conn))}
        sourceGene = conn$Gene[i]
        targetRepressedGenes = unlist(strsplit(conn$OutgoingRepressorLinkGenes[i],";"))
        targetActivatedGenes = unlist(strsplit(conn$OutgoingActivatorLinkGenes[i],";"))
        if(length(targetRepressedGenes)>0) {
            repressorConn[sourceGene,targetRepressedGenes] = repressorConn[sourceGene,targetRepressedGenes] + 1
        }
        if(length(targetActivatedGenes)>0) {
            activatorConn[sourceGene,targetActivatedGenes] = activatorConn[sourceGene,targetActivatedGenes] + 1
        }
    }
}
# 1-2 min

# connections in > 50% of NWs:
activatorConnBinary = ifelse(activatorConn>length(connectFiles)/2,1,0)
repressorConnBinary = ifelse(repressorConn>length(connectFiles)/2,1,0)

sum(activatorConnBinary) # number of consistent activators, 8658
sum(activatorConn)
sum(repressorConnBinary) # number of consistent repressors, 477
sum(repressorConn)

# convert to long format for plotting:

connectionsActive = NULL 
for(i in 1:nrow(activatorConnBinary)){
    wh = which(activatorConnBinary[i,]>0)
    if(length(wh)>=1) {
        connectionsActive =  
            c(connectionsActive,list(
                data.frame(stringsAsFactors = F, source = rownames(activatorConnBinary)[i], 
                             target=colnames(activatorConnBinary)[wh], type="activator")))
    }
}
length(connectionsActive)
connections = do.call(rbind,connectionsActive)

connectionsRepr = NULL 
for(i in 1:nrow(repressorConnBinary)){
    wh = which(repressorConnBinary[i,]>0)
    if(length(wh)>=1) {
        connectionsRepr =  
            c(connectionsRepr,list(
                data.frame(stringsAsFactors = F, source = rownames(repressorConnBinary)[i], 
                             target=colnames(repressorConnBinary)[wh], type="repressor")))
    }
}
length(connectionsRepr)
connections = rbind(connections,do.call(rbind,connectionsRepr))

head(connections,100)

# gene annotation from expression file
geneAnnot = read.csv(file = "Data/MelBrainSys-expression.csv", header = T,sep = "\t",stringsAsFactors = F)
geneAnnot = geneAnnot[,1:3]
rownames(geneAnnot) = geneAnnot$geneSymbol
head(geneAnnot,3)
unique(geneAnnot$chr)

# sort annotation by chromosomes:
o = order(as.numeric(gsub("chr","",geneAnnot$chr)), geneAnnot$pos)
geneAnnot = geneAnnot[o,]
unique(geneAnnot$chr)
# move chrY last
wh = which(geneAnnot$chr == "chrY")
geneAnnot = geneAnnot[ c(1:(min(wh)-1),(max(wh)+1):nrow(geneAnnot),wh),]
unique(geneAnnot$chr)

# get 1 number for each gene, so we can plot them normally
genesNum = 1:nrow(geneAnnot)
names(genesNum) = geneAnnot$geneSymbol
head(genesNum[ geneAnnot$geneSymbol[1:10]])
head(as.numeric(genesNum))

typeCol = c("activator"="red", "repressor"="blue")

chrBorders = c(0,sapply(unique(geneAnnot$chr),function(chr) max(which(geneAnnot$chr == chr)))+0.5)
chrMiddles = sapply(unique(geneAnnot$chr),function(chr) mean(range(which(geneAnnot$chr == chr))))

plotInteractions = function() {
    par(mar=c(4,4,1,1)) # slim border
    # empty plot to get the limit right:
    plot(NA, xlim=c(1,nrow(geneAnnot)), ylim=c(1,nrow(geneAnnot)),
         xlab = "source genes", ylab  = "target genes", xaxt='n', yaxt='n')
    # rectangles for chromosomes:
    for(i in 1:(length(chrBorders-1))) {
        rect(xleft = chrBorders[i], ybottom = chrBorders[i], xright = chrBorders[i+1], ytop = chrBorders[i+1])
    }
    # points for interactions;
    points(x = genesNum[connections$source], y = genesNum[ connections$target],
         bg = typeCol[connections$type], pch=22, col=NA,cex=0.3)
    axis(side = 1,at = chrBorders, labels=F, tick = T)
    axis(side = 1,at = chrMiddles, labels = gsub("chr","",unique(geneAnnot$chr)), tick = F)
    axis(side = 2,at = chrBorders, labels=F)
    axis(side = 2,at = chrMiddles, labels = gsub("chr","",unique(geneAnnot$chr)), tick = F)
    
}

plotInteractions()

png(width = 2400, height = 2400, res = 300, 
    filename = paste0(outDirectory,"/SupplFigure-S6-source-target-interactions.png"))
plotInteractions()
dev.off()

numOfNWs

linkStat3d = array(data=0,dim=c(nrow(conn),6,numOfNWs), 
      dimnames = c(list(conn$Gene),
                   list(c("IncomingRepressorLinks","IncomingActivatorLinks","IncomingLinks","OutgoingRepressorLinks",
                        "OutgoingActivatorLinks","OutgoingLinks")),
                      list(as.character(1:numOfNWs))))

nwIdx = 0
for( f in connectFiles) {
    nwIdx = nwIdx + 1
    conn = read.csv(file = f, header = T,sep = "\t",stringsAsFactors = F)
    rownames(conn) = conn$Gene
    linkStat3d[ conn$Gene,,nwIdx] = as.matrix(conn[,dimnames(linkStat3d)[[2]]])
}

rowMeans(linkStat3d[1:6,1,])

linkStatMeans = sapply(1:dim(linkStat3d)[2], function(i) rowMeans(linkStat3d[,i,]))
colnames(linkStatMeans) = dimnames(linkStat3d)[[2]]

head(linkStatMeans)

# merge activated and repressed target genes that appear in > 50% of NWs
wh = which(connections$type=="activator")
commonActivatedGenes = sapply( conn$Gene, function(gene) {
    wh2=which(connections$source==gene)
    paste0(connections[ intersect(wh,wh2),"target"],collapse=";")
})
wh = which(connections$type=="repressor")
commonRepressedGenes = sapply( conn$Gene, function(gene) {
    wh2=which(connections$source==gene)
    paste0(connections[ intersect(wh,wh2),"target"],collapse=";")
})

head(connections[ connections$type=="repressor",])

basePath

currGeneAnnot = geneAnnot
load(paste0(basePath,"annotation/Homo_sapiens.GRCh37.75-chrRename-noHaplo.RData"))
geneAnnotFull = geneAnnot
geneAnnot = currGeneAnnot

geneAnnotFull = geneAnnotFull[ geneAnnotFull$V3=="gene",c("V1","V2","V4","V5","V7","V9")]

colnames(geneAnnotFull) = c("chr","type","start","end","strand","V9")
geneAnnotFull$chr = as.character(geneAnnotFull$chr)
geneAnnotFull$strand = as.character(geneAnnotFull$strand)
geneAnnotFull$type = as.character(geneAnnotFull$type)
geneAnnotFull$V9 = as.character(geneAnnotFull$V9)

geneAnnotFull$EnsemblID = gsub(".*(ENSG[^;]+).*","\\1",geneAnnotFull$V9)
geneAnnotFull$gene = gsub(".*gene_name ([^;]+).*","\\1",geneAnnotFull$V9)

wh = which(duplicated(geneAnnotFull$gene))
length(wh)
geneAnnotFull = geneAnnotFull[-wh,]

rownames(geneAnnotFull) = geneAnnotFull$gene

head(geneAnnotFull)

connectStatsOut = data.frame(stringsAsFactors = F, source_gene = conn$Gene, 
                             geneAnnotFull[conn$Gene,c("EnsemblID","chr","type","start","end","strand")],
                             linkStatMeans[,c("OutgoingRepressorLinks","OutgoingActivatorLinks")], 
                            activated_target_genes=commonActivatedGenes, 
                            repressed_target_genes = commonRepressedGenes)
connectStatsOut

# load gene pathway memberships (from circos plot notbook)
geneList = readRDS(file = paste0(basePath,"annotation/gene-pathways-membership-lists.rds"))

names(geneList)
head(geneList$cancerCensus)

genePWannot = data.frame(stringsAsFactors = F,matrix(
    data = "no",nrow = nrow(connectStatsOut), ncol = length(geneList),
    dimnames = c(list(connectStatsOut$source_gene), list(colnames=names(geneList))) ))

for(i in 1:length(geneList)) {
    wh = which(connectStatsOut$source_gene %in% geneList[[i]])
    if(length(wh)>0) {
        genePWannot[wh,i] = "yes"
    }
}

colnames(genePWannot) = c("immune system", "metabolism","signaling pathway","cancer census","oncogene","tumor suppressor","transcription factor")
genePWannot

connectStatsOut = cbind(connectStatsOut,genePWannot[connectStatsOut$source_gene,])
connectStatsOut

write.xlsx(x = connectStatsOut, row.names = F,
           file = paste0(outDirectory,"SupplTable-S5-connectivity-stats.xls"))
