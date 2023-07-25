library(ggplot2)
library(gridExtra)
library(xlsx)
library(dplyr)
library(scales)
library(ggplotify)
library(survival)
library(survminer)

basePath = "/data/bcu_projects/MelBrainSys_PostdocProject_Gruetzmann/publications/2022-my-MelBrainSys-paper/scripts-etc-for-publication/"
setwd(basePath)

outDirectory = paste0(basePath,"FiguresTables/")

# for each samples -> which subgroup it belongs to 
subgroupPerSample = c(P04_Br = 'SG1', P04_Sk_GA = 'SG1', P08_Br = 'SG1 | SG2 | SG3', P08_St_BA = 'SG1', P16_Br = 'SG1', P16_Lu='SG1', 
  P42_Br_GA = 'SG1 | SG3', P42_Ly_GA = 'SG1', P08_St_GA = 'SG2', P18_Br = 'SG2 | SG3',
  P18_Lu_YA = 'SG2', P39_Br = 'SG2', P39_Lu = 'SG2', P03_Br = 'SG3', P03_Lu = 'SG3', 
  P08_St_YA = 'SG3', P18_Lu_GA = 'SG3', P42_Ly_YA ='SG3')
subgroupPerSample
subgroupPerSampleNA = c(P04_Br = 'SG1', P04_Sk_GA = 'SG1', P08_Br = "NA", P08_St_BA = 'SG1', P16_Br = 'SG1', P16_Lu='SG1', 
  P42_Br_GA = "NA", P42_Ly_GA = 'SG1', P08_St_GA = 'SG2', P18_Br = "NA",
  P18_Lu_YA = 'SG2', P39_Br = 'SG2', P39_Lu = 'SG2', P03_Br = 'SG3', P03_Lu = 'SG3', 
  P08_St_YA = 'SG3', P18_Lu_GA = 'SG3', P42_Ly_YA ='SG3')
subgroupPerSampleNA
melBrainSysSamples = names(subgroupPerSample)
melBrainSysSamples

TCGAinfo = read.xlsx(file = "annotation/NIHMS698912-supplement-3.xlsx", 
                     sheetName = "Supplemental Table S1D" , startRow = 2, as.data.frame = T, stringsAsFactors=F)
head(TCGAinfo,3)

melBrainSysExpr = read.csv(file = "data/MelBrainSys-expression.csv", header = T,sep = "\t", stringsAsFactors = F)
head(melBrainSysExpr,3)
melBrainSysExpr = melBrainSysExpr[ ,melBrainSysSamples ]
head(melBrainSysExpr,3)
nrow(melBrainSysExpr)

tcgaExpr = read.csv(file = "data/TCGA-expression.csv", header = T,sep = "\t", stringsAsFactors = F)
head(tcgaExpr,3)
colnames(tcgaExpr) = gsub("\\.","\\-",colnames(tcgaExpr))

sampleCor = matrix(data=NA, nrow=ncol(tcgaExpr)-3, ncol=length(melBrainSysSamples), 
                   dimnames = c(list(colnames(tcgaExpr)[-c(1:3)]), list(melBrainSysSamples)))
sampleCorPval = matrix(data=NA, nrow=ncol(tcgaExpr)-3, ncol=length(melBrainSysSamples), 
                   dimnames = c(list(colnames(tcgaExpr)[-c(1:3)]), list(melBrainSysSamples)))
head(sampleCor,2)

for(s1 in rownames(sampleCor)) {
    for(s2 in melBrainSysSamples) {
        sampleCor[ s1, s2] = cor(x = tcgaExpr[,s1], y = melBrainSysExpr[,s2] )
        sampleCorPval[ s1, s2] = cor.test(x = tcgaExpr[,s1], y = melBrainSysExpr[,s2] )$p.value
    }
}

quantile(sampleCorPval, probs = c(0,50,97,98,99,100)/100) # all cor. are highly signif
head(sampleCor,2)
quantile(sampleCor) # 0.21 .. 0.74, median 0.63

sampleCor[1:2,1:3]
sampleCorWhMax = apply(sampleCor,1,which.max)
table(sapply(1:length(sampleCorWhMax), function(i) length(which(sampleCor[i,] == sampleCor[i,sampleCorWhMax[i]]))))
# for each TCGA sample, there is exactly 1 MelBrainSys sample with max correl., not ambiguous
#  cbind(sampleCorWhMax,sampleCor)  
table(sampleCorWhMax); # 18 MelBrainSys samples, TCGA samples are assigned to 9 of them

TCGAassignSubgroups = unname(subgroupPerSample[melBrainSysSamples[sampleCorWhMax]])
table(TCGAassignSubgroups) 
TCGAassignSubgroups = unname(subgroupPerSampleNA[melBrainSysSamples[sampleCorWhMax]])
table(TCGAassignSubgroups) # 80 of 270 cannot be assigned unambiguously
# -> omit these TCGE samples

names(TCGAassignSubgroups) = rownames(sampleCor)
TCGAassignSubgroups = TCGAassignSubgroups[ which(!TCGAassignSubgroups=="NA")]
head(TCGAassignSubgroups)
table(TCGAassignSubgroups) # each subgroup is present, none marginally frequent

length(intersect(TCGAinfo$Name, names(TCGAassignSubgroups))) # 190 of 270 patients have subgroup
TCGAinfo$Subgroup = paste0("TCGA: ",TCGAassignSubgroups[TCGAinfo$Name]) 
table(TCGAinfo$Subgroup)
TCGAinfoClean = TCGAinfo[ !(TCGAinfo$Subgroup=="TCGA: NA"),]
table(TCGAinfoClean$Subgroup)

head(TCGAinfoClean$CURATED_T_STAGE_AT_DIAGNOSIS_SIMPLE)
TCGAinfoClean$T_STAGE = as.character(TCGAinfoClean$CURATED_T_STAGE_AT_DIAGNOSIS_SIMPLE)
TCGAinfoClean$T_STAGE[ grep("A", TCGAinfoClean$T_STAGE)] = NA
TCGAinfoClean$T_STAGE[ grep("\\-", TCGAinfoClean$T_STAGE)] = NA
table(TCGAinfoClean$T_STAGE)
TCGAinfoClean$T_STAGE_ord = factor(TCGAinfoClean$T_STAGE, ordered = T, 
                                       levels=c("T0", "T1", "T2", "T3", "T4" ))

tStageDistrib = TCGAinfoClean %>% group_by(Subgroup) %>% 
    summarise(T_stage_T0 = length(which(T_STAGE == "T0")),
              T_stage_T1 = length(which(T_STAGE == "T1")),
              T_stage_T2 = length(which(T_STAGE == "T2")),
              T_stage_T3 = length(which(T_STAGE == "T3")),
              T_stage_T4 = length(which(T_STAGE == "T4")))
tStageDistrib

# Kruskal Wallis test:
kruskal.test(T_STAGE_ord ~ Subgroup, data=TCGAinfoClean ) # 0.0465
kruskal.test(T_STAGE_ord ~ Subgroup, 
                 data=TCGAinfoClean[ TCGAinfoClean$Subgroup %in% c("TCGA: SG1","TCGA: SG2"),] ) # 0.01678
kruskal.test(T_STAGE_ord ~ Subgroup, 
                 data=TCGAinfoClean[ TCGAinfoClean$Subgroup %in% c("TCGA: SG1","TCGA: SG3"),] ) # 0.48
kruskal.test(T_STAGE_ord ~ Subgroup, 
                 data=TCGAinfoClean[ TCGAinfoClean$Subgroup %in% c("TCGA: SG2","TCGA: SG3"),] ) # 0.22

plotTstage = ggplot(TCGAinfoClean[ !is.na(TCGAinfoClean$T_STAGE_ord),], aes(x=T_STAGE_ord, fill=T_STAGE_ord)) + 
    geom_bar() + facet_wrap(~Subgroup, ncol = 1) + xlab("T stage") + ylab("Number of patients") + 
    theme(plot.margin = margin(2,0,0,1, "cm"))
    #+ theme(plot.background = element_rect(colour = "grey", fill=NA, size=5))
plotTstage  = plotTstage + guides(fill=guide_legend(title="T stage"))
plotTstage = plotTstage + labs(tag = "D") + theme(plot.tag.position = c(0.0, 1.08))
plotTstage
#png(filename = paste0(outDirectory,"Figure-tumor-stage-vs-subgroup.png"), width = 600,height = 400)
#plotTstage + guides(fill=guide_legend(title="T stage"))
#dev.off()

table(TCGAinfoClean$MUTATIONSUBTYPES)
TCGAinfoClean$MUTATIONSUBTYPES[ which(TCGAinfoClean$MUTATIONSUBTYPES=="-")] = NA

plotMutSubtype = ggplot(TCGAinfoClean[ !is.na(TCGAinfoClean$MUTATIONSUBTYPES),], aes(x=Subgroup,fill=MUTATIONSUBTYPES)) + 
    geom_bar(position = "fill") + ylab("Ratio of patients") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
         plot.margin = margin(1,0,0,1, "cm"))
#, 
         #plot.background = element_rect(colour = "grey", fill=NA, size=5))
plotMutSubtype = plotMutSubtype + guides(fill=guide_legend(title="Mutation subtype"))
plotMutSubtype = plotMutSubtype + labs(tag = "A") + theme(plot.tag.position = c(-0.05, 1.05))
plotMutSubtype
#grid.arrange(p1, p2, widths=c(1,1.5), ncol=2)
#png(filename = paste0(outDirectory,"mutationsubtypes-vs-subgroup.png"), width=700)
#grid.arrange(p1, p2, widths=c(1,1.55), ncol=2)
#dev.off()

m = table(TCGAinfoClean$MUTATIONSUBTYPES, TCGAinfoClean$Subgroup)
m
chisq.test(m) # p-value = 0.01103

table(TCGAinfoClean$MethTypes.201408)

plotMethSubtype = ggplot(TCGAinfoClean[ !is.na(TCGAinfoClean$MethTypes.201408),], aes(x=Subgroup,fill=MethTypes.201408)) + 
    geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
         plot.margin = margin(1,0,0,1, "cm")) + 
    ylab("Ratio of patients") #, 
         #plot.background = element_rect(colour = "grey", fill=NA, size=5))
plotMethSubtype = plotMethSubtype + guides(fill=guide_legend(title="Methylation subtype"))
plotMethSubtype = plotMethSubtype + labs(tag = "B") + theme(plot.tag.position = c(-0.05, 1.05))
plotMethSubtype
#grid.arrange(p1, p2, widths=c(1,1.5), ncol=2)
#png(filename = paste0(outDirectory,"meth-types-vs-subgroup-barplot.png"), width=700)
#grid.arrange(p1, p2, widths=c(1,1.55), ncol=2)
#dev.off()

m = table(TCGAinfoClean$MethTypes.201408, TCGAinfoClean$Subgroup)
m
chisq.test(m) # p-value =  0.02773

TCGAinfoClean$TOTAL.MUTATIONS = as.numeric(TCGAinfoClean$TOTAL.MUTATIONS)
quantile(TCGAinfoClean$TOTAL.MUTATIONS,na.rm = T)

p = ggplot(TCGAinfoClean, aes(x=TOTAL.MUTATIONS, y=Subgroup, fill=Subgroup)) + geom_violin(scale = "count")
p 

kruskal.test(TOTAL.MUTATIONS ~ Subgroup, data=TCGAinfoClean ) # 0.4392

TCGAinfoClean$RNASEQ.CLUSTER_CONSENHIER[ which(TCGAinfoClean$RNASEQ.CLUSTER_CONSENHIER=="-")] = NA
table(TCGAinfoClean$RNASEQ.CLUSTER_CONSENHIER)

m = table(TCGAinfoClean$RNASEQ.CLUSTER_CONSENHIER, TCGAinfoClean$Subgroup)
m
chisq.test(m) # p-value = 0.02078

plotExprCluster = ggplot(TCGAinfoClean[ !is.na(TCGAinfoClean$RNASEQ.CLUSTER_CONSENHIER),], aes(x=Subgroup,fill=RNASEQ.CLUSTER_CONSENHIER)) + 
    geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
         plot.margin = margin(1,0,0,1, "cm")) + 
    ylab("Ratio of patients") #, 
         #plot.background = element_rect(colour = "grey", fill=NA, size=5))
plotExprCluster = plotExprCluster + guides(fill=guide_legend(title="Expression cluster"))
plotExprCluster = plotExprCluster + labs(tag = "C") + theme(plot.tag.position = c(-0.05, 1.05))
plotExprCluster
# grid.arrange(p1, p2, widths=c(1,1.5), ncol=2)
# png(filename = paste0(outDirectory,"RNAseq-clusters-vs-subgroup-barplot.png"), width=700)
# grid.arrange(p1, p2, widths=c(1,1.55), ncol=2)
# dev.off()

# main text figure A-D
png(paste0(outDirectory,"Figure-6-TCGA-associations.png"), width = 2400, height = 2400, res=210)
margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))
pl = list(plotMutSubtype,plotMethSubtype, plotExprCluster,
          plotTstage)
grid.arrange(grobs = lapply(pl, "+", margin),  ncol=2)
dev.off()

TCGAinfoClean$vital_status = ifelse(TCGAinfoClean$CURATED_VITAL_STATUS=="Alive",0,
                                    ifelse(TCGAinfoClean$CURATED_VITAL_STATUS=="Dead",1,NA))
table(TCGAinfoClean$vital_status)

TCGAinfoClean$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU = as.numeric(TCGAinfoClean$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU)
quantile(TCGAinfoClean$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU, na.rm = T)

fit <- survfit(Surv(CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU, vital_status) ~ Subgroup,
               data = TCGAinfoClean)

plotSurvival = ggsurvplot(fit, data = TCGAinfoClean, risk.table = T, pval = T, 
               conf.int = F, break.time.by = 1000, ggtheme = theme_minimal(), 
               risk.table.y.text.col = T, risk.table.y.text = F)
#plotSurvival$plot =  plotSurvival$plot + theme(plot.margin = margin(1,0,0,1, "cm"))
#plotSurvival$plot = plotSurvival$plot + labs(tag = "E") + theme(plot.tag.position = c(0, 0.97))
plotSurvival
#png(paste0(outDirectory,"Figure-Kaplan-Meier-subgroup.png"), width=1000,height=700, res = 100)
#p
#dev.off()

# supplement figure E
plSurv = list(plotSurvival$plot, plotSurvival$table)
png(paste0(outDirectory,"SupplFigure-S8-TCGA-survival.png"), width = 3500, height = 1800, res=320)
grid.arrange(grobs = plSurv, ncol=1, heights=c(3,1))
dev.off()
