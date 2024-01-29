# retrieval of gene disease associations from DisGeNET via disgenet2r
# vignette: https://www.disgenet.org/static/disgenet2r/disgenet2r.html

# install disgenet2r:
# library(devtools)
# install_bitbucket("ibi_group/disgenet2r")

library(xlsx)
library(disgenet2r)

# register with DisGeNET and enter credentials here:
disgenet_api_key = get_disgenet_api_key(email = "your@email.com", password = "your password" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

setwd("..")

# load all target genes
# separate SG1 genes into log2 impact ration >/<0 ??
targetGeneFile = "FiguresTables/SupplTable-S10-target-gene-candidates.xls"
genesSG1 = read.xlsx(file = targetGeneFile, header = T, sheetIndex = 1)$gene
genesSG3 = read.xlsx(file = targetGeneFile, header = T, sheetIndex = 2)$gene
genesSG2 = read.xlsx(file = targetGeneFile, header = T, sheetIndex = 3)$gene
genesSG1
genesSG2
genesSG3

minScore = 0

res1 = extract(gene2disease( gene = genesSG1, vocabulary = "HGNC",
                     database = "ALL" , score =c(minScore, 1), warnings = F))
nrow(res1)
# filter disease_class_name on neoplasm (beign/malign), and optionally disease_name on melanoma, count
wh =  grep("neoplasm",res1$disease_class_name,ignore.case = T)
res1 = res1[wh,]

wh = grep("melanoma",res1$disease_name, ignore.case = T)
message(nrow(res1[wh,]), " SG1 melanoma interaction:") # 3 interactions, 3 genes: LDLR MST1R GGCT
cat(length(unique(res1$gene_symbol[wh])), "genes:", unique(res1$gene_symbol[wh]))
res1[wh,]
message(nrow(res1), " SG1 neoplasm interaction:")
cat(length(unique(res1$gene_symbol)), "genes:", unique(res1$gene_symbol))
res1

res2 = extract(gene2disease( gene = genesSG2, vocabulary = "HGNC",
                     database = "ALL" , score =c(minScore, 1), warnings = F))
nrow(res2)
# filter disease_class_name on neoplasm (beign/malign), and optionally disease_name on melanoma, count
wh =  grep("neoplasm",res2$disease_class_name,ignore.case = T)
res2 = res2[wh,]

wh = grep("melanoma",res2$disease_name, ignore.case = T)
message(nrow(res2[wh,]), " SG2 melanoma interaction:") # 4 interactions, 4 genes: AIM2 RGS1 FES GATA3
cat(length(unique(res2$gene_symbol[wh])), "genes:", unique(res2$gene_symbol[wh]))
res2[wh,]
message(nrow(res2), " SG2 neoplasm interaction:")
cat(length(unique(res2$gene_symbol)), "genes:", unique(res2$gene_symbol))
res2

res3 = extract(gene2disease( gene = genesSG3, vocabulary = "HGNC",
                     database = "ALL" , score =c(minScore, 1), warnings = F))
nrow(res3)
# filter disease_class_name on neoplasm (beign/malign), and optionally disease_name on melanoma, count
wh =  grep("neoplasm",res3$disease_class_name,ignore.case = T)
res3 = res3[wh,]

wh = grep("melanoma",res3$disease_name, ignore.case = T)
message(nrow(res3[wh,]), " SG3 melanoma interaction:") # 5 interactions, 3 genes: PYCARD GGCT SLAMF6
cat(length(unique(res3$gene_symbol[wh])), "genes:", unique(res3$gene_symbol[wh]))
res3[wh,]
message(nrow(res3), " SG3 neoplasm interaction:")
cat(length(unique(res3$gene_symbol)), "genes:", unique(res3$gene_symbol))
res3

outfile = 'FiguresTables/Suppl-Table-S12-DisGeNET-analysis.xls'

write.xlsx(x = res1, file = outfile, sheetName = "SG1 Higher", row.names = F)
write.xlsx(x = res3, file = outfile, sheetName = "SG3 Slightly lower", row.names = F, append=T)
write.xlsx(x = res2, file = outfile, sheetName = "SG2 Lower", row.names = F, append=T)
