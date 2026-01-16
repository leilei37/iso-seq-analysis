############normalization
#sort data
expr<-read.table('/Users/au760484/Downloads/legume/isoseq/expressionMatrix_nocotig/classification_all_nocotig_nonovel_aggregatesamegene.tsv',header = T)
rownames(expr)<-expr[,2]
expr<-expr[,-c(1:2)]
expr<-expr[,-c(1)]

colnames(expr)<-gsub('FL.F','Core',colnames(expr))
#hist
hist(rowSums(expr!=0),labels = T)
#filter genes
left_gene<-rownames(expr)[which(rowSums(expr!=0)>(19.8*3))]
#normalize
condition<-data.frame(cbind(colnames(expr),1))
colnames(condition)[2]<-'condition'
library(DESeq2)
dds <- DESeqDataSetFromMatrix(expr[left_gene,],
                              group,
                              design = ~ V3)
vst<-vst(dds,blind = T)
pdf(file = '~/Downloads/legume/isoseq/vst.pdf')
p<-plotPCA(vst,'condition')
print(p)
dev.off()
expr.nor<-vst@assays@data@listData[[1]]
#############################co-expression
library(WGCNA)
raw_data <- read.table('~/Downloads/legume/isoseq/osca_Mac_v0.45/expr_176_filter5re_v2.txt')
rownames(raw_data)<-raw_data[,1]
colnames(raw_data)<-raw_data[1,]
datfilt <- raw_data[-1,-1]
datfilt<-apply(datfilt,2,as.numeric)
rownames(datfilt)<-raw_data[-1,1]

#hclust
sampleTree = hclust(dist(datfilt), method = "average")
dev.off()
pdf('Downloads/Expression_Data/hclust.pdf',width = 30,height = 12)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex = 0.5)
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 155, minSize = 200)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
#datExpr0 = dataExpr[keepSamples, ]

nGenes = ncol(datfilt)
nSamples = nrow(datfilt)
type = "signed"
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datfilt, powerVector=powers,
                        networkType=type, verbose=5)
png('Downloads/Expression_Data/power.png')
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# R-square=0.85
abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
power = sft$powerEstimate
power = 20

dev.off()
#power<-18 for signed

type = "signed"
corType = "pearson"
exprMat <- "/Users/au760484/Downloads/legume/isoseq/faba_isoseqnor250610.txt"

maxPOutliers = ifelse(corType=="pearson",1,0.05)

net = blockwiseModules(datfilt, power = power, maxBlockSize = nGenes,
                       #                       deepSplit = 4,
                       TOMType = type, minModuleSize = 30,
                       #                       reassignThreshold = 0, 
                       #                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,networkType=type,
                       #                       maxPOutliers=maxPOutliers, 
                       #                       loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       #                       detectCutHeight = 0.995,
                       #                       maxCoreScatter = 1,maxAbsCoreScatter=1,
                       #                       pamStage=F,
 #                      mergeCutHeight = 0.1,
                       #                       trapErrors = T,
                       #                      stabilityCriterion='Common fraction',
                       #                      pamRespectsDendro = T,
                       #corFnc = cor,corOptions = list(use='p',method='spearman')
                       #corFnc = cor,corOptions = list(use='p',method='kendall'),
                       verbose = 3)
net.gene<-table(net$colors)
#write.table(net.gene,file = 'Downloads/Expression_Data/each_module_gene_filtermore_addprotein.txt')
# Convert labels to colors for plotting
moduleColors = net$colors

# Plot the dendrogram and the module colors underneath
png('Downloads/medicago/WGCNA/net_filtermore_addprotein.png')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors.df<-as.data.frame(moduleColors)
write.csv(moduleColors.df,file = '~/Downloads/legume/isoseq/module_gene_filterv2nor.csv')
