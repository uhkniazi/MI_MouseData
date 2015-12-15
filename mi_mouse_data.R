# Name: mi_mouse_data.R
# Date: 14/12/15
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Desc: summarize the mouse experiments


library(org.Mm.eg.db)
library(limma)

dfDat = read.csv('Data_external/input.csv', header = T)

# choose genes with refseq id nm_ present
i = grepl('NM_', dfDat$RefSeq)
dfDat = dfDat[i,]

# replace the refseq column with only the nm_ id and remove other ids
ref = gsub('.*(NM_\\d+).*', replacement = '\\1', dfDat$RefSeq)
# remove duplicate ids
i = duplicated(ref)
ref = ref[!i]
dfDat = dfDat[!i,]
rownames(dfDat) = ref

# get annotation for the data
dfAnnotation = select(org.Mm.eg.db, ref, columns = c('GENENAME', 'SYMBOL', 'UNIPROT'), keytype = 'REFSEQ')
# remove duplicated results, as some transcripts map to multiple alternates
dfAnnotation = dfAnnotation[!duplicated(as.character(dfAnnotation$REFSEQ)),]
rownames(dfAnnotation) = as.character(dfAnnotation$REFSEQ)
# some symbols and uniprot ids are duplicated, remove those
dfAnnotation = dfAnnotation[!duplicated(as.character(dfAnnotation$UNIPROT)),]
dfAnnotation = na.omit(dfAnnotation)
# remove these from the data table
dfDat = dfDat[rownames(dfAnnotation),]
rownames(dfDat) = as.character(dfAnnotation$UNIPROT)

dfDat.sub = dfDat[,1:288]

fGroups = rep(NA, times=288)
cn = colnames(dfDat.sub)
f = strsplit(cn[17:288], split = '\\.')
f2 = sapply(f, function(x) x[1])

# first 16 are controls and others are treatments
fGroups = rep(NA, times=288)
fGroups[1:16] = '0.Cont'
fGroups[17:288] = f2
fGroups[284:288] = 'SPS'
fGroups = factor(fGroups)

## perform DE analysis
mDat = as.matrix(dfDat.sub)

design = model.matrix(~fGroups)
colnames(design) = levels(fGroups)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# sanity check
nrow(dfAnnotation) == nrow(mDat)
fit$genes = dfAnnotation
topTable(fit, adjust='BH')
fSamples = fGroups

# look at top tables for each comparison
for (i in 2:length(levels(fSamples))){
  print(paste(levels(fSamples)[i], i))
  print(topTable(fit, coef=i, adjust='BH'))
}

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.adj = p.adjust(fit$p.value[,i], method = 'BH')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.01]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

######### Volcano plots
# plot volcano plots
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  #dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]],]
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < 0.1)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=paste(names(n[i])), xlim=c(-5, 5))
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

### group the genes
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)
mCommonGenes = sapply(seq_along(n), function(x) cvCommonGenes %in% lSigGenes.adj[[names(n[x])]])
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(n)


########### pathway analysis using CGraph library
source('../CGraphClust/CGraphClust.R')
library(downloader)
p.old = par()
# uniprot annotation for data
cvGenes = rownames(mCommonGenes)
dfGenes = dfAnnotation
rownames(dfGenes) = as.character(dfGenes$UNIPROT)

# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
# map reactome pathways
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfGenes$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfGenes$UNIPROT)
dfReactome.sub$SYMBOL = dfGenes$SYMBOL[i]
dfGraph = dfReactome.sub[,c('SYMBOL', 'V2')]
dfGraph = na.omit(dfGraph)

# get expression data
# separate the factor and the count matrix
mCounts = mDat[rownames(dfGenes),]
fGroups = as.character(fSamples)
# select only the groups with significant genes
n = (which(sapply(lSigGenes.adj, length) >= 20)) + 1
i = which(fSamples %in% levels(fSamples)[c(1, n)])
fGroups = factor(fGroups[i], levels = levels(fSamples)[c(1, n)])

# subset the count matrix
n = (which(sapply(lSigGenes.adj, length) >= 20)) + 1
i = which(fSamples %in% levels(fSamples)[c(1, n)])
mCounts = mCounts[,i]
colnames(mCounts) = fGroups
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)
# replace colnames from uniprot to symbols
cn = dfGenes[colnames(mCounts), 'SYMBOL']
colnames(mCounts) = cn
# select genes that have a reactome term attached
n = unique(dfGraph$SYMBOL)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
levels(fGroups)

# create a correlation matrix to decide cor cutoff
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# stabalize the data and check correlation again
# mCounts.bk = mCounts
# # stabalize the data
# mCounts.st = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
# rownames(mCounts.st) = fGroups
# 
# # create a correlation matrix
# mCor = cor(mCounts.st)
# # check distribution 
# hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
# axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)

## community structure
# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 10)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 10)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

## centrality diagnostics
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
par(p.old)

# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

### function to get gene annotation
f_dfGetGeneAnnotation = function(cvSymbol = NULL) {
  df = dfGenes[dfGenes$SYMBOL %in% cvSymbol,]
  return(df)
}

dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$SYMBOL),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

dir.create('Results', showWarnings = F)
write.csv(dfTopGenes.cent, file='Results/Top_Centrality_Genes.csv')

## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
m = plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
# principal component plots
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
#m = getClusterMarginal(oGr, t(mCounts))

csClust = rownames(m)
length(csClust)
i = 1
plot.cluster.variance(oGr, m[csClust[i:(i+1)],], fGroups, log = FALSE); i = i+2

i = 1
temp = t(as.matrix(m[csClust[i],]))
rownames(temp) = csClust[i]
plot.cluster.variance(oGr, temp, fGroups, log=FALSE); i = i+1

boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)


dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
sort(table(dfCluster$cluster))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))

df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$SYMBOL),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Results/Clusters.csv')


# plot the graphs at each contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=10)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}






