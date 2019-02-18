library(DESeq2)
library(BiocParallel)
library(edgeR)
library(limma)
library(factoextra)
#register(MulticoreParam(1))
### common functions
nsv_filter<-function(counts,verbose=FALSE) {
  if(verbose) {
    print(paste0("Input:",ncol(counts)))
  }
  nsv<-caret::nearZeroVar(counts)
  counts<-counts[,-nsv]
  if(verbose) {
    print(paste0("Output:",ncol(counts)))
  }
  return(counts)
}
### parameters
nCore<-32
dir<-"/home/imlay/storage/GTEx/data/"
### loading data
setwd(dir)
meta<-read.table("GTEx_v7_Annotations_SubjectPhenotypesDS.txt",sep="\t",header = TRUE,na.strings = "na")
meta2<-data.table::fread("GTEx_v7_Annotations_SampleAttributesDS.txt",sep="\t",nThread=nCore)
counts<-data.table::fread("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct",nThread=nCore)
### fix frames
get_counts<-function(raw_frame) {
  gene_names<-raw_frame$Name
  counts<-as.matrix(raw_frame[,-c(1:2)])
  dimnames(counts)[[1]]<-gene_names
  return(t(counts))
}
mat_counts<-get_counts(counts)
rm(counts)
#mat_counts<-nsv_filter(mat_counts,verbose=TRUE)
load("mat_counts.Rdata")
meta2<-meta2[meta2$SAMPID %in% rownames(mat_counts),]
meta2<-meta2[order(match(meta2$SAMPID,dimnames(mat_counts)[[1]])),] #orders runs in meta with the raw data
rownames(meta2)<-meta2$SAMPID
meta2<-meta2[,c("SAMPID","SMTS")]
meta2$SMTS<-factor(meta2$SMTS)
#SAMPID_SUBJID<-paste0(strsplit(meta2$SAMPID,"-")[1:2])

## EDA
plot_EDA<-function(mat_counts,meta) {
  par(las=2)
  par(mar=c(5,8,4,2))
  barplot(table(meta2$SMTS)[order(table(meta2$SMTS))],main="Tissues",horiz=TRUE)
}
## edgeR
x<-DGEList(t(mat_counts),samples=meta2)
x <- calcNormFactors(x, method = "TMM")
cpm<-cpm(x)
lcpm<-cpm(x,log=TRUE)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
ntop <- 500
rv <- rowVars(lcpm)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- lcpm[select, ]
res.pca<-prcomp(t(mat))
Tissue<-meta2$SMTS
res.df<-as.data.frame(res.pca$x)
percentage <- round(res.pca$sdev^2 / sum(res.pca$sdev^2) * 100, 2)
percentage <- paste( colnames(res.df), "(", paste( as.character(percentage), "%", ")", sep=""))
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p<-ggplot(res.df,aes(x=PC1,y=PC2,color=Tissue ))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])+labs(title = "PCA of all samples",subtitle = "lcpm;n=500")
print(p)

## DESeq2
#dds = DESeqDataSetFromMatrix(countData = t(mat_counts), colData = meta2,
#                            design = formula(~SMTS))
#dds<-DESeq(dds)
#res<-results(dds)
#vsd <- varianceStabilizingTransformation(dds,blind = FALSE) # blind is more appropriate since 

#stopCluster(cl)
