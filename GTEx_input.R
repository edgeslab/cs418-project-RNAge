library(edgeR)
library(limma)
library(foreach)
library(doSNOW)

### common functions
nzv_filter<-function(counts,verbose=FALSE) {
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
get_counts<-function(raw_frame) {
  gene_names<-raw_frame$Name
  counts<-as.matrix(raw_frame[,-c(1:2)])
  dimnames(counts)[[1]]<-gene_names
  return(t(counts))
}

plotmanual <- function(data, grp, title = NULL,alpha=.7,subtitle="lcpm;n=500") {
  res.pca <- prcomp(data, scale. = TRUE, center = TRUE)
  res.df <- as.data.frame(res.pca$x)
  percentage <-
    round(res.pca$sdev ^ 2 / sum(res.pca$sdev ^ 2) * 100, 2)
  percentage <-
    paste(colnames(res.df), "(", paste(as.character(percentage), "%", ")", sep =
                                         ""))
  theme <-
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      plot.margin = unit(c(1, 1, 1, 1), "line")
    )
  p <- ggplot(res.df, aes(
    x = PC1,
    y = PC2,
    color = grp
  ))
  p <-
    p + geom_point(alpha=alpha) + theme + xlab(percentage[1]) + ylab(percentage[2]) + labs(title = title, subtitle = subtitle)
  return(p)
}
metaBar<-function(meta,by_col,mar=c(5,8,4,2),main=FALSE) {
  if(!(main)) main=by_col
  par(las=2)
  par(mar=mar)
  barplot(table(meta[[by_col]])[order(table(meta[[by_col]]))],main=main,horiz=TRUE)
}
get_PCAmat<-function(lcpm) {
  #mat will need to be horizontal.
  ntop <- 500
  rv <- matrixStats::rowVars(lcpm)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- lcpm[select, ]
  return(mat)
}
tissuePCA<-function(lcpm,meta,tissue,grp) {
  sel<-meta[SMTS==tissue]
  mat<-get_PCAmat(lcpm[,sel$SAMPID])
  plotmanual(t(mat),sel[[grp]],title=tissue,subtitle=paste0("lcpm ; n=",nrow(sel)))
  ggsave(file.path(root_dir,"plots",paste0(tissue,"_",grp,"_PCA.png")),width = 8,height=8)
}
### parameters
nCore<-32
root_dir<-"/home/imlay/storage/cs418-project-RNAge/"
data_dir<-"/home/imlay/storage/cs418-project-RNAge/data"

### loading data
setwd(data_dir)
meta<-data.table::fread("merged_meta.tsv",sep="\t",header=TRUE,na.strings="na",nThread=nCore)
counts<-data.table::fread("All_Tissue_Site_Details.combined.reads.gct",nThread=nCore)
setwd(root_dir)

### process data and metadata
counts<-get_counts(counts)
rownames(meta)<-meta$SAMPID
keep<-meta$AGE
### edgeR pre-processing
counts<-DGEList(t(counts),samples=meta)
counts <- calcNormFactors(counts, method = "TMM")
cpm<-cpm(counts)
lcpm<-cpm(counts,log=TRUE)

## EDA
EDA<-function(){
png(file.path(root_dir,"plots","tissue_BAR.png"))
metaBar(meta,"SMTS")
dev.off()
png(file.path(root_dir,"plots","age_BAR.png"))
metaBar(meta,"AGE")
dev.off()

mat<-get_PCAmat(lcpm)
plotmanual(t(mat),meta$AGE,title="All samples")
ggsave(file.path(root_dir,"plots","All_age_PCA.png"))

for(t in unique(meta$SMTS)){
  tissuePCA(lcpm,meta,t,"SMTSD")
}

for(t in unique(meta$SMTS)){
  tissuePCA(lcpm,meta,t,"AGE")
}
}
#EDA()

# Writing per-tissue tsv

keep<-!meta$AGE==""
filtered_meta<-meta[keep,]
filtered_lcpm<-lcpm[,keep]
filtered_cpm<-cpm[,keep]
tissueSubset<-function() {
  foreach(t=unique(filtered_meta$SMTS)) %do%
    {
      sel<-filtered_meta$SMTS==t
      ldat<-filtered_cpm[,sel]
      dat<-filtered_lcpm[,sel]
      ldat<-data.table::data.table(t(ldat))
      dat<-data.table::data.table(t(dat))
      rownames(ldat)<-filtered_meta[filtered_meta$SMTS==t,]$SAMPID
      rownames(dat)<-filtered_meta[filtered_meta$SMTS==t,]$SAMPID
      data.table::fwrite(ldat,file = file.path(data_dir,"tissue-specific",paste0(t,"_lcpm.tsv")),sep = "\t",row.names = TRUE,nThread = 6)
      data.table::fwrite(dat,file = file.path(data_dir,"tissue-specific",paste0(t,"_cpm.tsv")),sep = "\t",row.names = TRUE,nThread = 6)
      
    }
}
tissueSubset()
## Optional saving of pre-processed matrices
#tlcpm<-data.table(t(lcpm))
#tcpm<-data.table(t(cpm))
#rownames(tlcpm)<-dimnames(lcpm)[[2]]
#rownames(tcpm)<-dimnames(tcpm)[[2]]
#data.table::fwrite(tlcpm,file = file.path(data_dir,"lcpm.tsv"),sep = "\t",row.names = TRUE,nThread = nCore)
#data.table::fwrite(tcpm,file = file.path(data_dir,"cpm.tsv"),sep = "\t",row.names = TRUE,nThread = nCore)
