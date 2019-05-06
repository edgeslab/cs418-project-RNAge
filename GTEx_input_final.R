
#Running this script requires the following:
# Plyr, edgeR, limma, Glimma, data.table, stringr, foreach, ggplot2
library(plyr)
library(edgeR)
library(limma)
library(Glimma)
library(ggplot2)

### Params
nCore<-6
data_dir<-"data"
plot_dir<-"Final_source_code"

### Functions
get_counts<-function(raw_frame) {
  gene_names<-raw_frame$Name
  counts<-as.matrix(raw_frame[,-c(1:2)])
  dimnames(counts)[[1]]<-gene_names
  return(t(counts))
}
get_PCAmat<-function(lcpm) {
  #mat will need to be horizontal.
  ntop <- 500
  rv <- matrixStats::rowVars(lcpm)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- lcpm[select, ]
  return(mat)
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
    p + geom_point(alpha=alpha,size=.1) + theme + xlab(percentage[1]) + ylab(percentage[2]) + labs(title = title, subtitle = subtitle)+ guides(colour = guide_legend(override.aes = list(size=2)))
  return(p)
}
### loading data
meta<-data.table::fread(file.path(data_dir,"merged_meta.tsv"),sep="\t",header=TRUE,na.strings="na",nThread=nCore)
counts<-data.table::fread(file.path(data_dir,"All_Tissue_Site_Details.combined.reads.gct"),nThread=nCore)

### process data and metadata
counts<-get_counts(counts)
colnames(counts)<-stringr::str_extract(colnames(counts),"[:alpha:]+[:digit:]+")
rownames(meta)<-meta$SAMPID

### edgeR pre-processing
counts<-DGEList(t(counts),samples=meta)
counts <- calcNormFactors(counts, method = "TMM")
cpm<-cpm(counts)
lcpm<-cpm(counts,log=TRUE)
tiss_table<-table(meta$SMTS)>200
sel_tissues<-names(tiss_table[tiss_table==TRUE])
sel<-meta$SMTS %in% sel_tissues
mat<-get_PCAmat(lcpm[,sel])
png(filename = file.path(plot_dir,"All_SMTS_PCA.png"),width=700,res = 100)
plotmanual(t(mat),meta$SMTS[sel],title="PCA of all samples by Tissue Origin")
dev.off()

### Creating per tissue data
dir.create(path = file.path(data_dir,"tissue-specific"),showWarnings = F)
keep<-!meta$AGE==""
meta<-meta[keep,]
counts<-counts$counts[,keep]
cpm<-cpm[,keep]
library(foreach)
tissueSubset<-function() {
  foreach(t=unique(filtered_meta$SMTS)) %do%
  {
    sel<-meta$SMTS==t
    cdat<-counts[,sel]
    dat<-cpm[,sel]
    cdat<-data.table::data.table(t(cdat))
    dat<-data.table::data.table(t(dat))
    rownames(cdat)<-filtered_meta[filtered_meta$SMTS==t,]$SAMPID
    rownames(dat)<-filtered_meta[filtered_meta$SMTS==t,]$SAMPID
    data.table::fwrite(cdat,file = file.path(data_dir,"tissue-specific",paste0(t,"_c.tsv")),sep = "\t",row.names = TRUE,nThread = 6)
    data.table::fwrite(dat,file = file.path(data_dir,"tissue-specific",paste0(t,"_cpm.tsv")),sep = "\t",row.names = TRUE,nThread = 6)
    
  }
}
tissueSubset()
system("zip -r data/tissue-specific.zip data/tissue-specific")