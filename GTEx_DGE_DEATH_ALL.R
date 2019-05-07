library(plyr)
library(edgeR)
library(limma)
library(Glimma)
library(stringr)

### Params
biomart_file_path=file.path("annotation","2018-04-12_biomart_ensemblGene_hgnc_uniprot.tab")
uniprot_meta_path=file.path("annotation","uniprot_meta.tsv")
nCore<-6
data_dir<-"data"

### common functions
filterBiomart <- function(infile) {
  biomart1=read.table(infile, header=T, sep='\t', stringsAsFactors=FALSE,na.strings = "")
  names(biomart1)=c("Gene","hgnc","uniprot")
  #remove duplicated ensembl
  biomart=biomart1[!duplicated(biomart1$Gene),]
  #biomart$desc=sub(" [[].*","",biomart$desc)
  return(biomart)
}
get_counts<-function(raw_frame) {
  ### Different than same name function in GTEx_input.R
  gene_names<-raw_frame$Name
  counts<-as.matrix(raw_frame[,-c(1:2)])
  dimnames(counts)[[1]]<-gene_names
  return(counts)
}

### loading data
meta<-data.table::fread(file.path(data_dir,"merged_meta.tsv"),sep="\t",header=TRUE,na.strings="na",nThread=nCore)
counts<-data.table::fread(file.path(data_dir,"All_Tissue_Site_Details.combined.reads.gct"),nThread=nCore)

### process data and metadata
counts<-get_counts(counts)
rownames(counts)<-stringr::str_extract(rownames(counts),"[:alpha:]+[:digit:]+") # removes trailing ensembl number.
rownames(meta)<-meta$SAMPID
meta$DTHHRDY<-factor(meta$DTHHRDY)
levels(meta$DTHHRDY)<-list(Vent=c("0"),Viol_Fast=c("1"),Fast_Nat=c("2"),Medium=c("3"),Slow=c("4"),None=c(NA))
### edgeR pre-processing
counts<-DGEList(counts,samples=meta)
counts <- calcNormFactors(counts, method = "TMM")

### DEG Analysis
# Remove tissues with too few samples.
tiss_table<-table(meta$SMTS)>200
DEG_tissues<-names(tiss_table[tiss_table==TRUE])
sel<-meta$SMTS %in% DEG_tissues
meta<-meta[sel,]
counts<-counts[,sel]
keep<-!is.na(meta$DTHHRDY)
DEG_meta<-meta[keep,]
DEG_counts<-counts[,keep]
DEG_meta$DTHHRDY<-droplevels(DEG_meta$DTHHRDY)
#rm(counts)

print(paste0("Input: ",nrow(DEG_counts)))
keep.exprs <- filterByExpr(DEG_counts$counts,group=factor(DEG_meta$DTHHRDY),min.count = 15) #first thing the function does is convert counts to matrix. MinSamples=10+(n-10)*.7
DEG_counts <- DEG_counts[keep.exprs,, keep.lib.sizes=FALSE]
print(paste0("Output: ",nrow(DEG_counts)))

DEG_lcpm<-cpm(DEG_counts,log = TRUE)

## DEG
design<-model.matrix(~0+DTHHRDY,data=DEG_meta)
colnames(design)<-str_remove(colnames(design),"DTHHRDY")
contr.matrix <- makeContrasts(
  FastvVent = Vent-Fast_Nat,
  levels = colnames(design))
png(filename = file.path("DGE_plots","DEATH","Init_SA.png"))
v <- voom(DEG_counts, design, plot=TRUE)
dev.off()
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
png(filename = file.path("DGE_plots","DEATH","Final_SA.png"))
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit,p.value = .05)
png(file.path("DGE_plots","DEATH","MDplot_DEATH_all.png"))
plotMD(tfit,status = dt)
dev.off()
## Positive FC is enrichment in vent deaths
## Negative FC is enrichment in natural, fast deaths


### Getting gene annotation
biomart<-filterBiomart(biomart_file_path)
## one time uniprot file creation for uniprot_topology_parser.py
#uniprot<-biomart[!is.na(biomart$uniprot),]$uniprot
#write.table(uniprot,"uniprot_ids.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
##
uniprot_meta<-data.table::fread(uniprot_meta_path,header=FALSE,sep = "\t")
names(uniprot_meta)<-c("uniprot","name","desc","sub","topology","func")
biomart<-join(biomart,uniprot_meta,by = "uniprot")
genes<-join(data.frame(Gene=rownames(DEG_counts),stringsAsFactors = FALSE),biomart)
genes$func<-NULL
genes$name<-NULL
genes$uniprot<-NULL
genes<-genes[!duplicated(genes$Gene),]
names(genes)[1]<-"ENSEMBL"

## Final plot and placement into directories
DGE_results<-topTreat(tfit, coef=1, n=Inf)
DGE_results$Gene<-rownames(DGE_results)
DGE_results<-join(DGE_results,biomart)
DGE_results<-DGE_results[!duplicated(DGE_results$Gene),]
rownames(DGE_results)<-DGE_results$Gene
data.table::fwrite(DGE_results,file=file.path("DGE_plots","DEATH","DGE_results.tsv"),sep="\t",row.names = TRUE)
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],folder=file.path("DGE_plots","DEATH"),
         side.main="hgnc", counts=DEG_lcpm, groups=DEG_meta$DTHHRDY,anno=genes,launch = FALSE)