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

### edgeR pre-processing
counts<-DGEList(counts,samples=meta)
counts <- calcNormFactors(counts, method = "TMM")

### DEG Analysis
# Remove tissues with too few samples.
tiss_table<-table(meta$SMTS)>200
DEG_tissues<-names(tiss_table[tiss_table==TRUE])
for(TISSUE in DEG_tissues) {
dir.create(path = file.path("DGE_plots",TISSUE),recursive = TRUE,showWarnings = FALSE)
keep<-meta$AGE!="" & meta$SMTS==TISSUE
DEG_meta<-meta[keep,]
DEG_counts<-counts[,keep]
#rm(counts)

### Changing age to numeric
DEG_meta$AGE<-paste0("A",as.numeric(factor(DEG_meta$AGE))) # convert to friendly group
DEG_meta$SMTS<-gsub(" ","",DEG_meta$SMTS) # still relevant
print(paste0("Input: ",nrow(DEG_counts)))
keep.exprs <- filterByExpr(DEG_counts$counts,group=factor(DEG_meta$AGE),min.count = 10) #first thing the function does is convert counts to matrix. MinSamples=10+(n-10)*.7
DEG_counts <- DEG_counts[keep.exprs,, keep.lib.sizes=FALSE]
print(paste0("Output: ",nrow(DEG_counts)))

DEG_lcpm<-cpm(DEG_counts,log = TRUE)

### DGE
design<-model.matrix(~0+AGE,data=DEG_meta)
colnames(design)<-str_remove(colnames(design),"AGE")
contr.matrix <- makeContrasts(
  LAvHA= (A1+A2+A3)/3-(A4+A5+A6)/3,
  levels = colnames(design))
png(filename = file.path("DGE_plots",TISSUE,"Init_SA.png"))
v <- voom(DEG_counts, design, plot=TRUE)
dev.off()
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
png(filename = file.path("DGE_plots",TISSUE,"Final_SA.png"))
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit,p.value = .05)
png(file.path("DGE_plots",TISSUE,paste0(TISSUE,"MDplot.png")))
plotMD(tfit,status = dt)
dev.off()

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

### Final Glimma plot, writing DGE results and placement into directories

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],folder=file.path("DGE_plots",TISSUE),
         side.main="hgnc", counts=DEG_lcpm, groups=DEG_meta$AGE,anno=genes,launch = FALSE)
glMDSPlot(DEG_lcpm, groups=DEG_meta[,c("SMTS","SMTSD","AGE","SEX","DTHHRDY")], main=paste0(TISSUE," MDS Plot"), folder=file.path("DGE_plots",TISSUE),launch=FALSE)
data.table::fwrite(topTreat(tfit,coef=1,n=Inf),file = file.path("DGE_plots",TISSUE,paste0(TISSUE,"DGE_results.tsv")),sep = "\t",row.names = TRUE)
}