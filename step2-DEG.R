## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2019-01-25 15:37:51
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2019-01-25 First version
### Update Log: 2019-01-25 second version
###
### ---------------

### https://github.com/jmzeng1314/GEO/blob/master/airway/DEG_rnsseq.R

rm(list=ls())
options(stringsAsFactors = F)
source('functions.R')
load(file = 'TCGA_AML_input.Rdata')
exprSet[1:4,1:4]
exprSet=2^exprSet-1
exprSet[1:4,1:4]

table(clin$status)
group_list=clin$status

### ---------------
###
### Firstly run DESeq2 
###
### ---------------

if(F){
  library(DESeq2)
  
  (colData <- data.frame(row.names=colnames(exprSet), 
                         group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  tmp_f=file.path(Rdata_dir,'TCGA-KIRC-miRNA-DESeq2-dds.Rdata')
  if(!file.exists(tmp_f)){
    dds <- DESeq(dds)
    save(dds,file = tmp_f)
  }
  load(file = tmp_f)
  res <- results(dds, 
                 contrast=c("group_list","tumor","normal"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG =as.data.frame(resOrdered)
  DESeq2_DEG = na.omit(DEG)
  
  nrDEG=DESeq2_DEG[,c(2,6)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')  
  draw_h_v(exprSet,nrDEG,'DEseq2',group_list,1)
}

### ---------------
###
### Then run edgeR 
###
### ---------------
if(T){
  library(edgeR)
  d <- DGEList(counts=exprSet,group=factor(group_list))
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  dge=d
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  edgeR_DEG =nrDEG 
  nrDEG=edgeR_DEG[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'edgeR',group_list,1)
  
}


### ---------------
###
### Lastly run voom from limma
###
### --------------- 
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('dead-alive'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='dead-alive', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'limma',group_list,1)
  
}

 

nrDEG1=DEG_limma_voom[,c(1,4)]
colnames(nrDEG1)=c('log2FoldChange','pvalue') 

nrDEG2=edgeR_DEG[,c(1,5)]
colnames(nrDEG2)=c('log2FoldChange','pvalue') 

nrDEG3=DESeq2_DEG[,c(2,6)]
colnames(nrDEG3)=c('log2FoldChange','pvalue')  

mi=unique(c(rownames(nrDEG1),rownames(nrDEG1),rownames(nrDEG1)))
lf=data.frame(lf1=nrDEG1[mi,1],
              lf2=nrDEG2[mi,1],
              lf3=nrDEG3[mi,1])




