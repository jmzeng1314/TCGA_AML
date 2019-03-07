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

# https://tcga.xenahubs.net/download/TCGA.LAML.sampleMap/LAML_clinicalMatrix.gz

# samples 200
# version 2016-04-27
# type of data  phenotype
# raw datahttps://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/bcr/
#   

rm(list = ls())
options(stringsAsFactors = F)
a=read.table('data/LAML_clinicalMatrix.gz',header = T,sep = '\t')
colnames(a)
grep('code',colnames(a))
# table(a[,45])
# a=a[,c("submitter_id.samples",
#        "leukemia_french_american_british_morphology_code",
#        "vital_status.diagnoses",
#        "days_to_death.diagnoses" 
#        )]
colnames(a)
table(a[,64])
a=a[,c(1,64,7,6)]
colnames(a)=c('id','fab','status','os')
a$id=gsub('-','.',a$id)
head(a)

# https://tcga.xenahubs.net/download/TCGA.LAML.sampleMap/miRNA_GA_gene.gz

# miRNA mature strand expression RNAseq

# samples 188
# version 2017-09-08
# type of data  miRNA mature strand expression RNAseq
# unit  log2(RPM+1)
# platform  IlluminaGA_miRNASeq
# author  British Columbia Cancer Agency TCGA genome characterization center
# raw data  https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/cgcc/bcgsc.ca/illuminaga_mirnaseq/mirnaseq/
#   

b=read.table('data/miRNA_GA_gene.gz',header = T,sep = '\t')
head(colnames(b))
b[1:4,1:4]
rownames(b)=b[,1]
b=b[,-1]
## too many NA values. 

substring(colnames(b),14,15)
clin=a[match(colnames(b),a$id),]
exprSet=b
table(clin$fab)
save(exprSet,clin,file = 'TCGA_AML_input.Rdata')









