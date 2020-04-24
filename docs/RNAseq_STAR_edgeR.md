# RNAseq analysis

RNAseq analysis of Illumina HiSeq 50 bp single-end reads.

### Prerequisites
Programs used
```
STAR(v2.6.0a)
Rsubread(v1.32.2)
edgeR(v3.24.3)

```
### Mapping RNAseq reads
Build STAR index
```
STAR --runThreadN 32 --runMode genomeGenerate --genomeFastaFiles genome.fsa --sjdbOverhang 49 --sjdbGTFfile genome.gff --sjdbGTFtagExonParentTranscript=Parent --genomeDir Nara
```
Run STAR to map reads to genome
```
STAR --genomeDir Nara --runThreadN 8 --readFilesCommand gunzip -c --readFilesIn /home/Pam/S1_braker/RNAseq/$1.fastq.gz --outFileNamePrefix $1.Nara.STAR --alignIntronMax 1000
```
$1 = read filename

### Obtaining read counts
Read counts with Rsubread in R
```
library(Rsubread)
#Read in STAR sam files from samlist.txt file
files <-readLines("samlist.txt")
counts<-featureCounts(files=files,annot.ext="genome.gff",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="Parent",nthreads=64,primaryOnly=TRUE,strandSpecific=2)
write.csv(counts$counts,file="Nara.STAR.subreadscount.csv")
```
### Identifying differentially expressed genes
Identifying differentially expressed genes in EdgeR

In R
```
raw <-read.csv("Nara.STAR.subreadscount.csv",sep=",",row.names=1)
x<-as.factor(substring(colnames(raw),1,2))
y<-DGEList(counts=raw,group=x)
keep <-filterByExpr(y)
y<-y[keep, , keep.lib.sizes=FALSE]

#TMM normalization

y<-calcNormFactors(y,method="TMM")

plotMDS(y,col=colors)
design <-model.matrix(~0+x)
y<-estimateDisp(y,design)
fit<-glmQLFit(y,design)
#Setup contrasts to run
my.contrasts<-makeContrasts(
Nara1PV=xD1-xVH,
Nara3PV=xD3-xVH,
Nara6PV=xD6-xVH,
NaraRPV=xRT-xVH,
levels=design)
contrastlist<-c("Nara1PV","Nara3PV","Nara6PV","NaraRPV")
lapply(contrastlist,function(df){
  qlf<-glmTreat(fit,contrast=my.contrasts[,df],lfc=1)
  #summary(decideTests(qlf))
  deonly<-sum(head(summary(decideTests(qlf)),1),tail(summary(decideTests(qlf)),1))
  write.csv(topTags(qlf,n=deonly),paste(df,".csv",sep=''))
})
```
[C. fructicola Nara gc5 reference GFF](../annotation/Nara.gff)

## References

*  STAR user manual:http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
*  Rsubread user manual:https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
*  EdgeR user manual:https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
