# RNAseq analysis

RNAseq analysis of Illumina HiSeq 50 bp single-end reads.

### Prerequisites
Programs used
```
STAR(v2.6.0a)
Rsubread(v1.32.2)
edgeR(v3.24.3)
```
Strawberry reference sequence
```
#Genome sequence
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_x_ananassa/Fragaria_x_ananassa_Camarosa_Genome_v1.0.a1/assembly/F_ana_Camarosa_6-28-17_hardmasked.fasta.gz

#GFF
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_x_ananassa/Fragaria_x_ananassa_Camarosa_Genome_v1.0.a1/genes/Fxa_v1.2_makerStandard_MakerGenes_woTposases.gff.gz
```
### Mapping RNAseq reads
Build STAR index
```
STAR --runThreadN 24 --runMode genomeGenerate --genomeFastaFiles F_ana_Camarosa_6-28-17_hardmasked.fasta --sjdbOverhang 49 --sjdbGTFfile Fxa_v1.2_makerStandard_MakerGenes_woTposases.gff --genomeDir Fxa_Camarosa --sjdbGTFtagExonParentTranscript=Parent
```
Run STAR to map reads to genome
```
STAR --genomeDir Fxa_Camarosa --runThreadN 8 --readFilesCommand gunzip -c --readFilesIn $1.fastq.gz --outFileNamePrefix $1.Fxa.STAR
```
$1 = read filename

### Obtaining read counts
Read counts with Rsubread in R
```
library(Rsubread)
#Read in STAR sam files from samlist.txt file
files <-readLines("samlist.txt")
df<-featureCounts(files=files,annot.ext="Fxa_v1.2_makerStandard_MakerGenes_woTposases.gff",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="Parent",nthreads=32)
write.csv(df$counts,file="20190423_Nara_413_Fxa_subreads_genes.csv")
```
### Identifying differentially expressed genes
Identifying differentially expressed genes in EdgeR
Parsed the output file from Rsubread so that gene IDs match new gene IDs (ftp://ftp.bioinfo.wsu.edu/species/Fragaria_x_ananassa/Fragaria_x_ananassa_Camarosa_Genome_v1.0.a1/functional/FxaC_newGeneID.txt.gz)
In R
```
raw <-read.csv("rawsubreadcounts.parsed.txt",sep="\t",row.names=1)
library(edgeR)

#Filter by expression

y<-DGEList(counts=raw,group=x)
keep <-filterByExpr(y)
y<-y[keep, , keep.lib.sizes=FALSE]

#TMM normalization

y<-calcNormFactors(y,method="TMM")

design <-model.matrix(~0+x)
y<-estimateDisp(y,design)
fit<-glmQLFit(y,design)
logcpm <-cpm(y,prior.count=2,log=TRUE)

#Set up contrasts
my.contrasts<-makeContrasts(
v1=xNara.1-xC413.1,
v3=xNara.3-xC413.3,
v6=xNara.6-xC413.6,
levels=design)

contrastlist<-c("v1",v3","v6")
lapply(contrastlist,function(df){
  qlf<-glmTreat(fit,contrast=my.contrasts[,df],lfc=1)
  deonly<-sum(head(summary(decideTests(qlf)),1),tail(summary(decideTests(qlf)),1))
  write.csv(topTags(qlf,n=deonly),paste("20190423-Fxa_STAR_edgeR_",df,".csv",sep=''))
})
```
For fungal genes

[C. fructicola Nara gc5 reference GFF](../annotation/Nara.gff)

## References

*  STAR user manual:http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
*  Rsubread user manual:https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
*  EdgeR user manual:https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
