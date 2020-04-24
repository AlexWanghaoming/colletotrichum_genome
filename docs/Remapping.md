# Remapping Reads

Remapping reads to PacBio assembled genomes

## Remapping

### Prerequisites

Programs used

```
bowtie2 v2.3.4.1
samtools v1.8
bedtools v2.29.2
ngmlr v0.2.7
```
### For paired-end illumina reads

Build bowtie2 index
```
bowtie2-build Nara.fsa Nara
```
Perform read mapping. For 100 bp paired-end reads of 500 bp insert-sized libraries pretrimmed with trimgalore
```
bowtie2 -x Nara --end-to-end --fr -p 32 -X 700 -I 400 --very-sensitive -1 $1_500_R1_val_1.fq.gz -2 $1_500_R2_val_2.fq.gz -S $1.500.Nararef.bt2.sam
samtools view -b -S -@ 32 $1.500.Nararef.bt2.sam| samtools sort -@ 32 -o $1.500.Nararef.bt2.sort.bam
samtools index $1.500.Nararef.bt2.sort.bam

#Remove duplicate reads
samtools sort -n -o $1.500.Nararef.bt2.nsort.bam $1.500.Nararef.bt2.sort.bam
samtools fixmate -@ 8 -m $1.500.Nararef.bt2.nsort.bam $1.500.Nararef.bt2.fm.bam
samtools sort -@ 8 -o $1.500.Nararef.bt2.sort.fm.bam $1.500.Nararef.bt2.fm.bam
samtools markdup -@ 8 -r $1.500.Nararef.bt2.sort.fm.bam $1.500.Nararef.bt2.sort.md.bam
samtools index $1.500.Nararef.bt2.sort.md.bam
```
Where $1 = strain name

### For PacBio reads
```
zcat $1.PB.fastq.gz|ngmlr -t 8 -r /home/Pam/submissions/repeats/Nara/Nara.fsa > $1.PB.Nara.ng.sam
samtools view -@ 8 -bS $1.PB.Nara.ng.sam > $1.PB.Nara.ng.bam
samtools sort -@ 8 -o $1.PB.Nara.ng.sort.bam $1.PB.Nara.ng.bam
rm $1.PB.Nara.ng.sam
rm $1.PB.Nara.ng.bam
```
### For read depths

Make windows of 10 kb in size throughout the whole genome
```

bedtools makewindows -g Nara.genome -w 10000 > Nara.10k.windows.bed
```
Count number of reads per 10 kb
```
bedtools coverage -a Nara.10k.windows.bed -b $1.sort.md.bam > $1.Nara.10k.bedtoolscoverage.bed
```
where $1 = strain name

For plots, see [Remapping jupyter notebook](../jupyter_notebooks/Plots_of_depths_coverage_densities.ipynb)
