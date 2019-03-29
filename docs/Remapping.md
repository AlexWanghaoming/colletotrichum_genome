# Remapping Reads

Remapping reads to PacBio assembled genomes

## Remapping

### Prerequisites

Programs used

```
bowtie2 v2.3.4.1
samtools v1.8
blasr
bamCoverage (from the deeptools 3.1.3 package:https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
unique-kmers.py (from khmer version 3.0.0a2)
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

### For RPKM normalized read depth values

Estimate Effective Genome Size
```
unique-kmers Nara.fsa
```
Calculate normalized read depths
```
samtools index $1.500.Nararef.bt2.sort.md.bam -@ 32
bamCoverage --bam $1.500.Nararef.bt2.sort.md.bam --outFileName $1.bcdt.rpkm.bedgraph --outFileFormat bedgraph -p 32 --binSize 10000 --normalizeUsing RPKM --effectiveGenomeSize 57859888
```
where $1 = strain name

To draw heatmaps, see [Remapping jupyter notebook](../jupyter_notebooks/merge_bcdt_rpkm_output_final.ipynb)

## Acknowledgments

The following resources were referred to while constructing this pipeline:
