# Remapping Reads

Remapping reads to PacBio assembled genomes

## Remapping

### Prerequisites

Programs used

```
bowtie2 v2.3.4.1
samtools v1.8
blasr (for pacbio reads)
bamCoverage (from deeptools package:https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html )
```
### For paired-end illumina reads

Build bowtie2 index
```
bowtie2-build /home/Pam/submissions/repeats/Nara/Nara.fsa Nara
```
Perform read mapping. For 100 bp paired-end reads of 500 bp insert sized libraries
```
bowtie2 -x Nara --end-to-end --fr -p 32 -X 700 -I 400 --very-sensitive -1 /home/Pam/submissions/remapping/fqlinks/$1_500_R1_val_1.fq.gz -2 /home/Pam/submissions/remapping/fqlinks/$1_500_R2_val_2.fq.gz -S $1.500.Nararef.bt2.sam
samtools view -b -S -@ 32 $1.500.Nararef.bt2.sam| samtools sort -@ 32 -o $1.500.Nararef.bt2.sort.bam
samtools index $1.500.Nararef.bt2.sort.bam

```
Where $1 = strain name

### For RPKM normalized coverage values

```
bamCoverage --bam $1.500.Nara.sort.bam --outFileName $1.bcdt.rpkm.bedgraph --outFileFormat bedgraph -p 32 --binSize 10000 --normalizeUsing RPKM --effectiveGenomeSize 59598044
```

where $1 = strain name

To draw heatmaps, seaborn was used. See [Remapping jupyter notebook](jupyter_notebooks/merge_bcdt_rpkm_output_final.ipynb)

### For pacbio reads

Run blasr
```
blasr /home/share/r54244_20180912_055659/2_B01/m54244_180912_164749.subreadset.xml /home/Pam/submissions/repeats/Nara/Nara.fsa --out Nara.Nararef.PB.bam --bam --nproc 32 --bestn 1
```
