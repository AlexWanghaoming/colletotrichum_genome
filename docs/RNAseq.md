# RNAseq analysis

RNAseq analysis of Illumina HiSeq 50 bp single-end reads.

## Remapping

### Prerequisites
Programs used
```
salmon(v0.13.1)
wasabi(v0.3)
sleuth(v0.30.0)
```
For more details see: [conda salmon environment](../conda_environment_yml_files/salmonenv.yml)

Input reference transcriptomes
```
f.ananassa_GDR_reftransV1.fasta (ftp://ftp.bioinfo.wsu.edu/species/Fragaria_x_ananassa/reftransV1/)

```
### RNAseq run
Build salmon index
```
Pam@thesium Strawberry]$ salmon index -t f.ananassa_GDR_reftransV1.fasta -i FANref_index
```
Run salmon to quantify reads
```
salmon quant -i FANref_index -l A -r $1.fastq.gz -p 8 -o quants/$1.quant --validateMappings
```
$1 = unique sample name

Quantify read counts
```
salmon quantmerge --quants quants/413-S-1DPI-1.quant quants/413-S-1DPI-2.quant quants/413-S-1DPI-3.quant quants/413-S-3DPI-1.quant quants/413-S-3DPI-2.quant quants/413-S-3DPI-3.quant quants/413-S-6DPI-1.quant quants/413-S-6DPI-2.quant quants/413-S-6DPI-3.quant quants/413-SH-1.quant quants/413-SH-2.quant quants/413-SH-3.quant quants/413-VH-1.quant quants/413-VH-2.quant quants/413-VH-3.quant quants/Nara-con-1.quant quants/Nara-con-2.quant quants/Nara-con-3.quant quants/Nara-S-1DPI-1.quant quants/Nara-S-1DPI-2.quant quants/Nara-S-1DPI-3.quant quants/Nara-S-3DPI-1.quant quants/Nara-S-3DPI-2.quant quants/Nara-S-3DPI-3.quant quants/Nara-S-6DPI-1.quant quants/Nara-S-6DPI-2.quant quants/Nara-S-6DPI-3.quant quants/Nara-SH-1.quant quants/Nara-SH-2.quant quants/Nara-SH-3.quant quants/Nara-VH-1.quant quants/Nara-VH-2.quant quants/Nara-VH-3.quant quants/Nb-Uninf-1.quant quants/Nb-Uninf-2.quant quants/Nb-Uninf-3.quant quants/Nyo-l-uninf-1.quant quants/Nyo-l-uninf-2.quant quants/Nyo-l-uninf-3.quant quants/Root-2DPI-1.quant quants/Root-2DPI-2.quant quants/Root-2DPI-3.quant quants/S1-N-1DPI-1.quant quants/S1-N-1DPI-2.quant quants/S1-N-1DPI-3.quant quants/S1-N-3DPI-1.quant quants/S1-N-3DPI-2.quant quants/S1-N-3DPI-3.quant quants/S1-N-6DPI-1.quant quants/S1-N-6DPI-2.quant quants/S1-N-6DPI-3.quant quants/S1-S-1DPI-1.quant quants/S1-S-1DPI-2.quant quants/S1-S-1DPI-3.quant quants/S1-S-3DPI-1.quant quants/S1-S-3DPI-2.quant quants/S1-S-3DPI-3.quant quants/S1-S-6DPI-1.quant quants/S1-S-6DPI-2.quant quants/S1-S-6DPI-3.quant quants/S1-SH-1.quant quants/S1-SH-2.quant quants/S1-SH-3.quant quants/S1-VH-1.quant quants/S1-VH-2.quant quants/S1-VH-3.quant quants/Sa-l-uninf-1.quant quants/Sa-l-uninf-2.quant quants/Sa-l-uninf-3.quant quants/Sa-r-uninf-1.quant quants/Sa-r-uninf-2.quant quants/Sa-r-uninf-3.quant --output FANref.salmon.quant.out
