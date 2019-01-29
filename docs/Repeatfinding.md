# Custom Repeat Library Finding Pipeline

Repeat identification for *Colletotrichum fructicola, Colletotrichum siamense, Colletotrichum aenigma* and *Colletotrichum higginsianum*.

## Running the pipeline

### Prerequisites

Programs used

```
RepeatModeler version open-1.0.11
RepeatMasker v 1.332
rmblastn 2.6.0+ (installed via conda)
TransposonPSI (http://transposonpsi.sourceforge.net/)
LTRpred (https://github.com/HajkD/LTRpred)
LTR_retriever (https://github.com/oushujun/LTR_retriever)
Genometools v1.5.7 (https://github.com/genometools)
samtools
bedtools
R v3.5.1
usearch v11.0.667_i86linux32

```
Repeat libraries used
```
RepBase23.12
```
Custom scripts used are included in the folder "custom_scripts". Prerequisites for the scripts are python3, biopython and (sometimes) pandas which can be installed via conda (https://anaconda.org/anaconda/)

Before running these custom scripts, run the following in the shell

```
export CUSTOM_SCRIPTS=<location of custom scripts>
```


## Finding repeats

### 1. RepeatModeler

run_repeatmodeler.sh

```
BuildDatabase -name $1 -engine ncbi $1.fsa
RepeatModeler -database $1
```
Example: To run in folder containing genome fasta file Nara.fsa
```
chmod 755 run_repeatmodeler.sh
./run_repeatmodeler Nara
```
### 2. Run ltr_retriever

run_ltr_retriever.sh
```
gt suffixerator -db $1.fsa -indexname $1.fsa -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index $1.fsa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 > $1.harvest.scn
ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 $1.fsa > $1.finder.scn
gt ltrharvest -index $1.fsa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 > $1.harvest.TGCA.scn
LTR_retriever -genome $1.fsa -inharvest $1.harvest.TGCA.scn -infinder $1.finder.scn -nonTGCA $1.harvest.scn
```
Example: To run in folder containing genome fasta file Nara.fsa
```
chmod 755 run_ltr_retriever.sh
./run_ltr_retriever.sh Nara
```
### 3. Run LTRpred

run_LTRpred.R
```
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("LTRpred")
genomename=paste(args[1],".fsa")
ltrpredresults <- paste(args[1],"_ltrpred", sep="")
ltrpreddatasheet <-paste(args[1],"_LTRpred_DataSheet.tsv",sep="")
strainfile<-file.path(ltrpredresults,ltrpreddatasheet)
strain_LTRpred<-read.ltrpred(strainfile)
functional<-quality.filter(strain_LTRpred,sim=97,n.orfs=1,strategy="stringent")
outfilename<-file.path(ltrpredresults,args[1],",sep="")

pred2fasta(LTRpred.tbl=functional,prediction.file=,program="LTRpred")
```
Example run on genome fasta file Nara.fsa
```
Rscript ./run_LTRpred.R Nara.fsa
```
### 4. Run TransposonPSI
```
perl transposonPSI.pl Nara.fsa nuc
```
Pull out fasta sequences related to besthits. Run custom script parse_gff_to_fasta.py

Example:
```
python parse_gff_to_fasta.py -fasta Nara.fsa -gff Nara.fsa.TPSI.allHits.chains.bestPerLocus.gff3 -ID 1 -feature translated_nucleotide_match
```

## Combining and filtering the predictions

Combine TransposonPSI, LTRpred and LTR_retriever outputs, filter out sequences that are less than 400 bp in length and have

Create BLAST databases for sequences to be used
```
makeblastdb -in Nara.fsa -out ~/db/genome/Nara -dbtype nucl
makeblastdb -in CF413.fsa -out ~/db/genome/CF413 -dbtype nucl
makeblastdb -in Cg363.fsa -out ~/db/genome/Cg363 -dbtype nucl
makeblastdb -in Cg56.fsa -out ~/db/genome/Cg56 -dbtype nucl
makeblastdb -in Colhig2_AssemblyScaffolds.fasta -out ~/db/genome/Chig_IMI_JGI -dbtype nucl

makeblastdb -dbtype prot -in ~/bin/RepeatMasker/RepeatMasker/Libraries/RepeatPeps.lib -out RepeatPeps.db

```
Run the following commands for all outputs for a given strain.

```
#!/bin/bash
#combine outputs
cat $1.ltretriever.fa  $1.ltrpred.fa  $1.TPSI.fasta > $1.cat3.fa
samtools faidx $1.cat3.fa

#Filter out sequences less than 400 bp in length and rename sequences
awk '$2>=400{print $1}' $1.cat3.fa.fai > $1.c3
python $CUSTOM_SCRIPTS/pull_fasta.py $1.cat3.fa $1.c3 > $1.cat3.400.fasta
awk '{print $1"\t"FILENAME"_"NR}' $1.c3 > $1.cat3.400.rn.list
python $CUSTOM_SCRIPTS/rename.py $1.cat3.400.rn.list $1.cat3.400.fasta > $1.cat3.400.rn.fasta

#Perform blastx. Change -db depending on location of RepeatPeps database
blastx -db /home/share/DB/RepeatPeps.db -num_threads 16 -evalue 1E-5 -outfmt 6 -out $1.blastx.repbase.out -query $1.cat3.400.rn.fasta

#Perform blastn against each genome
blastn -db ~/db/genome/Nara -num_threads 16 -evalue 1E-15 -outfmt 6 -out $1.blastn.Nara.out -query $1.cat3.400.rn.fasta
blastn -db ~/db/genome/CF413 -num_threads 16 -evalue 1E-15 -outfmt 6 -out $1.blastn.CF413.out -query $1.cat3.400.rn.fasta
blastn -db ~/db/genome/Cg363 -num_threads 16 -evalue 1E-15 -outfmt 6 -out $1.blastn.Cg363.out -query $1.cat3.400.rn.fasta
blastn -db ~/db/genome/Cg56 -num_threads 16 -evalue 1E-15 -outfmt 6 -out $1.blastn.Cg56.out -query $1.cat3.400.rn.fasta
blastn -db ~/db/genome/Chig_IMI_JGI -num_threads 16 -evalue 1E-15 -outfmt 6 -out $1.blastn.Chig_IMI_JGI.out -query $1.cat3.400.rn.fasta

#Save all outputs to a new folder
mkdir output_$1
mv $1* output_$1
```
## Filter outputs
### Filter outputs by length, number of occurrences in genomes and homology to repbase sequence

Using the BLAST outputs from the above step, repeat sequences were filtered to include only sequences with more than five BLASTn hits in any single genome (E-value cutoff 1E-15) and/or with significant homology to any sequence in the RepBase peptide database (BLASTx E-value cutoff 1E-5).

Create a fasta file of the filtered sequences where $1.compiled.out.list is the list of sequences that meet the above criteria.
```
python $CUSTOM_SCRIPTS/pull_fasta.py output_$1/$1.cat3.400.rn.fasta $1.compiled.out.list > $1.compiled.out.fasta
```
Combine sequences from all strains using usearch
```
cat CF413.compiled.out.fasta Cg363.compiled.out.fasta Cg56.compiled.out.fasta CH.compiled.out.fasta Nara.compiled.out.fasta > CGSChig.compiled.out.fasta
usearch -sortbylength CGSChig.compiled.out.fasta --fastaout CGSChig.merged.sorted.fasta --log usearch.log
usearch -cluster_fast CGSChig.merged.sorted.fasta --id 0.8 --centroids CGSChig.centroids.fa --uc CGSChig.result.uc -consout CGSChig.nr.cons.fa -msaout CGSChig.aligned.fa --log usearch2.log
```
###Annotate custom library
RepeatClassifer from the RepeatModeler package was used to annotate the custom library
```
RepeatClassifier -consensi CGSChig.nr.cons.fa -engine ncbi
```
In addition, transposonPSI was used to annotate the library.
```
perl transposonPSI.pl CGSChig.nr.cons.fa nuc
```
### Filter out protein coding sequences
Some of the sequences that are not classified (Unknown) could include protein sequences with repeats in them. To get rid of these, our repeats were filtered against coding sequences that are not related to transposable elements.

The predicted proteomes of our genomes of interest were used as queries against the RepeatPep database.

```
blastp -query $1 -db /home/share/DB/RepeatPeps.db -outfmt 6 -max_target_seqs 25 -culling_limit 2 -num_threads 8 -evalue 1e-5 -out $1.vs.RepeatPeps.25cul2.1e5.blastp.out
```
From the BLASTp output, a list of proteins in our predicted proteomes with homology to RepBase peptides can be obtained
```
awk '{print $1}' CGSChig.cat5.vs.RepeatPeps.25cul2.1e5.blastp.out |sort | uniq > CGSChig.cat5.exclude.list
```
Then, our custom library was searched for homology to CDS sequences corresponding to fungal proteins without homology to TEs.
```
#Get CDS sequences without homology to TEs
python $CUSTOM_SCRIPTS/exclude_fasta.py CGSChig.cat5.cds.fasta CGSChig.cat5.exclude.list > CGSChig.cat5.cds.noTEs.fasta
#BLASTn against custom library
blastn -task megablast -query CGSChig.nr.cons.fa.classified -db CGSChig.cat5.cds.noTEs.fasta -outfmt 6 -max_target_seqs 25 -culling_limit 2 -num_threads 32 -evalue 1e-10 -out CGSChig.cat5.repeatmodeller.lib.vs.noTES.25cu2.1e10.megablast.out
```
Sequences had a BLASTn hit against the CDS sequences and which were not annotated by RepeatClassifier (Unknown) or TransposonPSI were excluded from the final repeat library.

## Annotation of Genomes with Custom Repeat Library
RepeatMasker was used to annotate each genome
```
RepeatMasker -lib CGSChig.nr.cons.cds.filtered.fa Nara.fsa -pa 16 -dir ./ -gff
```
## Acknowledgments
The following resources were referred to while constructing this pipeline:

*  https://blaxter-lab-documentation.readthedocs.io/en/latest/repeat-masking.html
*  Coghlan et al. 2018 Protocol Exchange https://doi:10.1038/protex.2018.054
*  Castanera et al. 2016 PLOS Genetics https://doi.org/10.1371/journal.pgen.1006108
*  The very helpful Avrilomics blog by Avril Coghlan (http://avrilomics.blogspot.com)

## Authors

* **Pamela Gan** - *Initial work* - [Colletotrichum](https://)


## License
