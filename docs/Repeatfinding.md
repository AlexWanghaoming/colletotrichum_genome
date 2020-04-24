# Custom Repeat Library Building Pipeline

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
samtools v1.8 (using htslib 1.8)
bedtools v2.27.1
R v3.5.1
vsearch v2.10.4_linux_x86_64
fasta_formatter (fastx_toolkit 0.0.13:http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)

```
Repeat libraries used
```
RepBase23.12
```
Custom scripts used are included in the folder "custom_scripts". Prerequisites for the scripts are python3, biopython and (sometimes) pandas which can be installed via conda (https://anaconda.org/anaconda/), bedtools and fasta_formatter

Before running these custom scripts, run the following in the shell

```
export CUSTOM_SCRIPTS=<location of custom scripts>
```
## Overview
1. Finding repeats using different programs  
2. Merge repeats from different programs: Filter by length, frequency and homology to known repeat sequences  
3. Filter to exclude sequences with homology to predicted proteome

## 1. Finding repeats

### 1.1 RepeatModeler

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
### 1.2 Run ltr_retriever

run_ltr_retriever.sh
```
gt suffixerator -db $1.fsa -indexname $1.fsa -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index $1.fsa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 > $1.harvest.scn
ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 $1.fsa > $1.finder.scn
gt ltrharvest -index $1.fsa -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 > $1.harvest.TGCA.scn
LTR_retriever -genome $1.fsa -inharvest $1.harvest.TGCA.scn -infinder $1.finder.scn -nonTGCA $1.harvest.scn -noanno
```
Example: To run in folder containing genome fasta file Nara.fsa
```
chmod 755 run_ltr_retriever.sh
./run_ltr_retriever.sh Nara
```
### 1.3 Run LTRpred

run_LTRpred.R
```
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load LTRpred package
library("LTRpred")
genomename<-paste(args[1],"fsa",sep=".")

# de novo LTR transposon prediction
LTRpred(genome.file=genomename,cluster=TRUE,clust.sim=0.9,copy.number.est=TRUE)

ltrpredresults <- paste(args[1],"ltrpred", sep="_")
ltrpreddatasheet <-paste(args[1],"LTRpred_DataSheet.tsv",sep="_")
strainfile<-file.path(ltrpredresults,ltrpreddatasheet)

# import LTRpred prediction output
strain_LTRpred<-read.ltrpred(strainfile)

# apply filter to predicted LTRs
functional<-quality.filter(strain_LTRpred,sim=97,n.orfs=1,strategy="stringent")
ltrdigest<- paste(args[1],"ltrdigest",sep="_")
complete<- paste(args[1],"ltrdigest_complete.fas",sep="-")
outfilename<-file.path(paste(args[1],"ltrpred.97.fa",sep="."))
predictionfile<- file.path(ltrpredresults,ltrdigest,complete)

# export filtered predictions in fasta format
pred2fasta(LTRpred.tbl=functional,prediction.file=predictionfile,output=outfilename)
```
Example run on genome fasta file Nara.fsa
```
Rscript ./run_LTRpred.R Nara
```
Output is Nara.ltrpred.97.fa

### 1.4 Run TransposonPSI
```
perl transposonPSI.pl Nara.fsa nuc
```
Pull out fasta sequences related to besthits. Run custom script parse_gff_to_fasta2.py

Example:
```
$CUSTOM_SCRIPTS/parse_gff_to_fasta2.py -gff Nara.fsa.TPSI.allHits.chains.bestPerLocus.gff3 -num 1 -feature translated_nucleotide_match -s = -s \; -outfile Nara.TPSI.fsa -fasta Nara.fsa

```
parse_gff_to_fasta2.py requires bedtools (https://bedtools.readthedocs.io/en/latest/). Run: parse_gff_to_fasta2.py -h for more details on use.

## 2. Merge the predictions

Combine TransposonPSI, LTRpred and LTR_retriever outputs for each strain.
#### Select sequences that are at least 400 bp in length
```
cat $STRAIN-families.fa $STRAIN.fsa.LTRlib.fa  $STRAIN.ltrpred.97.fa  $STRAIN.TPSI.fsa > $STRAIN.cat4.fa
#Repeats were renamed
samtools faidx $STRAIN.cat4.fa
awk '$2>=400{print $1}' $STRAIN.cat4.fa.fai > $STRAIN.c4
python $CUSTOM_SCRIPTS/pull_fasta.py -f $STRAIN.cat4.fa -l $STRAIN.c4 > $STRAIN.cat4.400.fasta
awk '{print $1"\t"FILENAME"_"NR}' $STRAIN.c4 > $STRAIN.cat4.400.rn.list
python $CUSTOM_SCRIPTS/rename.py $STRAIN.cat4.400.rn.list $STRAIN.cat4.400.fasta > $STRAIN.cat4.400.rn.fasta
```
Combine all outputs
```
cat *cat4.400.rn.fasta > 400.rn.fasta
```
#### BLAST each genome to identify sequences with more than 5 hits

Create BLAST databases for sequences to be used
```
makeblastdb -in $STRAIN.genome.fsa -out ~/db/genome/$STRAIN -dbtype nucl

```
Run blast on each genome
```
blastn -db ~/db/genome/$STRAIN -num_threads 16 -evalue 1E-15 -outfmt 6 -out blastn.$STRAIN.out -query 400.rn.fasta
blastn -db ~/db/genome/CF413 -num_threads 16 -evalue 1E-15 -outfmt 6 -out blastn.CF413.out -query 400.rn.fasta
blastn -db ~/db/genome/Cg363 -num_threads 16 -evalue 1E-15 -outfmt 6 -out blastn.Cg363.out -query 400.rn.fasta
blastn -db ~/db/genome/Cg56 -num_threads 16 -evalue 1E-15 -outfmt 6 -out blastn.Cg56.out -query 400.rn.fasta
blastn -db ~/db/genome/Chig_IMI_JGI -num_threads 16 -evalue 1E-15 -outfmt 6 -out blastn.Chig_IMI_JGI.out -query 400.rn.fasta
```
Select sequences with more than five blast hits in a single genome.

#### Identify sequences with homology to RepBase sequences
Perform blastx on RepBase peptide sequences.

```
makeblastdb -dbtype prot -in ~/bin/RepeatMasker/RepeatMasker/Libraries/RepeatPeps.lib -out RepeatPeps.db
blastx -db /home/share/DB/RepeatPeps.db -num_threads 16 -evalue 1E-5 -outfmt 6 -out blastx.repbase.out -query 400.rn.fasta
```
#### Filter outputs
Using the BLAST outputs from the above steps, repeat sequences were filtered to include only sequences with:
- more than five BLASTn hits in any single genome (E-value cutoff 1E-15); and/or
- significant homology to any sequence in the RepBase peptide database (BLASTx E-value cutoff 1E-5).

Create a list of sequences that meet the above criteria (e.g. compiled.out.list).
```
$CUSTOM_SCRIPTS/pull_fasta.py -f 400.rn.fasta -l compiled.out.list > compiled.out.fasta
```
### Merge sequences from all strains using vsearch
```
vsearch -sortbylength compiled.out.fasta --output CGSChig.merged.sorted.fasta --log vsearch.log
vsearch -cluster_fast CGSChig.merged.sorted.fasta --id 0.8 --centroids CGSChig.centroids.fa --uc CGSChig.result.uc -consout CGSChig.nr.cons.fa -msaout CGSChig.aligned.fa --log vsearch2.log
```
The consensus sequence was used for the next steps.

### Annotate custom library
RepeatClassifer from the RepeatModeler package was used to annotate the custom library
```
RepeatClassifier -consensi CGSChig.nr.cons.fa -engine ncbi
```
In addition, transposonPSI was used to annotate the library.
```
perl transposonPSI.pl CGSChig.nr.cons.fa nuc
```
TransposonPSI annotations were used if repeats were not classified by RepeatClassifier. Sequences that were not annotated by both programs were annotated as "Unknown". Fastaheaders were adjusted to fit the format recognized by RepeatModeler (e.g. rc_nn#LTR/Copia)

## 3. Filter out protein coding sequences
Some of the sequences could include sequences from multigene families. To remove these, filter repeats against coding sequences that are not related to transposable elements.

To identify coding sequences that are not related to transposable elements, the predicted proteomes of genomes of interest were used as queries against the RepeatPep database.

```
blastp -query $STRAIN.protein.fa -db /home/share/DB/RepeatPeps.db -outfmt 6 -max_target_seqs 25 -culling_limit 2 -num_threads 8 -evalue 1e-5 -out $STRAIN.vs.RepeatPeps.25cul2.1e5.blastp.out
```
Combine blastp outputs.
```
cat *.vs.RepeatPeps.25cul2.1e5.blastp.out > CGSChig.cat5.vs.RepeatPeps.25cul2.1e5.blastp.out
```
From the BLASTp outputs, a list of proteins with homology to TEs can be obtained.  
Then, the custom repeat library was blasted against CDS sequences excluding those with homology to known TEs. Here, the fastaheaders in protein and CDS sequences are the same and all CDS sequences are in cat5.cds.fasta
```
awk '{print $1}' CGSChig.cat5.vs.RepeatPeps.25cul2.1e5.blastp.out |sort | uniq > CGSChig.cat5.exclude.list
python $CUSTOM_SCRIPTS/exclude_fasta.py cat5.cds.fasta CGSChig.cat5.exclude.list > CGSChig.cat5.cds.noTEs.fasta
#BLASTn against custom library
makeblastdb -in CGSChig.cat5.cds.noTEs.fasta -out ./CGSChig.cat5.cds.noTEs.fasta -dbtype nucl
blastn -task megablast -query CGSChig.nr.cons.fa -db CGSChig.cat5.cds.noTEs.fasta -outfmt 6 -max_target_seqs 25 -culling_limit 2 -num_threads 32 -evalue 1e-10 -out CGSChig.megablast.out
```
Sequences with a BLASTn hit against the CDS sequences were excluded from the final repeat library.
```
awk '{print $1}' CGSChig.megablast.out |sort|uniq > CGSChig.exclude
python $CUSTOM_SCRIPTS/exclude_fasta.py CGSChig.nr.cons.fa CGSChig.exclude > CGSChig.nr.cons.cds.filtered.fa
```


### Add taxon-specific repeats from RepBase
To add known Colletotrichum repeats from RepBase (excluding those that were annotated as artefacts, simple repeats and low complexity sequences).
```
~/RepeatMasker/RepeatMasker/util/queryRepeatDatabase.pl -species "Colletotrichum" > repeatmasker.Colletotrichum.taxon.fa
#Remove artefacts, simple repeats and low complexity sequences
awk '/>/' repeatmasker.Colletotrichum.taxon.fa |sed 's/>//g' | sed '/ARTEFACT/d' | sed '/Simple_repeat/d' |sed '/Low_complexity/d'|awk '{print $1}' > repeatmasker.Colletotrichum.taxon.fa.list2
tail -n+2 repeatmasker.Colletotrichum.taxon.fa > repeatmasker.Colletotrichum.taxon.parsed.fa
$CUSTOM_SCRIPTS/pull_fasta.py -l repeatmasker.Colletotrichum.taxon.fa.list2 -f repeatmasker.Colletotrichum.taxon.parsed.fa > repeatmasker.Colletotrichum.taxon.parsed2.fa
cat CGSChig.nr.cons.cds.filtered.fa repeatmasker.Colletotrichum.taxon.parsed2.fa > CGSChig.nr.cons.cds.filtered.rb.fa
```

## Annotation of genomes with custom repeat library
RepeatMasker was used to annotate each genome
```
RepeatMasker -lib CGSChig.nr.cons.cds.filtered.rb.fa Nara.fsa -pa 16 -dir ./ -gff
```
Final repeat library: [CGSChig.nr.cons.cds.filtered.rb.fa](../additional_data/CGSChig.nr.cons.cds.filtered.rb.fa)

## Annotation of telomeric repeats
Identify potential telomeric regions based on density of telomeric repeats across windows of specified intervals. Only intervals which are greater than the percentile specified will be reported (e.g. for -pc 95, only the top 5% intervals ranked by density of telomeric repeats will be given). Output will be a bed format file.
```
python identify_potential_telomeres.py -f fastafile -w windowsize -g bedtoolsgenomefile -o outputfile -pc percentile
```
## Acknowledgments

The following resources were referred to while constructing this pipeline:

*  https://blaxter-lab-documentation.readthedocs.io/en/latest/repeat-masking.html
*  Coghlan et al. 2018 Protocol Exchange https://protocolexchange.researchsquare.com/article/nprot-6761/v1
*  Castanera et al. 2016 PLOS Genetics https://doi.org/10.1371/journal.pgen.1006108
*  The very helpful Avrilomics blog by Avril Coghlan (http://avrilomics.blogspot.com)

## Contributors

* **Pamela Gan**
* **Ayako Tsushima** Pipeline testing
