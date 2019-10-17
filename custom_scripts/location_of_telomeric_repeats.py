from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser(description='Searches for TTAGGG subtelomeric repeat and prints location of repeats in gff format to shell')
parser.add_argument('-fasta', metavar='fasta', type=str,
                    help='fasta file')

args = parser.parse_args()

multifastafile  = args.fasta
records = list(SeqIO.parse((multifastafile),"fasta"))
def find_sequence(motif,sequence,sequencename,strand):
    from re import finditer
    for match in finditer(motif,sequence,flags=re.IGNORECASE):
        start,end = match.span()
        start = start+1
        end = end
        print(sequencename+"\tTelomericRepeat\tTelRepeat\t"+str(start)+"\t"+str(end)+"\t0\t"+strand+"\t.\tNote=(TTAGGG)")
for sequences in records:
    find_sequence("TTAGGG",str(sequences.seq),str(sequences.id),"+")
    find_sequence("CCCTAA",str(sequences.seq),str(sequences.id),"-")
