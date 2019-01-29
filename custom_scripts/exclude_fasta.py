# -*- coding: utf-8 -*-
#Sample usage: python script.py input.fasta sequence_list.txt > output.fasta
from Bio import SeqIO
import sys

notwanted = [line.strip() for line in open(sys.argv[2])]
seqiter = SeqIO.parse(open(sys.argv[1]), 'fasta')
SeqIO.write((seq for seq in seqiter if seq.id not in notwanted), sys.stdout, "fasta")
