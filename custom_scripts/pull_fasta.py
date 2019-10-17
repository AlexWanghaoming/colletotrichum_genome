#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Run as follows: pull_fasta.py -f fasta -l list > output.fasta')
parser.add_argument('-l', metavar='seqlist', type=str,
                    help='sequence list')
parser.add_argument('-f', metavar='fasta', type=str,
                    help='fasta file')

args = parser.parse_args()
seqlist = args.l
fasta = args.f
wanted = [line.strip() for line in open(seqlist)]
seqiter = SeqIO.parse(open(fasta), 'fasta')
SeqIO.write((seq for seq in seqiter if seq.id in wanted), sys.stdout, "fasta")
