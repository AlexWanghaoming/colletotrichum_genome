#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import find_seq as fn

parser = argparse.ArgumentParser(description='Finds TTAGGG sequences in genome and outputs gff file')
parser.add_argument('-f', metavar='fasta', type=str,
                     help='fasta file')
parser.add_argument('-o', metavar='outfilename', type=str,
                     help='output filename')
args = parser.parse_args()

fasta = args.f
outfilename = args.o


fn.make_telomere_gff(fasta,outfilename)
