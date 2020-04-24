#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 20:40:56 2020

@author: Pam
"""
import argparse,os
import custom_functions as fn

parser = argparse.ArgumentParser(description='Identify potential telomeres based on density of TTAGGG repeats')
parser.add_argument('-f', metavar='fasta', type=str,
                     help='fasta file')
parser.add_argument('-w', metavar='window', type=int,
                     help='window size')
parser.add_argument('-g', metavar='genomefile', type=str,
                     help='genomefile')
parser.add_argument('-o', metavar='outfilename', type=str,
                     help='output filename')
parser.add_argument('-pc', metavar='percentile', type=float,
                     help='Example: for 95, only top 5% of intervals ranked by density of TTAGGG repeats will be reported')
args = parser.parse_args()

fasta = args.f
window = args.w
genomefile = args.g
outfilename = args.o
percentile = args.pc

telomeres_gff = "tmp.telrep.gff"

fn.make_telomere_gff(fasta,telomeres_gff)
fn.find_telomeres(window,genomefile,telomeres_gff,percentile,outfilename)
os.remove("tmp.telrep.gff")
