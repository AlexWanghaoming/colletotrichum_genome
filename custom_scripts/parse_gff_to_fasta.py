#!/usr/bin/env python3
import re
import argparse

parser = argparse.ArgumentParser(description='Pull features of interest from a gff. Follow up with fasta_formatter -in inputfasta -out outputfasta -w 60 to get properly wrapped fastafile')
parser.add_argument('-gff', metavar='gff', type=str,
                    help='gff file')
parser.add_argument('-fasta', metavar='fasta', type=str,
                    help='fasta file')
parser.add_argument('-feature', metavar='feature', type=str,
                    help='feature to be pulled out (e.g. gene, exon, mRNA)')
parser.add_argument('-ID', metavar='feature ID', help='Number referring to string to use as fastaheader for each feature from the last column in a gff file. Follows python list conventions. Last column split by :;=')


def get_feature_intervals (gfffile,featureofinterest,plusorminus,fieldnumber):
    featurelist = []
    d = {}
    for lines in open(gfffile):
        try:
            line = lines.strip()
            line1 = line.split("\t")
            chrom = line1[0]
            start = line1[3]
            feature = line1[2]
            end = line1[4]
            tn = re.split(r"[=;]",line1[8])
            #tnscore = tn[5].split()
            name = tn[int(fieldnumber)]
            #tntype = tn[3]
            strand = line1[6]
            startstop = (chrom,start,end)
            if (feature!=featureofinterest):
                pass
            else:
                if (strand!=plusorminus):
                    pass
                else:
                    if name not in featurelist:
                        interval = [startstop]
                        #interval.append(startstop)
                        featurelist.append(name)
                        d[featurelist[-1]]=interval
                    else:
                        interval.append(startstop)
                        d[featurelist[-1]]=interval
        except:
            pass
    return d,featurelist

def pull_fasta_seqs(chrom,start,end,fastafile):
    from Bio import SeqIO, Seq
    fastafh = SeqIO.parse(open(fastafile),"fasta")
    for my_seq in fastafh:
        if my_seq.id in chrom:
            seq2 = my_seq.seq[int(start)-1:int(end)]
    return seq2

def pull_fasta_seqs_rc(chrom,start,end,fastafile):
    from Bio import SeqIO, Seq
    from Bio.Alphabet import generic_dna
    fastafh = SeqIO.parse(open(fastafile),"fasta")
    for my_seq in fastafh:
        if my_seq.id in chrom:
            seq2 = my_seq.seq[int(start)-1:int(end)]
            seq3=seq2.reverse_complement()
    return seq3

def get_fasta_from_plus_intervals(intervaldic,fastafile):
    for feature, value in intervaldic.items():
         print(">"+feature)
         for interval in value:
             chrom = interval[0]
             start = interval[1]
             end = interval[2]
             seq = pull_fasta_seqs(chrom,start,end,fastafile)
             print(seq)

def get_fasta_from_minus_intervals(intervaldic,fastafile):
    for feature, value in intervaldic.items():
        print(">"+feature)
        #
        value.sort()
        value.reverse()
        print(value)
        for interval in value:
             chrom = interval[0]
             start = interval[1]
             end = interval[2]
             seq = pull_fasta_seqs_rc(chrom,start,end,fastafile)
             print(seq)

args = parser.parse_args()

gfffile=args.gff
fastafile=args.fasta
featureofinterest=args.feature
fieldnumber=args.ID

plusintervals, plusfeaturelist = get_feature_intervals(gfffile,featureofinterest,"+",fieldnumber)
minusintervals, minusfeaturelist = get_feature_intervals(gfffile,featureofinterest,"-",fieldnumber)
get_fasta_from_plus_intervals(plusintervals,fastafile)
get_fasta_from_minus_intervals(minusintervals,fastafile)
