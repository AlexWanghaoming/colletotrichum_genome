#!/usr/bin/env python3

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="""To choose specific field from a fastaheader to be the new fasta header.
Default separator is "|"
Run as follows:
parse_fastaheaders.py -f CCHL.nt.fa -num 4 -s : > CCHL.rn.nt.fa""",formatter_class=RawTextHelpFormatter)
parser.add_argument('-f', metavar='fasta', type=str,
                    help='fasta file')
parser.add_argument('-num', metavar='number', type=int,
                    help='field number')
#parser.add_argument('-s',metavar="separators",nargs="+",type=str,help="separators")
parser.add_argument('-s',metavar="separators",action="append",type=str,help="separators")

args = parser.parse_args()
fastafile = args.f
fieldnum = args.num
seps = args.s


if seps is None:
    seps = []

seps.append(">")
seps.append("|")

#from stackoverflow:https://stackoverflow.com/questions/4697006/python-split-string-by-list-of-separators
def split(txt, seps):
    default_sep = seps[0]

    # we skip seps[0] because that's the default seperator
    for sep in seps[1:]:
        txt = txt.replace(sep, default_sep)
    return [i.strip() for i in txt.split(default_sep)]

def rename_fastaheaders(fastafile,seps,fieldnum):
    with open(fastafile) as fh:
        for lines in fh:
            lines = lines.strip()
            if lines.startswith(">"):
                line = split(lines,seps)
                field = line[fieldnum]
                print (">"+field)
            else:
                print (lines)
rename_fastaheaders(fastafile,seps,fieldnum)
