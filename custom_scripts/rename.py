# -*- coding: utf-8 -*-
#Input1: tab delimited file with original name in column 1 and new name in column 2
#Input2: fasta file to be modified
#Run: python3 rename.py input1 input2 > newfastafile
#Set up a dictionary holding the new and old fasta names
import sys
namelist = sys.argv[1]
modfasta = sys.argv[2]
d = {}
with open(namelist) as f:
    for lines in f:
        line = lines.strip()
        (key, val) = line.split("\t")
        d[key] = val
#Open new file
fh = open(modfasta)
f = open("missing.txt","a+")
for lines in fh:
    line2 = lines.strip()
#Look for fasta headers
    try:
        if lines.startswith(">"):
            line2 = lines.split()
            soi = line2[0]
#Look up new names in dictionary d
            newname = d[soi[1:]]
            print (">"+newname)
        else:
            print (line2)
    except:
        f.write(soi+"\n")
