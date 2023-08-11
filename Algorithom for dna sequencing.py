#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import re
from io import BytesIO
from xml.etree import ElementTree as ET
from collections import defaultdict
from io import BytesIO, StringIO
from pathlib import Path
from zipfile import ZipFile
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import (
    NNITreeSearcher,
    ParsimonyScorer,
    ParsimonyTreeConstructor,
)
from Bio.SeqRecord import SeqRecord
from ncbi.datasets import GeneApi
import sys, os, time
from Bio import Entrez, SeqIO

import random
random.choice('ACGTTA')
random.seed(7)
seq=''
for i in range(10): 
    seq+=random.choice('ACGTTA')
print(seq)
seq=''.join([random.choice('ACGTTA') for i in range(10)])
print(seq)
seq[1:3]
seq[ : 3]
seq[0:3]
seq[7:len(seq)]
#mManipulating two strings
def longestcommonprefix(s1,s2):
    i=0
    while i < len(s1) and i < len(s2) and s1[i]==s2[i]:        
        i+=1
    return s1[:i]
longestcommonprefix('ACCATGT','ACCAGAC')
#Mastch betwwen two strings
def match (s1,s2):
    if not len(s1)==len(s2):
        return False
    for i in range(len(s1)):
        if not s1[i]==s2[i]:
            return False 
    return True
    for i in range(len(s1)):
        if not s1[i]!=s2[i]:
            return False 
    return True
match('ATGCT','ATGCGAT')
print(match)
match('ATGCT','ATGCGAT')
match('ATGCT','ATGCT')
complemet={'A':'T','C':'G'}
complemet['C']

def  reversecomplement (s):
    complement={'A':'T','C':'G','G':'C','T':'A'}
    t=''
    for base in s:
        t=complement[base]+t
    return t
reversecomplement('AGTCCTACGATCGC')
#!wget -q -O- "https://www.ncbi.nlm.nih.gov/nuccore/KP857730.1?report=fasta"
#!wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=CAA37914&rettype=fasta" -O CAA37914.fa
get_ipython().system('wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_049243.1&rettype=fasta" -O NG_049243.1.fa ')
def readGenome(filename):
    genome = ''
    with open (filename,'r') as f:
        for line in f:
            if not line[0] =='>':
                genome+= line.rstrip()
                
    return genome
genome=readGenome ('NG_049243.1.fa')
genome[:100]
print(genome)
x=open('sequence.fasta','r')#downloaded from ncbi
x.read()
def readGenome(filename):
    genome = ''
    with open (filename,'r') as f:
        for line in f:
            if not line[0] =='>':
                genome+= line.rstrip()
                
    return genome
genome=readGenome ('sequence.fasta')
genome[:100]
print(genome)
len(genome)
#Count number of bases in genome fasta file.
counts={'A':0,'G':0,'C':0,'T':0}
for bases in genome:
        counts[bases] += 1
print(counts)
import collections
collections.Counter(genome)
y=open('ERR10747512_1.fastq','r')
# examine first few lines
y.read(10)
get_ipython().system('head ERR10747512_1.fastq')#downloaded from ENA
# examine first few lines 
get_ipython().system('tail ERR10747512_1.fastq')
# convert quality score character to numeric 
ord("?")
ord("F")
ord(":")
ord(",")
seq = y.read()
seq[:500]
def dnaGenome(filename):
    genome = ''
    with open (filename,'r') as f:
        for line in f:
            if not line[0] =='>':
                fastq+= line.rstrip()
                
    return genome
fastq=readGenome ('ERR10747512_1.fastq')
fastq[:100]
len(fastq)
seq[100]
fastq[:500]
fastq= fastq.replace("\n", "")
fastq = fastq.replace("\t", "")
fastq[:500]
#remove name line 
with open("ERR10747512_1.fastq", "r") as f: 
    fastq = f.readline()
    fastq = f.read()
fastq[:500]
fastq= fastq.replace("\n", "")
fastq[:500]
fastq[:1000]
# # FASTQ Format Handling Template
def readFastq(filename):
    """Reads FASTQ file and remove the special characters!"""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
seqs, quals = readFastq('ERR10747512_1.fastq')
seqs[:100]
quals[:100]
def phred33ToQ(qual):
    return ord(qual) - 33
phred33ToQ('F')
def createHist(qualities):
    # Create a histogram of quality scores
    hist = [0]*50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist
h = createHist(quals)
print(h)
# Plot the histogram
#%matplotlib inline turns on “inline plotting”, where plot graphics will appear in your notebook. This has important implications for interactivity: for inline plotting, commands in cells below the cell that outputs a plot will not affect the plot.
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.bar (range(len(h)), h)#plt.plot (range(len(h)), h)
plt.show()


# # Find GC by position in the read:

import matplotlib.pyplot as plt

def findGCByPos(reads):
    #Find the GC ratio at each position in the read
    # Initialize lists to store the G/C counts and total counts at each position
    max_read_length = max(len(read) for read in reads)
    gc = [0] * max_read_length
    totals = [0] * max_read_length
    
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    
    # Calculate GC ratios by position
    gc_ratios = [gc[i] / float(totals[i]) if totals[i] > 0 else 0 for i in range(max_read_length)]
    
    return gc_ratios

# Call the function to calculate GC ratios by position
gc = findGCByPos(seqs)

# Plot the GC ratios
plt.plot(range(len(gc)), gc)
plt.xlabel("Position")
plt.ylabel("GC Ratio")
plt.title("GC Ratio at Each Position")
plt.show()






