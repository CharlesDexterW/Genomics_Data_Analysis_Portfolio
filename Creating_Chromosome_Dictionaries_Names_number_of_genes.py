# -*- coding: utf-8 -*-
"""
Created on Thu Jan 4 19:23:40 2025

@author: Benjamin Garcés Cifuentes
email='agarces2381@gmail.com'

"""
#%%
'''
# Sequence Analysis
Importing necessary packages & Functions

csv will provide functionality for reading and writing tabular data in CSV (Comma Separated Values) format.
gc will provide an interface to the garbage collector, allowing manual control over memory management.
io will provide Python's main facilities for handling various types of I/O (input/output), including file I/O, string I/O, and in-memory streams.
zipfile will allow you to create, read, write, append, and list files in a ZIP archive.
'''

import csv
import gc
import io
import zipfile
import pyflakes

from LoadFASTA_Function import LoadFastaFile, LoadGene, TSSChroms

#%%
'''
First I will assign a file name to the variable to work smoothly.
'''

gene_file= 'mm9_sel_chroms_knownGene.txt'

#%% Loading genome into a variable

'''
I've been trying to use rstrip() to create a cleaner output. Still figuring it out.
This command creates a "file_open" variable containing the sequence as a list. we can check its lenght.
Meanwhile we use fileread
'''

file_open = open(gene_file).readlines() 

#%%
'''
Now to interpret the file we use "LoadGene" to create a dictionary from called  "gene_inf"
'''
gene_inf=LoadGene(gene_file)

#%%
'''
To access any needed information such as the chromosome and position of 
the gene on the chromosome, just use brackets like this:
'''
gene_inf['uc009auw.1']['chr']

#%% Evaluating the Data
'''
To count the number of chromosomes we first create a list:
'''

chroms=[]

#%% 
'''
For each gene, we look to see if the chromosome is in list.  If not, 
add it:
'''
for k in gene_inf.keys():
  chr=gene_inf[k]['chr']
  if chr not in chroms:
    chroms=chroms+[chr]

#%% Counting Genes
'''
The gene count can be accessed by navigating through the dictionary of information.
but first, this codeline creates a dictionary to store the relevant information:
'''
gene_counts={}
#%%
'''
Now this new dictionary must be filled using two for loops: one for each 
chromosome and one for each gene:
'''
for chr in chroms:
  chrom_count=0 
  for k in gene_inf.keys():
    if gene_inf[k]['chr']==chr:
      chrom_count+=1
  gene_counts[chr]=chrom_count
#%%
'''
for each chromosome, I determined the number of associated genes and reported it in the output with this code snippet.
'''
gene_counts

#%% Locating identified Genes in Genomic Data
"""
With the gene information at hand, the raw gene sequences can be read from the FASTA format. As always, it’s easier to work with these type of files by assigning them to a variable.
"""
fasta_file='selChroms_mm9.fa.zip'
#%%
'''
The preloaded module allows you to load the sequences into a dictionary. Now be patient, this might take a minute...or two depending on your computer’s RAM and processor.
'''
seq_dict=LoadFastaFile(fasta_file)
#%%
'''
In this dictionary, Chromosomes are individual sequences, which are shown if you type “seq_dict.keys()“
Now, you can start using the gene information within the chromosomes. 
For instance, gene Cntn4, which has was assigned the identifier of 
uc009dcr.2 by the genome browser:
'''
cntn4="uc009dcr.2"

#%%
'''
Now we can find out the chromosomal location of Cntn4 using our gene_info 
dictionary to get all the relevant data we collected:
To get the sequence for Cntn4, we need to first figure out what chromosome it 
is on:
'''
gene_inf[cntn4]['chr']

#%% Dictionary of Data
'''
Now the len() function can give you the length of any chromosome, like chromosome 16.
'''
len(seq_dict['chr6'])

#%% Extracting Genetic Information
'''
You can peek into the entire chromosomal sequence where we know Cntn4 to reside:
'''
hchr=gene_inf[cntn4]['chr']
#%%
hst=gene_inf[cntn4]['start']
#%%
hen=gene_inf[cntn4]['end']
#%%
'''
Now the following variable pulls out the sequence of interest.
'''
cntn4_seq=seq_dict[hchr][hst:hen]
#%%
'''
Type cntn4_seq at the prompt to see it for yourself (this may take a minute):
'''
cntn4_seq
#%%
'''
Any full chromosome sequence might roll off the screen, since most chromosome 
sequences are very long. However the python string processing function will 
calculate the length of the gene.
'''
len(cntn4_seq)
#%%
# Extracting Information about Genes
'''
You can query the genomic information for genes of interest and 
determine some of their information, like its length and location in the chromosome 
by using the command >>gene_inf[“gene_identifier”][‘chr’]

after having assigned variables like previous ones. Also, use >> len(“gene_identifier”_seq) 
to find the amount of basepairs in the desired gene. 
'''
'''
So far we can now use the raw sequence data peek into the genes themselves.
'''
#%% Reading the Genome
'''
Each sequence is a string of characters now, you can identify specific parts of that gene.  
'''
cntn4_seq[5:200]


#%%
'''
So far you still wouldn’t be able to tell the function of the gene. The following
steps are to translate the genetic sequence after identifying at which part 
exactly does the genetic sequence begin.

As you know, the translation of the protein occurs at the start codon, represented 
by the three  base pairs ATG.  We can identify where translation begins in the gene 
sequence using the index command, built into the Python string library:

'''
cntn4_seq.index('ATG')
#%% Translating the Genome
'''
Now let's do it with Matn2
'''
Matn2="uc007vll.1"
#%% 
gene_inf[Matn2]
#%%
hchr=gene_inf[Matn2]['chr']
#%%
hst=gene_inf[Matn2]['start']
#%%
hen=gene_inf[Matn2]['end']
#%%
Matn2_seq=seq_dict[hchr][hst:hen]
#%%
'''
Use the  gene identifier and the previously defined Matn2 sequence, to identify 
the first ten nucleotides in the Matn2 gene by typing:
'''
Matn2_seq[0:10]

#%%
'''
Also, to identify the index of the translation initiation codon ATG in the Matn2 gene, type:
'''
Matn2_seq.index('ATG')
#%% Genomic Statistics
"""

Here I create data structures. So far it works properly with python 3.12.7

Using the gene_inf dictionary, find out the length of every gene on the four chromosomes collected. 


"""
#%%
"""
Now, 
First, import the NumPy module to include more complicated math calculations:
"""
import numpy as np
#%%
"""
Now you can find the start and end of each gene, and store the difference:
"""

gene_lengths={}
for g in gene_inf.keys():
        st=gene_inf[g]['start']
        en=gene_inf[g]['end']
        gene_lengths[g]=np.absolute(en-st)
#%%
"""
Now Let’s use pyplot's basic plotting tools in matplotlib to create a histogram:
"""

import matplotlib.pyplot as plt
plt.hist(gene_lengths.values(), bins=50)

#%%
"""
You can change the color to green or any other by adding an additional color argument:
"""

plt.hist(gene_lengths.values(), bins=50, facecolor='green')

#%%
"""
For a better view let’s switch to a logarithmic scale:
"""

plt.hist(gene_lengths.values(), bins=50, log=True, facecolor='green')
#%%
"""
now if I want to use the max function here...
"""
max(gene_lengths,key=gene_lengths.get)
#%% basic mathematic operations with genes.
"""
To compare lengths of genes we can call out a gene by its identifier in the dictionary and add or substract it from 
another gene sequence such as cntn4.

"""
gene_lengths["uc012enb.1"]-len(cntn4_seq)

#%% Beyond the Coding Region
"""
Now let’s use the gene_inf data structure to count the base pairs in genes compared to those that aren’t in coding regions.

By using a set data structure, which stores unique items, you can remove duplicates. To make sure there are no double countings, let’s track the indices, of chromosome x so as not to add in a position in the list of coding regions if they’re already considered in another gene this way.
"""

ingene=set()
#%%
"""
By looping through every gene in the chromosome, every sequence index is assigned to the gene.
"""
for gene in gene_inf.keys():    
    if gene_inf[gene]['chr']=='chr6':
        gene_ind=range(gene_inf[gene]['start'],gene_inf[gene]['end']) #After running this line, update the set. Even if the indices were entered already as a different gene, the command will not yield duplicates.
        ingene.update(gene_ind)
len(ingene)
#%%
"""
As usual, it could take a while. Anyway, it can be done using a boolean array.
"""
chr6_len=len(seq_dict['chr6'])
#%%
ingene_numpy=np.zeros(chr6_len,dtype=bool)
#%%
for gene in gene_inf.keys():
    if gene_inf[gene]['chr']=='chr6':
        start_in=gene_inf[gene]['start']
        end_in=gene_inf[gene]['end']
        ingene_numpy[start_in:end_in]=True
#%% Counting the coding base pairs
"""
Now by summing this array, we will find out the number of coding index sites (lenght of coding sequence in chr6)
"""
sum_gene=ingene_numpy.sum()
print(sum_gene)
#%% Amount of non coding base-pairs of the chromosome 6
len(ingene_numpy)-sum_gene
#%% Fractions of non-coding DNA
sum_gene/len(ingene_numpy)
#%% Regulatory Regions
'''
Now we can see how much of the chromosome corresponds to encoding genes and what fraction of the DNA corresponds to non-consent regions
While non-coding regions do not encode RNA molecules, they may serve other purposes such as regulating DNA transcription. 

Let's check for instance the “TATA” Box protein binding single motif.
'''
# Finding Motifs
seq_dict['chr6'].upper().count('TATA')
#%%
'''
however not every TAT motif instance mean a sequence initiation. Other motifs have different functions.
Let's search for the binding sites more relevant for each gene using this function derived of a module and our dictionary.
'''
chr6_starts=TSSChroms(gene_inf,'chr6')
'''
This dictionary contains the starting positions of every gene in the chromosome.
By default, the TATA binding motif relevant to transcriptions is located in a small window of the starting sites.
Let's give it a try with a 40 base pair window upstream.
'''
#%%
'''
Now with this dictionary, we can identify which genes contain a TATA binding motif near the starting site, which means
these genes can be regulated with a TATA binding protein.
'''
tata_dis = {}
for g in chr6_starts.keys():
    e = chr6_starts[g]
    strand = gene_inf[g]['strand']
    if strand == '+':
        s = e - 40
        if 'TATA' in seq_dict['chr6'][s:e].upper():
            tata_dis[g] = seq_dict['chr6'][s:e].upper().rindex('TATA')
    else:
        s = e
        e = s + 40
        if 'TATA' in seq_dict['chr6'][s:e].upper():
            tata_dis[g] = seq_dict['chr6'][s:e].upper().index('TATA')
#%% A Simple analysis.
'''
Now how many genes have a motif within 40 bps of each transcription starting site?
'''
len(tata_dis)
#%%
'''
what's the mean distance of the transcription start sites and the TATA motif?
'''
sum(tata_dis.values())/len(tata_dis.values())
