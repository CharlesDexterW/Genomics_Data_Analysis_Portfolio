# Sequence Analysis Script
 
So far, the following codeline is able to work with data file formats of known gene chromosomes from UCSC Genome Browser, commonly used to describe genomic features like location and chromosome structures. I have adapted it from materials provided in the MITbiox online program Quantitative Biology Workshop from MIT. I am grateful for the educational resources.


## Description

This Python script performs a series of computational analyses on genomic sequence data. It processes gene information and raw sequence data to:

* Count the number of genes per chromosome.
* Extract gene sequences from chromosome sequences.
* Calculate gene lengths.
* Analyze the distribution of gene lengths.
* Determine the amount of coding and non-coding DNA.
* Identify potential regulatory regions by searching for Transcription Start Site (TSS) motifs (specifically, the TATA box).

The script is designed to handle data from a specific organism (mouse, *Mus musculus*), using gene annotations and genomic sequence files.  However, with modification to the input files, it could be adapted to other organisms.

## Dependencies

The script relies on the following Python libraries:

* **csv**: For working with comma-separated values files.
* **gc**: For garbage collection and memory management.
* **io**: For handling various types of input/output operations.
* **zipfile**: For working with ZIP archives.
* **numpy**: For numerical computations, especially for calculating gene lengths and handling genomic statistics.
* **matplotlib.pyplot**: For generating plots, such as the histogram of gene lengths.
* **LoadFASTA_Function**: This is a user-defined module (assumed to be in the same directory) that contains functions specific to this analysis:
    * `LoadFastaFile(fasta_file)`: Loads sequence data from a FASTA-formatted ZIP archive into a dictionary.
    * `LoadGene(gene_file)`:  Loads gene information from a text file into a dictionary.
    * `TSSChroms(gene_inf,chromosome)`: Returns a dictionary containing the Transcription Start Sites (TSS) for a given chromosome.

## Input Files

The script requires the following input files:

1.  **Gene Information File (`mm9_sel_chroms_knownGene.txt`)**: A text file containing gene annotations.  This file provides information such as gene identifiers, chromosome locations, and start/end positions.  The file format is specific to the UCSC Genome Browser's knownGene table.
2.  **FASTA Sequence File (`selChroms_mm9.fa.zip`)**: A ZIP archive containing the genomic sequence data in FASTA format.  The FASTA format is a standard text-based format for representing nucleotide or amino acid sequences.  The sequences are expected to be organized by chromosome.

## How to Use

The Genome Dictionary file tells what each step in the code does. As per the breaks used in the file, I've written this in Spyder. 
