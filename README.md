# Sequence Analysis Script
 
## Overview
This Python script performs a detailed analysis of genomic sequence data, specifically focusing on gene information and chromosomal sequences from the mouse genome (mm9 build). It demonstrates fundamental bioinformatics tasks such as loading genomic data, extracting gene-specific sequences, calculating gene lengths, identifying coding vs. non-coding regions, and searching for regulatory motifs like the TATA box.

This script is ideal to perform basic genomic data handling and analysis using Python.

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
* **numpy**: For numerical computations, especially for calculating gene lengths and handling genomic statistics. If by any chance you run into trouble trying to import these, try first running "pip install numpy" in the terminal.
* **matplotlib.pyplot**: For generating plots, such as the histogram of gene lengths. If by any chance you run into trouble trying to import these, try first running "pip install matplotlib in the terminal.
* **LoadFASTA_Function**: This is a user-defined module (assumed to be in the same directory) that contains functions specific to this analysis:
    * `LoadFastaFile(fasta_file)`: Loads sequence data from a FASTA-formatted ZIP archive into a dictionary.
    * `LoadGene(gene_file)`:  Loads gene information from a text file into a dictionary.
    * `TSSChroms(gene_inf,chromosome)`: Returns a dictionary containing the Transcription Start Sites (TSS) for a given chromosome.

## Input Files

The script requires the following input files:

1.  **Gene Information File (`mm9_sel_chroms_knownGene.txt`)**: A text file containing gene annotations.  This file provides information such as gene identifiers, chromosome locations, and start/end positions.  The file format is specific to the UCSC Genome Browser's knownGene table.
2.  **FASTA Sequence File (`selChroms_mm9.fa.zip`)**: A ZIP archive containing the genomic sequence data in FASTA format.  The FASTA format is a standard text-based format for representing nucleotide or amino acid sequences.  The sequences are expected to be organized by chromosome. You have to download it from https://drive.google.com/file/d/1jxOqX3W5Dsdhlr0Q3i2_oeujS7AJ3D2K/view?usp=drive_link

## How to Use

The Main_file tells what each step in the code does. As per the breaks used in the file, I've written this in Spyder. 


### Script Structure and Functionality

The script is organized into logical sections, each performing a specific part of the genomic analysis.

**Example Usage:**

To access the chromosome of a specific gene (e.g., `'uc009auw.1'`):

[insert screen capture]

#### 3\. Evaluating and Counting Genes

This section processes the `gene_inf` dictionary to identify all unique chromosomes present in the dataset and then counts the number of genes on each identified chromosome.

**Example Usage:**

After running this section, the `gene_counts` dictionary will hold the gene count for each chromosome:

[insert screen capture]

#### 4\. Loading Genomic Sequences

Reads raw genomic sequences from the zipped FASTA file (`selChroms_mm9.fa.zip`) into a dictionary called `seq_dict`. This process can be memory-intensive and may take some time depending on system resources.

**Example Usage:**

To check the length of a specific chromosome sequence (e.g., chromosome 6):

[insert screen capture]

To access the chromosomal location of a specific gene (e.g., `cntn4` which is 'uc009dcr.2'):

[insert screen capture]

#### 5\. Extracting Gene Sequences

Demonstrates how to extract the full sequence of a specific gene (e.g., `cntn4`) using its start, end, and chromosomal location from the `gene_inf` and `seq_dict` dictionaries. It also shows how to find the starting position of an 'ATG' codon within a gene sequence.

**Example Usage:**

After executing the `cntn4_seq` extraction:

[insert screen capture]

#### 6\. Genomic Statistics

Calculates the length of every gene in the `gene_inf` dictionary and visualizes the distribution of gene lengths using a histogram generated by `matplotlib`.

**Example Usage:**

To access the calculated length of a specific gene:

[insert screen capture]

To compare gene lengths:

[insert screen capture]

#### 7\. Beyond the Coding Region (Non-Coding DNA Analysis)

This section illustrates a method to identify and quantify coding vs. non-coding regions on a specific chromosome (e.g., chromosome 6). It uses a boolean NumPy array to mark all positions that fall within a gene's coding region.

**Example Usage:**

After this section runs, you can inspect the calculated values:

[insert screen capture]

#### 8\. Transcription Start Site (TSS) Analysis and TATA Motif Search

This advanced section uses the `TSSChroms` function to find all transcription start sites (TSS) for genes on a specified chromosome. It then searches for the 'TATA' binding motif within a 40 base pair window upstream (for positive strand genes) or downstream (for negative strand genes) of the TSS. This helps identify genes potentially regulated by TATA-binding proteins.

**Example Usage:**

To see how many genes contain a TATA motif within the specified window:

[insert screen capture]

To calculate the mean distance of the TATA motif from the TSS for the identified genes:

[insert screen capture]

