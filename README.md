# Scripts used for analyzing PAIso-Seq data
Welcome to PAIso-Seq.

## Software Dependencies
This analysis pipeline is implemented in Perl and R. Packages and software Versions used are listed below:
```
Perl v5.26.2
R 3.5.1
minimap2 2.15-r905
```

## Analysis of PAIso-Seq data
### Step 0. Known the library structure of PAIso-Seq data
5'-Adapter--->cDNA--->Barcode->Adapter-3'
5'-Adapter Sequence: 5'-AAGCAGTGGTATCAACGCAGAGTACATGGG-3' (30 nt)
3'-Adapter Sequence: 5'-GTACTCTGCGTTGATACCACTGCTT-3' (25 nt)
Barcode Sequence: GAGTGCTACTCTAGTA (16 nt)
### Step 1. Circular Consensus Sequence calling
`ccs movieX.subreads.bam movieX.ccs.bam --noPolish --minPasses 1 &>ccs.log`
However,We suggest the following command
`ccs movieX.subreads.bam movieX.ccs.bam --richQVs &>ccs.log`

### Step 2. Demultiplex and trim adapters

### Step 3. Calculate poly(A) tail length and call non-A residues
