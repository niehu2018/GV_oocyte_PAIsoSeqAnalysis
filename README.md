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
### Step 1. Known the library structure of PAIso-Seq data</br>
Library structure: 5'-Adapter--->cDNA--->Barcode->Adapter-3'</br>
5'-Adapter Sequence: 5'-AAGCAGTGGTATCAACGCAGAGTACATGGG-3' (30 nt)</br>
3'-Adapter Sequence: 5'-GTACTCTGCGTTGATACCACTGCTT-3' (25 nt)</br>
Barcode Sequence: 5'-GAGTGCTACTCTAGTA-3' (16 nt)</br>
### Step 2. Circular Consensus Sequence calling
`ccs movieX.subreads.bam movieX.ccs.bam --noPolish --minPasses 1 &>ccs.log`
</br>However,We now suggest the following command:</br>
`ccs movieX.subreads.bam movieX.ccs.bam --richQVs &>ccs.log`
</br>The output file is movieX.ccs.bam
### Step 3. Demultiplex and trim adapters
</br>Convert movieX.ccs.bam to fasta</br>
`bam2fasta -u -o ccs movieX.ccs.bam`
</br>The output file is movieX.ccs.fasta
</br>
</br>Extract data and trim adapters
``</br>
### Step 4. Calculate poly(A) tail length and call non-A residues
