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
Library structure: </br>
&nbsp;&nbsp;&nbsp;&nbsp;Forward style: 5'-TSO--->cDNA--->AAAAA...AAAAA->Barcode->Adapter-3'</br>
&nbsp;&nbsp;&nbsp;&nbsp;Reverse style: 5'-Adapter-Barcode->TTTTT...TTTTT->cDNA->TSO-3'</br>
5'-TSO Sequence: 5'-AAGCAGTGGTATCAACGCAGAGTACATGGG-3' (30 nt)</br>
3'-Adapter Sequence: 5'-GTACTCTGCGTTGATACCACTGCTT-3' (25 nt)</br>
Barcode Sequence: 5'-GAGTGCTACTCTAGTA-3' (16 nt)</br>
### Step 2. Circular Consensus Sequence calling
```
ccs movieX.subreads.bam movieX.ccs.bam --noPolish --minPasses 1 &>ccs.log
```
However,We now suggest the following command:</br>
```
ccs movieX.subreads.bam movieX.ccs.bam --richQVs &>ccs.log
```
The output file is movieX.ccs.bam</br>
Reference: https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.2.md

### Step 3. Demultiplex and trim adapters
</br>Convert movieX.ccs.bam to fasta:</br>
```
bam2fasta -u -o ccs movieX.ccs.bam
```
The output file is movieX.ccs.fasta</br>
</br>Extract data and trim adapters:</br>
```
./scripts/trim.py movieX.ccs.fasta sample GAGTGCTACTCTAGTAGTACTCTGCGTTGATACCACTGCTT 22 2 1>sample.out.fasta 2>sample.err.fasta
```
### Step 4. Calculate poly(A) tail length and call non-A residues
Create minimap2 index :
```
minimap2 -x splice -t 20 -d Mus_musculus.mmi Mus_musculus.GRCm38.dna.toplevel.chromosome.fa &>index.log
```
Collect CCS passes :
```
./scripts/GetCCSpass.pl movieX.ccs.bam >ccs.pass.txt
```
Run the pipeline:
```
perl run.pl --ccs sample.out.fasta --sample sampleName --species mm10 --minimap2_index Mus_musculus.mmi --minimap2_thread  10 --pass ccs.pass.txt &>run.log
```
