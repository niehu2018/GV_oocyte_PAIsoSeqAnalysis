#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
try:
	import regex
	import sys
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna
except ImportError:
	print("Python modules import error.")
	quit()

# functions
def usage():
	print("""
Description:
	This script is used to trim 5' barcoded primer of Pacbio Iso-Seq CCS reads
Usage:
	python trim_barcoded_adapter.py input barcoded_primer cut mismatch 1>out.fasta 2>err.fasta
Examples:
	python trim_adapter.py ccs.fasta CAGCAGTATAGACTGTGTACTCTGCGTTGATACCACTGCTT 22 2 1>ccs.out.fasta 2>ccs.err.fasta
	python trim_adapter.py ccs.fasta CAGCAGTATAGACTGTGTACTCTGCGTTGATACCACTGCTT 22 2 1>ccs.out.fasta 2>ccs.err.fasta
Authors:
	Hu Nie, niehu@genetics.ac.cn
Versions:
	Version 0.3.0, 2019-04-24
Notes:
	Step 1: Search barcoded adapter, orientate read and trim barcoded adapter and afterward nucleotides
	Step 2: Search 3' adapter, trim 3' adapter and the nucleotides before 3' adapter
	Step 3: Search 5' adapter, trim 5' adapter and the nucleotides before 5' adapter
""")

def tirm_chimera(cDNA, adapter, mismatch):
	"""
		This function is used to do step 2
	"""
	pattern = ""
	if mismatch != 0:
		pattern = r"(" + adapter + "){e<=" + str(mismatch) + "}"
	else:
		pattern = adapter
	
	matches = regex.finditer( pattern, cDNA, overlapped=False)
	clean_cDNA = ""
	for it in matches:
		match_end = it.end()
		clean_cDNA = cDNA[match_end:]
	if clean_cDNA != "":
		return clean_cDNA
	else:
		return cDNA

def trim_5p_adapter(cDNA, p5_adapter, mismatch):
	"""
		This function is used to do step 3
	"""
	pattern = ""
	if mismatch != 0:
		pattern = r"(" + p5_adapter + "){e<=" + str(mismatch) + "}"
	else:
		pattern = p5_adapter
	
	matches = regex.finditer( pattern, cDNA, overlapped=False)
	clean_cDNA = ""
	for it in matches:
		match_end = it.end()
		clean_cDNA = cDNA[match_end:]
	
	if clean_cDNA != "":
		return clean_cDNA
	else:
		return cDNA

# main program
if __name__ == "__main__":
	
	# Check arguments
	if len(sys.argv) != 5:
		usage()
		quit()
	
	# Read command line arguments
	target_file          =  sys.argv[1] 			   # CCS reads
	barcode_seq          =  sys.argv[2] 			   # Barcode sequence
	barcode_cut          =  int(sys.argv[3])
	mismatch             =  int(sys.argv[4])
	

	# Step 1
	barcode_cut_seq = barcode_seq[:barcode_cut]
	target_id = ""
	target_seq = ""
	# Scan CCS reads
	for record in SeqIO.parse(target_file, "fasta"):
		
		# read CCS
		target_id = record.id
		target_seq1 = str(record.seq)
		
		# CCS reverse complement
		target_seq2 = Seq(target_seq1)
		target_seq2 = target_seq2.reverse_complement()
		
		# Search pattern
		pattern = ""
		if mismatch != 0:
			pattern = r"("+ barcode_cut_seq + "){e<=" + str(mismatch) + "}"
		else:
			pattern = barcode_cut_seq
		
		# match CCS
		matches = regex.finditer( pattern, target_seq1, overlapped=False)
		cDNA = ""
		for it in matches:
			match_start = it.start()
			match_end = it.end()
			cDNA = target_seq1[:match_start]

		# match CCS RC
		if cDNA == "":
			matches = regex.finditer( pattern, str(target_seq2), overlapped=False)
			for it in matches:
				match_start = it.start()
				match_end = it.end()
				cDNA = target_seq2[:match_start]

			if cDNA == "":
				print(">" + target_id , file=sys.stderr)
				print(target_seq1, file=sys.stderr)
				continue

		# tirm chimera
		adapter = barcode_seq[16:]
		clean_cDNA = tirm_chimera(str(cDNA), str(adapter), 3)

		# trim 5p adapter
		p5_adapter = Seq(barcode_seq[16:])
		p5_adapter = p5_adapter.reverse_complement()
		p5_adapter = p5_adapter + "ATGGG"
		clean_cDNA = trim_5p_adapter(str(clean_cDNA), str(p5_adapter), 3)
		
		print(">"+ target_id)
		print(clean_cDNA)
