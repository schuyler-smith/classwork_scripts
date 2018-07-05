

def trim_fastq(fastq_file, trim_file):
	import os
	import sys
	from itertools import izip
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord

	trimmed_file = fastq_file.split('.')[0] + ".trimmed.fastq"
	with open(trimmed_file, 'w') as output_file:
		for sequence, trim in zip(SeqIO.parse(open(fastq_file, 'rU'), "fastq"), trim_file):
			new_sequence = SeqRecord(sequence.seq[trim:], id = sequence.id, name = sequence.id, description = "length=" + str(len(sequence.seq[trim:])), letter_annotations = {'phred_quality':sequence.letter_annotations.values()[0][trim:]})
			SeqIO.write(new_sequence, output_file, 'fastq')

	return str(trimmed_file + ".fastq")