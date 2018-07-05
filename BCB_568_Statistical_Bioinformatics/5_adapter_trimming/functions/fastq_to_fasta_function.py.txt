
def convert_fastq(fastq_file):
 	import sys
	from Bio import SeqIO

	file_name = fastq_file.split('.')[0]
	with open(str(file_name + ".fasta"), 'w') as output_file:
		for sequence in SeqIO.parse(open(fastq_file, 'rU'), "fastq"):
			SeqIO.write(sequence, output_file, 'fasta')
			# SeqIO.write(sequence, output_file, 'fasta')
	return str(file_name + ".fasta")
