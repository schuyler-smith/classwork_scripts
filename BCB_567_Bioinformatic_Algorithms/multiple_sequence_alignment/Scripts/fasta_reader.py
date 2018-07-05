def read_fasta(file):

	seqs_headers = []
	seqs = []
	for line in open(file):
		if ">" in line:
			seqs_headers.append(str(line.strip().strip('>')))
		else:
			seqs.append(line.strip())
	return seqs_headers, seqs

