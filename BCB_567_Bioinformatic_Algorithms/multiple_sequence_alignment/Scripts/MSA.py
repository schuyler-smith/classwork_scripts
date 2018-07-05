import sys
from fasta_reader import read_fasta
from word_code import wrd_code
from super_word import sw_code
from w_code import wcode
from id_seq import seq_id

if len(sys.argv) != 4:
	print "\n\nERROR: incorrect arguments given to script.\n\nUSAGE: python <MSA.py> <sequence.fasta> <word_model.txt> wlcut\n\n"
	sys.exit(1)

with open(sys.argv[2], 'r') as wm:
	word_model = wm.readline().strip()
wlcut = int(sys.argv[3])
w = len(word_model)
if w < 1:
	print "\n\nERROR: word_model must be of legnth >= 1"
	sys.exit(1)

seqs_headers, seqs_inputs = read_fasta(sys.argv[1])
C = ''
k = len(seqs_inputs)
st = []
length = []
p = -1
for seq in seqs_inputs:
	C += (seq + '#')
	st.append(p+1)
	p += len(seq)+1
	length.append(len(seq))

word_code = wrd_code(C, word_model)
positions = range(0, len(C))
sup_wrd = sw_code(positions, word_code, w, wlcut)

n_block = 1
blocks = []
block_seq_pos = []
prev = []
pos = []
for i in range(0, len(sup_wrd)+1):
	try:
		current = wcode(word_code, w, wlcut, sup_wrd[i])
	except:
		current = -1
	pos.append(sup_wrd[i-1]+1)
	if current == prev and -1 not in current:
		n_block += 1
	else:
		if n_block == k:
			check_seq_id = []
			block_start_pos = []
			for index in pos:
				number, loc = seq_id(index, st, length)
				check_seq_id.append(number)
				block_start_pos.append(loc)
			if check_seq_id == range(1,k+1):
				blocks.append(prev)
				block_seq_pos.append(block_start_pos)
		pos = []
		n_block = 1
	prev = current

for i in range(k):
	print '>' + seqs_headers[i]
	print seqs_inputs[i]
print "\nWORD_MODEL:  %s" %(word_model)
print "WLCUT:  %s" %(wlcut)
chain = []
for i in range(len(blocks)):
	for j in range(w*wlcut):
		chain.append(block_seq_pos[i][0]+j)
print "Longest superword block chain: %s" %(len(set(chain)))
print "Number of blocks in chain: %s" %(len(blocks))
for i in range(len(blocks)):
	print "\nBlock %s:" %(i+1)
	for j in range(k):
		nucleotides = range(block_seq_pos[i][j]-1, (block_seq_pos[i][j]) + w*wlcut)
		print "%-10s %-9s %s" %(seqs_headers[j], block_seq_pos[i][j], seqs_inputs[j][min(nucleotides):max(nucleotides)])


