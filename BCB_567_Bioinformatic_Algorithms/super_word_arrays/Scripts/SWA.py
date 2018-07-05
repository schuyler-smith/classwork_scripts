import sys

from fasta_reader import read_fasta
from word_code import wrd_code
from super_word import sw_code
from w_code import wcode

if len(sys.argv) != 4:
	print "\n\nERROR: incorrect arguments given to script.\n\nUSAGE: python <SWA.py> <sequence.fasta> <word_model.txt> wlcut\n\nEXAMPLE: python SWA.py Sample_Data/seq_1.fasta Word_Code_1.txt 3\n\n"
	sys.exit(1)

seqs_headers, seqs_inputs = read_fasta(sys.argv[1])
seqs = seqs_inputs[0]
with open(sys.argv[2], 'r') as wm:
    word_model = wm.readline().strip()
wlcut = int(sys.argv[3])
w = len(word_model)
if w < 1:
	print "\n\nERROR: word_model must be of legnth >= 1"
	sys.exit(1)

word_code = wrd_code(seqs, word_model)
positions = range(0, len(seqs))
sup_wrd = sw_code(positions, word_code, w, wlcut)

start = []
max_start = []
prev = []
pos = []
for i in range(0, len(sup_wrd)):
	current = wcode(word_code, w, wlcut, sup_wrd[i])
	start.append(i)
	pos.append(sup_wrd[i-1]+1)
	if current == prev and -1 not in current:
		if len(start) > len(max_start):
			max_start = start
			max_pos = pos
			code_max = current
	else:
		start = []
		pos = []
	prev = current


print "\nWORD_MODEL:  %s" %(word_model)
print "WLCUT:  %s" %(wlcut)
print "LENGTH_LARGEST_BLOCK:  %s" %(len(max_start))
print "POSITION_BLOCK_IN_SEQUENCE:  %s" %(max_pos)
print "POSITION_BLOCK_IN_SW_ARRAY:  %s"%(max_start)
print "SUPERWORD_OF_BLOCK:  %s" %(code_max)



