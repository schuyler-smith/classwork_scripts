import sys

from fasta_reader import read_fasta
from populate_matrices import populate_matrices, delta_calc
from trace_back import traceback

if len(sys.argv) != 6:
    print "\n\nERROR: incorrect arguments given to script.\n\nUSAGE: sh SDS_BLAST.sh <sequence_1.fasta> <sequence_2.fasta> mismatch_score gap_open_Penalty gap_extension_penalty\n\nEXAMPLE: sh SDS_BLAST.sh Data/seq_11.fasta Data/seq_22.fasta -20 40 2\n\n"
    sys.exit(1)

seqs1_inputs = str(sys.argv[1])
seqs2_inputs = str(sys.argv[2])
match = int(10)
mismatch = int(sys.argv[3])
gap_open = int(sys.argv[4])
gap_extend = int(sys.argv[5])

gap_init_penalty = gap_open + gap_extend

seqs1_headers, seqs1 = read_fasta(seqs1_inputs)
seqs2_headers, seqs2 = read_fasta(seqs2_inputs)

S, D, I, score, start_A, start_B  = populate_matrices(seqs1[0], seqs2[0], match, mismatch, gap_open, gap_extend)
align_1, align_2, end_A, end_B, num_matches, num_mismatches, num_gaps, middle_array = traceback(S, D, I, seqs1[0], seqs2[0], gap_init_penalty, start_A, start_B)


print "\n" + "SUCCESSFULLY ALIGNED SEQUENCES:".center(40) + "\n\n(A) %s  by  (B) %s \n" %(seqs1_headers[0], seqs2_headers[0])
print "LENGTH_OF_SEQUENCE_A = %s\nLENGTH_OF_SEQUENCE_B = %s\nLENGTH_OF_ALIGNMENT = %s" %(len(seqs1[0]), len(seqs2[0]), str(len(align_1)))
print "START_POS_A = %s\nSTART_POS_B = %s\nEND_POS_A = %s\nEND_POS_B = %s\n" %(start_A, start_B, end_A, end_B)
a_1 = []
a_2 = []
mid = []
s_a = start_A
s_b = start_B
num = 0
for i in range(0, int(len(align_1))):
	a_1.append(align_1[i])
	a_2.append(align_2[i])
	mid.append(middle_array[i])
	num += 1
	if num == 70 or i == int(len(align_1))-1:
		print "%4s %4s %s" %("(A)", s_a, ''.join(list(a_1)))
		print "%3s %5s %s" %("|", "", ''.join(list(mid)))
		print "%4s %4s %s\n" %("(B)", s_b, ''.join(list(a_2)))
		s_a = s_a + num - a_1.count('-')
		s_b = s_b + num - a_2.count('-')
		a_1 = []
		a_2 = []
		mid = []
		num = 0
print "SIMILARITY_SCORE = %s\nPERCENT_IDENTITY = %s" %(score, ((num_matches * 100)/int(len(align_1))))
print "NUM_OF_MATCHES = %s\nNUM_OF_MISMATCHES = %s\nLENGTH_OF_GAPS = %s\n" %(num_matches, num_mismatches, num_gaps)
print "### ALIGNMENT_PARAMETERS ###\n".center(40)
print "MATCH_SCORE = %s\nMISMATCH_SCORE = %s\nGAP_OPEN_PENALTY = %s\nGAP_EXTENSION_PENALTY = %s\n" %(match, mismatch, gap_open, gap_extend)

