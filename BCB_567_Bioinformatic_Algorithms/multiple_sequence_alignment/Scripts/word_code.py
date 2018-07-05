

def wrd_code(seq, wrd_model):
	w = len(wrd_model)
	word_code = []
	for i in range(0, len(seq)):
		word = []
		try:
			for j in range(0, w):
				if int(wrd_model[j]) == 1:
					word.append(seq[i+j])
			word = int(''.join(word).replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3'),4)
		except:
			word = -1
		word_code.append(word)
	return word_code
