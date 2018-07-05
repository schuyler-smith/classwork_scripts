
def sw_code(position, word_code, w, wlcut):
	flatten = lambda l: [i for sl in l for i in sl]
	sw = position
	n = len(word_code)
	for wlev in range(wlcut):
		bucket = [[] for _ in range(0, (4**w)+1)]
		for i in range(n):
			pos = sw[i]-w
			if pos >= 0:
				code = word_code[pos]
				bucket[code+1].append(pos)
		for i in range(n-w, n):
			bucket[0].append(i)
		sw = flatten(bucket)
	return sw