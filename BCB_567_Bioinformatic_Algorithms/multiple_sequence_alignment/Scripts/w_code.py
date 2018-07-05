def wcode(word_code, w, wlcut, i):
	wcode = []
	for x in range(wlcut):
		try:
			wcode.append(word_code[i+(x*w)])
		except:
			wcode.append(-1)
	return wcode

