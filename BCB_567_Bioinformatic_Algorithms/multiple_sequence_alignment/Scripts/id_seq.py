

def seq_id(index, start, length):

	for i in range(len(start)):
		if index in range(start[i], start[i]+length[i]):
			return (i+1, index-start[i])
	return(0)
