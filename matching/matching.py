#!/usr/bin/python

import urllib

def scoring_matrix(filename):
	scoring = urllib.urlopen(filename).readlines()
	scoring = [i.decode("utf-8").strip('\n') for i in scoring[1:]]
	keys = [i[0] for i in scoring]
	scoring = [i.split()[1:] for i in scoring]
	scoring_dict = {}
	
	for ii,i in enumerate(keys):
		scoring_dict[i] = {}
		for ji,j in enumerate(keys):
			scoring_dict[i][j] = int(scoring[ii][ji])
	return scoring_dict
	

def global_alignment(v,w,penalty,matrix_file):
	score_matrix = scoring_matrix(matrix_file)
	s = [[0]*(len(w)+1) for i in range(len(v)+1)]
	s[0] = list(range(0,penalty*len(s[0]),penalty))
	
	for i in range(1,len(s)):
		s[i][0] = penalty + s[i-1][0]
	
	path = [['||'] + ['--']*(len(w)) for i in range(len(v)+1)]
	path[0][0] = '\\'
	
	for i in range(1,len(v)+1):
		for j in range(1,len(w)+1):
			score = score_matrix[v[i-1]][w[j-1]]
			s[i][j] = max(s[i-1][j-1] + score, s[i][j-1] + penalty, s[i-1][j] + penalty)
			
			if s[i][j] == s[i-1][j] + penalty:
				path[i][j] = "||"
			if s[i][j] == s[i][j-1] + penalty:
				path[i][j] = "--"
			if s[i][j] == s[i-1][j-1] + score:
				path[i][j] = "\\"
				
	score = 0
	align1 = ''
	align2 = ''
	
	while i >= 1 or j >= 1:
		if path[i][j] == "||":
			align1 += v[i-1]
			align2 += '-'
			score += penalty
			i -= 1
		elif path[i][j] == "--":
			align1 += '-'
			align2 += w[j-1]
			score += penalty
			j -= 1
        else:
			align1 += v[i-1]
			align2 += w[j-1]
			score += score_matrix[w[j-1]][v[i-1]]
			i -= 1
			j -= 1
	align1 = align1[::-1]
	align2 = align2[::-1]
	print('\n'.join([str(score),align1,align2]))
	return [score, align1,align2]
	
blosum62 = 'https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/34c4a34ce49d58b595c3fb2dc77b89a5b6a9b0af/blosum62.dat'
x = 'ISTHISALL'
y = 'ALIGNED'
penalty = -5
global_alignment(x,y,penalty,blosum62)
