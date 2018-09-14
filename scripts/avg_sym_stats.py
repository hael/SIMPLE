import sys
import numpy as np

class sym_entity:

	def __init__(self,file):
		self.n      = 0
		self.rank   = []
		self.pgrp   = []
		self.score  = []
		self.corr   = []
		self.zscore = []
		lines = open(file,'r').readlines()
		for line in lines:
			if line.startswith('RANK'):
				self.n += 1
				vals = line.strip().split()
				self.rank.append(float(vals[1]))
				self.pgrp.append(vals[3])
				self.score.append(float(vals[5]))
				self.corr.append(float(vals[7]))
				self.zscore.append(float(vals[9]))
		self.rank   = np.array(self.rank, dtype=np.float)
		self.pgrp   = np.array(self.pgrp)
		self.score  = np.array(self.score, dtype=np.float)
		self.corr   = np.array(self.corr, dtype=np.float)
		self.zscore = np.array(self.zscore, dtype=np.float)

class sym_avger:
	def __init__(self):
		self.collection = []

	def add(self, sym_entity):
		if len(self.collection) == 0:
			self.collection.append(sym_entity)
		else:
			if sym_entity.n != self.collection[0].n:
				print 'Inconsitent stats between inputs'
				sys.exit(-1)
			else:
				self.collection.append(sym_entity)

	def calc(self):
		print "PGRP     RANK    SCORE   CORR   ZSCORE"
		n       = self.collection[0].n
		nobjs   = len(self.collection)
		ranks   = np.zeros(n)
		scores  = np.zeros(n)
		zscores = np.zeros(n)
		corrs   = np.zeros(n)
		for ipgrp in range(n):
			pgrp = self.collection[0].pgrp[ipgrp]
			ranks   = np.zeros(nobjs)
			scores  = np.zeros(nobjs)
			zscores = np.zeros(nobjs)
			corrs   = np.zeros(nobjs)
			for ientity in range(nobjs):
				for jpgrp in range(self.collection[ientity].n):
					if pgrp == self.collection[ientity].pgrp[jpgrp]:
						ranks[ientity]   = self.collection[ientity].rank[jpgrp]
						scores[ientity]  = self.collection[ientity].score[jpgrp]
						corrs[ientity]   = self.collection[ientity].corr[jpgrp]
						zscores[ientity] = self.collection[ientity].zscore[jpgrp]
			print "%5s%8.2f%8.2f%8.2f%8.2f" % (pgrp, ranks.mean(), scores.mean(), corrs.mean(), zscores.mean())


files = []
nfiles = len(sys.argv)-1
if nfiles==1:
	print 'you do not need this script for on file'
	sys.exit(0)

avger = sym_avger()
for i in range(1,nfiles+1):
	file   = sys.argv[i]
	symobj = sym_entity(file)
	avger.add(symobj)
avger.calc()


