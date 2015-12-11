#!/usr/bin/python 
# -*- coding: utf-8 -*-
import csv
import math
import pandas as pn
import numpy as np
from sklearn import preprocessing
from scipy.stats.stats import pearsonr




def Kernel():
	cell_len = 85
	kernel_cor = np.zeros((cell_len,cell_len),float)
	kernel_dot = np.zeros((cell_len,cell_len),float)
	kernel_ang = np.zeros((cell_len,cell_len),float)
	isle = pn.read_csv('methyl_islands.csv')
	isle_len = len(isle['C32'])
	isle = isle.sort_index(axis=1)
	feature = np.zeros((cell_len,isle_len),float)
	x = isle.values
	for i in range(0,isle_len):
		for j in range(0,cell_len):
			feature[j,i]=x[i,j]
	isle_nor = preprocessing.normalize(feature, norm='l2')
	for m in range(0,cell_len):
		for n in range(0,cell_len):
			t = isle_nor[m,].dot(isle_nor[n,])
			kernel_dot[m,n] = t
			if m == n: 
				t = 1
			else:
				t = 1-np.arccos(t)/np.pi
			kernel_ang[m,n] = t
			t = pearsonr(isle_nor[m,],isle_nor[n,])
			kernel_cor[m,n] = t[0]
	np.savetxt('../../data/round2/kernels/dot_product_methylation_islands.txt',kernel_dot, delimiter = '\t')
	np.savetxt('../../data/round2/kernels/corr_methylation_islands.txt',kernel_cor, delimiter = '\t')	
	np.savetxt('../../data/round2/kernels/angular_similarity_methylation_islands.txt',kernel_ang, delimiter = '\t')			
	
Kernel()