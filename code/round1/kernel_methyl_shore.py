#!/usr/bin/python 
# -*- coding: utf-8 -*-
import csv
import math
import pandas as pn
import numpy as np
from sklearn import preprocessing
from scipy.stats.stats import pearsonr
count_locus = 0
count_celllines = 0
with open('cell_line_order.csv','rb') as csvfile1:
	reader = csv.DictReader(csvfile1, dialect='excel')
	celllines = [row['x'] for row in reader]

def Kernel():
	dot = 0
	number=[]
	cell_len = len(celllines)
	kernel_cor = np.zeros((cell_len,cell_len),float)
	kernel_dot = np.zeros((cell_len,cell_len),float)
	shore = pn.read_csv('probe_info_shore_1.csv')
	probe = pn.read_csv('methyl_probe_m.csv')
	probe_len = len(probe['C32'])
	shore_len = len(shore['Name'])
	pro = probe.values
	sho = shore.values
	start = 0
	for j in range(0,shore_len):
		for i in range(start,probe_len):
			if pro[i,0] == sho[j,0]:
				start = i
				number.append(i)
				break
	print len(number)
	k = 0
	for i in range(0,probe_len):
		for j in range(0,len(number)):
			if i == number[j]:
				k = 1
		if k != 1:
			probe = probe.drop(i)
		k = 0
	del probe['Name']
	probe = probe.sort_index(axis=1)
	probe_len = len(probe['C32'])
	feature = np.zeros((cell_len,probe_len),float)
	x = probe.values
	for i in range(0,probe_len):
		for j in range(0,cell_len):
			feature[j,i]=x[i,j]
	probe_nor = preprocessing.normalize(feature, norm='l2')
	for m in range(0,cell_len):
		for n in range(0,cell_len):
			t = pearsonr(probe_nor[m,],probe_nor[n,])
			kernel_cor[m,n] = t[0]
	csvfile_output = file('methyl_shore_kernel_cor.csv', 'wb')
	writer = csv.writer(csvfile_output)
	data = kernel_cor
	writer.writerows(data)
	csvfile_output.close()


Kernel()