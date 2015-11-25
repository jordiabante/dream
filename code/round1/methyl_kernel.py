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

def ML():
	dot = 0
	length = len(celllines)
	kernel_cor = np.zeros((length,length),float)
	kernel_dot = np.zeros((length,length),float)
	df = pn.read_csv('methyl_ilse_m.csv')
	df = df.sort_index(axis=1)
	lon = len(df['C32'])
	g = np.zeros((length,lon),float)
	print df
	x = df.values
	for i in range(0,lon):
		for j in range(0,length):
			g[j,i]=x[i,j]
	df_nor = preprocessing.normalize(g, norm='l2')
	print df_nor
	for m in range(0,length):
		for n in range(0,length):
			t = pearsonr(df_nor[m,],df_nor[n,])
			kernel_cor[m,n] = t[0]
	csvfile_output = file('methyl_kernel_cor.csv', 'wb')
	writer = csv.writer(csvfile_output)
	data = kernel_cor
	writer.writerows(data)
	csvfile_output.close()
ML()