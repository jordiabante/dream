#!/usr/bin/python 
# -*- coding: utf-8 -*-
import csv
import sys
import math
from scipy.stats.stats import pearsonr 
import numpy as np
with open('ch1_train_combination_and_monoTherapy.csv','rb') as csvfile1:
	reader = csv.DictReader(csvfile1, dialect='excel')
	combination_id = [row['COMBINATION_ID'] for row in reader]
with open('ch1_train_combination_and_monoTherapy.csv','rb') as csvfile2:
	reader = csv.DictReader(csvfile2, dialect='excel')
	compound_a = [row['COMPOUND_A'] for row in reader]
with open('ch1_train_combination_and_monoTherapy.csv','rb') as csvfile3:
	reader = csv.DictReader(csvfile3, dialect='excel')
	compound_b = [row['COMPOUND_B'] for row in reader]	
with open('Drug_info_release.csv','rb') as csvfile4:
	reader = csv.DictReader(csvfile4, dialect='excel')
	target = [row['Target(Official Symbol)'] for row in reader]
with open('Drug_info_release.csv','rb') as csvfile5:
	reader = csv.DictReader(csvfile5, dialect='excel')
	drug = [row['ChallengeName'] for row in reader]

targetline = []
target_count = 0




def targets():
	global targetline
	global target_count
	for i in range(0,len(drug)):
		s = target[i].split(',')
		length = len(s)
		if length == 1:
			targetline.append(s[0])
		else:
			for k in range(0,length):
				targetline.append(s[k])
	l = len(targetline)
	for i in range(0,l-1):
		for j in range(i+1,l):
			if targetline[i] == targetline[j]:
				targetline[j] = '0'
	j = 0
	for i in range(0,l):
		if targetline[i] != '0':
			targetline[j] = targetline[i]
			j = j+1
	while len(targetline) != j:
		del(targetline[j])
	for o in range(0,j):
		t = targetline[o]
		star_position = t.find('*')
		if star_position != -1:
			t = t[0:star_position]
			for p in range(0,j):
					if targetline[p].find(t) != -1 and o != p:
						targetline[o] = '0'

	k = 0
	index = 0
	while 1==1:
		if targetline[index] == '0':
			del(targetline[index])
			k = k+1
		else:
			index = index+1
			if index == j-k:
				break
	target_count = len(targetline)
targets()

def combination():
	global compound_a
	global compound_b
	for i in range(0,len(compound_a)-1):
		p = 1
		while compound_a[i] == compound_a[i+p]:
			if compound_b[i] == 0:
				break
			if compound_b[i] == compound_b[i+p]:
				compound_b[i+p] = 0 
			p = p+1
			if i+p == len(compound_a):
				break
	j = 0
	for i in range(0,len(compound_a)):
		if compound_b[i] != 0:
			compound_a[j] = compound_a[i]
			compound_b[j] = compound_b[i]
			j = j+1
	while len(compound_b) != j:
		del(compound_b[j])
		del(compound_a[j])
	return j


	

def ML():
	global compound_a
	global compound_b
	add = 0
	drug_target = np.zeros((len(drug),target_count),float)
	drug_combination = np.zeros((combination(),target_count),float)
	kernel_1 = np.zeros((combination(),combination()),float)
	kernel_2 = np.zeros((combination(),combination()),float)
	t=[]
	for i in range(0,len(drug)):
		v = target[i].split(',')
		for m in range(0,len(v)):
			if v[m].find('*') == -1:
				for n in range(0,target_count):
					if targetline[n] == v[m]:
						drug_target[i,n] = drug_target[i,n]+1
			else:
				group_name = v[m]
				group_name = group_name[0:v[m].find('*')]
				for n in range(0,target_count):
					if targetline[n].find(group_name) != -1:
						add = add+1
				for n in range(0,target_count):
					if targetline[n].find(group_name) != -1:
						drug_target[i,n] = drug_target[i,n]+(math.pow(2,add-1)/(math.pow(2,add)-1))
				add = 0
	for i in range(0,combination()):
		for m in range(0,len(drug)):
			if drug[m] == compound_a[i]:
				address_a = m
			elif drug[m] == compound_b[i]:
				address_b = m
		for n in range(0,target_count):
			drug_combination[i,n] = drug_target[address_a,n]+drug_target[address_b,n]
	csvfile_output = file('drug_comb_feature.csv', 'wb')
	writer = csv.writer(csvfile_output)
	data = drug_combination
	writer.writerows(data)
	csvfile_output.close()
	norm = 0
	for i in range(0,combination()):
		for m in range(0,target_count):
			norm = norm+math.pow(drug_combination[i,m],2)
		norm = math.sqrt(norm)

		for n in range(0,target_count):
			drug_combination[i,n] = drug_combination[i,n]/norm
		norm = 0
	summation = 0
	for m in range(0,combination()):
		for n in range(0,combination()):
			for i in range(0,target_count):
				summation = summation+drug_combination[m,i]*drug_combination[n,i]
			kernel_1[m,n] = summation
			summation = 0
	for m in range(0,combination()):
		for n in range(0,combination()):
			t = pearsonr(drug_combination[m,],drug_combination[n,])
			kernel_2[m,n] = t[0]
	csvfile_output = file('drug_kernel_dot.csv', 'wb')
	writer = csv.writer(csvfile_output)
	data = kernel_1
	writer.writerows(data)
	csvfile_output.close()
	csvfile_output = file('drug_kernel_cor.csv', 'wb')
	writer = csv.writer(csvfile_output)
	data = kernel_2
	writer.writerows(data)
	csvfile_output.close()
	csvfile_output = file('target_line.csv', 'wb')
	writer = csv.writer(csvfile_output)
	data = [(targetline)]
	writer.writerows(data)
	print len(targetline)
	





ML()
















