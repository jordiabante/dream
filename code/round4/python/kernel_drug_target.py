#!/home/mostafa_karimi/anaconda2/bin/python 
import math
import csv
from scipy.stats.stats import pearsonr 
import numpy as np

with open('../../data/originals/ch1_train_combination_and_monoTherapy.csv','rb') as csvfile1: 
        reader = csv.DictReader(csvfile1, dialect='excel')
        combination_id = [row['COMBINATION_ID'] for row in reader]
with open('../../data/originals/ch1_train_combination_and_monoTherapy.csv','rb') as csvfile2:
        reader = csv.DictReader(csvfile2, dialect='excel')
        compound_a = [row['COMPOUND_A'] for row in reader]
with open('../../data/originals/ch1_train_combination_and_monoTherapy.csv','rb') as csvfile3:
        reader = csv.DictReader(csvfile3, dialect='excel')
        compound_b = [row['COMPOUND_B'] for row in reader]
with open('../../data/originals/drug_info_release.csv','rU') as csvfile4:
        reader = csv.DictReader(csvfile4, dialect='excel')
        target = [row['Target(Official Symbol)'] for row in reader]
with open('../../data/originals/drug_info_release.csv','rU') as csvfile5:
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
	kernel_3 = np.zeros((combination(),combination()),float)
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
	np.savetxt('../../data/round4/features/drug_comb_feature.txt',drug_combination, delimiter = '\t')
	norm = 0
	for i in range(0,combination()):
		for m in range(0,target_count):
			norm = norm+math.pow(drug_combination[i,m],2)
		norm = math.sqrt(norm)

		for n in range(0,target_count):
			drug_combination[i,n] = drug_combination[i,n]/norm
		norm = 0
	for m in range(0,combination()):
		for n in range(0,combination()):
			dot = drug_combination[m,].dot(drug_combination[n,])
			kernel_1[m,n] = dot
			if m == n: 
				dot = 1
			else:
				if dot > 1:
					dot = 1
				dot = 1-np.arccos(dot)/np.pi
			kernel_3[m,n] = dot
	for m in range(0,combination()):
		for n in range(0,combination()):
			t = pearsonr(drug_combination[m,],drug_combination[n,])
			kernel_2[m,n] = t[0]
	np.savetxt('../../data/round4/kernels/dot_product_drug_target.txt',kernel_1, delimiter = '\t')
	np.savetxt('../../data/round4/kernels/corr_drug_target.txt',kernel_1, delimiter = '\t')
	np.savetxt('../../data/round4/kernels/angular_similarity_drug_target.txt',kernel_1, delimiter = '\t')
	

print("kernel drug target completed")	





ML()
















