#!/usr/bin/env Rscript

# Load client
require(synapseClient)
synapseLogin()

# Retrieve drug_info_release.csv
file_entity <- synGet('syn4925534') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/drug_info_release.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve ch1_leaderBoard_monoTherapy.csv
file_entity <- synGet('syn4925538') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/ch1_leaderBoard_monoTherapy.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve ch1_test_monoTherapy.csv
file_entity <- synGet('syn4925540') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/ch1_test_monoTherapy.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve ch1_train_combination_and_monoTherapy.csv
file_entity <- synGet('syn4925542') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/ch1_train_combination_and_monoTherapy.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve ch2_leaderBoard_monoTherapy.csv
file_entity <- synGet('syn4925544') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/ch2_leaderBoard_monoTherapy.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve ch2_test_monoTherapy.csv
file_entity <- synGet('syn4925546') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/ch2_test_monoTherapy.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve gex.csv
file_entity <- synGet('syn4925562') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/gex.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve mutations.csv
file_entity <- synGet('syn4925564') 
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/mutations.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve cnv_segment.csv
file_entity <- synGet('syn4925568')
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/cnv_segment.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve methyl_ilse_m.csv
file_entity <- synGet('syn4925573')
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity[,2:ncol(read_file_entity)],file='../../data/originals/methyl_ilse_m.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve methyl_probe_m.csv
file_entity <- synGet('syn4925577')
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/methyl_probe_m.csv',
	row.names=F,sep=",",quote=FALSE)

# Retrieve probe_info.csv
file_entity <- synGet('syn4925579')
read_file_entity <- read.csv(getFileLocation(file_entity), header=T, na.strings=".")
write.table(read_file_entity,file='../../data/originals/probe_info.csv',
	row.names=F,sep=",",quote=FALSE)

