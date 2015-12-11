#!/usr/bin/env Rscript

# Arguments
args=commandArgs(trailingOnly = TRUE)
file=args[1]
outfile=args[2]

# Read in data
x=read.table(file,sep=",",header=T)
x=x[,-1]
drug_comb=levels(x[,2])

# Loop through all drug combinations
median_dc=0
conf_scores=0
out=c()
for( i in 1:length(drug_comb))
{
    scores=x[x[,2]==drug_comb[i],3]
    median_dc=as.numeric(median(scores))
    z=density(scores)
    f=approxfun(z$x, z$y, yleft = 0, yright = 0)
    conf_score=integrate(f, median_dc-2,median_dc+2)
    out=rbind(out,c(drug_comb[i],median_dc,conf_score$value))
}

write.table(out,file=outfile,sep=",",quote=F,row.names=F)
