#!/usr/bin/env Rscript

# Arguments
args=commandArgs(trailingOnly = TRUE)
file=args[1]
outfile=args[2]

# Read in data
x=read.table(file,sep="\t",header=T)
drug_comb=levels(x[,2])

# Loop through all drug combinations
median_dc=0
conf_scores=0
out=c()
for( i in 1:length(drug_comb))
{
    scores=as.numeric(as.vector(x[x[,2]==drug_comb[i],3]))
    median_dc=as.numeric(median(scores))
    z=density(scores)
    f=approxfun(z$x, z$y, yleft = 0, yright = 0)
    conf_score=integrate(f, median_dc-1,median_dc+1)
    if(conf_score$value>1){conf_score$value=1}
    out=rbind(out,c(drug_comb[i],median_dc,conf_score$value))
}

write.table(out,file=outfile,sep=",",quote=F,row.names=F)
