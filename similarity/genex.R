## Based on:
## http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
## No gene expression available for the cell lines MDA-MB-175-VII and NCI-H1437
## Matrix 83 cell lines times 17,419 genes
## HUGO gene names
## Normalised cell line names
## RMA normalised

# Libraries
library('limma')

# Read in the data
x=read.table('gex.csv',sep=",",row.names=1,header=T)

# Generate design matrix
samples=as.factor(colnames(x))
design=model.matrix(~0 + samples)
colnames(design)=colnames(x)

# Fit the linear model to the filtered expression set
fit=lmFit(x, design)

# Set up a contrast matrix (comparisons we want to do)
contrast.matrix=makeContrasts("X22RV1-A549","X22RV1-C32",levels=design)

# Combine fit and contrast matrix
fit2=contrasts.fit(fit, contrast.matrix)

# Bayes
ebFit=eBayes(fit2)

# Top 10 genes
top_genes=topTable(ebFit, number=10, coef=1,, lfc=5)
