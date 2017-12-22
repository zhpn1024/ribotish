args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
n1 <- as.numeric(args[2])
n2 <- as.numeric(args[3])
outfile <- args[4]
sep <- "\t"
if(length(args)>4){
  sep <- args[5]
}
library(DESeq2)

x <- read.delim(infile, row.names=1, sep=sep)
condition <- c(rep("c1",n1), rep("c2",n2))
coldata <- as.data.frame(condition)
dds <- DESeqDataSetFromMatrix(countData=x, colData=coldata, design= ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","c2","c1"))

write.table(as.data.frame(res), file=outfile, quote=F, sep='\t')

