args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
n1 <- as.numeric(args[2])
n2 <- as.numeric(args[3])
outfile <- args[4]

library(edgeR)

x <- read.delim(infile, row.names=1)
group <- factor(cbind(rep(1,n1),rep(2,n2)))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
#topTags(et)
write.table(et$table, file=outfile, quote=F, sep='\t')

