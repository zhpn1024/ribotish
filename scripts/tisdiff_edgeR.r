args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
n1 <- as.numeric(args[2])
n2 <- as.numeric(args[3])
outfile <- args[4]
sep <- "\t"
if(length(args)>4){
  sep <- args[5]
}

library(edgeR)

x <- read.delim(infile, row.names=1, sep=sep)
group <- factor(cbind(rep(1,n1),rep(2,n2)))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
#topTags(et)
table <- et$table
table$FDR <- p.adjust(table$PValue, method = "BH")
write.table(table, file=outfile, quote=F, sep='\t')

