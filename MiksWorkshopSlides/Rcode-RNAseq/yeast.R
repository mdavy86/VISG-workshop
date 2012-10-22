yy<-read.csv('yeast-counts-chrom1.csv',row.names=1)
colnames(yy)<-c(paste("WT",1:3,sep=''),paste("MT",1:3,sep=''))

GType<-rep(c("WT","MT"),c(3,3))
Counts<-yy

# Limma analysis of count data

library(limma)
library(edgeR)

# All genes have a count of at least one
sort(apply(yy,1,sum),dec=T)

nf <- calcNormFactors(Counts)

design <- model.matrix(~GType)

y <- voom(Counts,design,plot=TRUE,lib.size=colSums(Counts)*nf)

plotMDS(y,top=10,labels=substring(GType,1,1),col=ifelse(GType=="WT","blue","red"),gene.selection="common")

fit <- lmFit(y,design)
fit <- eBayes(fit)
options(digits=3)
topTable(fit,coef=2,n=16)#[,-c(3,5,6,7,8)]
tt<-topTable(fit,coef=2,n=nrow(Counts),adjust.method='bonferroni')

sum(tt$adj.P.Val<0.05)

volcanoplot(fit,coef=2)
sig<-tt$adj.P.Val<0.05
points(tt$logFC[sig],tt$B[sig],pch=20,col='red',cex=0.5)

##########################################
# DESeq analysis of count data

library(DESeq)

cds <- newCountDataSet( Counts, GType )

cds <- estimateSizeFactors( cds )
sizeFactors( cds )

cds <- estimateDispersions( cds )

str( fitInfo(cds) )

plotDispEsts <- function( cds ) {
  plot(
       rowMeans( counts( cds, normalized=TRUE ) ),
       fitInfo(cds)$perGeneDispEsts,
       pch = '.', log="xy" )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}

plotDispEsts( cds )

res <- nbinomTest( cds, "WT", "MT" )

plotDE <- function( res ) {
  plot(
       res$baseMean,
       res$log2FoldChange,
       log="x", pch=20, cex=.3,
       col = ifelse( res$padj < .1, "red", "black" ) )
}

plotDE( res )

hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")

# DESeq "volcano plot"
plot(log2(res$foldChange),-log(res$pval),pch=20,cex=0.4)
sig<-res$padj<0.05
points(log2(res$foldChange[sig]),-log(res$pval[sig]),pch=20,cex=0.4,col='red')

#############################################################
# Compare limma and DESeq
limma.sig<-ifelse(p.adjust(fit$p.value[,2],"BH")<0.05,1,0)
deseq.sig<-ifelse(p.adjust(res$pval,"BH")<0.05,1,0)

vennDiagram(cbind(limma.sig,deseq.sig))

plot(-log(fit$p.value[,2]),-log(res$pval),col=rgb(limma.sig,rep(0,nrow(Counts)),deseq.sig),pch=16,cex=0.6)

