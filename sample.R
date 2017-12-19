#sample R code to generate the top  panel of Fig. 5
#paper: https://www.nature.com/articles/srep44016
#download in advance, GSE51808_series_matrix.txt from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51808/matrix/ 
x <- read.csv("GSE51808_series_matrix.txt",sep="\t",comment.char="!")
pca <- prcomp(scale(x[,-1]))
P <- pchisq(rowSums(scale(pca$x[,2:3])^2),2,lower.tail=F)
#downlaod in advance, GPL13158-5065.txt from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13158
y <- read.csv("GPL13158-5065.txt",sep="\t",comment.char="#")
index <- p.adjust(P,"BH")<0.01
gene <- y[y[,1] %in% x[index,1],11]
gene <-gene[gene!=""]
gene <- unique(unlist(strsplit(as.character(gene),"//")))
gene <- gsub(" ","",gsub("/ ","",gene))
write.table(file="gene_3.csv",gene,sep="\t",quote=F,col.names=F,row.names=F)

#download in advance, GSE13052_series_matrix.txt from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13052/matrix/
x <- read.csv("GSE13052_series_matrix.txt",sep="\t",comment.char="!")
pca <- prcomp(scale(x[,-1]))
P <- pchisq(rowSums(scale(pca$x[,2:3])^2),2,lower.tail=F)
index <- p.adjust(P,"BH")<0.01
y <- read.csv("GPL2700.soft",sep="\t",comment.char="!")
refseq <- y[y[,1] %in% x[index,1],3]
 refseq2 <- unlist(lapply(strsplit(as.character(refseq),".",fixed=T),"[",1))
require(biomaRt)
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene <- getBM(attributes=c("hgnc_symbol", "refseq_mrna"), filters ="refseq_mrna" , values = refseq2, mart=ensembl)
write.table(file="gene_4.csv",gene,sep="\t",col.names=F,row.names=F,quote=F)

gene3 <- read.csv("gene_3.csv",header=F)
gene4 <- read.csv("gene_4.csv",header=F,sep="\t")
write.table(file="gene_3_4.csv",intersect(gene3[,1],gene4[,1]),col.names=F,row.names=F,quote=F)

#download in advance, GSE25001_series_matrix.txt from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25001/matrix/
x <- read.csv("GSE25001_series_matrix.txt",sep="\t",comment.char="!")
#downlaod in advance, GPL6104-11576.txt from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6104
y <- read.csv("GPL6104-11576.txt",sep="\t",comment.char="#")
gene <- read.csv("gene_3_4.csv",header=F)
ii <- x[,1] %in% y[y[,6] %in% gene[,1],1]

pca2 <- prcomp(scale(x[ii,-1]))

class <- read.csv("DENG1.csv",sep="\t",header=F)
class_lst <- strsplit(as.character(class[,2]),"_")
class2 <- toupper(unlist(lapply(class_lst,"[",2)))

TP <- toupper(unlist(lapply(strsplit(as.character(class[,2]),"_"),"[",2)))
DIS <- unlist(lapply(strsplit(as.character(class[,2]),"_"),"[",1))

PCC <- NULL
for (DIS0 in c("DSS", "UC dengue"))
{
for (TP0 in c( "ACUTE", "[0-1]", "DIS",    "FOLLOWUP"))
{
PCC <- rbind(PCC,colMeans(pca2$rotation[ TP==TP0 & DIS == DIS0,2:3]))
}
}


pdf(file="gene_5pp.pdf")
rs <- 200
plot(pca2$rotation[,2:3]*rs,col=factor(class2),pch=3,lwd=2)
points(pca2$x[,2:3],col=6,pch=4+i1+i2)
legend(-25,-1,levels(factor(class2))[c(1,4,2,3)],col=c(1,4,2,3),pch=3,title="Patients")
legend(8,-20,c("Interferon","heme biosynthesis 1","heme biosynthesis 2"),col=6,pch=c(4,6,5),title="Gene annotation")
segments(0,-100,0,100,lty=2)
segments(-100,0,100,0,lty=2)
arrows(-4,4.5,-10,9.5)
text(-3,7,"DSS")
arrows(1,-2.8,-2,-13)
par(mai=c(2.5,0.5,0.5,0.5))
lines(PCC[c(1:4),]*rs,lwd=2,col=5,lty=2)
points(PCC[c(1:4),]*rs,lwd=2,col=5,type="p",pch=16,lty=2,cex=2)
lines(PCC[c(5:8),]*rs,lwd=2,col=5,lty=1)
points(PCC[c(5:8),]*rs,lwd=2,col=5,type="p",pch=16,lty=1,cex=2)
points(PCC[c(1:4),]*rs,lwd=2,col=c(1,4,2,3),type="p",pch=16,lty=2)
points(PCC[c(5:8),]*rs,lwd=2,col=c(1,4,2,3),type="p",pch=16,lty=2)
legend(-25,-20,c("DSS","uncomplicated"),lty=2:1,lwd=2,col=5,title="Disease progression")
text(2,-14,"Uncomplicated")
dev.off()

