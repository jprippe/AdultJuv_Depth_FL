setwd("~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/mcav/analysis_fst_new/")

genes=read.table("mcav_genes_jp.txt")
names(genes)=c("contig","start","end","gene")
head(genes)
gnames=read.table("mcav_gnames_jp.txt", sep = "\t", header = F, stringsAsFactors = F)
names(gnames)=c("gene","protein")
genes=merge(genes,gnames,by="gene",all.x=T)
genes$protein=as.character(genes$protein)
genes$protein[is.na(genes$protein)]="unknown"
nrow(genes[genes$protein!="unknown",])

fst12=read.table("p12.fst")
names(fst12)=c("contig","pos","a","b")
fst13=read.table("p13.fst")
names(fst13)=c("contig","pos","a","b")
fst14=read.table("p14.fst")
names(fst14)=c("contig","pos","a","b")
fst23=read.table("p23.fst")
names(fst23)=c("contig","pos","a","b")
fst24=read.table("p24.fst")
names(fst24)=c("contig","pos","a","b")
fst34=read.table("p34.fst")
names(fst34)=c("contig","pos","a","b")

fst12$contig=as.character(fst12$contig)
fst13$contig=as.character(fst13$contig)
fst14$contig=as.character(fst14$contig)
fst23$contig=as.character(fst23$contig)
fst24$contig=as.character(fst24$contig)
fst34$contig=as.character(fst34$contig)

# removing zero-only bases
fst12[,3:4]=round(fst12[,3:4],3)
ch12=apply(fst12[,3:4],1,sum)
fst13[,3:4]=round(fst13[,3:4],3)
ch13=apply(fst13[,3:4],1,sum)
fst14[,3:4]=round(fst14[,3:4],3)
ch14=apply(fst14[,3:4],1,sum)
fst23[,3:4]=round(fst23[,3:4],3)
ch23=apply(fst23[,3:4],1,sum)
fst24[,3:4]=round(fst24[,3:4],3)
ch24=apply(fst24[,3:4],1,sum)
fst34[,3:4]=round(fst34[,3:4],3)
ch34=apply(fst34[,3:4],1,sum)

chh=(ch12>0 | ch13>0 | ch14>0 | ch23>0 | ch24>0 | ch34>0)
table(chh)
fst12=fst12[chh,]
fst13=fst13[chh,]
fst14=fst14[chh,]
fst23=fst23[chh,]
fst24=fst24[chh,]
fst34=fst34[chh,]

# computing weighted Fst per gene
i=1
gfst12=c();gfst13=c();gfst14=c();gfst23=c();gfst24=c();gfst34=c()
ns12=c();ns13=c();ns14=c();ns23=c();ns24=c();ns34=c()
pb=txtProgressBar(0,nrow(genes))
for (i in 1:nrow(genes)) {
	setTxtProgressBar(pb,i)
	sub=subset(fst12,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
	  gfst12=append(gfst12,NA)
	} else {
	  gfst12=append(gfst12,sum(sub$a)/sum(sub$b))
	  ns12=ns12+nrow(sub)
	}
	sub=subset(fst13,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
	  gfst13=append(gfst13,NA)
	} else {
	  gfst13=append(gfst13,sum(sub$a)/sum(sub$b))
	  ns13=ns13+nrow(sub)
	}
	sub=subset(fst14,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
	  gfst14=append(gfst14,NA)
	} else {
	  gfst14=append(gfst14,sum(sub$a)/sum(sub$b))
	  ns14=ns14+nrow(sub)
	}
	sub=subset(fst23,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
	  gfst23=append(gfst23,NA)
	} else {
	  gfst23=append(gfst23,sum(sub$a)/sum(sub$b))
	  ns23=ns23+nrow(sub)
	}
	sub=subset(fst24,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
	  gfst24=append(gfst24,NA)
	} else {
	  gfst24=append(gfst24,sum(sub$a)/sum(sub$b))
	  ns24=ns24+nrow(sub)
	}
	sub=subset(fst34,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
	if (is.null(sub[1,1]) | sum(sub$b)==0) { 
	  gfst34=append(gfst34,NA)
	} else {
	  gfst34=append(gfst34,sum(sub$a)/sum(sub$b))
	  ns34=ns34+nrow(sub)
	}
}

plot(density(na.omit(gfst13)),xlim=c(-0.1,0.5))
lines(density(na.omit(gfst12)),col="red")
lines(density(na.omit(gfst14)),col="red")
lines(density(na.omit(gfst23)),col="red")
lines(density(na.omit(gfst24)),col="red")
lines(density(na.omit(gfst34)),col="red")
abline(v=0, lty=2)
genes$fst12=gfst12
genes$fst13=gfst13
genes$fst14=gfst14
genes$fst23=gfst23
genes$fst24=gfst24
genes$fst34=gfst34
save(genes,file="PerGeneFst.RData")

require(reshape2)
require(ggplot2)
genes.melt <- melt(genes, id.vars = c("gene","contig","start","end","protein"), variable.name = "pair", value.name = "fst")
ggplot(data=genes.melt, aes(fst, color = pair, fill = pair))+
  theme_bw()+
  geom_density(alpha = 0.1)+
  coord_cartesian(xlim = c(0,0.3), ylim = c(0,70), expand = FALSE)

plot(gfst14~gfst34,pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.1))
outliers=identify(gfst14, gfst34,n=10)
genes[outliers,]

f12=genes[,c("gene","fst12")]
f13=genes[,c("gene","fst13")]
f14=genes[,c("gene","fst14")]
f23=genes[,c("gene","fst23")]
f24=genes[,c("gene","fst24")]
f34=genes[,c("gene","fst34")]
write.csv(f12,file="f12.csv",quote=F,row.names=F)
write.csv(f13,file="f13.csv",quote=F,row.names=F)
write.csv(f14,file="f14.csv",quote=F,row.names=F)
write.csv(f23,file="f23.csv",quote=F,row.names=F)
write.csv(f24,file="f24.csv",quote=F,row.names=F)
write.csv(f34,file="f34.csv",quote=F,row.names=F)

gos=read.table("../../GO_MWU/mcav_go_jp.txt", sep="\t", stringsAsFactors = F)
gos$V2[gos$V2==""]="unknown"
write.table(gos, file = "../../GO_MWU/mcav_go_jp.tab", sep = "\t", quote = F, row.names = F)

head(gos)

table(f12$gene %in% gos$V1)
