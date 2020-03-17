#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("geno file, ngsLD file, bams file and species (ssid or mcav) must be supplied in that order.\n", call.=FALSE)
}

require(tools)
require(data.table)

# ------------ genotypes

geno=read.table(args[1], sep="\t")
geno=geno[,-ncol(geno)]
row.names(geno)=paste(geno[,1],geno[,2],sep=":")
bams=read.table(args[3])[,1]
bams=sub(".fastq.bt2.bam","",bams)
# remove leader columns from geno file
cols.rm=ncol(geno) %% length(bams)
geno=geno[,-c(1:cols.rm)]

# depending on geno file format, return called genotypes or calculate weighted genotypes based on likelihoods
# split data frame into grouped columns by sample
indCols=ncol(geno)/length(bams)
index=split(1:ncol(geno), rep(1:ncol(geno), each=indCols, length=ncol(geno)))
# calculate weighted genotypes if file matches expected angsd output formats
  # -doGeno 6 or 7
if(indCols==2){
  geno=geno[,c(T,F)]
  # -doGeno 8 or 9
} else if(indCols==3){
  geno=data.frame(do.call(cbind, lapply(index, function(i) apply(geno[,i], 1, function(x) sum(x*c(0,1,2))))))
  # -doGeno 10, 11, 12, 13 or 15
} else if(indCols==4 | indCols==5){
  extra=indCols-3
  geno=data.frame(do.call(cbind, lapply(index, function(i) apply(geno[,i[-c(1:extra)]], 1, function(x) sum(x*c(0,1,2))))))
} else if(!(indCols %in% c(1:5))){
  stop('Check format of geno file! Incorrect number of columns detected.')
}
names(geno)=bams
datt=t(geno)
#datt=datt+0.01
datt[1:30,1:4]
dim(datt)

# ------------- traits

if (args[4] == "ssid"){
  c1=read.table("cluster1")[,1]
  c2=read.table("cluster2")[,1]
  c3=read.table("cluster3")[,1]
  c4=read.table("cluster4")[,1]
  bams=read.table(args[3])[,1]
  k=rep(1,length(bams))
  k[bams %in% c2]=2
  k[bams %in% c3]=3
  k[bams %in% c4]=4
  bams=sub(".nosymbio.fastq.bam","",bams)
  site=gsub("SS|[0-9]|[JA]|c-|_filt","",bams)
  age=gsub("SS|[0-9]|[NOD]|c-|_filt","",bams)
} else if (args[4] == "mcav"){
  c1=read.table("cluster1")[,1]
  c2=read.table("cluster2")[,1]
  c3=read.table("cluster3")[,1]
  c4=read.table("cluster4")[,1]
  bams=read.table(args[3])[,1]
  k=rep(1,length(bams))
  k[bams %in% c2]=2
  k[bams %in% c3]=3
  k[bams %in% c4]=4
  bams=sub(".fastq.bt2.bam","",bams)
  site=gsub("MC|[0-9]|[JA]|c|-","",bams)
  age=gsub("MC|[0-9]|[NOD]|c|-","",bams)
}
traits=data.frame(cbind(site,age,k))
row.names(traits)=bams
str(traits)
traits$is1=as.numeric(traits$k==1)
traits$is2=as.numeric(traits$k==2)
traits$is3=as.numeric(traits$k==3)
traits$is4=as.numeric(traits$k==4)
traits$isA=as.numeric(traits$age=="A")
traits$isD=as.numeric(traits$site=="D")
traits$isN=as.numeric(traits$site=="N")
traits$isO=as.numeric(traits$site=="O")

# ---------------- Reformatting ngsLD output into square matrix
LD.long <- fread(args[2], nThread = 20, data.table = F, showProgress = F)

for(i in c('V5','V7')){
  if(i == 'V5'){stat <- 'DEM'} else if(i == 'V7'){stat <- 'rEM'}
  ldtab <- dcast(LD.long, V2~V1, value.var = i)
  
  # Add empty row and column for the two sites that only occur in V1 or V2 in long-format table
  tempdf <- data.frame(matrix(rep(NA, ncol(ldtab)), nrow = 1)); colnames(tempdf) <- colnames(ldtab)
  ldtab2 <- rbind(tempdf, ldtab)
  ldtab2[ldtab$V2[which(!(ldtab$V2 %in% colnames(ldtab)[-1]))]] <- rep(NA, nrow(ldtab2))
  
  # Sort matrix into lower triangle of the square matrix
  ldmat <- as.matrix(ldtab2[,-1])
  ldmat.rsort <- ldmat[rev(order(rowSums(is.na(ldmat)))),]
  ldmat.csort <- ldmat.rsort[,order(colSums(is.na(ldmat.rsort)))]
  
  # Transpose and copy upper triangle
  tm <- t(ldmat.csort)
  ldmat.csort[upper.tri(ldmat.csort)] <- tm[upper.tri(tm)]
  diag(ldmat.csort) <- 1
  rownames(ldmat.csort) <- colnames(ldmat.csort)
  assign(paste0('ldmat.sq.', stat), ldmat.csort)
}

# ---------------- extracting minor allele freq from ngsLD output
ald <- LD.long
names(ald)=c("s1","s2","dist","r2","D1","D2","rEM","n","af1","af2","h1","h2","h3","h4","chi2")
ald2 <- ald
s2 <- ald2$s2
ald2$s2 <- ald2$s1
ald2$s1 <- s2
ald3 <- rbind(ald,ald2)
ald3$s1 <- as.character(ald3$s1)
ald3$s2 <- as.character(ald3$s2)
mafs <- ald3[match(unique(ald3$s1), ald3$s1),]$af1
names(mafs) <- sites <- unique(ald3$s1)

save(ldmat.sq.rEM, ldmat.sq.DEM, mafs, sites, datt, traits, file = 'LDsquare_datt_traits.RData')
