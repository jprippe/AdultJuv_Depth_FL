require(WGCNA)
require(data.table)
require(flashClust)
require(vegan)
require(qqman)
require(plyr)
require(dplyr)

#############
spp <- 'mcav' #ssid or mcav
ld.stat <- 'rEM' #rEM or DEM (r2 or D calculated with EM algorithm)
#############

setwd(paste0('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/',spp))
load('LDsquare_datt_traits.RData')

if(ld.stat == 'rEM'){rm(ldmat.sq.DEM)} else if(ld.stat == 'DEM'){rm(ldmat.sq.rEM)}

#ngsLD is removing sites for some reason...
datt <- datt[,-which(!(colnames(datt) %in% names(mafs)))]

# ---------------- pruning physically linked sites from LD table; retaining ones with higher allele freq

mindist=1000

#folding mafs
#here "minor" is defined as the lower frequency allele (as opposed to the non-reference allele)
#if the non-reference allele frequency is greater than 0.5 in the mafs file, we subtract from 1 to get the frequency of the "minor" allele
mafs.fold=mafs
mafs.fold[mafs.fold>0.5]=1-mafs.fold[mafs.fold>0.5]

chrom.tm=factor(sub("\\:.+","",sites))
pos.tm=as.numeric(as.character(sub(".+\\:","",sites)))
sites.order <- sites[order(as.numeric(chrom.tm), pos.tm)]
chrom=factor(sub("\\:.+","",sites.order))
pos=as.numeric(as.character(sub(".+\\:","",sites.order)))

#Create data frame of loci ordered by position within each chromosome
sites.df <- data.frame(sites = sites.order, 
                       chrom = sub("\\:.+","",sites.order), 
                       pos = as.numeric(as.character(sub(".+\\:","",sites.order))),
                       maf = mafs.fold[match(sites.order, names(mafs.fold))]) #which mafs dataset you choose here is consequential
sites.df$sites <- as.character(sites.df$sites)

#For de novo 2bRAD genomes, linkage filtering must account for SNPs on the same tag but there is no way of knowing the genomic distance between each tag
if(spp == 'ssid'){
  tags.df <- data.frame(beg = seq(1,max(sites.df$pos), by=36), end = seq(36,max(sites.df$pos)+36, by=36))
  tags <- c()
  for(i in 1:nrow(sites.df)){tags <- c(tags, which(tags.df$beg <= sites.df$pos[i] & tags.df$end >= sites.df$pos[i]))}
  sites.df$tag <- as.factor(tags)
  sites.df <- sites.df %>% group_by(chrom,tag) %>% mutate(max.maf = max(maf))
  prune.link <- colnames(datt)[which(sites.df$maf != sites.df$max.maf)]
} else if(spp == 'mcav'){
  #Add column of distance between each site
  sites.df <- sites.df %>% group_by(chrom) %>% mutate(diff = c(Inf, diff(pos)))
  links <- which(sites.df$diff < mindist)
  
  linkGrp <- function(x){
    if(!is.numeric(x)) x <- as.numeric(x)
    n <- length(x)
    y <- x[-1] != x[-n] + 1
    i <- c(which(y|is.na(y)),n)
    list(
      nsites = diff(c(0,i)), #number of sites in each linkage group
      begin = x[head(c(0,i)+1,-1)] - 1 #first index of each linkage group
    )
  }
  link.groups <- linkGrp(links)
  
  #Retain the site within each "linkage" group with the highest allele frequency
  #All others are filtered out of the dataset
  prune.link <- c()
  for(i in 1:length(link.groups$nsites)){
    prune.link <- prune.link
    link.grp <- seq(0,link.groups$nsites[i]) + link.groups$begin[i]
    sites.grp <- sites.df[link.grp,]
    max.maf <- which.max(sites.grp$maf)
    prune.link <- c(prune.link, as.character(sites.grp$sites[-max.maf]))
  }
}

#remove monomorphic sites (this should be NULL, all majority-heterozygote sites should be removed using HetMajorityProb.py script)
prune.mono <- names(which(apply(datt, 2, var, na.rm=TRUE) == 0))

prune.all <- c(prune.link, prune.mono)
message(length(prune.all)) #1465 (mcav), 2872 (ssid)
if(ld.stat == 'DEM'){
  ldmat.filt <- ldmat.sq.DEM[-which(colnames(ldmat.sq.DEM) %in% prune.all),-which(colnames(ldmat.sq.DEM) %in% prune.all)]
} else if (ld.stat == 'rEM'){
  ldmat.filt <- ldmat.sq.rEM[-which(colnames(ldmat.sq.rEM) %in% prune.all),-which(colnames(ldmat.sq.rEM) %in% prune.all)]
}
mafs.fold=mafs.fold[-which(names(mafs.fold) %in% prune.all)]
mafs=mafs[-which(names(mafs) %in% prune.all)]
datt=datt[,-which(colnames(datt) %in% prune.all)]

#--------------- WGCNA

# # TOM and initial clustering
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5)
# 
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Lowest power for which the the scale-free topology fit index reaches 0.90
power=1
ldmat.filt[ldmat.filt==Inf]=1
TOM = TOMsimilarity(ldmat.filt^power, TOMType="unsigned");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

minModuleSize = 20; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
data.frame(table(dynamicColors))
# Label 0 (grey) represents unassigned genes

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0,
                    addGuide = FALSE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes (first principal component of geno matrix subset by each gene module)
# eg <- datt[,dynamicMods==2]
# eg.pca <- rda(eg, scale = T)

MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
MEs$MEgrey=NULL
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average")

MEDissThres = 0.5 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
#pdf("clusterModuleEigen.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
# plotting the fabulous ridiculogram
pdf(paste0("./wgcna_",ld.stat,"/",spp,"_geneModuleTree.pdf"), width = 7.5, height = 5.5)
plotDendroAndColors(geneTree,mergedColors,
                    dendroLabels = FALSE, hang = 0.0,
                    addGuide = F, guideHang = 0.05,lwd=0.6,main=paste0("merged at r = ", MEDissThres),autoColorHeight=F,colorHeight=0.15)
dev.off()
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEs$MEgrey=NULL

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
#	sizeGrWindow(7, 6)
plot(METree, main = "mc_all",
     xlab = "", sub = "")
dev.off()

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs$MEgrey=NULL

# correlations of genes with eigengenes (same as signedkME output between MEs and datt, aka "module membership")
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

# ------------ module-trait correlations

Traits=traits[,-c(1:3)]
Traits$isDA <- as.numeric(traits$site=="D" & traits$age=="A")
Traits$isDJ <- as.numeric(traits$site=="D" & traits$age=="J")
moduleTraitCor = cor(MEs, Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
head(MEs)

pdf(paste0("./wgcna_",ld.stat,"/",spp,"_MEheatmap.pdf"), width = 8, height = 4)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3));
# reorder table according to previous order
#	moduleTraitCor= moduleTraitCor[mtrows,]

# Display the correlation values within a heatmap plot
# Positive correlation indicates elevated frequency of derived (minor) alleles
# Negative correlation indicates elevated frequency of ancestral (major/reference) alleles
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(Traits),
               #    yLabels = names(MEs),
               yLabels = sub("ME","",row.names(moduleTraitCor)),
               ySymbols = row.names(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               #		cex.text = 0.0001,
               zlim = c(-1,1),
               main=paste0(spp,"_all"))
dev.off()
#	mtrows=row.names(moduleTraitCor)
print(data.frame(table(moduleColors))) # gives numbers of genes in each module
# MCAV
# moduleColors Freq
#         blue 1832
#        brown  287
#    turquoise 3410
#       yellow  266
# SSID
# moduleColors Freq
#         blue 4282
#        brown  868
#        green   94
#    turquoise 5135
#       yellow  145

row.names(MEs)=row.names(traits)
colnames(MEs)=sub("ME","",colnames(MEs))
head(MEs)

geneModuleMembership = as.data.frame(signedKME(datt, MEs));
colnames(geneModuleMembership) = colnames(MEs)
save(datt,traits,Traits,MEs,METree,geneTree,moduleColors,moduleLabels,nSamples, nGenes, moduleGeneCor, moduleGenePvalue, moduleTraitPvalue, moduleTraitCor,geneModuleMembership,file=paste0("./wgcna_",ld.stat,"/",spp,"_wgcna_p1.RData"))

# Discriminant analysis of individuals based on module eigengenes
rr=rda(MEs,scale=F)
#biplot(rr, col=traits$k)
plot(rr$CA$u,pch=as.numeric(traits$k)-1, col="grey50", xlim=c(-0.25,0.2),ylim=c(-0.25,0.4))
arrowScale=6
mod.cols <- unique(moduleColors)
nmods=length(mod.cols)
arrows(rep(0,nmods),rep(0,nmods),rr$CA$v[c(1:nmods),1]/arrowScale,rr$CA$v[c(1:nmods),2]/arrowScale,length=0.1,lwd=2,col=mod.cols)
legend("bottomleft",pch=unique(as.numeric(traits$k)-1),legend=unique(traits$k),bty="n",title="cluster",cex=0.9)
dev.off()

conds2=cbind(traits,MEs)
conds2$k=factor(conds2$k)

pdf(paste0("./wgcna_",ld.stat,"/",spp,"_eigPopCluster.pdf"), width = 8, height = 6)
par(mfrow=c(2,3))
for(i in 1:nmods){plot(get(mod.cols[i])~k, conds2, main=mod.cols[i], mgp=c(2.1,1,0), ylab = "Eigengene")}
dev.off()

# Discriminant analysis of SNPs based on module membership/correlation
par(mfrow=c(1,1))
rr=rda(t(moduleGeneCor))
plot(rr$CA$eig)
arrowScale=20
axes2plot=c(1,2)
plot(rr$CA$u[,axes2plot], pch=16, cex=0.5, col=rgb(0,0,0,alpha=0.1), mgp=c(2.3,1,0), ylim = c(-0.05, 0.05))
arrows(rep(0,nmods), rep(0,nmods), rr$CA$v[c(1:nmods), axes2plot[1]]/arrowScale, rr$CA$v[c(1:nmods), axes2plot[2]]/arrowScale, length=0.1, lwd=2, col=mod.cols)

# -------------- Manhattan plot
require(qqman)
chrom.num <- as.factor(sub("\\:.+", "", rownames(geneModuleMembership)))
levels(chrom.num) <- seq(1:length(levels(chrom.num)))
chrom.num <- as.numeric(chrom.num)
par(mfrow=c(2,3))
#for(i in 1:nrow(moduleGenePvalue)){
#  mandf <- data.frame(SNP = seq(1, nGenes),
#                      CHR = chrom.num,
#                      BP = as.numeric(sub(".+\\:", "", rownames(geneModuleMembership))),
#                      P = -log(moduleGenePvalue[i,], 10))
#  manhattan(mandf, main = rownames(moduleGenePvalue)[i])
#}
dev.off()

# -------------- by-trait scatters

setwd(paste0('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/',spp))
load(paste0("./wgcna_",ld.stat,"/",spp,"_wgcna_p1.RData"))
library(WGCNA)
library(plotrix) # for color.scale
length(moduleColors)
data.frame(table(moduleColors))

nGenes = ncol(datt);
nSamples = nrow(datt);
# names (colors) of the modules
modNames = names(MEs)
geneModuleMembership = as.data.frame(signedKME(datt, MEs)); # same as moduleGeneCor
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); # same as moduleGenePvalue
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

for(whichTrait in c("is1","is2","is3","is4","isA","isD","isO","isN")){
  selTrait = as.data.frame(traits[,whichTrait]);
  names(selTrait) = whichTrait
  
  # correlation of genes to traits
  # (positive value suggests membership to trait category is correlated with more derived alleles at that site)
  geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
  names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
  
  pdf(paste0("./wgcna_",ld.stat,"/",spp, "_geneSigModuleCorr_", whichTrait, ".pdf"), width = 7.5, height = 5.5)
  par(mfrow=c(2,3))
  counter=0
  for(module in modNames){
    counter=counter+1
    if (counter>9) {
      quartz()
      par(mfrow=c(3,3))
      counter=1
    }
    column = match(module, modNames);
    moduleGenes = which(moduleColors==module );
    moduleMafs = mafs[moduleColors==module]; #use mafs here instead of mafs.fold to represent frequency of derived allele
    moduleGenes=moduleGenes[order(moduleMafs)]
    length(moduleGenes)
    moduleMafs=moduleMafs[order(moduleMafs)]
    
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       ylim=c(-1,1),xlim=c(0,1),
                       geneTraitSignificance[moduleGenes, 1],
                       #	abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste(module,"kME"),
                       ylab = paste("GS for", whichTrait),
                       #	col = rgb(0,0,0,alpha=0.1),
                       col=color.scale(moduleMafs,c(0,1),c(0.8,0),c(1,0),alpha=0.5),
                       pch=16,mgp=c(2.3,1,0)
    )
    abline(0,1,lty=3,col="grey50")
    abline(0,-1,lty=3,col="grey50")
  }
  
  # plotting color scale bar
  colscale=data.frame(cbind("color"=1,"maf"=seq(0,1,0.01)))
  plot(color~maf,colscale,pch=16,cex=3,col=color.scale(seq(0,1,0.01),c(0,1),c(0.8,0),c(1,0)))
  
  dev.off()
}

# Plot interpretation:
# isD trait:
# For loci in the yellow module (n = 298) with a HIGH minor allele frequency across all samples (redder)...
# more "yellow-ness" = stronger POSITIVE correlation between # of derived alleles and being a deep individual
# For loci in the yellow module (n = 298) with a LOW minor allele frequency across all samples (bluer)...
# more "yellow-ness" = stronger NEGATIVE correlation between # of derived alleles and being a deep individual

# For ancestral loci (LOW minor allele frequency) in the yellow module, deep samples are more ancestral
# For derived loci (HIGH minor allele frequency) in the yellow module, deep samples are more derived

#------------------ density plots of derived/ancestral alleles in modules

#setwd("~/Dropbox/keysMcav/mcav_dec2018/")
#ll=load("mc_wgcna_p1_sep11_2019.RData")

for(module in modNames){
  pdf(paste0("./wgcna_",ld.stat,"/",spp, "_alleleDensity_", module, ".pdf"), width = 7.5, height = 5.5)
  par(mfrow=c(3,3))
  kmeCutoff=0.25
  tt=traits[,c("is1","is2","is3","is4","isA","isD","isO","isN")]
  modNames = names(MEs)
  for(i in 1:length(tt)) {
    selTrait = as.data.frame(tt[,i]);
    whichTrait=names(tt)[i]
    geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
    moduleGenes = which(moduleColors==module);
    #	plot(density(geneTraitSignificance[moduleGenes,1]),main=module,xlim=c(-1,1))
    
    #retain top percentage of genes in specified module based on absolute value of module membership
    gm=abs(geneModuleMembership[moduleGenes, paste("MM",module,sep="")])
    moduleGenesTop=moduleGenes[which(gm>quantile(gm,1-kmeCutoff))]
    
    #subset geno data by NON-module genes and samples associated with the trait
    moduledatt0=datt[as.logical(tt[,whichTrait]),which(moduleColors!=module)]
    #or subset geno data by module genes and samples NOT associated with the trait
    #moduledatt0=datt[!as.logical(tt[,whichTrait]),moduleGenes]
    
    #and calculate average number of derived alleles for each gene
    meangt0=apply(moduledatt0,2,mean)
    d0=density(meangt0)
    d0$y=d0$y/max(d0$y)
    plot(d0,col="cyan3",main=paste(module,":",whichTrait),xlab="average number of derived alleles",mgp=c(2.3,1,0),bty="n",yaxt="n",ylab="")
    #	mtext(side=2,"Normalized density",cex=0.5)
    
    #subset geno data by module genes and samples that ARE associated with the trait
    moduledatt=datt[as.logical(tt[,whichTrait]),moduleGenes]
    meangt=apply(moduledatt,2,mean)
    d1=density(meangt)
    d1$y=d1$y/max(d1$y)
    lines(d1,col="goldenrod")
    
    # subset geno data by module genes with top 25% module membership and samples that ARE associated with the trait
    moduledatt=datt[as.logical(tt[,whichTrait]),moduleGenesTop]
    meangt=apply(moduledatt,2,mean)
    d1=density(meangt)
    d1$y=d1$y/max(d1$y)
    lines(d1,col="red")
    
    #	abline(v=0,lty=3)
  }
  plot(d0,col="grey50",type="n",xlab="",mgp=c(2.3,1,0),bty="n",yaxt="n",xaxt="n",ylab="",main="")
  legend("topright",lty=1,col=c("red","goldenrod","cyan3"),lwd=2,legend=c("top 25% kME genes","all genes in module","genes not in module"),bty="n")
  
  dev.off()
}

# ---- module, kme per gene

setwd(paste0('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/',spp))
ll=load(paste0("./wgcna_",ld.stat,"/",spp,"_wgcna_p1.RData"))
genes=read.table("mcav_genes_clean.txt")
names(genes)=c("contig","start","end","gene")
genes$contig=as.character(genes$contig)
modNames = names(MEs)

# Generate list of each module's SNPs with associated contig, position and log10 p-value of module membership
coords=c()
for(module in modNames){
  moduleGenes = names(moduleGenePvalue[paste("ME",module,sep=""), moduleColors==module])
  coords[[module]]=ldply(strsplit(moduleGenes,":"))
  names(coords[[module]])=c("contig","pos")
  coords[[module]]$pos=as.numeric(coords[[module]]$pos)
  coords[[module]]$lpv=-log(moduleGenePvalue[paste("ME",module,sep=""),moduleColors==module],10)
  coords[[module]]$lpv=as.numeric(coords[[module]]$lpv)
}

# Generate data frame of module membership for functional genes (max membership if multiple SNPs in gene)
ngenes=nrow(genes)
pb=txtProgressBar(0,ngenes)
Gmax=c()
for (i in 1:ngenes) {
  setTxtProgressBar(pb,i)
  gmax=c()
  for(module in modNames){
    s=subset(coords[[module]],contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
    if (nrow(s)>0) { gmax=c(gmax,max(s$lpv)) } else { gmax=c(gmax,0)} #if there are multiple SNPs within a gene, keep the max log10 p-value of module membership
  }
  Gmax=rbind(Gmax,gmax)
}
Gmax=data.frame(Gmax)
names(Gmax)=modNames
row.names(Gmax)=genes$gene
table(Gmax[,1]>0) # 428 blue
table(Gmax[,2]>0) # 87 brown
table(Gmax[,3]>0) # 919 turquoise
table(Gmax[,4]>0) # 84 yellow
save(Gmax,file=paste0("./wgcna_",ld.stat,"/",spp,"_wgcna_modmemByGene.RData"))
head(Gmax)
dim(Gmax)
for(module in modNames){
  togo=data.frame(cbind("gene"=as.character(genes$gene),"moduleLPV"=Gmax[,module]))
  write.csv(togo,row.names=F,quote=F,file=paste("./wgcna_",ld.stat,"/",module,".csv",sep=""))
}

# --------------- examining sites in modules

source("~/Dropbox/Documents/ecogeno2018/manhattanPlot.R")

setwd(paste0('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/',spp))
ll=load(paste0("./wgcna_",ld.stat,"/",spp,"_wgcna_p1.RData"))

chrom=sub("\\:.+","",colnames(datt))
chrlen=table(chrom)
goods=names(chrlen)[chrlen>12]

module="turquoise";tops=0.75
ldmods=data.frame(cbind(chrom))
ldmods$pos=as.numeric(sub(".+\\:","",colnames(datt)))
moduleGenes = which(moduleColors==module)
ldmods$lpv=-log(moduleGenePvalue[paste("ME",module,sep=""),],10)
ldmods$gs=geneModuleMembership[,paste("kME",module,sep="")]
#ldmods$gs=geneModuleMembership[,'MMtur']
ldmods$gs[moduleColors!=module]=0
ldmods$sig=as.numeric(moduleColors==module)
gm=abs(geneModuleMembership[,paste("kME",module,sep="")])
topkme=(gm>quantile(gm[moduleGenes],tops) & ldmods$sig>0)
ldmods$sig2=as.numeric(topkme)
col=rgb(0.5,0.5,0.5,alpha=0.5)
dim(ldmods)

#quartz(height=3.5,width=6)
manhattanPlot(ldmods[ldmods$chrom %in% goods,],chr="chrom",colors=c(col,col,module),measure="gs",bp="pos",pch=16,hlite="sig2",main=paste(module,"mc_all",tops,sum(ldmods$sig2)),cex=1,smoothSpan=0)
abline(h=0,lty=3,lwd=2)

#----------------- saving data for GO analysis

setwd(paste0('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/',spp))
ll=load(paste0("./wgcna_",ld.stat,"/",spp,"_wgcna_p1.RData"))

# bads=which(is.na(apply(geneModuleMembership,1,mean)))
# geneModuleMembership= geneModuleMembership[-bads,]
# moduleColors= moduleColors[-bads]

#names(geneModuleMembership)=sub("MM","", names(geneModuleMembership))

coords=row.names(geneModuleMembership)
coo=c();i=1
for (i in 1:length(coords)) {
  co=strsplit(coords[i],"\\.")
  coo=rbind(coo,c(paste(co[[1]][1:3],collapse="|"),co[[1]][4]))
}

coo=data.frame(coo)
names(coo)=c("contig","pos")
coo$contig=as.character(coo$contig)
coo$pos=as.numeric(coo$pos)
str(coo)

blue=data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"blue"])))
blue[moduleColors!="blue",3]=0
brown=data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"brown"])))
brown[moduleColors!="brown",3]=0
turquoise =data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"turquoise"])))
turquoise[moduleColors!="turquoise",3]=0
yellow =data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"yellow"])))
yellow[moduleColors!="yellow",3]=0

head(blue)

plot(density(blue[,3]),ylim=c(0,2))
plot(density(turquoise[,3]),ylim=c(0,2))
plot(density(yellow[,3]),ylim=c(0,2))
plot(density(brown[,3]),ylim=c(0,2))
table(moduleColors)
# blue     brown turquoise    yellow 
# 790       570      1524       178 

write.table(blue,sep="\t",quote=F,file="blue_kme.tab")
write.table(brown,sep="\t",quote=F,file="brown_kme.tab")
write.table(turquoise,sep="\t",quote=F,file="turquoise_kme.tab")
write.table(yellow,sep="\t",quote=F,file="yellow_kme.tab")

# kme per gene
genes=read.table("mcav_genes.txt")
names(genes)=c("contig","start","end","gene")
genes$contig=as.character(genes$contig)
genes=genes[genes$contig %in% blue$contig,]

dat=brown
i=1;kme=c();ns01=0
pb=txtProgressBar(0,nrow(genes))
for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  sub=subset(dat,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
  if (nrow(sub)==0) { 
    kme=append(kme,NA)
  } else {
    kme=append(kme,max(sub$kme))
    ns01=ns01+nrow(sub)
  }
}
table(is.na(kme))
plot(density(na.omit(kme)))

res=data.frame(cbind("gene"=as.character(genes$gene),kme))
res=na.omit(res)
res$kme=as.numeric(as.character(res$kme))
res

turquoise.bygene=res
blue.bygene=res
yellow.bygene=res
brown.bygene=res


bygene=data.frame(cbind("blue"=blue.bygene$kme,"yellow"=yellow.bygene$kme,"turquoise"= turquoise.bygene$kme,"brown"=brown.bygene$kme))
row.names(bygene)=blue.bygene$gene

pairs(bygene,pch=16,col=rgb(0,0,0,alpha=0.5),cex=0.7)

save(bygene,file="bygene_kme.RData")

write.csv(turquoise.bygene,file="turquoise.bygene.csv",row.names=F, quote=F)
write.csv(blue.bygene,file="blue.bygene.csv",row.names=F, quote=F)
write.csv(yellow.bygene,file="yellow.bygene.csv",row.names=F, quote=F)
write.csv(brown.bygene,file="brown.bygene.csv",row.names=F, quote=F)

#------------

setwd("~/Dropbox/keysMcav/mcav_dec2018/")
ll=load("bygene_kme.RData")
head(bygene)
bygene$gene=row.names(bygene)

pcn=read.table("protComNumber.tab",sep="\t")
names(pcn)=c("hgene","interactions")
g2g=read.table("mcav_hgenes.tab",sep="\t")
names(g2g)=c("gene","hgene")


bygene2=merge(bygene,g2g,by="gene") 
dim(bygene2)
bygene3=merge(bygene2,pcn,by="hgene",all.x=T)
bygene3$interactions[is.na(bygene3$interactions)]=0.1

plot(brown~interactions,bygene3)



#############################################################
########## Old function for filtering linked sites ##########
#############################################################

#NOTE: for this loop to work correctly, "sites" must be ordered by position within each chromosome
prune=c();d=1
for (i in 2:length(sites.order)) {
  tp=0
  #start the counter over when we hit a new chromosome
  if(chrom[i]!=chrom[i-d]) { d=1 }
  #if the loci are on the same chromosome and are separated by less than 1000bp...
  if (chrom[i]==chrom[i-d] & abs(pos[i]-pos[i-d]) < mindist) {
    #print the "sites" index followed by the scaffold ID, positions and allele frequencies of the two loci
    message(paste(i,":",chrom[i],pos[i-d],mafs.fold[sites.order[i-d]]," - ",pos[i],mafs.fold[sites.order[i]],d))
    #retain the locus with the higher allele frequency
    if (round(mafs.fold[sites.order[i]],2)<=round(mafs.fold[sites.order[i-d]],2)) {
      tp=i
      d=d+1 
    } else { 
      tp=i-d
    }
    #print which locus is being removed
    message(paste("     pruning",sites.order[tp]))
    prune=append(prune,tp)
  }
}

