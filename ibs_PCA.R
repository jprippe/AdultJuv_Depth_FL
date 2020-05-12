require(vegan)
require(ggplot2)

#---------------
spp <- 'mcav' #choose species
npops <- 4 #choose number of populations
#---------------

dir <- paste0('set1/', spp, '/') # change this to where your scp'd files are
bams <- read.table(paste0(dir, "bams_noclones"))[,1] # list of bam files
if(spp=='ssid'){
  prefix <- 'ss'
  bams <- sub(".nosymbio.fastq.bam","",bams)
} else if (spp=='mcav'){
  prefix <- 'mc'
  bams <- sub(".fastq.bt2.bam","",bams)
  }

ibsMat <- as.matrix(read.table(paste0(dir,prefix,"1.ibsMat")))
dimnames(ibsMat) <- list(bams,bams)

hc <- hclust(as.dist(ibsMat),"ave")
plot(hc,cex=0.7) # clustering of samples by IBS (great to detect clones or closely related individuals)

if(spp=='ssid'){
  clones <- c('SSNJ20c-1','SSNJ20c-2','SSNJ20c-3','SSDA20c-1','SSDA20c-2','SSDJ20c-1','SSDJ19c-2','SSDJ19c-3','SSOJ20c-2','SSOJ20c-3','SSNA1')
} else if (spp=='mcav'){
  clones <- c("MCNA16","MCNA17","MCDA20c-1","MCDA20c-2","MCOA10c","MCDJ19c-2","MCDJ19c-3","MCNJ5","MCNA10c","MCOJ5","MCNJ15","MCOJ15c")
}
goods <- which(!(bams %in% clones))
clonn <- which((bams %in% clones))
ibsMat <- ibsMat[goods,goods]
bams <- bams[goods]
if(spp=='ssid'){
  bamlist <- data.frame(paste(bams,".nosymbio.fastq.bam",sep=""))
} else if (spp=='mcav'){
  bamlist <- data.frame(paste(bams,".fastq.bt2.bam",sep=""))
}
write.table(bamlist,file=paste0(dir,"bams_noclones"),row.names=F,col.names=F,quote=F) # use this bam list to rerun ANGSD for ADMIXTURE

#==================================

# subsetters:
adults=grep("A",bams)
juveniles=grep("J",bams)
near=grep("N",bams)
off=grep("O",bams)
deep=grep("D",bams)

# colors:
cols=rep("royalblue",length(bams))
cols[off]="gold"
cols[near]="green3"

# factors:
site=rep("deep",length(bams))
site[off]="offshore"
site[near]="inshore"
aj=rep("adult",length(bams))
aj[juveniles]="juvenile"

# loading admixture clusters
load(paste0(dir, prefix, "_k", npops, ".qopt_clusters.RData"))

# significance of by-site divergence
conds=data.frame(cbind(site))
adonis(ibsMat~site,conds)

# performing PCoA and CAP
conds=data.frame(cbind(site))
pp0=capscale(ibsMat~1)
pp=capscale(ibsMat~site,conds)

# eigenvectors
plot(pp0$CA$eig)

axes2plot=c(1,2)  
#plot(pp0$CA$u[adults,axes2plot],pch=site[adults],col=cols[adults],cex=1.5)
#points(pp0$CA$u[juveniles,axes2plot],pch=site[juveniles],col=cols[juveniles])
#ordispider(pp0$CA$u[,axes2plot],groups=cluster.admix,col="grey80",label=T)

pp0.admix <- data.frame(id=bams, age=substr(bams,4,4), habitat=substr(bams,3,3), admix=cluster.admix, pp0$CA$u)
admix.spider <- merge(pp0.admix, aggregate(cbind(mean.x=pp0$CA$u[,axes2plot[1]], mean.y=pp0$CA$u[,axes2plot[2]])~admix, pp0.admix, mean), by='admix')
ggplot(data=admix.spider)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.title = element_blank())+
  geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(pp0$CA$u)[axes2plot[1]], yend=colnames(pp0$CA$u)[axes2plot[2]]), col='grey90')+
  geom_point(aes_string(x=colnames(pp0$CA$u)[axes2plot[1]], y=colnames(pp0$CA$u)[axes2plot[2]], shape='age', col='habitat'), size = 2)+
  geom_label(aes(x=mean.x, y=mean.y, label=admix))+
  scale_color_manual(values = c('skyblue1', 'palegreen1', 'plum1'))
ggsave(paste0(dir, prefix, 'PCA_k', npops, '_admix_fig1.png'), height = 3, width = 4)

#--------------
# k-means clustering (instead of admixture-based)
kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
# This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
# Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering

ks <- kmeans(pp0$CA$u[,1:kmeans.axes], npops, nstart = 100, iter.max = 1000)
pp0.admix.ks <- cbind(pp0.admix, kmeans = ks$cluster)
ks.spider <- merge(pp0.admix.ks, aggregate(cbind(mean.x=pp0$CA$u[,axes2plot[1]], mean.y=pp0$CA$u[,axes2plot[2]])~kmeans, pp0.admix.ks, mean), by='kmeans')
ggplot(data=ks.spider)+
  theme_bw()+
  geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(pp0$CA$u)[axes2plot[1]], yend=colnames(pp0$CA$u)[axes2plot[2]]), col='grey90')+
  geom_point(aes_string(x=colnames(pp0$CA$u)[axes2plot[1]], y=colnames(pp0$CA$u)[axes2plot[2]], shape='age', col='habitat'), size = 2)+
  geom_label(aes(x=mean.x, y=mean.y, label=kmeans))+
  scale_color_manual(values = c('skyblue1', 'palegreen1', 'plum1'))
ggsave(paste0(dir, prefix, 'PCA_k', npops, '_kmeans.png'), height = 6, width = 8)

#plot(pp0$CA$u[,axes2plot],pch=aj,col=cols)
#ordispider(pp0$CA$u[,axes2plot],groups=ks$cluster,col="grey80")
#
#plot(pp0$CA$u[,axes2plot],pch=site,col=cols)
#ordihull(pp0$CA$u[,axes2plot],groups=ks$cluster,label=T)
#
#plot(pp0$CA$u[,axes2plot],pch=aj,col=cols)
#ordiellipse(pp0$CA$u[,axes2plot],groups=ks$cluster,label=T,draw="polygon",col=1:4)

newClusters=data.frame(cbind(bams,"kmeans"=ks$cluster,"admix"=cluster.admix))
write.table(newClusters,quote=F,row.names=F,file=paste0(dir,"newClusters_k", npops, "_", spp, ".tab"))

# 3D PCA in plot.ly -------------------------------------------------------

library(plotly)

fig <- plot_ly(pp0.admix, x = ~MDS1, y = ~MDS2, z = ~MDS3, mode = 'markers', color = ~habitat, symbol = ~age, 
               colors = c('skyblue1', 'palegreen1', 'plum1'), marker = list(size=4)) %>% 
  add_markers() %>% 
  layout(scene = list(xaxis = list(title = 'MDS1'),
                                   yaxis = list(title = 'MDS2'),
                                   zaxis = list(title = 'MDS3'),
                                   camera = list(eye = list(x = 1, y = 1, z = 1))))
p <- plotly_build(fig)
p$x$data[[2]]$marker$symbol <- 'diamond'
p

# Output a rotational image series to create a gif
for(i in seq(0,6.3,by=0.1)){
  outfile <- paste("PCA",round(i,digits=2), sep = "_")
  cam.zoom = 2
  ver.angle = 0
  fig <- plot_ly(pp0.admix, x = ~MDS1, y = ~MDS2, z = ~MDS3, 
                 mode = 'markers', color = ~habitat, symbol = ~age, 
                 colors = c('skyblue1', 'palegreen1', 'plum1'), marker = list(size=4)) %>% 
    add_markers() %>% 
    layout(scene = list(xaxis = list(title = 'MDS1'),
                        yaxis = list(title = 'MDS2'),
                        zaxis = list(title = 'MDS3'),
                        camera = list(eye = list(x = cos(i)*cam.zoom,y = sin(i)*cam.zoom, z=0.2),
                                      center = list(x = 0, y = 0, z = 0))))
  p <- plotly_build(fig)
  p$x$data[[2]]$marker$symbol <- 'diamond'
  p
  
  cat("Now rendering iteration:", i,"\n")
  orca(fig, paste(outfile,"png", sep="."), width = 1200, height = 1050)
}

