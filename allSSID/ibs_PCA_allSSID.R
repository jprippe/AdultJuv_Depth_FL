require(vegan)
require(ggplot2)

#---------------
npops <- 3 #choose number of populations
#---------------

setwd(paste0('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/allSSID/', spp)) # change this to where your scp'd files are
bams <- read.table("bams_all")[,1] # list of bam files
bams_prefix <- sub(".nosymbio.fastq.bam","",bams)
samples <- read.csv('sampleNames.csv', header = T)
samples$sampleID <- as.character(samples$sampleID)
samples$sampleName <- as.character(samples$sampleName)
bams2 <- as.character(samples$sampleName[match(bams_prefix, samples$sampleID)])
bams2 <- c(bams2[!is.na(bams2)], bams_prefix[is.na(bams2)])
bams2 <- factor(bams2)

ibsMat <- as.matrix(read.table(paste0(prefix,"nc.ibsMat")))
dimnames(ibsMat) <- list(bams2,bams2)

hc <- hclust(as.dist(ibsMat),"ave")
plot(hc,cex=0.7) # clustering of samples by IBS (great to detect clones or closely related individuals)

#==================================

# subsetters:
adults=grep("A",bams2)
juveniles=grep("J",bams2)
near=grep("SSN",bams2)
off=grep("SSO",bams2)
deep=grep("SSD",bams2)

dg=grep("DG",bams2)
pg=grep("PG",bams2)
pl=grep("PL",bams2)
bc=grep("BC",bams2)
sp=grep("SP",bams2)

flk=c(near,off,deep)
bel=c(dg,pg,pl,bc,sp)

north=c(dg,bc,sp)
south=c(pl,pg)

ir=grep("IR",bams2)
or=grep("OR",bams2)

# colors:
cols=rep("red3",length(bams2))
cols[or]="dodgerblue3"

cols=rep("royalblue",length(bams2))
cols[off]="gold"
cols[near]="green3"

# factors:
site=rep("FL",length(bams2))
site[dg]="DG"
site[pg]="PG"
site[pl]="PL"
site[bc]="BC"
site[sp]="SP"
ns=rep("FL",length(bams2))
ns[north]="N"
ns[south]="S"
region=rep("FL",length(bams2))
region[bel]="Belize"

rz=rep("deep",length(bams2))
rz[c(off,or)]="offshore"
rz[c(near,ir)]="inshore"
aj=rep("unknown",length(bams2))
aj[adults]="adult"
aj[juveniles]="juvenile"

# loading admixture clusters
load(paste0(prefix, "_k", npops, ".qopt_clusters.RData"))

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

pp0.admix <- data.frame(id=bams2, region=region, rz=rz, age=aj, ns=ns, site=site, admix=cluster.admix, pp0$CA$u)
admix.spider <- merge(pp0.admix, aggregate(cbind(mean.x=pp0$CA$u[,axes2plot[1]], mean.y=pp0$CA$u[,axes2plot[2]])~admix, pp0.admix, mean), by='admix')
ggplot(data=admix.spider)+
  theme_bw()+
  #geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(pp0$CA$u)[axes2plot[1]], yend=colnames(pp0$CA$u)[axes2plot[2]]), col='grey90')+
  geom_point(aes_string(x=colnames(pp0$CA$u)[axes2plot[1]], y=colnames(pp0$CA$u)[axes2plot[2]], shape='rz', col='region', alpha='site'), size = 2)+
  #geom_label(aes(x=mean.x, y=mean.y, label=admix))+
  scale_alpha_manual(values = c(0.2,1,0.2,0.2,0.2,0.2))
  #scale_color_manual(values = c('skyblue1', 'palegreen1', 'plum1'))
ggsave(paste0(prefix, 'PCA_k', npops, '_admix.png'), height = 6, width = 8)

#--------------
# k-means clustering (instead of admixture-based)
kmeans.axes <- sum(diff(summary(eigenvals(pp0))[1,])*-1 > 0.01) #based on eigen plot

ks <- kmeans(pp0$CA$u[,seq(1:kmeans.axes)], npops, nstart = 100, iter.max = 1000)
pp0.admix.ks <- cbind(pp0.admix, kmeans = ks$cluster)
ks.spider <- merge(pp0.admix.ks, aggregate(cbind(mean.x=pp0$CA$u[,axes2plot[1]], mean.y=pp0$CA$u[,axes2plot[2]])~kmeans, pp0.admix.ks, mean), by='kmeans')
ggplot(data=ks.spider)+
  theme_bw()+
  geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(pp0$CA$u)[axes2plot[1]], yend=colnames(pp0$CA$u)[axes2plot[2]]), col='grey90')+
  geom_point(aes_string(x=colnames(pp0$CA$u)[axes2plot[1]], y=colnames(pp0$CA$u)[axes2plot[2]], shape='age', col='habitat'), size = 2)+
  geom_label(aes(x=mean.x, y=mean.y, label=kmeans))+
  scale_color_manual(values = c('skyblue1', 'palegreen1', 'plum1'))
ggsave(paste0(prefix, 'PCA_k', npops, '_kmeans.png'), height = 6, width = 8)

#plot(pp0$CA$u[,axes2plot],pch=aj,col=cols)
#ordispider(pp0$CA$u[,axes2plot],groups=ks$cluster,col="grey80")
#
#plot(pp0$CA$u[,axes2plot],pch=site,col=cols)
#ordihull(pp0$CA$u[,axes2plot],groups=ks$cluster,label=T)
#
#plot(pp0$CA$u[,axes2plot],pch=aj,col=cols)
#ordiellipse(pp0$CA$u[,axes2plot],groups=ks$cluster,label=T,draw="polygon",col=1:4)

newClusters=data.frame(cbind(bams,"kmeans"=ks$cluster,"admix"=cluster.admix))
write.table(newClusters,quote=F,row.names=F,file=paste0("newClusters_k", npops, "_", spp, ".tab"))

# 3D PCA in plot.ly -------------------------------------------------------

library(plotly)

fig <- plot_ly(admix.spider, x = ~MDS1, y = ~MDS2, z = ~MDS3, 
               mode = 'markers', color = ~region, symbol = ~rz, 
               colors = c('#F8766D','#00BFC4'), marker = list(size=4)) %>% 
  add_markers() %>% 
  layout(scene = list(xaxis = list(title = 'MDS1'),
                      yaxis = list(title = 'MDS2'),
                      zaxis = list(title = 'MDS3'),
                      camera = list(eye = list(x = cos(i)*cam.zoom,y = sin(i)*cam.zoom, z=0.2),
                                    center = list(x = 0, y = 0, z = 0))))
p <- plotly_build(fig)
p$x$data[[2]]$marker$symbol <- 'diamond'
p

#Output a rotational image series to create a gif
for(i in seq(0,6.3,by=0.1)){
  outfile <- paste("PCA",round(i,digits=2), sep = "_")
  cam.zoom = 2
  ver.angle = 0
  fig <- plot_ly(admix.spider, x = ~MDS1, y = ~MDS2, z = ~MDS3, 
                 mode = 'markers', color = ~region, symbol = ~rz, 
                 colors = c('#F8766D','#00BFC4'), marker = list(size=4)) %>% 
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


