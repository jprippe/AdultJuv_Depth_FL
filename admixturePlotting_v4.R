#~~~~~~
# Code for plotting ngsAdmix and PCAngsd results (Figure 1, S2)
#~~~~~~

source("plot_admixture_v4_function.R")

#---------------
spp <- 'ssid' #choose species
npops <- 4 #choose number of populations
#---------------

# assembling the input table
dir=paste0("set2/2brad/", spp, "/") # path to input files
if(spp=='ssid'){
  prefix <- 'ss'
  inName <- paste0('ss_k', npops, '.qopt')
} else if (spp=='mcav'){
  prefix <- 'mc'
  inName <- paste0('mc_k', npops, '.qopt')
}
pops="inds2pops" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

#------------

npops=as.numeric(sub("\\D+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
##########
# OR WITH PCAngsd (Figure S2)
tbl=as.data.frame(RcppCNPy::npyLoad('set1/mcav/pcangsd/pcangsd.admix.4.Q.npy'))
##########
i2p=read.table(paste(dir,pops,sep=""),header=F)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
rownames(tbl) <- tbl$ind <- sub("(.*?)\\..*$", "\\1", tbl$ind)
head(tbl,20) # this is how the resulting dataset must look

# putting populations in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))

plotAdmixture(data=tbl,npops=npops,grouping.method="distance")
ggsave(paste0(dir,prefix, '_k', npops, '_admixPlot.png'), height = 2, width = 10)

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
save(cluster.admix,file=paste0(dir,inName,"_clusters.RData"))

