library(lostruct)


# Statistics by chromosome ------------------------------------------------

chrom.lens <- c(Sc0000000=1870516, Sc0000001=1787856, Sc0000002=1766234, Sc0000003=1700957, Sc0000004=1572259)
chrom.stats <- sapply(paste0("Sc000000",0:4), function (chrom) {
  bcf.file <- file.path(sprintf("%s.vcf",chrom))
  sites <- vcf_positions(bcf.file)
  spacings <- diff(sites[[chrom]])
  out <- c( nsnps = length(sites[[chrom]]),
            nbp = chrom.lens[chrom],
            spacing.mean = mean(spacings),
            spacing = quantile(spacings,.05),
            spacing = quantile(spacings,.25),
            spacing = quantile(spacings,.50),
            spacing = quantile(spacings,.75),
            spacing = quantile(spacings,.95),
            spacing.max = max(spacings,.95) )
  return(out)
} )
floor(chrom.stats)

##               Sc0000000 Sc0000001 Sc0000002 Sc0000003 Sc0000004
## nsnps                24        38        32        40        34
## nbp.Sc0000000   1870516   1787856   1766234   1700957   1572259
## spacing.mean      69841     43282     53389     41288     45053
## spacing.5%          128         7         3         2         4
## spacing.25%       12375      2068      5577      2298        27
## spacing.50%       45956     38111     28482     16494     19854
## spacing.75%      101314     64057     82026     66626     51725
## spacing.95%      184704    116261    205382    137018    180844
## spacing.max      345622    230095    254066    151989    281578


# Recoding ----------------------------------------------------------------

get.indivs <- c("MCDA19.fastq.bt2.bam","MCNA12.fastq.bt2.bam","MCNA14.fastq.bt2.bam","MCNA18.fastq.bt2.bam")
bcf.file <- "Sc0000000.vcf"
# first read in the info about sites
bcf.con <- pipe(paste("bcftools query -f '%POS\\n'",bcf.file),open="r")
bcf.sites <- list(Sc0000000=scan(bcf.con))
close(bcf.con)
# now use bcftools to pull out sites we want
win.fn <- function (n, win.snps=100) {
  sites <- list(Sc0000000=bcf.sites[[1]][seq(1+(n-1)*win.snps,length.out=win.snps)])
  regions <- enclosing_region(sites)
  query_genotypes(bcf.file, regions, get.indivs)
}
win.fn(20,10)  # get the 20th window of 10 SNPs


# Test run on Sc0000001 ---------------------------------------------------

bcf.file <- "Sc0000001.vcf.gz"
sites <- vcf_positions(bcf.file)
do.indivs <- system(paste0("bcftools query -l ", bcf.file), intern = T)
win.fn.snp <- vcf_windower(bcf.file, size=8, type="snp", sites=sites, samples=do.indivs) 
snp.pca <- eigen_windows(win.fn.snp,k=2)
Sc0000001.pcdist <- pc_dist(snp.pca)

# there may be windows with missing data
na.inds <- is.na(snp.pca$values[1,])
Sc0000001.mds <- cmdscale(Sc0000001.pcdist, eig=TRUE, k=2)
mds.coords <- Sc0000001.mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))









