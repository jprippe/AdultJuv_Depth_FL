bams <- read.table('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/mcav/bams_noclones')
bams[,1] <- gsub("\\..*","",bams[,1])
bams$id <- seq(0,104)
bams$id <- paste0('ind',bams$id)

setwd('~/Google Drive/Work/~Research/~Projects_Active/AdultJuv_Depth_FL/work/mcav/BEAST_tree/')
fasta <- read.table('mc1_tree.fasta', header = F)
fasta[,1] <- as.character(fasta[,1])
for(i in seq(1, nrow(fasta), by=2)) {fasta[i,1] <- paste0(">",bams$V1[which(bams$id == gsub(">","",as.character(fasta[i,1])))])}

write.table(fasta, 'mc1_tree_relabel.fasta', row.names = F, col.names = F, quote = F)
