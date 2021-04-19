#~~~~~~
# Which MCAV lineage is most closely related to the reference genome? (Figure S10)
#~~~~~~

load('set1/mcav/datt.RData')
pops <- read.table('set1/mcav/newClusters_k4_mcav.tab', header = T)

datt.summ <- read.table('set1/mcav/bams_noclones') %>%
  rename(bams = V1) %>%
  mutate(bams = gsub(".fastq.bt2.bam", "", bams),
         sum = rowSums(datt), avg = rowMeans(datt)) %>%
  left_join(pops, by = 'bams') %>%
  filter(admix %in% c(1,2,3,4))
datt.summ$admix.names <- 'Near'
datt.summ$admix.names[datt.summ$admix == 2] <- 'Off'
datt.summ$admix.names[datt.summ$admix == 3] <- 'Deep1'
datt.summ$admix.names[datt.summ$admix == 4] <- 'Deep2'

scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
}

geno.plot <- ggplot(datt.summ)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0.2,0.1,0.2,0.7),"cm"))+
  geom_boxplot(aes(x = admix.names, y = avg, group = admix.names))+
  labs(y = expression(paste("Derived alleles per locus \n (variant and invariant)")))+
  scale_y_log10(label = scientific_10)
geno.plot

cowplot::save_plot("set1/mcav/reference_similarity.pdf", plot = geno.plot, base_height = 3, base_width = 4)

pop1 <- read.table('set1/mcav/sfsSites_1.mafs.gz', header = T)
pop2 <- read.table('set1/mcav/sfsSites_2.mafs.gz', header = T)
pop3 <- read.table('set1/mcav/sfsSites_3.mafs.gz', header = T)
pop4 <- read.table('set1/mcav/sfsSites_4.mafs.gz', header = T)
pop1$pop <- 'p1'
pop2$pop <- 'p2'
pop3$pop <- 'p3'
pop4$pop <- 'p4'
allpops <- rbind(pop1, pop2, pop3, pop4)
allpops2 <- allpops[allpops$knownEM > 0,]

maf.plot <- ggplot(allpops2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_boxplot(aes(x = pop, y = knownEM))+
  labs(x = "Minor allele frequency")
maf.plot

cowplot::save_plot("reference_similarity_maf.pdf", plot = maf.plot, base_height = 3, base_width = 4)
