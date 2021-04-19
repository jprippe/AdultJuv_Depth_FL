#~~~~~~
# Fst tree (Figure 2A)
#~~~~~~

library(tidyverse)
library(ggtree)
library(ape)
library(cowplot)

fst.mcav <- as.matrix(read.table('set1/mcav/populationFst_mcav.txt'))
fst.ssid <- as.matrix(read.table('set1/ssid/~denovo/populationFst_ssid.txt'))

mcav.tree <- nj(fst.mcav)
mcav.tree.plot <- ggtree(mcav.tree, size = 2, layout = 'daylight')+
#  geom_tiplab(size = 8, hjust = 0.5) +
  ggplot2::xlim(-0.06, 0.11)+
  ggplot2::ylim(-0.09, 0.1)+
  geom_treescale(width = 0.05, x = 0.06, y = -0.02, linesize = 2)
mcav.tree.plot

ssid.tree <- nj(fst.ssid)
ssid.tree.plot <- ggtree(ssid.tree, size = 2, layout = 'daylight')+
#  geom_tiplab(size = 8, hjust = 0.5) +
  ggplot2::xlim(-0.14, 0.13)+
  ggplot2::ylim(-0.2, 0.07)+
  geom_treescale(width = 0.05, x = 0.018, y = -0.2, linesize = 2)
ssid.tree.plot

save_plot(paste0('set1/mcav/popFst_tree_mcav.pdf'), mcav.tree.plot, base_width = 4, base_height = 3)
#save_plot(paste0('set1/ssid/~denovo/popFst_tree_ssid.pdf'), ssid.tree.plot, base_width = 4, base_height = 3)
