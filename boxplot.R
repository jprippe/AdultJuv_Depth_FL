align_bel = read.table('set1/ssid/alignmentRates_belRef')
align_denovo = read.table('set1/ssid/alignmentRates_denovo')
align = rbind(align_bel[order(align_bel$V1),], align_denovo[order(align_denovo$V1),])
align$xpos = c(rep(0.2,128), rep(0.3,128))
align$dataset = c(rep('bel',128), rep('denovo',128))

ggplot()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.ticks.x = element_blank())+
  geom_point(data=align, aes(x=xpos, y=V2), alpha = 0.1, position = position_jitter(width = 0.04))+
  geom_boxplot(data=align, aes(x=xpos, y=V2, group=dataset), width = 0.07, alpha = 0.6, outlier.alpha = 0)+
  coord_cartesian(ylim = c(0,100))
