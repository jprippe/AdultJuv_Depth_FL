library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

spp <- 'mcav'
dir <- paste0('set1/', spp, '/analysis_sfs/moments/')
files <- list.files(dir)

param.ls <- list()
if(spp == 'mcav'){
  for(direc in files[grep('_uf', files)]){
    temp.dir <- paste0(dir,direc)
    if(length(list.files(temp.dir, pattern = "params")) > 1){
      mod.choice <- readline(prompt = paste0("Specify the model ID (", direc, "):"))
      param.ls[[direc]] <- read.table(paste0(temp.dir, '/', mod.choice, "_params.txt"), header = T, stringsAsFactors = F)
    } else if(length(grep("params", list.files(paste0(dir,direc)))) == 1){
      param.file <- list.files(temp.dir, pattern = "params")
      param.ls[[direc]] <- read.table(paste0(temp.dir, '/', param.file), header = T, stringsAsFactors = F)
    }
  }
} else if(spp == 'ssid'){
  for(direc in files[grep('c[1-4]', files)]){
    temp.dir <- paste0(dir,direc)
    if(length(list.files(temp.dir, pattern = "params")) > 1){
      mod.choice <- readline(prompt = paste0("Specify the model ID (", direc, "):"))
      param.ls[[direc]] <- read.table(paste0(temp.dir, '/', mod.choice, "_params.txt"), header = T, stringsAsFactors = F)
    } else if(length(grep("params", list.files(paste0(dir,direc)))) == 1){
      param.file <- list.files(temp.dir, pattern = "params")
      param.ls[[direc]] <- read.table(paste0(temp.dir, '/', param.file), header = T, stringsAsFactors = F)
    }
  }
}

prop.gen <- lapply(param.ls, function(x) {
  i <- which(x$Parameter == 'P')
  if(length(i) > 0) {x[i, ]} else {NA} 
})
pgen.df <- do.call(rbind, prop.gen)

pgen.pairs <- data.frame(Var1 = paste0('Pop', substr(rownames(pgen.df), 2, 2)),
                         Var2 = paste0('Pop', substr(rownames(pgen.df), 3, 3)),
                         value = pgen.df$Value, uncert = pgen.df$Uncertainty)

times.migs <- lapply(param.ls, function(x) {
  i <- c(grep("T", x$Parameter), grep("^m", x$Parameter))
  if(length(i) > 0) {x[i, ]} else {NA} 
})

ratio.df <- data.frame(mod = names(times.migs), m12 = NA, m21 = NA)
for(mod in names(times.migs)){
  temp.df <- times.migs[[mod]]
  times <- temp.df[grep("^T", temp.df$Parameter),]
  times$scaleValue <- times$Value / sum(times$Value)
  migs <- temp.df[grep("^m", temp.df$Parameter),]
  isl.migs <- migs$Parameter[endsWith(migs$Parameter, 'i')]
  non.isl.migs <- migs$Parameter[!endsWith(migs$Parameter, 'i')]
  
  migs$scaleValue <- NA
  migs$Ratio <- NA
  for(i in seq(length(non.isl.migs))){
    mig.par <- non.isl.migs[i]
    if(length(strsplit(non.isl.migs, '_')[[i]]) == 1){
      mig.time <- times$Parameter[nrow(times)-1]
    } else if (length(strsplit(non.isl.migs, '_')[[i]]) == 2){
      mig.time <- times$Parameter[as.numeric(strsplit(non.isl.migs, '_')[[i]][2])]
    }
    migs$scaleValue[migs$Parameter == mig.par] <- migs$Value[migs$Parameter == mig.par] * times$scaleValue[times$Parameter == mig.time]
    migs$scaleValue[migs$Parameter == paste0(mig.par,'i')] <- migs$Value[migs$Parameter == paste0(mig.par,'i')] * times$scaleValue[times$Parameter == mig.time]
    
    migs$Ratio[migs$Parameter == mig.par] <- migs$scaleValue[migs$Parameter == mig.par]/migs$scaleValue[migs$Parameter == paste0(mig.par,'i')]
  }
  m12 <- migs[c(grep("^m12", migs$Parameter), grep("_", migs$Parameter, invert = T)),]
  ratio.df$m12[ratio.df$mod == mod] <- sum(m12$scaleValue[!endsWith(m12$Parameter,"i")])/sum(m12$scaleValue[endsWith(m12$Parameter,"i")])
  m21 <- migs[c(grep("^m21", migs$Parameter), grep("_", migs$Parameter, invert = T)),]
  ratio.df$m21[ratio.df$mod == mod] <- sum(m21$scaleValue[!endsWith(m21$Parameter,"i")])/sum(m21$scaleValue[endsWith(m21$Parameter,"i")])
}

pair.tri <- matrix(nrow = 4, ncol = 4)
rownames(pair.tri) <- paste0('Pop', 1:4); colnames(pair.tri) <- paste0('Pop', 1:4)
for(i in 1:nrow(ratio.df)){
  pair.tri[as.numeric(substr(ratio.df[i,'mod'],3,3)), as.numeric(substr(ratio.df[i,'mod'],2,2))] <- log(ratio.df$m12[i], 2)
  pair.tri[as.numeric(substr(ratio.df[i,'mod'],2,2)), as.numeric(substr(ratio.df[i,'mod'],3,3))] <- log(ratio.df$m21[i], 2)
}

pair.melt <- melt(pair.tri, na.rm = TRUE)
pair.melt$value <- round(pair.melt$value, digits = 3)
pair.melt$prop <- NA
for(i in 1:nrow(pgen.df)){
  pops <- unlist(strsplit(substr(rownames(pgen.df)[i],2,3), ""))
  pops.index <- c(which(with(pair.melt, grepl(pops[1], Var1) & grepl(pops[2], Var2))), 
                  which(with(pair.melt, grepl(pops[2], Var1) & grepl(pops[1], Var2))))
  pair.melt$prop[pops.index] <- pgen.df$Value[i]
}
pair.melt$prop.rescale <- 6-rescale(pair.melt$prop, c(0,5))
prop.rescale.key <- c(0.1, 0.5, 1)

pgen.plot <- ggplot(data = pgen.pairs, aes(Var1, Var2)) +
  geom_tile(color = "black", fill = "white") +
  coord_fixed(ratio = 0.35)+
  #  geom_text(aes(Var1, Var2, 
  #                label = paste0(format(round(value, 2), nsmall = 1), "\n(", 
  #                               format(round(value-uncert, 2), nsmall = 2), " - ",
  #                               format(round(value+uncert, 2), nsmall = 2), ")")),
  #            color = "black", size = 3.5) +
  geom_text(aes(Var1, Var2, label = paste0(format(round(value, 2), nsmall = 2), " \u00B1 ", 
                                           format(round(uncert, 2), nsmall = 2))), color = "black", size = 3.5) +
  theme(axis.text = element_text(vjust = 0.5, size = 10, hjust = 0.5),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold")) +
  scale_y_discrete(limits = rev(levels(pair.melt$Var1)))+
  scale_x_discrete(limits = levels(pair.melt$Var1))+
  #  scale_size_identity()+
  labs(title = paste0('Genomic islands of differentiation - ', toupper(spp)), subtitle = "Proportion of genome with reduced migration", x = "\nSOURCE", y = "SINK\n")

mig.plot <- ggplot(data = pair.melt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "grey40", high = "grey40", mid = "white", midpoint = 0, limit = c(0,12), space = "Lab", name=expression("log"[2]*"(Ratio)")) +
  coord_fixed(ratio=0.35) +
  geom_text(aes(Var2, Var1, label = format(round(value, 1), nsmall = 1)), color = "black", size = 3.5) +
  theme(axis.text = element_text(vjust = 0.5, size = 10, hjust = 0.5),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  scale_y_discrete(limits = rev(levels(pair.melt$Var1)))+
#  scale_size_identity()+
  labs(title = paste0('Genomic islands of differentiation - ', toupper(spp)), subtitle = expression("Relative reduction in migration rates (log"[2]*")"), x = "\nSOURCE", y = "SINK\n")

mig.plot2 <- mig.plot +
  theme(legend.position = "right",
        legend.direction = "vertical") +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, title.position = "top", title.hjust = 0.5))

pgen.plot
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_prop_", spp, ".pdf"), pgen.plot, base_width = 7, base_height = 4)
mig.plot
mig.plot2
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_mig_", spp, "_nolegend.pdf"), mig.plot, base_width = 7, base_height = 4)
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_mig_", spp, "_legend.pdf"), mig.plot2, base_width = 7, base_height = 4)

multi_plot <- plot_grid(plotlist=list(pgen.plot, mig.plot), nrow = 1)
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_", spp, ".pdf"), multi_plot, base_width = 12, base_height = 6)

