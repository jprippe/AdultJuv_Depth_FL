library(tidyverse)
library(reshape2)
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

ratio.df <- data.frame(mod = names(times.migs), m12 = NA, m21 = NA, whichIntrogress = NA)

# Define a function to determine which portion of the genome (island or non-island as specified in the model) is experiencing higher migration rates
whichIsIntrogress <- function(mod){
  temp.df <- times.migs[[mod]]
  migs.df <- temp.df[grep("^m", temp.df$Parameter),]
  migs.df$isl <- ifelse(endsWith(migs.df$Parameter, 'i'), 1, 0)
  isl.is.introgress <- migs.df$Value[migs.df$isl == 1] > migs.df$Value[migs.df$isl == 0]
  if (sum(isl.is.introgress) > length(isl.is.introgress)/2) {
    which.is.introgress <- "isl"
  } else if (sum(isl.is.introgress) < length(isl.is.introgress)/2) {
    which.is.introgress <- "non-isl"
  } else if (sum(isl.is.introgress) == length(isl.is.introgress)/2) {
    which.is.introgress <- "split"
  }
  return(which.is.introgress)
}

for(mod in names(times.migs)){
  temp.df <- times.migs[[mod]]
  
  # Calculate proportional duration of each epoch relative to others
  times <- temp.df[grep("^T", temp.df$Parameter),]
  times$scaleValue <- times$Value / sum(times$Value)
  
  # Subset island and non-island migration parameters
  migs <- temp.df[grep("^m", temp.df$Parameter),]
  isl.migs <- migs$Parameter[endsWith(migs$Parameter, 'i')]
  non.isl.migs <- migs$Parameter[!endsWith(migs$Parameter, 'i')]
  
  migs$scaleValue <- NA
  migs$Ratio <- NA
  
  # For each migration parameter (paired island and non-island), specify the corresponding epoch
  # If there is no suffix, the migration rate corresponds to the penultimate epoch
  # If there is a suffix (_i), the migration rate corresponds to the ith epoch
  for (mig.par in non.isl.migs){
    if (nrow(times)>1){
      if (length(strsplit(mig.par, '_')[[1]]) == 1){
        mig.time <- times$Parameter[nrow(times)-1]
      } else if (length(strsplit(mig.par, '_')[[1]]) == 2){
        mig.time <- times$Parameter[as.numeric(strsplit(mig.par, '_')[[1]][2])]
      }
    } else {mig.time <- times$Parameter}
    
    # Scale the migration rate by the proportional duration of the epoch relative to other epochs
    migs$scaleValue[migs$Parameter == mig.par] <- migs$Value[migs$Parameter == mig.par] * times$scaleValue[times$Parameter == mig.time]
    migs$scaleValue[migs$Parameter == paste0(mig.par,'i')] <- migs$Value[migs$Parameter == paste0(mig.par,'i')] * times$scaleValue[times$Parameter == mig.time]
    
    # Calculate the non-island:island ratio of each pair of migration parameters
    migs$Ratio[migs$Parameter == mig.par] <- migs$scaleValue[migs$Parameter == mig.par]/migs$scaleValue[migs$Parameter == paste0(mig.par,'i')]
    migs$Ratio[migs$Parameter == paste0(mig.par,'i')] <- migs$scaleValue[migs$Parameter == paste0(mig.par,'i')]/migs$scaleValue[migs$Parameter == mig.par]
  }
  
  # Subset asymmetrical migration rates in each direction (symmetrical migration included in both directions)
  # Sum of scaled non-island migration rates divided by sum of scaled island migration rates
  m12 <- migs[c(grep("^m12", migs$Parameter), grep("_", migs$Parameter, invert = T)),]
  m21 <- migs[c(grep("^m21", migs$Parameter), grep("_", migs$Parameter, invert = T)),]
  if (whichIsIntrogress(mod) == "non-isl"){
    ratio.df$m12[ratio.df$mod == mod] <- sum(m12$scaleValue[!endsWith(m12$Parameter,"i")])/sum(m12$scaleValue[endsWith(m12$Parameter,"i")])
    ratio.df$m21[ratio.df$mod == mod] <- sum(m21$scaleValue[!endsWith(m21$Parameter,"i")])/sum(m21$scaleValue[endsWith(m21$Parameter,"i")])
  } else if (whichIsIntrogress(mod) == "isl"){
    ratio.df$m12[ratio.df$mod == mod] <- sum(m12$scaleValue[endsWith(m12$Parameter,"i")])/sum(m12$scaleValue[!endsWith(m12$Parameter,"i")])
    ratio.df$m21[ratio.df$mod == mod] <- sum(m21$scaleValue[endsWith(m21$Parameter,"i")])/sum(m21$scaleValue[!endsWith(m21$Parameter,"i")])
  } else {
    ratio.df$m12[ratio.df$mod == mod] <- ratio.df$m21[ratio.df$mod == mod] <- "split"
  }
  ratio.df$whichIntrogress[ratio.df$mod == mod] <- whichIsIntrogress(mod)
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

pgen.pairs2 <- pgen.pairs
pgen.pairs2$value <- 1-pgen.pairs$value
pgen.pairs$isl <- "Island"; pgen.pairs2$isl <- "Non-Island"
pgen.pairs.long <- rbind(pgen.pairs, pgen.pairs2)
pgen.pairs.long$pair <- factor(paste(gsub("Pop","",pgen.pairs.long$Var1), gsub("Pop","",pgen.pairs.long$Var2), sep = "_"))

pie.list <- list()
for(i in levels(pgen.pairs.long$pair)){
  pgen.pie <- ggplot(data = subset(pgen.pairs.long, pair == i), aes(x = "", y = value, fill = isl))+
    geom_bar(width = 1, stat = "identity")+
    coord_polar("y", start=0)+
    scale_fill_grey()+
    ylab(paste0(as.character(sprintf("%.1f", round(subset(pgen.pairs.long, pair == i & isl == "Island")$value, 3)*100)), "%"))+
    theme_minimal()+
    theme(
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(), 
      legend.position = "none"
    )
  pie.list[[i]] <- pgen.pie
}
if (spp == 'mcav'){
  #Rename factor levels to match modified pop ids
  pair.melt$Var1 <- fct_recode(pair.melt$Var1, Nearshore='Pop1', Offshore='Pop2', Deep1='Pop3', Deep2='Pop4')
  pair.melt$Var2 <- fct_recode(pair.melt$Var2, Nearshore='Pop1', Offshore='Pop2', Deep1='Pop3', Deep2='Pop4')
  pgen.pairs$Var1 <- fct_expand(pgen.pairs$Var1, 'Deep2') %>%
    fct_recode(Nearshore='Pop1', Offshore='Pop2', Deep1='Pop3') %>%
    fct_relevel('Nearshore', 'Offshore', 'Deep1', 'Deep2')
  pgen.pairs$Var2 <- fct_expand(pgen.pairs$Var2, 'Nearshore') %>%
    fct_recode(Offshore='Pop2', Deep1='Pop3', Deep2='Pop4') %>%
    fct_relevel('Nearshore', 'Offshore', 'Deep1', 'Deep2')
  
  #Piecharts displaying proportion of genome with reduced migration rates
  multi.pie <- plot_grid(pie.list$`1_2`, NULL, NULL,
                         pie.list$`1_3`, pie.list$`2_3`, NULL,
                         pie.list$`1_4`, pie.list$`2_4`, pie.list$`3_4`,
                         ncol = 3)
} else if (spp == 'ssid'){
  #Rename and reorder factor levels to match modified pop ids
  pair.melt$Var1 <- fct_recode(pair.melt$Var1, Shallow1='Pop1', Shallow2='Pop2', Deep2='Pop3', Deep1='Pop4') %>%
    fct_relevel('Shallow1', 'Shallow2', 'Deep1', 'Deep2')
  pair.melt$Var2 <- fct_recode(pair.melt$Var2, Shallow1='Pop1', Shallow2='Pop2', Deep2='Pop3', Deep1='Pop4') %>%
    fct_relevel('Shallow1', 'Shallow2', 'Deep1', 'Deep2')
  pgen.pairs$Var1 <- fct_expand(pgen.pairs$Var1, 'Deep1') %>%
    fct_recode(Shallow1='Pop1', Shallow2='Pop2', Deep2='Pop3') %>%
    fct_relevel('Shallow1', 'Shallow2', 'Deep1', 'Deep2')
  pgen.pairs$Var2 <- fct_expand(pgen.pairs$Var2, 'Shallow1') %>%
    fct_recode(Shallow2='Pop2', Deep2='Pop3', Deep1='Pop4') %>%
    fct_relevel('Shallow1', 'Shallow2', 'Deep1', 'Deep2')
  #Flip last row indices for plotting
  pgen.pairs$Var1[6] <- factor('Deep1', levels = levels(pgen.pairs$Var1))
  pgen.pairs$Var2[6] <- factor('Deep2', levels = levels(pgen.pairs$Var2))
  
  #Piecharts displaying proportion of genome with reduced migration rates
  #(Plot layout rearranged to match modified pop ids)
  multi.pie <- plot_grid(pie.list$`1_2`, NULL, NULL,
                         pie.list$`1_4`, pie.list$`2_4`, NULL,
                         pie.list$`1_3`, pie.list$`2_3`, pie.list$`3_4`,
                         ncol = 3)
}
multi.pie
save_plot(paste0('set1/', spp, "/~denovo/analysis_sfs/pairwiseIOD_pies_", spp, "_modified.pdf"), multi.pie, base_width = 4, base_height = 3.3)

pgen.plot <- ggplot(data = pgen.pairs, aes(Var1, Var2)) +
  geom_tile(color = "black", fill = "white") +
  coord_fixed(ratio = 0.35)+
  #  geom_text(aes(Var1, Var2, 
  #                label = paste0(format(round(value, 2), nsmall = 1), "\n(", 
  #                               format(round(value-uncert, 2), nsmall = 2), " - ",
  #                               format(round(value+uncert, 2), nsmall = 2), ")")),
  #            color = "black", size = 3.5) +
  geom_text(aes(label = paste0(format(round(value, 2), nsmall = 2), " \u00B1 ", 
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
  scale_fill_gradient2(low = "grey40", high = "grey40", mid = "white", midpoint = 0, limit = c(0,8), space = "Lab", name=expression("log"[2]*"(Ratio)")) +
  coord_fixed(ratio=0.35) +
  geom_text(aes(Var2, Var1, label = format(round(value, 1), nsmall = 1)), color = "black", size = 6) +
  theme(axis.text = element_text(vjust = 0.5, size = 14, hjust = 0.5),
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
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_prop_", spp, "_modified.pdf"), pgen.plot, base_width = 7, base_height = 4)
mig.plot
mig.plot2
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_mig_", spp, "_nolegend_modified.pdf"), mig.plot, base_width = 7, base_height = 4)
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_mig_", spp, "_legend_modified.pdf"), mig.plot2, base_width = 7, base_height = 4)

multi_plot <- plot_grid(plotlist=list(pgen.plot, mig.plot), nrow = 1)
save_plot(paste0('set1/', spp, "/analysis_sfs/pairwiseIOD_", spp, "_modified.pdf"), multi_plot, base_width = 12, base_height = 6)

