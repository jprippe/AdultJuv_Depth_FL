#~~~~~~~
# Analysis to estimate the proportion of the genome exhibiting reduced introgression
# and the strength of this reduction (Figure 2A,B)
#~~~~~~~

library(tidyverse)
library(reshape2)
library(scales)
library(cowplot)

spp <- 'mcav'
dir <- paste0('set1/', spp, '/analysis_sfs/bootstrap/')
files <- list.files(dir)

param.ls <- list()
for(direc in grep('p[1-4]', files, value = T)){
  temp.dir <- paste0(dir,direc)
  load(paste0(temp.dir, '/', direc, '.winboots.res_bootres.RData'))
  temp.list <- list(maxlike = maxlike, mres = mres)
  param.ls[[direc]] <- temp.list
}

prop.gen <- lapply(param.ls, function(x) {
  i <- which(rownames(x$mres) == 'P')
  if(length(i) > 0) {x$mres[i, ]} else {NA} 
})
pgen.df <- do.call(rbind, prop.gen)

pgen.pairs <- data.frame(Var1 = paste0('Pop', substr(rownames(pgen.df), 2, 2)),
                         Var2 = paste0('Pop', substr(rownames(pgen.df), 3, 3)),
                         value = pgen.df$medians, q25 = pgen.df$q25, q75 = pgen.df$q75)

times.migs <- lapply(param.ls, function(x) {
  i <- c(grep("T", rownames(x$mres)), grep("^m", rownames(x$mres)))
  out <- x$mres[i, ]
  # Convert mig parameters from log10 back to raw estimates
  out[grep("^m", rownames(out)),] <- 10^(out[grep("^m", rownames(out)),])
  if(length(i) > 0) {out} else {NA} 
})

ratio.df <- data.frame(mod = names(times.migs), m12 = NA, m21 = NA, whichIntrogress = NA)

# Define a function to determine which portion of the genome (island or non-island as specified in the model) is experiencing higher migration rates
whichIsIntrogress <- function(mod){
  temp.df <- times.migs[[mod]]
  migs.df <- temp.df[grep("^m", rownames(temp.df)),]
  migs.df$isl <- ifelse(endsWith(rownames(migs.df), 'i'), 1, 0)
  isl.is.introgress <- migs.df$medians[migs.df$isl == 1] > migs.df$medians[migs.df$isl == 0]
  if (sum(isl.is.introgress) > length(isl.is.introgress)/2) {
    which.is.introgress <- "isl"
  } else if (sum(isl.is.introgress) < length(isl.is.introgress)/2) {
    which.is.introgress <- "non-isl"
  } else if (sum(isl.is.introgress) == length(isl.is.introgress)/2) {
    which.is.introgress <- "split"
  }
  return(which.is.introgress)
}

# For each migration parameter (paired island and non-island), specify the corresponding epoch
# MCAV: p12: 3, 2, 2; p13: 1, 2, 2; p14: 1, 2, 2; p23: 2, 3, 3; p24: 2, 1, 1; p34: 2, 1, 1
# SSID: p12: 1 2 3, 1 2 3; p13: 1, 2, 2; p14: 2, 3, 3; p23: 2 3, 2 3; p24: 2, 2; p34: 3, 3
migs.ls <- list()
for(mod in names(times.migs)){
  temp.df <- times.migs[[mod]]
  temp.df2 <- param.ls[[mod]][['mres']]
  
  # Calculate proportional duration of each epoch relative to others
  times <- temp.df[grep("^T", rownames(temp.df)),]
  times$scaleValue <- times$medians / sum(times$medians)
  
  # Subset island and non-island migration parameters
  migs <- temp.df[grep("^m", rownames(temp.df)),]
  isl.migs <- rownames(migs)[endsWith(rownames(migs), 'i')]
  non.isl.migs <- rownames(migs)[!endsWith(rownames(migs), 'i')]
  
  # Subset population size parameters
  popsize.1 <- temp.df2[grep("^nu1", rownames(temp.df2)),]
  popsize.2 <- temp.df2[grep("^nu2", rownames(temp.df2)),]
  
  # Determine the proper order of time and population size parameters
  all.pops.underscore <- all(grepl('_', rownames(popsize.1)))
  all.times.index <- all(!is.na(as.numeric(gsub('T', '', rownames(times)))))
  if(all.pops.underscore & all.times.index){
    time.param.order <- order(gsub('T', '', rownames(times)))
    pop.param.order <- order(as.numeric(gsub('.+_', '', rownames(popsize.1))))
  } else if(all.pops.underscore & !all.times.index){
    time.param.order <- rev(order(rownames(times)))
    pop.param.order <- order(as.numeric(gsub('.+_', '', rownames(popsize.1))))
  } else if(!all.pops.underscore & !all.times.index){
    time.param.order <- rev(order(rownames(times)))
    pop.param.order <- rev(order(rownames(popsize.1)))
  }
  
  # Identify the time parameter and population size estimates associated with the current epoch
  n.epochs <- nrow(times)
  migs$current.epoch <- rownames(times)[time.param.order][n.epochs]
  migs$nu1_current <- popsize.1$medians[pop.param.order][n.epochs]
  migs$nu2_current <- popsize.2$medians[pop.param.order][n.epochs]
  
  migs$theta <- param.ls[[mod]][['mres']]['theta','medians']
  migs$timeParams <- NA
  migs$scaleValue <- NA
  migs$Ratio <- NA
  migs$popParams <- NA
  
  for (mig.par in non.isl.migs){
    which.epoch <- readline(prompt = paste0("Which epoch(s) is ",  mig.par, " associated with (contrast: ", mod, ")? Separate with spaces. \nFull list of mig parameters: ", paste(non.isl.migs, sep = "", collapse = ", ")))
    
    # Identify time and population size parameters associated with asymmetric migration parameters
    # Names of asymmetric migration parameters specify the sink pop followed by source pop (i.e., m12 = migration from pop2 to pop1)
    mig.time <- rownames(times)[time.param.order[as.numeric(unlist(strsplit(which.epoch, " ")))]]
    if(grepl('m12', mig.par)){
      mig.popsize <- rownames(popsize.2)[pop.param.order[as.numeric(unlist(strsplit(which.epoch, " ")))]]
    } else if(grepl('m21', mig.par)){
      mig.popsize <- rownames(popsize.1)[pop.param.order[as.numeric(unlist(strsplit(which.epoch, " ")))]]
    } else {
      mig.popsize <- c(rownames(popsize.1)[pop.param.order[as.numeric(unlist(strsplit(which.epoch, " ")))]],
                       rownames(popsize.2)[pop.param.order[as.numeric(unlist(strsplit(which.epoch, " ")))]])
    }
    
    # Record time & population size parameters relevant to the migration parameter
    migs$timeParams[rownames(migs) == mig.par | rownames(migs) == paste0(mig.par,'i')] <- paste(mig.time, collapse = ',')
    migs$popParams[rownames(migs) == mig.par | rownames(migs) == paste0(mig.par,'i')] <- paste(mig.popsize, collapse = ',')

    # Scale the migration rate by the proportional duration of the epoch(s) relative to the full demographic history
    migs$scaleValue[rownames(migs) == mig.par] <- sum(migs$medians[rownames(migs) == mig.par] * times$scaleValue[rownames(times) %in% mig.time])
    migs$scaleValue[rownames(migs) == paste0(mig.par,'i')] <- sum(migs$medians[rownames(migs) == paste0(mig.par,'i')] * times$scaleValue[rownames(times) %in% mig.time])
    
    # Calculate the non-island:island ratio of each pair of migration parameters
    migs$Ratio[rownames(migs) == mig.par] <- migs$scaleValue[rownames(migs) == mig.par]/migs$scaleValue[rownames(migs) == paste0(mig.par,'i')]
    migs$Ratio[rownames(migs) == paste0(mig.par,'i')] <- migs$scaleValue[rownames(migs) == paste0(mig.par,'i')]/migs$scaleValue[rownames(migs) == mig.par]
    
    # Record the population parameter estimates if they are associated with the current epoch
    
  }
  migs.ls[[mod]] <- migs
  
  # Subset asymmetrical migration rates in each direction (symmetrical migration included in both directions)
  # Sum of scaled non-island migration rates divided by sum of scaled island migration rates
  m12 <- migs[c(grep("^m12", rownames(migs)), grep("_", rownames(migs), invert = T)),]
  m21 <- migs[c(grep("^m21", rownames(migs)), grep("_", rownames(migs), invert = T)),]
  if (whichIsIntrogress(mod) == "non-isl"){
    ratio.df$m12[ratio.df$mod == mod] <- sum(m12$scaleValue[!endsWith(rownames(m12),"i")])/sum(m12$scaleValue[endsWith(rownames(m12),"i")])
    ratio.df$m21[ratio.df$mod == mod] <- sum(m21$scaleValue[!endsWith(rownames(m21),"i")])/sum(m21$scaleValue[endsWith(rownames(m21),"i")])
  } else if (whichIsIntrogress(mod) == "isl"){
    ratio.df$m12[ratio.df$mod == mod] <- sum(m12$scaleValue[endsWith(rownames(m12),"i")])/sum(m12$scaleValue[!endsWith(rownames(m12),"i")])
    ratio.df$m21[ratio.df$mod == mod] <- sum(m21$scaleValue[endsWith(rownames(m21),"i")])/sum(m21$scaleValue[!endsWith(rownames(m21),"i")])
  } else {
    ratio.df$m12[ratio.df$mod == mod] <- ratio.df$m21[ratio.df$mod == mod] <- "split"
  }
  ratio.df$whichIntrogress[ratio.df$mod == mod] <- whichIsIntrogress(mod)
}

# Convert ratio.df to pairwise table
pair.tri <- matrix(nrow = 4, ncol = 4)
rownames(pair.tri) <- paste0('Pop', 1:4); colnames(pair.tri) <- paste0('Pop', 1:4)
for(i in 1:nrow(ratio.df)){
  pair.tri[as.numeric(substr(ratio.df[i,'mod'],3,3)), as.numeric(substr(ratio.df[i,'mod'],2,2))] <- ratio.df$m12[i]
  pair.tri[as.numeric(substr(ratio.df[i,'mod'],2,2)), as.numeric(substr(ratio.df[i,'mod'],3,3))] <- ratio.df$m21[i]
}

# Melt pairwise ratio.df table to long-form and add column from popgen.df (proportion of genome experiencing reduced migration)
pair.melt <- melt(pair.tri, na.rm = TRUE)
pair.melt$value <- round(pair.melt$value, digits = 3)
pair.melt$prop <- NA
for(i in 1:nrow(pgen.df)){
  pops <- unlist(strsplit(substr(rownames(pgen.df)[i],2,3), ""))
  pops.index <- c(which(with(pair.melt, grepl(pops[1], Var1) & grepl(pops[2], Var2))), 
                  which(with(pair.melt, grepl(pops[2], Var1) & grepl(pops[1], Var2))))
  pair.melt$prop[pops.index] <- pgen.df$medians[i]
}

# Determine if P represents the region of restricted or unrestricted migration
# In the moments models, P is associated with "i" parameters
# If non-"i" migration estimates are higher than "i" estimates (labeled "non-isl" in ratio.df)...
# ...then P reflects the proportion of the genome experiencing reduced migration
# If "i" migration estimates are higher than non-"i" estimates (labeled "isl" in ratio.df)...
# ...then P reflects the unrestricted proportion of the genome
pgen.pairs2 <- pgen.pairs
pgen.pairs2$value <- 1-pgen.pairs$value
pgen.pairs$isl <- ifelse(ratio.df$whichIntrogress == 'non-isl', 'Island', ifelse(ratio.df$whichIntrogress == 'isl', 'Non-Island', 'Split'))
pgen.pairs2$isl <- ifelse(pgen.pairs$isl == 'Non-Island', 'Island', 'Non-Island')
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
save_plot(paste0('set1/', spp, "/analysis_sfs/bootstrap/pairwiseIOD_pies_", spp, "_modified.pdf"), multi.pie, base_width = 4, base_height = 3.3)

pgen.plot <- ggplot(data = pgen.pairs, aes(Var1, Var2)) +
  geom_tile(color = "black", fill = "white") +
  coord_fixed(ratio = 0.35)+
  #  geom_text(aes(Var1, Var2, 
  #                label = paste0(format(round(value, 2), nsmall = 1), "\n(", 
  #                               format(round(value-uncert, 2), nsmall = 2), " - ",
  #                               format(round(value+uncert, 2), nsmall = 2), ")")),
  #            color = "black", size = 3.5) +
  #geom_text(aes(label = paste0(format(round(value, 2), nsmall = 2), " \u00B1 ", 
  #                             format(round(uncert, 2), nsmall = 2))), color = "black", size = 3.5) +
  geom_text(aes(label = paste0(format(round(value, 2), nsmall = 2), " (", format(round(q25, 2), nsmall = 2), ", ", format(round(q75, 2), nsmall = 2), ")")), color = "black", size = 3.5) +
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
  scale_fill_gradient(low = "white", high = "grey40", limit = c(1, 100), space = "Lab", name="Ratio", trans = "log", breaks = c(1,10,100)) +
  coord_fixed(ratio=0.35) +
  geom_text(aes(Var2, Var1, label = format(round(value, 1), nsmall = 1)), color = "black", size = 6) +
  theme(axis.text = element_text(vjust = 0.5, size = 14, hjust = 0.5),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none") +
  scale_y_discrete(limits = rev(levels(pair.melt$Var1)))+
#  scale_size_identity()+
  labs(title = paste0('Genomic islands of differentiation - ', toupper(spp)), subtitle = expression("Fold reduction in migration rates"), x = "\nSOURCE", y = "SINK\n")

mig.plot2 <- mig.plot +
  theme(legend.position = "right",
        legend.direction = "vertical") +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, title.position = "top", title.hjust = 0.5))

pgen.plot
save_plot(paste0('set1/', spp, "/analysis_sfs/bootstrap/pairwiseIOD_prop_", spp, "_modified.pdf"), pgen.plot, base_width = 7, base_height = 4)
mig.plot
mig.plot2
save_plot(paste0('set1/', spp, "/analysis_sfs/bootstrap/pairwiseIOD_mig_", spp, "_nolegend_modified.pdf"), mig.plot, base_width = 7, base_height = 4)
save_plot(paste0('set1/', spp, "/analysis_sfs/bootstrap/pairwiseIOD_mig_", spp, "_legend_modified.pdf"), mig.plot2, base_width = 7, base_height = 4)

multi_plot <- plot_grid(plotlist=list(pgen.plot, mig.plot), nrow = 1)
save_plot(paste0('set1/', spp, "/analysis_sfs/bootstrap/pairwiseIOD_", spp, "_modified.pdf"), multi_plot, base_width = 12, base_height = 6)


# Direction of asymmetrical migration -------------------------------------

allmigs <- do.call(rbind, migs.ls) %>%
  # Retain only estimates that apply to the current epoch
  filter(apply(., 1, function(x) x[which(colnames(.)=='current.epoch')] %in% unlist(strsplit(x[which(colnames(.)=='timeParams')], split = ',')))) %>%
  # Create columns for population contrasts and parameter ids
  mutate(pop.param = rownames(.)) %>%
  separate(pop.param, c('pop', 'param'), sep = '\\.') %>%
  # Retain only cross-depth scenarios
  filter(!(substr(pop, 2, 3) %in% c('12', '34'))) %>%
  # Define direction of migration estimate (m12: from pop2 to pop1, m21: from pop1 to pop2)
  # For "symmetrical" migration estimates, we input a list of both directions for unnest() to split into duplicate rows
  mutate(direction = ifelse(substr(param, 2, 3) == '12', 'toShallow', ifelse(substr(param, 2, 3) == '21', 'toDeep', strsplit('toShallow,toDeep', split = ',')))) %>%
  unnest(direction)%>%
  # Define which estimates are associated with restricted or unrestricted portion of the genome
  left_join(rename(ratio.df, pop = 'mod') %>% select(pop, whichIntrogress), by = 'pop') %>%
  mutate(isl = ifelse(endsWith(param, 'i') & whichIntrogress == 'non-isl', 'Restricted',
                      ifelse(!endsWith(param, 'i') & whichIntrogress == 'isl', 'Restricted', 'Unrestricted'))) %>%
  # Define column for error bar grouping
  mutate(direc.isl = paste0(direction, isl)) %>%
  mutate(scaledMig = ifelse(direction == 'toShallow', medians/(theta*nu2_current), ifelse(direction == 'toDeep', medians/(theta*nu1_current), NA)),
         q25_scaled = ifelse(direction == 'toShallow', q25/(theta*nu2_current), ifelse(direction == 'toDeep', q25/(theta*nu1_current), NA)),
         q75_scaled = ifelse(direction == 'toShallow', q75/(theta*nu2_current), ifelse(direction == 'toDeep', q75/(theta*nu1_current), NA)))

if(spp == 'mcav') {
  allmigs$pair <- ifelse(allmigs$pop == 'p13', 'Near-Deep1', ifelse(allmigs$pop == 'p14', 'Near-Deep2', ifelse(allmigs$pop == 'p23', 'Off-Deep1', 'Off-Deep2')))
} else if(spp == 'ssid'){
  allmigs$pair <- ifelse(allmigs$pop == 'p13', 'Shallow1-Deep2', ifelse(allmigs$pop == 'p14', 'Shallow1-Deep1', ifelse(allmigs$pop == 'p23', 'Shallow2-Deep2', 'Shallow2-Deep1')))
}

allmigs$isl <- factor(allmigs$isl, levels = c('Unrestricted', 'Restricted'))
allmigs$direc.isl <- factor(allmigs$direc.isl, levels = c(allmigs$direc.isl[2], allmigs$direc.isl[4], allmigs$direc.isl[1], allmigs$direc.isl[3]))

allmigs.plot <- ggplot(allmigs)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        axis.text.y = element_text(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 3, r = 0, b = 2, l = 8, 'mm'))+
  ylab('') +
  geom_col(aes(x = pair, y = scaledMig, fill = direction, alpha = isl), position = 'dodge') +
  geom_linerange(aes(x = pair, ymin = q25_scaled, ymax = q75_scaled, group = direc.isl), position = position_dodge(width = 0.9))+
#  scale_y_continuous(limits = c(0, 5e-04), expand = c(0,0), labels = scales::scientific)+
  scale_y_continuous(limits = c(0, 2e-04), expand = c(0,0), labels = scales::scientific, breaks = c(0, 1e-4, 2e-4))+
  scale_fill_manual(values = c('dodgerblue4', 'goldenrod1'), name = 'Direction') +
  scale_alpha_manual(values = c(1, 0.4), name = 'Genomic Region')

ylab1 <- 'Migration Probability'
ylab2 <- expression(paste("(M / (", theta, " * nu))"))

if(spp == 'mcav') {
  allmigs.plot.label <- ggdraw(allmigs.plot) + 
    draw_label(ylab1, x = 0.025, y = 0.60, angle = 90, size = 11) + 
    draw_label(ylab2, x = 0.075, y = 0.60, angle = 90, size = 11)
} else if(spp == 'ssid'){
  allmigs.plot.label <- ggdraw(allmigs.plot) + 
    draw_label(ylab1, x = 0.025, y = 0.62, angle = 90, size = 11) + 
    draw_label(ylab2, x = 0.075, y = 0.62, angle = 90, size = 11)
}

allmigs.plot.label

save_plot(paste0('set1/', spp, "/analysis_sfs/bootstrap/migDirection_", spp, "_scaled2.pdf"), allmigs.plot.label, base_width = 5, base_height = 4)


