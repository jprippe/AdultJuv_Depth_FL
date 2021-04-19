#~~~~~~~~~~~~
# Analysis to determine the strength of evidence in support of a secondary contact model vs.
# a non-secondary contact model with respect to each pairwise demographic history
#~~~~~~~~~~~~

library(tidyverse)

spp <- 'mcav'
dir <- paste0('set1/', spp, '/analysis_sfs/bootstrap/')
files <- list.files(dir)

aic.list <- list()
for(direc in grep('p[1-4]', files, value = T)){
  temp.dir <- paste0(dir,direc)
  modselResult <- paste0(temp.dir, '/', direc, '.modsel')
  system(paste("grep RESULT ", modselResult," -A 4 | grep -v Launcher | grep -E \"[0-9]|\\]\" | perl -pe 's/^100.+\\.o\\d+\\S//' | perl -pe 's/\\n//' | perl -pe 's/[\\[\\]]//g' | perl -pe 's/RESULT/\\nRESULT/g' | grep RESULT | perl -pe 's/array//g' | perl -pe 's/Generation.+$//' | perl -pe 's/#.+//' | perl -pe \"s/[\\(\\)\\']//g\"  >", modselResult,".res",sep=""))
  
  all.mods <- readLines(paste0(temp.dir, '/', direc, '.modsel.res'))
  all.mods2 <- t(matrix(unlist(lapply(strsplit(all.mods, " "), function(x) {x[c(2:6)]})), nrow = 5))
  all.mods3 <- data.frame(all.mods2[complete.cases(all.mods2),])
  npl <- mutate(all.mods3, across(everything(), as.character)) %>%
    mutate(X3 = as.numeric(X3), X4 = as.numeric(X4))
  colnames(npl) <- c("model","id","npara","ll","boot")
  
  contrast=sub("_.+","", npl$boot[1])
  npl$boot=factor(npl$boot)
  
  aics=list()
  for (b in 1:length(levels(npl$boot))) {
    bb=levels(npl$boot)[b]
    nplb=subset(npl,boot==bb)
    maxlike=c();nmod=c()
    for (m in unique(nplb$model)) {
      sub=subset(nplb,model==m)
      nmod=c(nmod,nrow(sub))
      maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
    }
    npara=maxlike$npara
    likes=maxlike$ll
    aic=2*npara-2*likes
    aicc=data.frame(cbind(model=as.character(unique(nplb$model))))
    aicc$aic=aic
    aicc$nmod=nmod
    aicc$boot=bb
    aics[[b]]=aicc
  }
  awt=data.frame(do.call(rbind,aics))
  models=unique(awt$model)
  med=c()
  for (m in models) {
    ss=subset(awt,model==m)
    med=c(med,median(ss$aic))
  }
  modmed=data.frame(cbind(mod=as.character(models)))
  modmed$med=med
  modmed=modmed[order(med),]
  modmed$mod=factor(modmed$mod,levels=modmed$mod)
  awt$model=factor(awt$model,levels=modmed$mod)
  
  aic.list[[direc]] <- modmed
}

##### MCAV
aics12 <- aic.list$p12$med
wtsc12.1 <- exp(-(aics12[1]-aics12[2])/2)/(exp(-(aics12[1]-aics12[2])/2)+exp(-(aics12[2]-aics12[1])/2))
woe12 <- aic.list$p12 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc12.2 <- woe12$woe[1]/woe12$woe[2]
wtsc12.3 <- woe12$woe[1]/(woe12$woe[1] + woe12$woe[2])

aics13 <- aic.list$p13$med
wtsc13.1 <- exp(-(aics13[5]-aics13[1])/2)/(exp(-(aics13[5]-aics13[1])/2)+exp(-(aics13[1]-aics13[5])/2))
woe13 <- aic.list$p13 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc13.2 <- woe13$woe[5]/woe13$woe[1]
wtsc13.3 <- woe13$woe[5]/(woe13$woe[5] + woe13$woe[1])

aics14 <- aic.list$p14$med
wtsc14.1 <- exp(-(aics14[3]-aics14[1])/2)/(exp(-(aics14[3]-aics14[1])/2)+exp(-(aics14[1]-aics14[3])/2))
woe14 <- aic.list$p14 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc14.2 <- woe14$woe[3]/woe14$woe[1]
wtsc14.3 <- woe14$woe[3]/(woe14$woe[3] + woe14$woe[1])

aics23 <- aic.list$p23$med
wtsc23.1 <- exp(-(aics23[2]-aics23[1])/2)/(exp(-(aics23[2]-aics23[1])/2)+exp(-(aics23[1]-aics23[2])/2))
woe23 <- aic.list$p23 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc23.2 <- woe23$woe[2]/woe23$woe[1]
wtsc23.3 <- woe23$woe[2]/(woe23$woe[2] + woe23$woe[1])

aics24 <- aic.list$p24$med
wtsc24.1 <- exp(-(aics24[3]-aics24[1])/2)/(exp(-(aics24[3]-aics24[1])/2)+exp(-(aics24[1]-aics24[3])/2))
woe24 <- aic.list$p24 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc24.2 <- woe24$woe[3]/woe24$woe[1]
wtsc24.3 <- woe24$woe[3]/(woe24$woe[3] + woe24$woe[1])

aics34 <- aic.list$p34$med
wtsc34.1 <- exp(-(aics34[1]-aics34[4])/2)/(exp(-(aics34[1]-aics34[4])/2)+exp(-(aics34[4]-aics34[1])/2))
woe34 <- aic.list$p34 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc34.2 <- woe34$woe[1]/woe34$woe[4]
wtsc34.3 <- woe34$woe[1]/(woe34$woe[1] + woe34$woe[4])

##### SSID
aics12 <- aic.list$p12$med
wtsc12.1 <- exp(-(aics12[1]-aics12[2])/2)/(exp(-(aics12[1]-aics12[2])/2)+exp(-(aics12[2]-aics12[1])/2))
woe12 <- aic.list$p12 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc12.2 <- woe12$woe[1]/woe12$woe[2]
wtsc12.3 <- woe12$woe[1]/(woe12$woe[1] + woe12$woe[2])

aics13 <- aic.list$p13$med
wtsc13.1 <- exp(-(aics13[3]-aics13[1])/2)/(exp(-(aics13[3]-aics13[1])/2)+exp(-(aics13[1]-aics13[3])/2))
woe13 <- aic.list$p13 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc13.2 <- woe13$woe[3]/woe13$woe[1]
wtsc13.3 <- woe13$woe[3]/(woe13$woe[3] + woe13$woe[1])

aics14 <- aic.list$p14$med
wtsc14.1 <- exp(-(aics14[1]-aics14[3])/2)/(exp(-(aics14[1]-aics14[3])/2)+exp(-(aics14[3]-aics14[1])/2))
woe14 <- aic.list$p14 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc14.2 <- woe14$woe[1]/woe14$woe[3]
wtsc14.3 <- woe14$woe[1]/(woe14$woe[1] + woe14$woe[3])

aics23 <- aic.list$p23$med
wtsc23.1 <- exp(-(aics23[2]-aics23[1])/2)/(exp(-(aics23[2]-aics23[1])/2)+exp(-(aics23[1]-aics23[2])/2))
woe23 <- aic.list$p23 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc23.2 <- woe23$woe[2]/woe23$woe[1]
wtsc23.3 <- woe23$woe[2]/(woe23$woe[2] + woe23$woe[1])

aics24 <- aic.list$p24$med
wtsc24.1 <- exp(-(aics24[1]-aics24[5])/2)/(exp(-(aics24[1]-aics24[5])/2)+exp(-(aics24[5]-aics24[1])/2))
woe24 <- aic.list$p24 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc24.2 <- woe24$woe[1]/woe24$woe[5]
wtsc24.3 <- woe24$woe[1]/(woe24$woe[1] + woe24$woe[5])

aics34 <- aic.list$p34$med
wtsc34.1 <- exp(-(aics34[1]-aics34[3])/2)/(exp(-(aics34[1]-aics34[3])/2)+exp(-(aics34[3]-aics34[1])/2))
woe34 <- aic.list$p34 %>%
  mutate(deltaAIC = med - min(med, na.rm = T)) %>%
  mutate(exps = exp(-deltaAIC/2)) %>%
  mutate(woe = exps/sum(exps, na.rm = T))
wtsc34.2 <- woe34$woe[1]/woe34$woe[3]
wtsc34.3 <- woe34$woe[1]/(woe34$woe[1] + woe34$woe[3])

##### Concatenation and plotting

wtsc <- data.frame(contrast = c('p12','p13','p14','p23','p24','p34'), 
                   wtsc1 = c(wtsc12.1, wtsc13.1, wtsc14.1, wtsc23.1, wtsc24.1, wtsc34.1),
                   wtsc2 = c(wtsc12.2, wtsc13.2, wtsc14.2, wtsc23.2, wtsc24.2, wtsc34.2),
                   wtsc3 = c(wtsc12.3, wtsc13.3, wtsc14.3, wtsc23.3, wtsc24.3, wtsc34.3))
write.csv(wtsc, file = paste0(dir, 'weight_sc_', spp,'.csv'), row.names = F, quote = F)


scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
}

wtsc.plot <- ggplot(wtsc)+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y = "Evidence Ratio")+
  geom_hline(yintercept = 1)+
  geom_col(aes(x=contrast, y=wtsc2))+
  scale_y_log10(label = scientific_10, breaks = c(1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6, 1e8), limits = c(1e-9, 1e9), expand = c(0,0))
wtsc.plot

cowplot::save_plot(filename = paste0(dir, 'wtsc_plot_', spp, '.pdf'), plot = wtsc.plot, base_height = 3, base_width = 4)
