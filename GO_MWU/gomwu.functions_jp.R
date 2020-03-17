clusteringGOs=function(gen2go,div,cutHeight) {
	inname=paste("dissim0_",div,"_",gen2go,sep="")
	outname=paste("cl_",inname,sep="")
	if (!file.exists(outname)) {
		diss=read.table(inname,sep="\t",header=T,check.names=F)
		row.names(diss)=names(diss)
		hc=hclust(as.dist(diss),method="complete")
		cc=cutree(hc,h=cutHeight)
		write.csv(cc,file=outname,quote=F)
	}
}



#---------------
gomwuStats=function(input,goDatabase,goAnnotations, goDivision, Module=FALSE, Alternative="t", adjust.multcomp="BH", clusterCutHeight=0.25,largest=0.1,smallest=5,perlPath="perl", shuffle.reps=10){

	extraOptions=paste("largest=",largest," smallest=",smallest," cutHeight=",clusterCutHeight,sep="")
	if (Module==TRUE) { adjust.multcomp="shuffle" }
	system(paste(perlPath,"./gomwu_a.pl",goDatabase,goAnnotations,input,goDivision,extraOptions))
	clusteringGOs(goAnnotations,goDivision,clusterCutHeight)
	system(paste(perlPath,"./gomwu_b.pl",goAnnotations,input,goDivision))

	inname=paste(goDivision,"_",input,sep="")	
	rsq=read.table(inname,sep="\t",header=T)
	rsq$term=as.factor(rsq$term)
	# 'value' here is the chosen variable of interest for each gene (-log10 p-value of module membership or log-fold change of gene expression)

	mwut.t=TRUE
	if (length(levels(as.factor(rsq$value)))==2) {
		cat("Binary classification detected; will perform Fisher's test\n");
		mwut.t=F
		rr=fisherTest(rsq)
	} else {
		if (Module==TRUE) {
			rsq.f=rsq
			rsq.f$value=as.numeric(rsq.f$value>0)
			rr=fisherTest(rsq.f)
			rsq.m=rsq[rsq$value>0,]
			rsq.m$co=1
			st=aggregate(rsq.m[,"co"],list(rsq.m$term),sum)
			goodgos=st$Group.1[st$x>=smallest]
			rsq.m=rsq.m[rsq.m$term %in% goodgos,]
			rsq.m$term=factor(rsq.m$term,levels=unique(rsq.m$term))
			rrm=mwuTest(rsq.m,"g")
			rr0=rr[rr$term %in% rrm$term,]
			rr1=rr[!(rr$term %in% rrm$term),]
			rr0=rr0[order(rr0$term),]
			rrm=rrm[order(rrm$term),]
			rr0$pval=rr0$pval*rrm$pval
			rr=rbind(rr0,rr1)
		} else {
			cat("Continuous measure of interest: will perform MWU test\n");		
			rr=mwuTest(rsq,Alternative)
		  }
	}
	
	if (adjust.multcomp=="shuffle"){
	 cat("shuffling values to calculate FDR, ",shuffle.reps," reps\n")
	 reps=shuffle.reps
	 spv=c()
	 for (i in 1:reps) {
		 print(paste("replicate",i))
		 rsqq=rsq
		 rsqq$value=sample(rsq$value)
		 if (mwut.t==TRUE) { rs=mwuTest(rsqq,Alternative) } else { rs=fisherTest(rsqq) }
		 spv=append(spv,rs$pval)
	 }
	 fdr=c()
	 for (p in rr$pval){
		 fdr=append(fdr,(sum(spv<=p)/reps)/sum(rr$pval<=p))
	 }
	 fdr[fdr>1]=1
	} else {
	  fdr=p.adjust(rr$pval,method=adjust.multcomp)
	}
	 
	cat(paste(sum(fdr<0.1)," GO terms at 10% FDR\n"))
	rr$p.adj=fdr
	fname=paste("MWU_",inname,sep="")
	write.table(rr,fname,row.names=F)
}

#---------------------
mwuTest=function(gotable,Alternative) {
	gos=gotable #rsq dataframe from gomwuStats
	terms=levels(gos$term)
	gos$seq=as.character(gos$seq)
	nrg=gos[!duplicated(gos$seq),5] #stores -log10 p-value of module membership for only the first instance of each gene, removes all duplicates that follow
	names(nrg)=gos[!duplicated(gos$seq),4]
#	nrg=nrg+rnorm(nrg,sd=0.01) # to be able to do exact wilcox test
	rnk=rank(nrg) #ranks genes by -log10 p-value of module membership (min(rnk)*2-1 should equal the number of zeros in nrg)
	names(rnk)=names(nrg)
	pvals=c();drs=c();nams=c();levs=c();nseqs=c()
	for (t in terms){ #for each GO term...
		got=gos[gos$term==t,] #subset 'gos' dataframe to rows with matching GO term
		got=got[!duplicated(got$seq),] #remove gene duplicates
		ngot=gos[gos$term!=t,] #subset 'gos' dataframe to rows that do NOT match GO term
		ngot=ngot[!duplicated(ngot$seq),] #remove gene duplicates
		ngot=ngot[!(ngot$seq %in% got$seq),] #ensure that 'ngot' dataframe does not include any genes from 'got' dataframe
		sgo.yes=got$seq #vector of genes associated with GO term
		n1=length(sgo.yes) #number of genes associated with GO term
		sgo.no=ngot$seq #vector of genes NOT associated with GO term
		n2=length(sgo.no) #number of genes NOT associated with GO term
		
		# Wilcoxon test determines whether the -log10 p-value of module membership for GO genes is significantly different than that of non-GO genes
		wi=wilcox.test(nrg[sgo.yes],nrg[sgo.no],alternative=Alternative)	# removed correct=FALSE -Misha
		
		r1=sum(rnk[sgo.yes])/n1 #average rank of GO genes with respect to module membership
		r0=sum(rnk[sgo.no])/n2 #average rank of non-GO genes with respect to module membership
		dr=r1-r0 #delta rank: difference between rank of GO genes and non-GO genes with respect to module membership
		drs=append(drs,round(dr,0))
		# positive delta rank indicates that GO genes are ranked higher with respect to module membership on average than non-GO genes
		levs=append(levs,got$lev[1])
		nams=append(nams,as.character(got$name[1]))
		pvals=append(pvals,wi$p.value)
		nseqs=append(nseqs,n1)	
	}
	res=data.frame(cbind("delta.rank"=drs,"pval"=pvals,"level"=levs,nseqs))
	res=cbind(res,"term"=as.character(terms),"name"=nams)
	res$pval=as.numeric(as.character(res$pval))
	res$delta.rank=as.numeric(as.character(res$delta.rank))
	res$level=as.numeric(as.character(res$level))
	res$nseqs=as.numeric(as.character(res$nseqs))
	return(res)
}
#------------------------
fisherTest=function(gotable) {
	gos=gotable
	terms=levels(gos$term)
	gos$seq=as.character(gos$seq)
	pft=c();nam=c();lev=c();nseqs=c()
	for (t in terms) {
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		go.sig=sum(got$value)
		go.ns=length(got[,1])-go.sig
		ngo.sig=sum(ngot$value)
		ngo.ns=length(ngot[,1])-ngo.sig
		sig=c(go.sig,ngo.sig) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(go.ns,ngo.ns) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
		ff=fisher.test(mm,alternative="greater")
		pft=append(pft,ff$p.value)
		nam=append(nam,as.character(got$name[1]))
		lev=append(lev,got$lev[1])
		nseqs=append(nseqs,length(got[,1]))
	}
	res=data.frame(cbind("delta.rank"=rep(0),"pval"=pft,"level"=lev,nseqs,"term"=terms,"name"=nam))
	res[,1]=as.numeric(as.character(res[,1]))
	res[,2]=as.numeric(as.character(res[,2]))
	res[,3]=as.numeric(as.character(res[,3]))
	res$nseqs=as.numeric(as.character(res$nseqs))
	return(res)
}

#-------------------------
gomwuPlot=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL) {
	require(ape)
  require(ggplot2)
  require(ggdendro)
	
	input=inFile
	in.mwu=paste("MWU",goDivision,input,sep="_")
	in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")
	
	cutoff=-log(level1,10)
	pv=read.table(in.mwu,header=T)
	row.names(pv)=pv$term
	in.raw=paste(goDivision,input,sep="_")
	rsq=read.table(in.raw,sep="\t",header=T)
	rsq$term=as.factor(rsq$term)
	
	if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
	heat=data.frame(cbind("pval"=pvals)) 
	row.names(heat)=pv$term
	heat$pval=-log(heat$pval+1e-15,10)
	heat$direction=0
	heat$direction[pv$delta.rank>0]=1
	if (cutoff>0) { 
		goods=subset(heat,pval>=cutoff) 
	} else {
		goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
		goods=heat[row.names(heat) %in% goods.names,]
	}
	
	if (is.null(colors) | length(colors)<4 ) {
		colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
		if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
			colors=c("black","black","grey50","grey50")
		}
	}
	goods.names=row.names(goods)
	
	# reading and subsetting dissimilarity matrix
	diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
	row.names(diss)=names(diss)
	diss.goods=diss[goods.names,goods.names]
	
	# how many genes out of what we started with we account for with our best categories?
	good.len=c();good.genes=c()
	for (g in goods.names) {
		sel=rsq[rsq$term==g,]	
		pass=abs(sel$value)>=absValue
		sel=sel[pass,]
		good.genes=append(good.genes,as.character(sel$seq))
		good.len=append(good.len,nrow(sel))
	}
	ngenes=length(unique(good.genes))
	
	#hist(rsq$value)
	totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
	row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
	row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
	row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
	
	# clustering terms better than cutoff
	GO.categories=as.dist(diss.goods)
	cl.goods=hclust(GO.categories,method="average")
	labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
	goods=goods[labs,]
	rownames(goods)=sub(" activity","",rownames(goods))
	ddata <- dendro_data(as.dendrogram(cl.goods), type = "rectangle")
	goods$x <- label(ddata)$x
	goods$y <- label(ddata)$y
	goods$color=1
	goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
	goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
	goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
	goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
	goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
	goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
	goods$fontface='italic'; goods$fontsize=3.5*txtsize
	goods$fontface[goods$pval > -log(level2,10)]='plain'
	goods$fontface[goods$pval > -log(level3,10)]='bold'
	goods$fontsize[goods$pval > -log(level3,10)]=4.5*txtsize
	
	gomwu_plot <- ggplot() + 
	  geom_segment(data = segment(ddata), aes(x = x, y = y, xend = xend, yend = yend)) +
	  geom_text(data = goods, aes(x = x, y = y, label = rownames(goods), size = fontsize, fontface = fontface, color = color), 
	            hjust = 'left', nudge_y = 0.3, family = font.family) +
	  theme_dendro() +
	  coord_flip(ylim = c(-20,5)) + 
	  scale_y_reverse() +
	  scale_color_identity() +
	  scale_size_identity() +
	  annotate("text", x = nrow(goods)-c(0,1,2), y = c(5,5,5), 
	           label = c(paste0("p < ", level3), paste0("p < ", level2), paste0("p < ", as.character(10^(-cutoff)))), 
	           parse = T, size = c(4.5*txtsize, 3.5*txtsize, 3.5*txtsize), color = c('black', 'black', 'grey50'),
	           hjust = 'left', family = font.family, fontface = c("bold", "plain", "italic"))
	print(gomwu_plot)
	
	cat(paste("GO terms dispayed: ",length(goods.names)),"\n")
	cat(paste("\"Good genes\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))
	goods$pval=10^(-1*goods$pval)
	return(goods[,1:2])
}

#------------------
# returns non-overlapping GO categories based on dissimilarity table
indepGO=function(dissim.table,min.similarity=1) {
	tt=read.table(dissim.table,sep="\t", header=TRUE)
	tt=as.matrix(tt)
	diag(tt)=1
	for (i in 1:ncol(tt)) {	
		mins=apply(tt,2,min)
		if (min(mins)>=min.similarity) { break }
		sums=apply(tt,2,sum)
		worst=which(sums==min(sums))[1]
		# cat("\n",worsts,"\n")
		# gw=c()
		# for(w in worsts) { gw=append(gw,sum(tt[,w]==1)) }
		# cat(gw)
		# worst=worsts[gw==min(gw)]
#		cat("\n",i," removing",worst,"; sum =",sums[worst])
		tt=tt[-worst,-worst]		
		mins=mins[-worst]		
#		cat(" new min =",min(mins))
	}
	goods=colnames(tt)
	goods=gsub("GO\\.","GO:",goods)
	goods=gsub("\\.GO",";GO",goods)
}

