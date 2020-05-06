#to create input file: run quick.py > uniq > grep -wf this against contig_len4R.txt
#vcf must be bgzip and tabix 

library(PopGenome)
WbL3All.div.neut=data.frame()

	#Add input file name
input_vector <- read.delim("R_input.fb20.hwe.HCinds", header=FALSE)

	#choose one
#input_vector=subset(input_vector,V2>=5000) #width=5000,jump=1000
input_vector=subset(input_vector,V2>=10000) #width=10000, jump=2000

	#this is two populations that only read into the F_ST.stats function, the other by set.populations before slide will give mutliple columns
M<-c("WbPNG17B","WbPNG17D")
S<-c("WbPNG48_53","WbPNG48_73")
window=10000

popgenome_stats <- function(contig,length){
  	#Add input file name
  pngvcf<-readVCF("Wb.AA.outgroup.vcf.gz", tid=as.character(contig),frompos=1,topos=length, include.unknown=TRUE,numcols=3000)
  
  	#must have outgroup sample name
  pngvcf<-set.outgroup(pngvcf,"Bmal")
  	
  	#set population, but then every statistic will generate #POP columns
  #pngvcf<-set.populations(pngvcf,list(c("WbPNG17B","WbPNG51","WbPNG17D","WbPNG17E"),c("WbPNG48_53","WbPNG48_51","WbPNG48_73","WbPNG48B")),diploid=T)
  
    #choose to correspond with above
  #pngvcf<-sliding.window.transform(pngvcf,width=2000,jump=500,type=2,whole.data=TRUE) #input_vector
  #pngvcf<-sliding.window.transform(pngvcf,width=5000,jump=1000,type=2,whole.data=TRUE) #input_vector
  pngvcf<-sliding.window.transform(pngvcf,width=10000,jump=10000,type=2,whole.data=TRUE) #input_vector
  
  pngvcf=F_ST.stats(pngvcf,list(c(M),c(S)),mode="ALL")
  pngvcf=F_ST.stats.2(pngvcf,list(c(M),c(S)),snn=T,Phi_ST=T)
  
  Div_w=pngvcf@nuc.diversity.within/window
  
  pngvcf=diversity.stats(pngvcf,pi=T) #calculates diversity stats
  pngvcf=neutrality.stats(pngvcf,detail=F,do.R2=T)
  pngvcf=sweeps.stats(pngvcf)
  pngvcf=Achaz.stats(pngvcf)
  
  positions=gsub("\\s","",unlist(strsplit(pngvcf@region.names,":")))
  pi_nei=pngvcf@Pi/window #per site
  SegSites=pngvcf@n.segregating.sites #per window
  thetaW=pngvcf@theta_Watterson/window #per site
  TajD=pngvcf@Tajima.D
  FuLiF=pngvcf@Fu.Li.F
  FuLiD=pngvcf@Fu.Li.D
  FayWuH=pngvcf@Fay.Wu.H
  ZengE=pngvcf@Zeng.E
  FST=pngvcf@nucleotide.F_ST
  Snn=pngvcf@Hudson.Snn #between all pops
  #PhiST=pngvcf@Phi_ST #between all pops
  #FuFs=pngvcf@Fu.F_S #if add need: detail=T in neutrality.stats, add FuFs in tempd and colnames
  #Div_w=pngvcf@nuc.diversity.within/window
  Div_b=t(pngvcf@nuc.diversity.between/window) #pairwise
  RozasR2=pngvcf@Rozas.R_2
  CL=pngvcf@CL
  CLR=pngvcf@CLR
  CLmax=pngvcf@CLmax
  AchazY=pngvcf@Yach
  tempd=data.frame(rep(contig,length(TajD)),positions,SegSites,thetaW,pi_nei,TajD,FuLiF,FuLiD,FayWuH,ZengE,Snn,FST,Div_w[,1],Div_w[,2],Div_b,RozasR2,CL,CLR,AchazY)
  WbL3All.div.neut=rbind(WbL3All.div.neut,tempd)
  return(WbL3All.div.neut)
}

for (i in seq(1:length(input_vector$V1))){
  i
  WbL3All.div.neut=popgenome_stats(input_vector[i,1],input_vector[i,2])
  }

colnames(WbL3All.div.neut)=c("CHROM","POS","Seg","thetaW","PI_nei","TajD","FuLiF","FuLiD","FayWuH","ZengE","Snn","FST","Div_w1","Div_w2","Div_b","RozasR2","CL","CLR","AchazY")
	#Add output file name
write.table(WbL3All.div.neut,file="WbPNGL3-stats.fb20.hwe.HCinds.10000-Fst2.out",sep="\t",quote=FALSE,row.names=F,col.names=T)


######################## code for significance

attach(WbL3All.div.neut)
#FuLiF
quantile(FuLiF,probs=(seq(0,1,.025)),na.rm=T) #3 and 39 are 5%
two=quantile(FuLiF,probs=(seq(0,1,.025)),na.rm=T)[2]
ninety_seven=quantile(FuLiF,probs=(seq(0,1,.025)),na.rm=T)[40]
fulif.Neg=which(FuLiF < (two))
fulif.Pos=which(FuLiF > (ninety_seven))
f.pos=CHROM[fulif.Pos]
f.neg=CHROM[fulif.Neg]
#FuLiD
quantile(FuLiD,probs=(seq(0,1,.025)),na.rm=T) #3 and 39 are 5%
two=quantile(FuLiD,probs=(seq(0,1,.025)),na.rm=T)[2]
ninety_seven=quantile(FuLiD,probs=(seq(0,1,.025)),na.rm=T)[40]
fulid.Neg=which(FuLiD < (two))
fulid.Pos=which(FuLiD > (ninety_seven))
d.pos=CHROM[fulid.Pos]
d.neg=CHROM[fulid.Neg]
#TajD
quantile(TajD,probs=(seq(0,1,.025)),na.rm=T) #3 and 39 are 5%
two=quantile(TajD,probs=(seq(0,1,.025)),na.rm=T)[2]
ninety_seven=quantile(TajD,probs=(seq(0,1,.025)),na.rm=T)[40]
tajd.Neg=which(TajD < (two))
tajd.Pos=which(TajD > (ninety_seven))
t.pos=CHROM[tajd.Pos]
t.neg=CHROM[tajd.Neg]
#ZengE
quantile(ZengE,probs=(seq(0,1,.025)),na.rm=T) #3 and 39 are 5%
two=quantile(ZengE,probs=(seq(0,1,.025)),na.rm=T)[2]
ninety_seven=quantile(ZengE,probs=(seq(0,1,.025)),na.rm=T)[40]
zenge.Neg=which(ZengE < (two))
zenge.Pos=which(ZengE > (ninety_seven))
z.pos=CHROM[zenge.Pos]
z.neg=CHROM[zenge.Neg]

#use merge to combine/intersect lists significant Neg and Pos enteries: merge(WbL3All.div.neut[fulid.Neg,], WbL3All.div.neut[zenge.Neg,], by=c("CHROM","POS"))
#take this table to python and using selection2bed.py produce a bed file: "CHROM START STOP TAJD ZENGE"
	#possibly merge overlapping regions
#now use bedtools getfasta against the reference >> blast the resulting file


#intersect of Chromosomes, deal with this list for canidates of selection
x=intersect(t.neg,z.neg) #28, 36
y=intersect(f.neg,z.neg) #29, 35
z=intersect(d.neg,z.neg) #28, 36

i=intersect(t.pos,z.pos) #46, 36
j=intersect(f.pos,z.pos) #42, 40
k=intersect(d.pos,z.pos) #46, 43

####ploting PI_nei and snps positions
ggplot(subset(WbPNGL3.stats.NoMiss.10000,CHROM=="PairedContig_262"), aes(x=as.numeric(POS2),y=PI_nei, group=T)) + geom_line() + geom_abline(intercept=1.582418, colour="red",slope=0) + geom_bar(data=snps,aes(x=pdog,y=.3),stat="identity", width=1000)
WbPNGL3.stats.NoMiss.10000$POS2=gsub("-.+",'',WbPNGL3.stats.NoMiss.10000$POS2)

#remove Na
WbPopGenome=WbPNGL3.stats.fb20.hwe.HCinds.10000[complete.cases(WbPNGL3.stats.fb20.hwe.HCinds.10000[,1:11]),] #removes NAs

#to calculate joint DH from Esteve-Codin 2013
WbPopGenome=WbPopGenome[ order(WbPopGenome[,6]),] #orders on col6 here TajD
r=rank(WbPopGenome$TajD) #ranks TajD
d=WbPopGenome$TajD * r/length(WbPopGenome$TajD) #weights based on Ramos-Onins
norm_d=scale(d)
WbPopGenome$Norm_D=norm_d
WbPopGenome=WbPopGenome[ order(WbPopGenome[,9]),] #orders on col9 here H
r=rank(WbPopGenome$FayWuH) #ranks FayWuH
h=WbPopGenome$FayWuH * r/length(WbPopGenome$FayWuH) #weights based on Ramos-Onins
norm_h=scale(h)
WbPopGenome$Norm_H=norm_h
WbPopGenome=WbPopGenome[ order(WbPopGenome[,4]),] #orders on col4 here thetaW
r=rank(WbPopGenome$thetaW) #ranks thetaW
w=WbPopGenome$thetaW * r/length(WbPopGenome$thetaW) #weights based on Ramos-Onins
norm_w=scale(w)
WbPopGenome$Norm_W=norm_w

fn=ecdf(WbPopGenome$Norm_H)
h=sapply(WbPopGenome$Norm_H,function(x) fn(x))

fn=ecdf(WbPopGenome$Norm_D)
d=sapply(WbPopGenome$Norm_D,function(x) fn(x))

fn=ecdf(WbPopGenome$Norm_W)
w=sapply(WbPopGenome$Norm_W,function(x) fn(x))

WbPopGenome$Hquantile = h
WbPopGenome$Dquantile = d
WbPopGenome$Wquantile = w #some idication if it is balancing or directional selection

neg=subset(WbPopGenome,Hquantile<=.05 & Dquantile<=.05)
pos=subset(WbPopGenome,Hquantile>=.95 & Dquantile>=.95)
dpos=subset(WbPopGenome,Dquantile>=.95 & Hquantile<=.05)
dneg=subset(WbPopGenome,Dquantile<=.05 & Hquantile>=.95)

selection_jointDH=rbind(neg,pos,dpos,dneg)

#MKT
input_vector <- read.delim("R_input", header=FALSE)
d=data.frame()

for (i in seq(1:length(input_vector$V1))){
  i
pngvcf <- readVCF("Wb.AA.outgroup.vcf.gz",tid=as.character(input_vector[i,1]),frompos=1,topos=input_vector[i,2],include.unknown=TRUE,numcols=10000, gffpath=paste0("geneGFF/", as.character(input_vector[i,1])))
pngvcf <- set.synnonsyn(pngvcf, ref.chr=paste0("geneFASTA/",as.character(input_vector[i,1])),save.codons=TRUE)
pngvcf<-set.outgroup(pngvcf,"Bmal")
M<-c("WbPNG17B","WbPNG51","WbPNG17D","WbPNG17E","WbPNG48_53","WbPNG48_51","WbPNG48_73","WbPNG48B")
S<-"Bmal"
pngvcf<-set.populations(pngvcf,list(M,S),diploid=T)
#nonsyn
nsyn[i]=sum(pngvcf@region.data@synonymous[[1]]==0, na.rm=TRUE)
#syn
syn[i]=sum(pngvcf@region.data@synonymous[[1]]==1, na.rm=TRUE)
pngvcf <- MKT(pngvcf)
mkt=get.MKT(pngvcf)
d=rbind(d,mkt)
}





