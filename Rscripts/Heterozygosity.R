library(ggplot2)

#Load table and output order of the samples
Inbreeding_table<-read.table("Heterozygosity.txt")
order<-readLines("Samples")
Inbreeding_table$sample<-factor(Inbreeding_table$sample,levels=order,ordered=TRUE)

#Compute subspecies mean heterozygosity, excluding Donald (central-western hybrid)
mean_hets<-by(Inbreeding_table[Inbreeding_table$sample!="9730_Donald",]$Het, Inbreeding_table[Inbreeding_table$sample!="9730_Donald",]$Subspecies, mean)
means=sort(as.numeric(mean_hets))

#Also standard deviations. These are not really used in the plot
sd_hets<-by(Inbreeding_table[Inbreeding_table$sample!="9730_Donald",]$Het, Inbreeding_table[Inbreeding_table$sample!="9730_Donald",]$Subspecies, sd)
lower=as.numeric(mean_hets)-as.numeric(sd_hets)
upper=as.numeric(mean_hets)+as.numeric(sd_hets)

#Range of the different subspecies
n=c(18,19,10,11,10)
n=n[5:1]

#Colots and labels
colors2=c("purple","#E31A1C","#FF7F00","#33A02C","#1F78B4")
labelcols=unlist(lapply(1:5, function(i) rep(colors2[i],each=n[i])))

#Plot
h <- ggplot(Inbreeding_table[Inbreeding_table$sample!="9730_Donald",],aes(x=Het,y=sample)) + geom_point(aes(color=Subspecies),size=5) + theme_bw()
h <- h + scale_color_manual(name="Population",labels=c("Bonobo","Nigeria-Cameroon","Eastern","Central","Western"),values=colors2) + theme(axis.text.y=element_text(color=labelcols))
h <- h + geom_point(aes(color=Subspecies),size=0.8,color="black") + annotate("segment",x=means,xend=means,y=c(0,10,22,32,51),yend=c(10,22,31.3,50.4,69),linetype=2,color=c("purple","#1F78B4","#E31A1C","#FF7F00","#33A02C"))
h <- h + ylab("") + xlab("\nHeterozygosity") + theme(panel.background = element_rect(colour = "black",size=0.7)) + theme(legend.key = element_blank())
h + theme(legend.title=element_text(face="bold",size=14),axis.title.x=element_text(vjust=-5,size=14),legend.text=element_text(size=14),axis.text.x=element_text(size=12),axis.text.y=element_text(size=11))
