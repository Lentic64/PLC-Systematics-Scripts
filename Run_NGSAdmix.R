library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
setwd('/Users/dangates/Desktop/Isoetes/NGSAdmix/')

#I run out ngsadmix with this bash script
#for j in {1..10}; do for i in {1..8}; do /Users/dangates/Downloads/angsd/misc/NGSadmix -likes ../Aligns/popFinal.beagle -K $i -outfiles "NGSA${i}_${j}.out"; done; done


#Lets looks at NGAdmix logs
#based on: https://github.com/alexkrohn/AmargosaVoleTutorials/blob/master/ngsAdmix_tutorial.md
fils<-list.files()[grep('.log',list.files())]
dat<-t(sapply(fils,function(x) {
  fil<-readLines(x)
  pop<-as.numeric(strsplit(strsplit(fil[1],'nPop=')[[1]][2],', ')[[1]][1])
  lik<-as.numeric(strsplit(strsplit(fil[9],'like=')[[1]][2],' after')[[1]][1])
  return(c(pop,lik))
}))

write.table(dat,file='../NGSAdmixForClumpack.txt',quote=FALSE,row.names = FALSE,sep='\t',col.names = FALSE)
  
#This goes to clumpak at :http://clumpak.tau.ac.il/bestK.html and it's the log probability file

#it thinks that k=4 is best for us


#A simple script to look at ngsadmix output:

anc<-fread('./NGSA4_1.out.qopt',data.table=FALSE)
colnames(anc)<-paste('grp',1:ncol(anc),sep="")
oN<-fread('../Aligns/popFinal.beagle',data.table=FALSE)
inds<-colnames(oN[,c(4:ncol(oN))[seq(from=1,to=length(4:ncol(oN)),by=3)]])

anc$inds<-inds

anc %>% gather(key=key,value=value,-inds) -> gplot

#Now the best way to order this?
acc<-readLines('../bamfileList4.txt')
pop<-sapply(acc,function(x) strsplit(x,'_')[[1]][1])
sp<-sapply(acc,function(x) strsplit(strsplit(x,'Isoetes_')[[1]][2],'.bam')[[1]][1])

trgt<-as.numeric(sapply(anc$inds,substr,4,10))+1

newp<-pop[trgt]
newsp<-sp[trgt]

gplot$pop<-newp
gplot$sp<-newsp

gplot<-gplot[with(gplot,order(sp,pop)),]

gplot$new_ind<-paste(gplot$sp,gplot$pop,gplot$inds,sep="_")

png("admixPlotNew.png",width=4000,height=1000,res=250)
ggplot(gplot,aes(fill=key,y=value,x=new_ind))+
  geom_bar(position="stack",stat="identity")+
  theme(axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1))
dev.off()




