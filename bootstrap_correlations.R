rm(list=ls())
################
#######  CHRONOLOGIES SIBERIA TREES
library(dplR)
library(graphics)
library(utils)
library(sciplot)
library(dplR) # Dendrochronology Program Library in R
library(bootRes) # Bootstrapped response and correlation functions
library(matrixStats) 
library(dplyr)
library(reshape2)
library(ggplot2)
library(AICcmodavg)
##########
##########
############
##########R###########
setwd("C:\\Users\\xjw5\\Documents\\NAU\\Data\\Siberia\\Siberia_TR_MB\\OldStuff")
############

############# NOTE that 20 trees only have onle measurement - ALL other trees have two measurements

ALL2<-read.csv("All_raw_ring_widths.csv",header=T) 
names(ALL2)

ALL<-dplyr::rename(ALL2, "year"="year")
names(ALL)
rownames<-ALL[,1]

ALL<-ALL[,c(2:151)]
names(ALL)

###convert to BAI
#ALL.bai<-bai.in(ALL)
#head(ALL.bai)


########################## now to the same with the spline (ie the detrended data)
############## 
####################################
require(bootRes)
head(ALL.bai)

tr.dat<-ALL.bai
dim(tr.dat)
head(tr.dat)
names(tr.dat)
tr.dat
ntrees<-ncol(tr.dat)
treeids<-colnames(tr.dat)

clim.dat<-read.csv("clim_1.csv")
names(clim.dat)
head(clim.dat)
climnames<-read.csv("clim_names.csv")
timespan<-c(1981,2011)

cor.dat<-matrix(nrow=34, ncol=ntrees)
colnames(cor.dat)<-treeids
#rownames(cor.dat)<-climnames
sig.dat<-matrix(nrow=34, ncol=ntrees)
colnames(sig.dat)<-treeids
#rownames(sig.dat)<-climnames
for (i in 2:ntrees){	
  tree<-matrix(nrow=nrow(tr.dat), ncol=1)
  row.names(tree)<-rownames
  treeid<-treeids[i]
  colnames(tree)<-treeid
  tree[,1]<-tr.dat[,i]
  dc<-dendroclim(tree, clim.dat, method="corr", timespan=timespan, start=-4,end=8)
  corrdat<-dc[,1]
  sigdat<-dc[,2]
  cor.dat[,i]<-corrdat
  sig.dat[,i]<-sigdat
}


write.csv(cor.dat, file="corr_Siberia_bai.csv")	
write.csv(sig.dat, file="sigcorr_Siberia_bai.csv")

######## check to make sure those  files look OK

############ here is some code to make a figure:
##########  

bai<-read.csv("corr_Siberia_bai.csv")
head(bai)
bai=dplyr::rename(bai, clim=variable)

bai=melt(bai, id="clim")
names(bai)
bai=dplyr::rename(bai, corr=value)
bai=dplyr::rename(bai, tree=variable)
head(bai)

bai=mutate(bai, temp.precip=ifelse(grepl("temp",bai$clim),'Temp','Precip'))
head(bai)

bai=mutate(bai, month=ifelse(grepl("prev.apr",bai$clim),'apr',
                             ifelse(grepl("prev.may",bai$clim),'may',
                                    ifelse(grepl("prev.jun",bai$clim),'jun',
                                           ifelse(grepl("prev.jul",bai$clim),'jul',
                                                  ifelse(grepl("prev.aug",bai$clim),'aug',
                                                         ifelse(grepl("prev.sep",bai$clim),'sep',
                                                                ifelse(grepl("prev.oct",bai$clim),'oct',
                                                                       ifelse(grepl("prev.nov",bai$clim),'nov',
                                                                              ifelse(grepl("prev.dec",bai$clim),'dec',
                                                                                     ifelse(grepl("curr.jan",bai$clim),'JAN',
                                                                                            ifelse(grepl("curr.feb",bai$clim),'FEB',
                                                                                                   ifelse(grepl("curr.mar",bai$clim),'MAR',
                                                                                                          ifelse(grepl("curr.apr",bai$clim),'APR',
                                                                                                                 ifelse(grepl("curr.may",bai$clim),'MAY',
                                                                                                                        ifelse(grepl("curr.jun",bai$clim),'JUN',
                                                                                                                               ifelse(grepl("curr.jul",bai$clim),'JUL',
                                                                                                                                      ifelse(grepl("curr.aug",bai$clim),'AUG', "stop"))))))))))))))))))
                                                                                                                                             
head(bai)

corr<-read.csv("sigcorr_Siberia_bai.csv")
head(corr)

library(reshape)
library(dplyr)
corr=melt(corr, id="clim")
names(corr)

corr=dplyr::rename(corr, sig.corr=value)
corr=dplyr::rename(corr, tree=variable)
head(corr)

corr2=merge(bai, corr, by=c("clim", "tree"))

head(corr2)


########## this is just my site level data
bai.site.clim<-read.csv("bai.site.clim.csv")
head(bai.site.clim)


bai.site.clim$density.names<-factor(bai.site.clim$density.names, levels=c("Low", "Med", "High", "Very High"))

bai.site.clim<-mutate(bai.site.clim, density.names2= ifelse(density.names=="Very High"|density.names=="High", "High", 
                                                            ifelse(density.names=="Med", "Medium", "Low")))

bai.site.clim$density.names2<-factor(bai.site.clim$density.names2, levels=c("Low", "Medium", "High"))

site=filter(bai.site.clim, year=="1980")

corr.site=merge(site, corr2, by="tree")
head(corr.site)

head(corr.site)
corr.site$temp.precip <- factor(corr.site$temp.precip, levels=c("Temp", "Precip"))
corr.site<-mutate(corr.site, temp.precip2= ifelse(temp.precip=="Temp", "Temperature", "Precipitation"))

corr.site$temp.precip2 <-factor(corr.site$temp.precip2,  levels=c("Temperature", "Precipitation"))
levels(corr.site$temp.precip2)
corr.site$month<-factor(corr.site$month, levels=c("apr", "may", "jun", "jul" ,"aug", "sep", "oct", "nov", "dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG"))
levels(corr.site$month)

corr.site$month <- fct_explicit_na(corr.site$month)
fct_count(corr.site$month)


levels(corr.site$density.names2)


names(corr.site)
corr.site=mutate(corr.site, sig.dens=interaction(sig.corr, density.names2))
head(corr.site)
levels(corr.site$sig.dens)


corr.site$corr<-as.numeric(corr.site$corr)
library(dplyr)


bai.site.sum <- corr.site %>%
  dplyr::group_by(temp.precip2, month, density.names2) %>%
  dplyr::summarise(corr_mean = mean(corr, na.rm=TRUE),
            corr_sd = sd(corr, na.rm=TRUE))

head(bai.site.sum)

levels(bai.site.sum$month)
levels(bai.site.sum$temp.precip)
levels(bai.site.sum$density.names2)
library(forcats)


names(bai.site.sum)
bai.site.sum<-mutate(bai.site.sum, y.min=corr_mean-corr_sd)
bai.site.sum<-mutate(bai.site.sum, y.max=corr_mean+corr_sd)
bai.site.sum <- mutate(bai.site.sum, sig=ifelse(y.min>0 |y.max<0, "TRUE", "FALSE"))
names(bai.site.sum)
head(bai.site.sum)

bai.site.sum=mutate(bai.site.sum, sig.dens=interaction(sig, density.names2))

levels(bai.site.sum$sig.dens)

bai.site.sum$sig.dens <- factor(bai.site.sum$sig.dens, levels=c("FALSE.Low", "FALSE.Medium", "FALSE.High", "TRUE.Low","TRUE.Medium" ,"TRUE.High"  ))


corr.site$sig.dens <- factor(corr.site$sig.dens, levels=c("FALSE.Low", "FALSE.Medium", "FALSE.High", "TRUE.Low","TRUE.Medium" ,"TRUE.High"  ))
levels(corr.site$sig.dens)
levels(bai.site.sum$sig.dens)

theme_set(theme_bw())


p=ggplot() 
p = p + theme(panel.border = element_rect(size=0.1))+ theme(panel.background = element_blank())
p
p=p+ geom_jitter(data=corr.site, aes(x=month, y=corr, fill=sig.dens), col="grey41",width=0.15, pch=21, size=1.5, stroke=0.01) 
p=p + theme(legend.position="none") 
p=p+scale_fill_manual(values=c("white" ,"white","white", "#A1D99B","#238B45","#00441B"))
p=p + geom_hline(yintercept = 0, size=0.1) 
p
#p=p+geom_errorbar(aes(x=month, ymin = y.min , ymax = y.max), width = 0.8, data=bai.site.sum, color="black")
p=p+ geom_point(data=bai.site.sum, aes(y=corr_mean, x=month, fill=sig.dens), size=4, pch=21)
p=p + facet_grid(density.names2~temp.precip2) +ylab("Correlation") +xlab("")
p=p+theme(strip.background = element_rect(colour = "black", fill = "white", size=0.1)) + theme(panel.spacing = unit(0.5, "lines"))
p=p + theme(axis.title= element_text(size=8), axis.text=element_text(size=8))
p=p+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(strip.text= element_text(size = 8))
p=p + theme(axis.line=element_line(size=0), axis.ticks=element_line(size=0.2))
p = p + theme(panel.grid  = element_line(size=0.1, color="gray90"))

p

