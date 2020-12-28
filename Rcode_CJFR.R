#######################################################################
rm(list=ls())
################
library(scales)
show_col(viridis_pal()(10))

setwd("C:\\Users\\xjw5\\Dropbox\\NAU\\Data\\Siberia\\Final files CJFR")
library(MuMIn)
library(dplR)
library(graphics)
library(utils)
library(sciplot)
library(dplR) # Dendrochronology Program Library in R
library(bootRes) # Bootstrapped response and correlation functions
library(matrixStats) 
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
#####
#Annual trends in climate, growth, productivity, and NDVI (OBJECTIVE 1)

#climate

climate<-read.csv("climate.csv")

climate=filter(climate, year<2012)
climate=filter(climate, year>1979)

names(climate)
x1=lm(data=climate, Pfall.temp~year)
summary(x1)
mean(x1$residuals)

plot(residuals(x1))
plot(acf(resid(x1)))

library(VGAM)
kendall.tau(climate$Pfall.temp, climate$year, exact = FALSE, max.n = 3000)
cor.test(climate$Pfall.temp, climate$year, method=c("pearson", "kendall", "spearman"))


x1=lm(data=climate, winter.temp~year)
summary(x1)
plot(residuals(x1))
plot(acf(resid(x1)))
cor.test(climate$winter.temp, climate$year, method=c("pearson", "kendall", "spearman"))


x1=lm(data=climate, sping.temp~year)
summary(x1)
plot(residuals(x1))
cor.test(climate$sping.temp, climate$year, method=c("pearson", "kendall", "spearman"))


x1=lm(data=climate, gs.temp~year)
summary(x1)
plot(residuals(x1))
cor.test(climate$gs.temp, climate$year, method=c("pearson", "kendall", "spearman"))


head(climate)
x1=lm(data=climate, Pfall.precip~year)
summary(x1)
plot(residuals(x1))
cor.test(climate$Pfall.precip, climate$year, method=c("pearson", "kendall", "spearman"))

x1=lm(data=climate, winter.precip~year)
summary(x1)
plot(residuals(x1))
cor.test(climate$winter.precip, climate$year, method=c("pearson", "kendall", "spearman"))


x1=lm(data=climate, gs.precip~year)
summary(x1)
plot(residuals(x1))
cor.test(climate$gs.precip, climate$year, method=c("pearson", "kendall", "spearman"))


x1=lm(data=climate, spring.precip~year)
summary(x1)
plot(residuals(x1))
cor.test(climate$spring.precip, climate$year, method=c("pearson", "kendall", "spearman"))


#tree growth
all<-read.csv("bai.csv")
names(all)
levels(all$site)
all$density.names<-factor(all$density.names, levels=c("Low", "Med", "High", "Very High"))

all<-mutate(all, density.names2= ifelse(density.names=="Very High"|density.names=="High", "High", 
                                        ifelse(density.names=="Med", "Med", "Low")))

all$density.names2<-factor(all$density.names2, levels=c("Low", "Med", "High"))


names(all)
str(all)

all.1980=filter(all, year>1979&year<2012)
names(all.1980)

all.19801=dplyr::select(all.1980, tree, site, year, bai, dens, density.names2)

all.19801=na.omit(all.19801)

lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000)
control=lmeControl(opt = "optim")

tree.bai.ar1 <-lme(log(bai) ~ year*density.names2, random=~1|site/year, correlation=corAR1(), data = all.19801, na.action=na.omit, method="ML")
plot(tree.bai.ar1)
summary(tree.bai.ar1)
r.squaredGLMM(tree.bai.ar1)

acf(residuals(tree.bai.ar1, type="normalized"))

tree.bai <-lme(log(bai) ~ year*density.names2, random=~1|site/tree, data = all.19801, na.action=na.omit, method="ML")
plot(tree.bai)
r.squaredGLMM(tree.bai)
acf(residuals(tree.bai, type="normalized"))

anova(tree.bai, tree.bai.ar1)

tree.bai.add<-lme(log(bai) ~ year+density.names2, random=~1|site/tree, data = all.19801, na.action=na.omit, method="ML")
anova(tree.bai.add, tree.bai)

tree.bai.final <-lme(log(bai) ~ year*density.names2, random=~1|site/tree, data = all.19801, na.action=na.omit, method="REML")
plot(tree.bai.final)
summary(tree.bai.final)
r.squaredGLMM(tree.bai.final)
acf(residuals(tree.bai.final, type="normalized"))

E <- residuals(tree.bai.final, type = "normalized")
plot(E~all.19801$year)
plot(E~all.19801$density.names2)
plot(E~all.19801$site)
plot(E~all.19801$tree)


dens<- emtrends(tree.bai.final, ~ density.names2|year, var = "year",adjust="tukey", transform="response")
summary(dens,infer=c(TRUE,TRUE))
pairs(dens)

densar<- emtrends(tree.bai.ar1, ~ density.names2|year, var = "year",adjust="tukey", transform="response")
summary(densar,infer=c(TRUE,TRUE))
pairs(densar)


# stand productivity
names(all.1980)
library(dplyr)
detach(package:plyr)

names(all.1980)

all.1980=mutate(all.1980, dens.ha=dens*10000)
all.1980=mutate(all.1980, bai.m2=bai/1000000)

all.sum.BAI=all.1980%>%
  dplyr::group_by(year,  site, density.names2)%>%
  dplyr:: summarise(mean_bai=mean(bai.m2),
                    mean_dens=mean(dens.ha))
head(all.sum.BAI)
all.sum.BAI<-mutate(all.sum.BAI, bai.site=(mean_bai*mean_dens))
library(glmmTMB)
library(lme4)
all.sum.BAI=as.data.frame(all.sum.BAI)
all.sum.BAI=na.omit(all.sum.BAI)

site.bai <-lme(log(bai.site) ~ year*density.names2, random=~1|site,  data = all.sum.BAI, na.action=na.omit, method="ML")
plot(site.bai)
summary(site.bai)
r.squaredGLMM(site.bai)


site.bai.ar <-lme(log(bai.site) ~ year*density.names2, random=~1|year, correlation=corAR1(), data = all.sum.BAI, na.action=na.omit, method="ML")
summary(site.bai.ar)
anova(site.bai.ar, site.bai)


site.bai.add <-lme(log(bai.site) ~ year+density.names2, random=~1|site,  data = all.sum.BAI, na.action=na.omit, method="ML")
plot(site.bai.add)
anova(site.bai.add, site.bai)

site.bai.final <-lme(log(bai.site) ~ year*density.names2,random=~1|site,  data = all.sum.BAI, na.action=na.omit, method="REML")
plot(site.bai.final)
AIC(site.bai.final)

summary(site.bai.final)
r.squaredGLMM(site.bai.final)
acf(residuals(site.bai.final, type="normalized"))

E <- residuals(site.bai.final, type = "normalized")
plot(E~all.sum.BAI$year)
plot(E~all.sum.BAI$density.names2)
plot(E~all.sum.BAI$site)

library(emmeans)
dens<- emtrends(site.bai.final, ~ density.names2|year, var = "year",adjust="tukey", transform="response")
summary(dens,infer=c(TRUE,TRUE))
pairs(dens)


densar<- emtrends(site.bai.ar, ~ density.names2|year, var = "year",adjust="tukey", transform="response")
summary(densar,infer=c(TRUE,TRUE))
pairs(densar)

#ndvi

ndvi<-read.csv("ndvi.csv")


all.ndvi<-merge(ndvi, all.19801, by=c("site", "year"), all=TRUE)
head(all.ndvi)
names(all.ndvi)
levels(all.ndvi$site)
levels(all.ndvi$density.names2)


all.ndvi1<-dplyr::select(all.ndvi, site, tree, year, ndvi.max, bai, dens, density.names2)
all.ndvi1<-na.omit(all.ndvi1)
levels(all.ndvi1$site)

names(all.ndvi1)
all.ndvi1=filter(all.ndvi1, year<2012)
all.ndvi1=filter(all.ndvi1, year>1998)
all.ndvi1$density.names2

names(all.ndvi1)
head(all.ndvi1)


all.ndvi1=mutate(all.ndvi1, dens.ha=dens*10000)
all.ndvi1=mutate(all.ndvi1, bai.m2=bai/1000000)

all.ndvi.site=all.ndvi1%>%
  dplyr::group_by(year,  site, density.names2)%>%
  dplyr:: summarise(mean_band=mean(ndvi.max), 
                    mean_bai=mean(bai.m2),
                    mean_dens=mean(dens.ha))


all.ndvi.site=mutate(all.ndvi.site, bai.site=mean_dens*mean_bai)
head(all.ndvi.site)
all.ndvi.site=na.omit(all.ndvi.site)

##
levels(all.ndvi.site$site)

ndvi.ar<-lme(mean_band~year*density.names2,  random=~1|year, correlation=corAR1(),data=all.ndvi.site, na.action=na.omit, method = "ML")
summary(ndvi.ar)
plot(ndvi.ar)
AIC(ndvi.ar)
plot(acf(residuals(ndvi.ar)))


ndvi.ar1<-lme(mean_band~year+density.names2,  random=~1|year, correlation=corAR1(),data=all.ndvi.site, na.action=na.omit, method = "ML")
summary(ndvi.ar1)
plot(ndvi.ar1)
AIC(ndvi.ar1)
anova(ndvi.ar1, ndvi.ar)


ndvi.ar2<-lme(mean_band~year,  random=~1|year, correlation=corAR1(),data=all.ndvi.site, na.action=na.omit, method = "ML")
anova(ndvi.ar1,ndvi.ar2)
summary(ndvi.ar2)


ndvi<-lme(mean_band~year*density.names2,  random=~1|site, data=all.ndvi.site, na.action=na.omit, method = "ML")
summary(ndvi)
AIC(ndvi)
anova(ndvi, ndvi.ar1)


ndvi1<-lme(mean_band~year+density.names2, random=~1|site,data=all.ndvi.site, na.action=na.omit, method = "ML")
summary(ndvi1)
anova(ndvi1, ndvi)

ndvi2<-lme(mean_band~year, random=~1|site,data=all.ndvi.site, na.action=na.omit, method = "ML")
summary(ndvi2)
anova(ndvi2, ndvi.ar2)


ndvi3<-lme(mean_band~1, random=~1|site,data=all.ndvi.site, na.action=na.omit, method = "ML")
summary(ndvi3)
anova(ndvi2, ndvi3)


ndvi.final<-lme(mean_band~year,  random=~1|site, data=all.ndvi.site, na.action=na.omit, method = "REML")
summary(ndvi.final)
r.squaredGLMM(ndvi.final)
plot(ndvi.final)
plot(acf(residuals(ndvi.final)))
E <- residuals(ndvi.final, type = "normalized")
plot(E~all.ndvi.site$year)
plot(E~all.ndvi.site$density.names2)
plot(E~all.ndvi.site$site)



#####
### Impacts of density on stand structure, tree growth, stand productivity, and NDVI (OBJECTIVE 2)
site1<-read.csv("site_data.csv")

y=aov(data=site1, age.mean~density.names2)
anova(y)
summary(y)
plot(y)
TukeyHSD(y)

levels(site1$site)

site1=mutate(site1, dens.ha=dens*10000)

site.ha=site1%>%
  group_by(density.names2)%>%
  summarise(max=max(dens.ha), 
            min=min(dens.ha))


site1=mutate(site1, larch.bio.ha=larch.bio/100)
site1=mutate(site1, tot.under.ha=tot.under/100)
site1=mutate(site1, total.above.ha=total.above/100)


y=aov(data=site1, dens.ha~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, diam~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, ba~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, td~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, sol~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, larch.bio.ha~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, tot.under.ha~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, total.above.ha~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

y=aov(data=site1, can.cov~density.names2)
summary(y)
plot(y)
TukeyHSD(y)

x=site1%>%
  group_by(density.names2)%>%
  summarise(y=max(dens), 
            x=min(dens))
x=site1%>%
  group_by(density.names2)%>%
  summarise(y=max(can.cov), 
            x=min(can.cov))

names(site1)
site2=dplyr::select(site1, age, dens.ha, dens, diam, ba, td, sol, larch.bio, larch.bio.ha,can.cov, tot.under, total.above, tot.under.ha, total.above.ha, snag, density.names2)
library(plotrix)

#Table 1 - dofferences in stand structure
all.sum<-site2%>%
  group_by(density.names2)%>%
  select_if(is.numeric) %>% 
  summarise_all(funs(mean, std.error))

head(all.sum)

all.sumt=t(all.sum)
head(all.sumt)

# table 1 - difference in tree growth
all.19801=na.omit(all.19801)

all.19801.sum=all.19801%>%
  group_by(density.names2)%>%
  summarize(mean.tree.bai=mean(bai), 
            se.tree.bai=std.error(bai))
all.19801.sum

tree.dens<-lme(log(bai) ~ density.names2, random=~1|site/year, data = all.19801, na.action=na.omit, method="ML")
summary(tree.dens)

tree.dens1<-lme(log(bai) ~ 1, random=~1|site/year, data = all.19801, na.action=na.omit, method="ML")
summary(tree.dens1)
anova(tree.dens, tree.dens1)

tree.dens.final <-lme(log(bai) ~ density.names2, random=~1|site/year, data = all.19801, na.action=na.omit, method="REML")
plot(tree.dens.final)
summary(tree.dens.final)
r.squaredGLMM(tree.dens.final)
acf(residuals(tree.dens.final, type="normalized"))
E <- residuals(tree.dens.final, type = "normalized")
plot(E~all.19801$year)
plot(E~all.19801$density.names2)
plot(E~all.19801$site)


t <- emmeans(tree.dens.final, ~ density.names2, var = "density.names2", adjust="tukey", transform="response")
summary(t, infer=c(TRUE, TRUE))
pairs(t)

summary(final.tree1)

# table 1 - difference in site growth
head(all.sum.BAI)
all.sum.BAI=na.omit(all.sum.BAI)

all.sum.BAI.1=all.sum.BAI%>%
  group_by(density.names2)%>%
  summarize(mean.site.bai=mean(bai.site), 
            se.site.bai=std.error(bai.site))

all.sum.BAI.1

names(all.sum.BAI)
site.dens<-lme(log(bai.site) ~ density.names2, random=~1|year, data = all.sum.BAI, na.action=na.omit, method="ML")
summary(site.dens)

site.dens1<-lme(log(bai.site) ~ 1, random=~1|year, data = all.sum.BAI, na.action=na.omit, method="ML")
summary(site.dens1)
anova(site.dens, site.dens1)

site.dens.final <-lme(log(bai.site) ~ density.names2, random=~1|year, data = all.sum.BAI, na.action=na.omit, method="REML")
plot(site.dens.final)
summary(site.dens.final)
r.squaredGLMM(site.dens.final)
acf(residuals(site.dens.final, type="normalized"))
E <- residuals(site.dens.final, type = "normalized")
plot(E~all.sum.BAI$year)
plot(E~all.sum.BAI$density.names2)


t <- emmeans(site.dens.final, ~ density.names2, var = "density.names2", adjust="tukey", transform="response")
summary(t, infer=c(TRUE, TRUE))
pairs(t)

summary(final.site1)



#Table 1 - differences in ndvi
all.ndvi1=na.omit(all.ndvi1)

all.ndvi1.sum=all.ndvi1%>%
  group_by(density.names2)%>%
  summarize(mean.ndvi.max=mean(ndvi.max), 
            se.ndvi.max=std.error(ndvi.max))

y1<-lme(ndvi.max~density.names2,  random=~1|year, data=all.ndvi1, na.action=na.omit, method = "ML")
summary(y1)
y<-lme(ndvi.max~1,  random=~1|year, data=all.ndvi1, na.action=na.omit, method = "ML")
anova(y, y1)

y1<-lme(ndvi.max~density.names2,  random=~1|year, data=all.ndvi1, na.action=na.omit, method = "REML")
summary(y1)
r.squaredGLMM(y1)
plot(y1)
acf(residuals(y1))
E <- residuals(y1, type = "normalized")
plot(E~all.ndvi1$year)
plot(E~all.ndvi1$density.names2)


library(emmeans)
NDVI<- emmeans(y1, ~ density.names2,  adjust="tukey", type="response")
summary(NDVI,infer=c(TRUE,TRUE))
pairs(NDVI)



######
##Density dependent climate growth analyses(OBJECTIVE 3) 

bai.climate=merge(all.19801, climate, all=TRUE, by="year")
head(bai.climate)

bai.site.clim=merge(bai.climate, site, all=TRUE, by=c("site", "density.names2"))
head(bai.site.clim)

names(bai.site.clim)
##
head(bai.site.clim)
full <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
             scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
             scale(winter.temp)*density.names2 + 
             scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
           +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
           data = bai.site.clim, na.action=na.omit, method="ML")
summary(full)

AIC(full)
acf(residuals(full))

full.REML <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
             scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
             scale(winter.temp)*density.names2 + 
             scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
           +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
           data = bai.site.clim, na.action=na.omit, method="REML")
summary(full.REML) # data in table S5

# Pspring temp 

a <-lme(log(bai) ~   scale(Pspring.temp)+density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree, 
        data = bai.site.clim, na.action=na.omit, method="ML")
summary(a)

anova(a, full)

###Ps pring precip 

b <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)+density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree, 
        data = bai.site.clim, na.action=na.omit, method="ML")

anova(full, b)

##Pgs temp 
c <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)+density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2,random=~1|site/tree, 
        data = bai.site.clim, na.action=na.omit, method="ML")

anova(full, c)


## Pgs precip 
d <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)+density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")

anova(full, d)


##Pfall.temp 
e <-lme(log(bai) ~  scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)+density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")

anova(full, e)


### Pfall.precip 
f <-lme(log(bai) ~  scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)+density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2,random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")

anova(full, f)

###winter temp 
g <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)+density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")
anova(full, g)

###winter precip 
h <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)+density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2,random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")
anova(full, h)

###spring temp 
i <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)+density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")
anova(full, i)

##spring precip 
j <-lme(log(bai) ~  scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)+density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")
anova(full, j)

##gs temp 
k <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)+density.names2 + scale(gs.precip)*density.names2, random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")
anova(full, k)

#gs precip 
l <-lme(log(bai) ~   scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)*density.names2 +
          scale(Pgs.precip)*density.names2+ scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+ 
          scale(winter.temp)*density.names2 + 
          scale(winter.precip)*density.names2 +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
        +scale(gs.temp)*density.names2 + scale(gs.precip)+density.names2,random=~1|site/tree,
        data = bai.site.clim, na.action=na.omit, method="ML")
anova(full, l)



newfull <-lme(log(bai) ~    scale(Pspring.temp)*density.names2 + scale(Pspring.precip)*density.names2+ scale(Pgs.temp)+scale(Pgs.precip)*density.names2
              + scale(Pfall.temp)*density.names2 + scale(Pfall.precip)*density.names2+scale(winter.temp)*density.names2 + 
                scale(winter.precip) +scale(sping.temp)*density.names2+scale(spring.precip)*density.names2
              +scale(gs.temp)*density.names2 + scale(gs.precip), random=~1|site/tree,
              data = bai.site.clim, na.action=na.omit, method="REML")

summary(newfull)
head(bai.site.clim)

AIC(newfull)
plot(newfull)
E <- residuals(newfull, type = "normalized")
I1 <- !is.na(bai.site.clim$bai)
Efull <- vector(length = length(bai.site.clim$bai))
Efull <- NA
Efull[I1] <- E
acf(Efull, na.action = na.pass,
    main = "Auto-correlation plot for residuals")


r.squaredGLMM(newfull)

emm_options(opt.digits = FALSE)
emm_options(pbkrtest.limit = 4936)
library(emmeans)

Pspring.slope <- emtrends(newfull, ~ density.names2|Pspring.temp, var = "Pspring.temp", adjust="tukey", transform="response")
summary(Pspring.slope,infer=c(TRUE,TRUE))
pairs(Pspring.slope)


Pgs.slope <- emtrends(newfull, ~ Pgs.temp, var = "Pgs.temp",adjust="tukey", transform="response")
summary(Pgs.slope,infer=c(TRUE,TRUE))
pairs(Pgs.slope)


Pfall.slope <- emtrends(newfull, ~ density.names2|Pfall.temp, var = "Pfall.temp",adjust="tukey", transform="response")
summary(Pfall.slope,infer=c(TRUE,TRUE))
pairs(Pfall.slope)


winter.slope <- emtrends(newfull, ~ density.names2|winter.temp, var = "winter.temp",adjust="tukey", transform="response")
summary(winter.slope,infer=c(TRUE,TRUE))
pairs(winter.slope)


sping.slope <- emtrends(newfull, ~ density.names2|sping.temp, var = "sping.temp",adjust="tukey", transform="response")
summary(sping.slope,infer=c(TRUE,TRUE))
pairs(sping.slope)


gs.slope <- emtrends(newfull, ~ density.names2|gs.temp, var = "gs.temp",adjust="tukey", transform="response")
summary(gs.slope,infer=c(TRUE,TRUE))
pairs(gs.slope)

precip.Pspring.slope <- emtrends(newfull, ~ density.names2|Pspring.precip, var = "Pspring.precip",adjust="tukey", transform="response")
summary(precip.Pspring.slope,infer=c(TRUE,TRUE))
pairs(precip.Pspring.slope)


precip.Pgs.slope <- emtrends(newfull, ~ density.names2|Pgs.precip, var = "Pgs.precip",adjust="tukey", transform="response")
summary(precip.Pgs.slope,infer=c(TRUE,TRUE))
pairs(precip.Pgs.slope)

precip.Pfall.slope <- emtrends(newfull, ~ density.names2|Pfall.precip, var = "Pfall.precip",adjust="tukey", transform="response")
summary(precip.Pfall.slope,infer=c(TRUE,TRUE))
pairs(precip.Pfall.slope)


precip.winter.slope <- emtrends(newfull, ~ winter.precip, var = "winter.precip",adjust="tukey", transform="response")
summary(precip.winter.slope,infer=c(TRUE,TRUE))
pairs(precip.winter.slope)

precip.spring.slope <- emtrends(newfull, ~ density.names2|spring.precip, var = "spring.precip",adjust="tukey", transform="response")
summary(precip.spring.slope,infer=c(TRUE,TRUE))
pairs(precip.spring.slope)

precip.gs.slope <- emtrends(newfull, ~ gs.precip, var = "gs.precip",adjust="tukey", transform="response")
summary(precip.gs.slope,infer=c(TRUE,TRUE))
pairs(precip.gs.slope)





#####
#NDVI and stand productivity analysis
names(all.ndvi.site)
head(all.ndvi.site)

low1<-lme(mean_band~bai.site, method="ML",random=~1|site,data =subset(all.ndvi.site, density.names2=="Low"), na.action=na.omit)
summary(low1)
low=lme(mean_band~1, method="ML",random=~1|site,data =subset(all.ndvi.site, density.names2=="Low"), na.action=na.omit)
anova(low, low1)

med1<-lme(mean_band~bai.site, method="ML",random=~1|site,data =subset(all.ndvi.site, density.names2=="Medium"), na.action=na.omit)
summary(med1)
med<-lme(mean_band~1, method="ML",random=~1|site,data =subset(all.ndvi.site, density.names2=="Medium"), na.action=na.omit)
anova(med1, med)

high1<-lme(mean_band~bai.site, method="ML",random=~1|site,data =subset(all.ndvi.site, density.names2=="High"), na.action=na.omit)
summary(high1)
high<-lme(mean_band~1, method="ML",random=~1|site,data =subset(all.ndvi.site, density.names2=="High"), na.action=na.omit)
anova(high, high1)

high.final<-lme(mean_band~bai.site, method="REML",random=~1|site,data =subset(all.ndvi.site, density.names2=="High"), na.action=na.omit)
summary(high.final)
r.squaredGLMM(high.final)
plot(high.final)
E=residuals(high.final)
acf(E)
high.data=subset(all.ndvi.site, density.names2=="High")
plot(E~high.data$bai.site)
plot(E~high.data$site)

med.final<-lme(mean_band~bai.site, method="REML",random=~1|site, data =subset(all.ndvi.site, density.names2=="Medium"), na.action=na.omit)
summary(med.final)
r.squaredGLMM(med.final)

plot(med.final)
E=residuals(med.final)
acf(E)
med.data=subset(all.ndvi.site, density.names2=="Medium")
plot(E~med.data$bai.site)
plot(E~med.data$site)

low.final<-lme(mean_band~bai.site, method="REML",random=~1|site,data =subset(all.ndvi.site, density.names2=="Low"), na.action=na.omit)
summary(low.final)
plot(low.final)
r.squaredGLMM(low.final)

plot(low.final)
E=residuals(low.final)
acf(E)
low.data=subset(all.ndvi.site, density.names2=="Low")
plot(E~low.data$bai.site)
plot(E~low.data$site)


