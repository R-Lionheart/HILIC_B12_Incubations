## NMDS Analysis testing
source("src/biostats.R")
library(vegan)
library(pastecs) # masks dplyr + tidyr

## Start analysis

## Lecture 2: Data Screening

envdata <- read.csv('data_extras/MAHA_environment.csv', header=TRUE, row.names=1)
spedata <- read.csv('data_extras/MAHA_speciesabu.csv', header=TRUE, row.names=1)

str(envdata)

stat.desc(envdata)
stat.desc(spedata) # Column summary

testdata <- replace.missing(spedata)
testdata <- replace.missing(spedata,method='mean')

testdata <- drop.var(spedata, pct.missing=5)

foa.plots(spedata)

testdata <- drop.var(spedata,min.po=5)
testdata <- drop.var(spedata,min.fo=5)
str(spedata)
str(testdata)

testdata <- drop.var(spedata,max.po=95)
testdata <- drop.var(spedata,min.cv=5)

ecdf.plots(spedata)
hist.plots(spedata)
box.plots(spedata)
qqnorm.plots(spedata)
uv.plots(spedata)
#
data.trans(spedata,method='log')
data.trans(spedata,method='power',exp=.5)
data.trans(spedata,method='power',exp=0)
data.trans(spedata,method='asin')
?decostand

tradata <- data.trans(spedata,method='log',plot=F)
tradata <- log10(spedata+1)

data.stand(spedata,method='total',margin='row')

uv.outliers(envdata, id='Sinuosity:BasinAre', var='Elev', sd.limit=1)

mv.outliers(envdata,method='euclidean', sd.limit=1)
mv.outliers(envdata, method='mahalanobis', sd.limit=1)



####################
## Lecture 8
# Non metric dimensional scaling (NMDS)
speabu <- read.csv("data_extras/MAHA_speciesabu.csv", header = TRUE, row.names = 1)
speabu.log<-log(speabu+1)

spe.nmds<-metaMDS(speabu.log, distance="bray", k=2, autotransform=FALSE,
                  trymax=100)
spe.nmds
names(spe.nmds)

spe.nmds2 <- metaMDS(speabu.log, distance='bray', k=3, autotransform=FALSE,
                     trymax=100)

nmds.scree(speabu.log, distance='bray', k=10, autotransform=FALSE, trymax=20)

nmds.monte(speabu.log, distance = "bray", k = 3, autotransform = FALSE, trymax = 20)

stressplot(spe.nmds2)

plot(spe.nmds,type='n')
text(spe.nmds,labels=row.names(speabu))

plot(spe.nmds,type='n')
points(spe.nmds,cex=speabu.log$ROSYDACE)

vec.sp<-envfit(spe.nmds$points, speabu.log, perm=1000)
vec.sp

ordiplot(spe.nmds, choices = c(1, 2), type="text", display='sites',
         xlab='Axis 1', ylab='Axis 2')
plot(vec.sp, p.max=.01,col = "blue")


speabu.d<-vegdist(speabu.log, "bray")
sitecl.ward<-hclust(speabu.d,method='ward.D')
sitecl.class<-cutree(sitecl.ward,k=4)
groups<-levels(factor(sitecl.class))

site.sc <- scores(spe.nmds)
p <- ordiplot(site.sc, type="n", main="NMDS combined with clustering")
for (i in 1:length(groups))
{
  points(site.sc[sitecl.class==i,], pch=(14+i), cex=2, col=i+1)
}
text(site.sc, row.names(speabu), pos=4, cex=0.7)

legend(locator(1), paste("Group",c(1:length(groups))),
       pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)

ordicluster(p, sitecl.ward, col="dark grey")

legend(locator(1), paste("Group", c(1:length(groups))),
       pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)