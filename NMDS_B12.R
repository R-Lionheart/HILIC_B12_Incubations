## NMDS with B12 Data
source("src/biostats.R")
library(vegan)
library(pastecs) # masks dplyr + tidyr
library(tidyverse) 
#library(stringr)

# DSW700 needs to be removed- not relevant
# Split the dataset into separate eddies

# Uploads and sampID filtering --------------------------------------------
HILIC_data <- read.csv("data_processed/BMIS_Output_2019-12-09_duplicatesremoved.csv", stringsAsFactors = FALSE) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  filter(!SampID %in% c("Sept29QC", "TruePooWeek1", "TruePooWeek2", "TruePooWeek3", "TruePooWeek4")) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ",")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Area.Ave = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, Area.Ave) %>%
  unique()

HILIC_wide_mid <- HILIC_data %>%
  ungroup() %>%
  tidyr::spread(SampID, Area.Ave) %>%
  as.data.frame()

HILIC_wide <- HILIC_wide_mid[,-1]
rownames(HILIC_wide) <- HILIC_wide_mid[,1]

HILIC_wide <- data.frame(HILIC_wide)
HILIC_wide[is.na(HILIC_wide)] <- NA

HILIC_noNA <- na.omit(HILIC_wide)

rm("HILIC_data")
rm("HILIC_wide_mid")


# Structure exploration + data screening --------------------------------------------
class(HILIC_wide)
stat.desc(HILIC_wide)

# Drops the variable (ie SampID) that is missing > 5, rather than the compound. Will stick with tallying method.
# testdata <- drop.var(HILIC_wide, pct.missing=5)

# foa.plots(HILIC_wide) # Not great visualization overall

# The empirical cumulative distribution function (ECDF) 
# is a simple rank order distribution of increasing values of the variable.
ecdf.plots(HILIC_wide)

hist.plots(HILIC_wide)
# box.plots(HILIC_wide)
qqnorm.plots(HILIC_wide)
#uv.plots(HILIC_wide)


data.trans(HILIC_wide, method='log')
#data.trans(HILIC_wide, method='power',exp=.5)
#data.trans(HILIC_wide, method='power',exp=0)
#data.trans(HILIC_wide, method='asin')


testdata <- decostand(HILIC_wide, method = 'normalize', na.rm = TRUE)
hist.plots(testdata)

data.stand(HILIC_wide, method='standardize', na.rm = TRUE) # method = 'log'


uv.outliers(HILIC_wide, id = 'IL1Control:IL1DMB', var='IT0', sd.limit=3) # sd.limit=1
mv.outliers(HILIC_wide, method='euclidean', sd.limit=3) # sd.limit=1



# NMDS Visualization --------------------------------------------
HILIC_log <- log(HILIC_wide)
HILIC_log <- na.omit(HILIC_log) # Stopgap measure for now
HILIC.nmds<-metaMDS(HILIC_log, distance='bray', k=2, autotransform=FALSE, trymax=100)
HILIC.nmds
names(HILIC.nmds)

# In the current HILIC_log configuration, k>2 will not reach a stress test solution.
# Take another look at this section once the initial data munging has occurred.
# HILIC.nmds2 <- metaMDS(HILIC_log, distance='bray', k=3, autotransform=FALSE, trymax=100)

nmds.scree(HILIC_log, distance='bray', k=10, autotransform=FALSE, trymax=20)
nmds.monte(HILIC_log, distance='bray', k=3, autotransform=FALSE, trymax=20)

stressplot(HILIC.nmds)
plot(HILIC.nmds,type='n')
text(HILIC.nmds,labels=row.names(HILIC_noNA))

plot(HILIC.nmds,type='n')
points(HILIC.nmds, cex=HILIC_log$IL1Control)

HILIC.vec <- envfit(HILIC.nmds$points, HILIC_log, perm=1000)
HILIC.vec

ordiplot(HILIC.nmds, choices = c(1, 2), type='text', display='sites', xlab='Axis 1', ylab='Axis 2')
plot(HILIC.vec, p.max=.01, col='red')

# Combine clustering and ordination --------------------------------------------
HILIC.d<-vegdist(HILIC_log, "bray")
HILIC.sitecl.ward<-hclust(HILIC.d, method='ward.D')
HILIC.sitecl.class<-cutree(HILIC.sitecl.ward, k=4)
groups<-levels(factor(HILIC.sitecl.class))

HILIC.site.sc <- scores(HILIC.nmds)
p <- ordiplot(HILIC.site.sc, type="n", main="NMDS combined with clustering")
for (i in 1:length(groups))
{
  points(HILIC.site.sc[HILIC.sitecl.class==i,], pch=(14+i), cex=2, col=i+1)
}
text(HILIC.site.sc, row.names(HILIC_noNA), pos=4, cex=0.7)

legend(locator(1), paste("Group",c(1:length(groups))),
       pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)

ordicluster(p, HILIC.sitecl.ward, col="dark grey")

legend(locator(1), paste("Group", c(1:length(groups))),
       pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)
