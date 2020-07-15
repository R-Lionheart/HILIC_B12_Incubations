## Originally in NMDS_figs.R after the scree plot visualization


# Experiment with NMDS visualizations ------------------------------------------------------------------------
Iso_wide_nmds <- vegan::metaMDS(df_wide_normalizedT, distance = "euclidean", 
                                k = 3, autotransform = FALSE, trymax = 100, wascores = FALSE)

# Visualize scree plot of potential ordination axes
dimcheckMDS(df_wide_normalizedT, distance="euclidean", k=6, autotransform=FALSE, trymax=20) 
vegan::stressplot(Iso_wide_nmds, main = paste("Stressplot, Eddy", EddyInformation, sep = " "))

# Quick vectors
myNMDS = data.frame(MDS1 = Iso_wide_nmds$points[,1], MDS2 = Iso_wide_nmds$points[,2], 
                    MDS3 = Iso_wide_nmds$points[,3]) #originally had just 2

myvec.sp <- envfit(Iso_wide_nmds$points, df_wide_normlizedT, perm=1000)
myvec.sp.df <- as.data.frame(myvec.sp$vectors$arrows*sqrt(myvec.sp$vectors$r))
myvec.sp.df$species<-rownames(myvec.sp.df)

myNMDS2 <- myNMDS %>%
  rownames_to_column()



# Check stressplots, scree diagrams
Iso_wide_nmds$stress # Add a flag if this is high?
nmds.monte(df_wide_normlizedT, distance="euclidean", k=3, autotransform=FALSE, trymax=20)

#Samples in ordinate space -----------------------------------------------
Plot 2 dimensional NMDS configuration.
plot(Iso_wide_nmds$points, type="n") # plotting the scores(iso_wide_nmds)
text(Iso_wide_nmds,labels=row.names(Iso_wideT), cex = 1)
title(paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))

# See how particular compound changes with location
iso_df <- as.data.frame(df_wide_normlizedT)
plot(scores(Iso_wide_nmds), type = "p")
points(Iso_wide_nmds, cex = iso_df$Ectoine, col = "red")
title(paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))


# Vectors and ordiplots - messy due to scale differences
ggplot(data = myNMDS, aes(MDS1, MDS2)) +
  #geom_point() + # comment out for neatness
  geom_segment(data=myvec.sp.df, aes(x=0, xend=MDS1, y=0, yend=MDS2),
               #arrow = arrow(),
               colour="grey") +
  #geom_text(data=myNMDS2, aes(x=MDS1, y=MDS2, label=rowname), size=3) +
  geom_text(data=myvec.sp.df, aes(x=MDS1, y=MDS2, label=species), size=3) +
  ggtitle(paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))


# NMDS combined with clustering
mysol.d <- vegdist(df_wide_normlizedT, "euclidean")
mysitecl.ward <- hclust(mysol.d,method='ward.D')
mysitecl.class <- cutree(mysitecl.ward, k=7) # customize k groups here
groups <- levels(factor(mysitecl.class))
mysite.sc <- scores(myNMDS)

my.p <- ordiplot(mysite.sc, type="n",
                 main=paste("Incubation Experiments: Eddy", EddyInformation, sep = " "))
for (i in 1:length(groups))
{
  points(mysite.sc[mysitecl.class==i,], pch=(14+i), cex=2, col=i+1)
}
text(mysite.sc, row.names(Iso_wideT), pos=4, cex=0.7)
ordicluster(my.p, mysitecl.ward, col="dark grey")

legend("bottomleft", paste("Group", c(1:length(groups))),
       pch=14+c(1:length(groups)), col=1+c(1:length(groups)), pt.cex=2)

Ordiplot with hulls CURRENTLY NOT WORKING
ordiplot(Iso_wide_nmds, type="n")
ordihull(Iso_wide_nmds, groups=Treatment$Treatment.Status, draw="polygon",col="grey90",label=T)
orditorp(Iso_wide_nmds, display="sites",
         air=0.01, cex=0.5)

# # Plot convex hulls with colors based on treatment NOT WORKING FOR SAME REASON as above
colors=c(rep("red",5), rep("blue",5))
ordiplot(Iso_wide_nmds, type="n")
for(i in unique(Treatment$Treatment.Status)) {
  ordihull(Iso_wide_nmds$point[grep(i,Treatment$Treatment.Status),], draw="polygon",
           groups=Treatment$Treatment.Status[Treatment$Treatment.Status==i],
           col=colors[grep(i, Treatment$Treatment.Status)], label=T)}
#orditorp(Iso_wide_nmds, display="species", col="red", air=0.01)
orditorp(Iso_wide_nmds, display="sites", air=0.01, cex=1.25)

# Experiments with species scores, not very useful
envfit(Iso_wide_nmds, df_wide_normlizedT)
test <- envfit(df_wide_normlizedT ~ Iso_pointlocation$Supergroup, data = iso_df, perm=1000) #???

scores(Iso_wide_nmds) # this is the same as the $points call, species scores do not appear to exist