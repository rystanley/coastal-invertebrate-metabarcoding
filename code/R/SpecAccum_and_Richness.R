
# NMDS and Species Diversity Plots Script ---------------------------------

## using the ASV table to create a matrix of sites by species to conduct NMDS, spec.accum curves, and plots of species/taxa found


# Load libraries ----------------------------------------------------------

library(vegan)
library(dplyr)
library(janitor)
library(ggplot2)
library(ggvegan)
library(vroom)
library(ggrepel)
library(RColorBrewer)
library(tibble)
library(stringr)
library(tidyr)
library(ggsci)


# Load the data -----------------------------------------------------------


load("data/Trask_CommunityStats.RData")

# Or start here from scratch

taxtable <- read.csv("data/Final_ASVTable_withTaxonomy_Filtered.csv", header = T)
head(taxtable) #use "head" to look at our data briefly 
#select some columns for this specific analysis 

#Update December 10 2024 - removed RSE.3 from this due to suspected contamination 
commat <- taxtable %>% select(CAB.1, CAB.2, CAB.3, CON.1, CON.2, CON.3, FAI.1, FAI.2, FAI.3, FRK.1, FRK.2, FRK.3, MOS.1, MOS.2, MOS.3, PLS.1, PLS.2, PLS.3, RSE.1, RSE.2, TAE.1, TAE.2, TAE.3, TAW.1, TAW.2, TAW.3, WRK02.1, WRK02.3, ASB.1, ASB.2, ASB.3, CHB.1, CHB.2, CHB.3, LIN.1, LIN.2, LIN.3, MIR.1, MIR.2, MIR.3, CHT.1, CHT.2, CHT.3, NOR.1, NOR.2, NOR.3, ESB.1, ESB.2, ESB.3)

groups<- c(rep("East",6), rep("South",3), rep("East", 9), rep("South",2), rep("East",8), rep("CB", 3),rep("South",3),rep("CB",15))


# Create stacked barplot of taxa by site ----------------------------------

tax_long <- taxtable %>% 
  gather(Site_Rep, Count,-c(Phylum, Class, Species, Probability)) %>% 
  separate(Site_Rep, into=c("Site", "Replicate"), sep="\\.")

spec_summary<- tax_long %>%
  group_by(Site, Phylum) %>% 
  summarise(TotalCount =sum(Count)) %>%
  ungroup()

#could probably combine this with the above piping 
spec_summary <- spec_summary %>% 
  group_by(Site) %>%
  mutate(RelativeCount = TotalCount/sum(TotalCount)) %>%
  ungroup() %>%
  mutate(Site=str_replace(Site, "WRK02","WRK"))

level_order <- c("CHB","RSE","FAI","FRK","WRK","PLS","CON","CAB","TAE","TAW","MOS","ESB","MIR","LIN","NOR","CHT","ASB")
phylum_order <- c("Porifera", "Cnidaria","Platyhelminthes", "Annelida", "Mollusca", "Nemertea", "Bryozoa", "Arthropoda", "Echinodermata", "Hemichordata", "Chordata")

#colour Site names by their colours in the map
axis_pal <- c(rep("#d1495b",3), rep("#edae49",8), rep("#00798c",6))
colourCount = length(unique(spec_summary$Phylum))
#getPalette = colorRampPalette(brewer.pal(9, "Set1"))
new_palette = pal_jco()(10)
my_palette <- c("#0072B2", "#D55E00", "#CC79A7","#009E73",  "#F0E442", "#9467BD", "#8E382C", "#E69F00", "#3182BD", "#766060")


p1 <- ggplot(spec_summary, aes(x = factor(Site, level=level_order), y = RelativeCount, fill = factor(Phylum, level = phylum_order))) +
  geom_bar(stat = "identity", colour="black") +
  scale_fill_manual(values=my_palette) +
  labs(x = "Site",
       y = "Relative Count") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(text=element_text(size=16),
        axis.text.x=element_text(colour=axis_pal,face="bold"))

p1 + guides(fill=guide_legend(title="Phylum"))

ggsave(filename = "Site_ClassBarPlots_AprilUpdate.png", plot = last_plot(),device = "png", path = "output/", width = 12, height = 8, units = "in", dpi = 300)


# Non-metric Multidimensional Scaling -------------------------------------

#Transpose the table for vegan
commat2<-t(commat)
colnames(commat2)<-taxtable[,3] #this selects the 3rd column from our taxon table and inserts the species IDs into this matrix

coi.nmds <-metaMDS(commat2, distance="jaccard", k=14, trymax = 300, maxit=500)
plot(coi.nmds) #this is not very informative without labels!



##plot NMDS with ggplot to look nicer
data.scores <- as.data.frame(scores(coi.nmds,"sites")) %>%
 mutate(ID=rownames(.), group=groups)
species.scores <- as.data.frame(scores(coi.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data <- data.scores %>%
  as.data.frame() %>%
  group_by(group) %>%
  slice(chull(x=NMDS1,y=NMDS2))

nmdstheme <- theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=18))


  # regular ggplot
  p1 <- ggplot() +
    #geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,color=group, fill=group),alpha=0.20) + #with filled polygons
    geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=group),color="black",alpha=0.8, linewidth=1.2) + # add the hulls
    # add the hulls
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels or remove   this if too messy
    geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=group,colour=group),size=3) +
    #scale_colour_manual()
    geom_text_repel(data=data.scores,aes(x=NMDS1,y=NMDS2, label=ID))+
    #scale_color_gradient(low="blue", high="black")+
    scale_colour_brewer(palette = "Accent") +
    coord_equal()+
    #scale_shape_manual(values=c(15,7,18,16,10,8,17))+
    geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(coi.nmds$stress,3),"k =",coi.nmds$ndim)))+
    nmdstheme
  p1

ggsave(filename = "COI_NMDS_3groups.png", plot = p1, device = "png", path = "output/", width = 12, height = 10, units = "in", dpi = 300)



# Species accumulation and rarefaction curves -----------------------------

gg<- specaccum(commat2, method="random", ci.type="polygon")
gg<-data.frame(Sites=gg$sites, Richness=gg$richness, sd=gg$sd)

ggplot(gg,aes(Sites, Richness))+
  geom_line(linewidth=1.5,color="#00798c")+
 geom_linerange(aes(x=Sites, ymin=Richness -2*sd, ymax=Richness + 2*sd),colour="#00798c")+
  labs(x="Environmental sample")+
  theme_bw()+
  theme(text=element_text(size=16))

ggsave("SpecAccum.png",plot=last_plot(), device = "png", path = "output/", height=6,width=8,units="in",dpi=320)
#rarefaction
coi.rarefaction <-rarecurve(commat2, step=10)


#this gives us Shannon index values per sample 
shan <- data.frame(diversity(commat2, index="shannon"))
shan
#we can aggregate these so we get a site level index
mean.shan <- aggregate(. ~ substr(rownames(shan), 1, 3), shan, mean, na.rm = TRUE)
colnames(mean.shan)<-c("Code","ShannonDiversity")
#get species richness as n species too
rich_df<-as.data.frame(commat2)

df_rich_rep <- rich_df %>%
  rownames_to_column(var="site") %>%
  mutate(Code=substr(site, start=1, stop=3)) %>%
  group_by(Code) %>%
  summarize(avg_richness = mean(rowSums(across(where(is.numeric))>0))) %>%
  ungroup() %>%
  data.frame()



# Regression of diversity metrics with latitude ---------------------------

coords <- read.csv("data/Map.csv")
lat.rich <- full_join(mean.shan,coords,by="Code")
lat.rich.2 <- left_join(lat.rich, df_rich_rep, by="Code") %>%
  pivot_longer(cols=c(ShannonDiversity, avg_richness), names_to="richness",values_to = "value") %>% 
  mutate(richness=str_replace(richness, "avg_richness","Taxonomic Richness")) %>%
  mutate(richness=str_replace(richness,"ShannonDiversity","Shannon Diversity"))
#run the linear regression model for stats
summary(lm(ShannonDiversity~Latitude, data=lat.rich))
summary(lm(value~Latitude, data=lat.rich.2 %>% filter(richness=="Taxonomic Richness") %>% data.frame))


ggplot(data=lat.rich.2,aes(x=Latitude,y=value))+
  geom_point(pch=21, size=3, colour="black", fill="black")+
  facet_wrap(~richness, scales="free_y")+
  geom_smooth(method = "lm",colour="firebrick", se=TRUE, linewidth=0.75)+
  labs(x="Latitude °N", y="Diversity value")+
  theme_bw()+
  theme(text=element_text(size=16),strip.background = element_rect(fill="white"))

ggsave(filename = "SpeciesRichness_by_Latitude.png",plot = last_plot(),device = "png",path = "output/", width = 12, height=8, dpi = 400,units = "in")



# Save our data -----------------------------------------------------------
save.image("data/NSeDNA_CommunityStats.RData")
