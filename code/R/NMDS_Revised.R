## NMDS plot revised

#load libraries ----
library(tidyr)
library(ggplot2)
library(vegan)
library(dplyr)
library(tibble)
library(ggnewscale)
library(ggrepel)
library(dendextend)
library(ggdendro)

#create a custom fill palette for the regions
fill_pal <- c("#00798c","#edae49","#d1495b")

#set a theme for plotting
nmdstheme <- theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=18))

#load in the new eDNA assignments and the cleaned taxonomy table

tax_table <- read.csv("data/Final_ASVTable_withTaxonomy_Filtered.csv")

comp_df <- tax_table%>%
            mutate(sp=gsub("_"," ",Species),
                   sp=gsub(" BOLD:.*", "", sp),
                   sp=gsub(" cf.","",sp),
                   sp=gsub(" CMC02","",sp),
                   sp=gsub("sp1","",sp),
                   sp=gsub(" sp.","",sp),
                   sp=gsub("\\s+sp$","",sp),
                   sp=trimws(sp),
                   sp=gsub("Amphitrite figulus","Neoamphitrite figulus",sp))%>%
              select(Species,sp)%>%
              rename(Species = 2, call = 1)

tax_df <- read.csv("data/taxonomy_comparison_formatted.csv")

#set up data for the nmds plot
site_grouping <- data.frame(site=colnames(tax_table[5:length(colnames(tax_table))]))%>%
                  mutate(site_simple = gsub("\\..*", "", site))%>%
                  mutate(group=case_when(site_simple %in% c("CHB","FAI","RSE") ~ "South",
                                         site_simple %in% c("LIN","NOR","ESB","MIR",
                                                            "CHT","ASB") ~ "CB",
                                         TRUE ~ "East"))

#data prep for nmds analyses
nmds_df <- tax_table %>% 
              select(5:length(tax_table))%>%
              mutate(Species = tax_table$Species)%>%
              pivot_longer(cols= -Species,
                           names_to = "site",
                           values_to = "value")%>%
              pivot_wider(names_from = Species,
                          values_from = value)%>%
              column_to_rownames('site')%>%
              mutate(tot = rowSums(.,na.rm=T))%>%
              filter(tot>0)%>% #filter stations with 0 assignments -- WRK2.2 didn't have any
              select(-tot)%>%
              decostand(method = "total")

#conduct nmds 
nmds_bray <- metaMDS(nmds_df,k=8,distance="bray")

pa_data <- nmds_df%>%
  mutate_if(is.numeric, ~1 * (. > 0))%>%
  vegdist(., method = "jaccard", binary = TRUE)

nmds_pa <- pa_data%>%metaMDS(.,k=8)


#extract nmds scores
data.scores <- as.data.frame(scores(nmds_bray,"sites"))%>%
  mutate(site=rownames(.),
         method="bray",
         stress=round(nmds_bray$stress,3))%>%
  rbind(.,
        as.data.frame(scores(nmds_pa,"sites"))%>%
          mutate(site=rownames(.),
                 method="pa",
                 stress=round(nmds_pa$stress,3)))%>%
  left_join(.,site_grouping)%>%
  select(NMDS1,NMDS2,site_simple,group,method,stress)

species.scores <- as.data.frame(scores(nmds_bray,"species"))%>%
                mutate(method="bray")%>%
                rbind(.,
                      as.data.frame(scores(nmds_pa,"species"))%>%
                        mutate(method="pa"))%>%
                rownames_to_column(var="call")%>%
                left_join(comp_df)%>%
                select(-call)

#calculate contex hull
hull_data <- data.scores%>%
  group_by(group,method)%>%
  slice(chull(x=NMDS1,y=NMDS2))

#extract the stress for each  
stresses <- hull_data%>%
            group_by(method)%>%
            summarise(stress=round(unique(stress),3))%>%
            ungroup()
            
#assemble plots
p1_pa <- ggplot()+
  geom_polygon(data=hull_data%>%filter(method=="pa"),aes(x=NMDS1,y=NMDS2,fill=group),alpha=0.30,col="grey30",lwd=0.2)+
  geom_point(data=data.scores%>%filter(method=="pa"),aes(x=NMDS1,y=NMDS2,shape=group,fill=group),col="black",size=3)+
  geom_text_repel(data=data.scores%>%filter(method=="pa"),aes(x=NMDS1,y=NMDS2, label=site_simple))+
  scale_shape_manual(values=c(21,22,23))+
  nmdstheme+
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.85),
        legend.title=element_blank(),
        legend.background = element_blank())+
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5, 
           label = paste("Stress =", stresses%>%filter(method=="pa")%>%pull(stress), "k = 10")) +
  labs(shape="",fill="")+
  scale_fill_manual(values=fill_pal)

p1_bray <- ggplot()+
  geom_polygon(data=hull_data%>%filter(method=="bray"),aes(x=NMDS1,y=NMDS2,fill=group),alpha=0.30,col="grey30",lwd=0.2)+
  geom_point(data=data.scores%>%filter(method=="bray"),aes(x=NMDS1,y=NMDS2,shape=group,fill=group),col="black",size=3)+
  geom_text_repel(data=data.scores%>%filter(method=="bray"),aes(x=NMDS1,y=NMDS2, label=site_simple))+
  scale_shape_manual(values=c(21,22,23))+
  nmdstheme+
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.9),
        legend.title=element_blank(),
        legend.background = element_blank())+
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5, 
           label = paste("Stress =", stresses%>%filter(method=="bray")%>%pull(stress), "k = 10")) +
  labs(shape="",fill="")+
  scale_fill_manual(values=fill_pal)

#save plots
ggsave("output/nmds_jaccard_updated.png",p1_pa,height=8,width=8,units="in",dpi=300)
ggsave("output/nmds_braycurtis_updated.png",p1_bray,height=8,width=8,units="in",dpi=300)

#now do the combined site analysis

#data prep for nmds - add the species counts per station 
nmds_group_df <- tax_table %>% 
                  select(5:length(tax_table))%>%
                  mutate(Species = tax_table$Species)%>%
                  pivot_longer(cols= -Species,
                               names_to = "site",
                               values_to = "value")%>%
                  left_join(site_grouping)%>%
                  group_by(site_simple,Species)%>%
                  summarise(value=sum(value))%>%
                  ungroup()%>%
                  rename(site=site_simple)%>%
                  pivot_wider(names_from = Species,
                              values_from = value)%>%
                  column_to_rownames('site')%>%
                  mutate(tot = rowSums(.,na.rm=T))%>%
                  filter(tot>0)%>% #filter stations with 0 assignments -- WRK2.2 didn't have any
                  select(-tot)%>%
                  decostand(method = "total")

#conduct nmds analyses
nmds_group_bray <- metaMDS(nmds_group_df,k=8,distance="bray")

pa_group_data <- nmds_group_df%>%
                    mutate_if(is.numeric, ~1 * (. > 0))%>%
                    vegdist(., method = "jaccard", binary = TRUE)

nmds_group_pa <- pa_group_data%>%metaMDS(.,k=8)

#if you wanted to play around with clustering

hc <- hclust(pa_data, method = "average")
cutree(hc,4)

dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)

cutree(pa_clusters,4)%>%
        data.frame()%>%
        rename(cluster=1)%>%
        rownames_to_column(var="site_simple")%>%
        select(site_simple,cluster)%>%
        left_join(.,site_grouping%>%distinct(site_simple,.keep_all = TRUE)%>%select(site_simple,group))%>%
        arrange(cluster,group)

plot(hc, hang = -1, main = "Hierarchical Clustering Dendrogram", sub = "", xlab = "")
rect.hclust(hc, k = 4, border = "red")

#extract scores from nmds analyses
data.scores.grouped <- as.data.frame(scores(nmds_group_bray,"sites"))%>%
                      mutate(site=rownames(.),
                             method="bray",
                             stress=round(nmds_group_bray$stress,3))%>%
                      rbind(.,
                            as.data.frame(scores(nmds_group_pa,"sites"))%>%
                              mutate(site=rownames(.),
                                     method="pa",
                                     stress=round(nmds_group_pa$stress,3)))%>%
                      rename(site_simple = site)%>%
                      left_join(.,site_grouping%>%distinct(site_simple,.keep_all=TRUE)%>%select(site_simple,group))%>%
                      select(NMDS1,NMDS2,site_simple,group,method,stress)%>%
                      rename(site=site_simple)

species.scores.grouped <- as.data.frame(scores(nmds_group_bray,"species"))%>%
                        mutate(method="bray")%>%
                        rbind(.,
                              as.data.frame(scores(nmds_group_pa,"species"))%>%
                                mutate(method="pa"))%>%
                        rownames_to_column(var="call")%>%
                        left_join(comp_df)%>%
                        select(-call)

#calculate convex hull
hull_group_data <- data.scores.grouped%>%
  group_by(group,method)%>%
  slice(chull(x=NMDS1,y=NMDS2))

#extract the stress for each  
stresses_grouped <- hull_group_data%>%
  group_by(method)%>%
  summarise(stress=round(unique(stress),3))%>%
  ungroup()

#assemble plots
p1_grouped_bray <- ggplot()+
  geom_polygon(data=hull_group_data%>%filter(method=="bray"),aes(x=NMDS1,y=NMDS2,fill=group),alpha=0.30,col="grey30",lwd=0.2)+
  geom_point(data=data.scores.grouped%>%filter(method=="bray"),aes(x=NMDS1,y=NMDS2,shape=group,fill=group),col="black",size=3)+
  geom_text_repel(data=data.scores.grouped%>%filter(method=="bray"),aes(x=NMDS1,y=NMDS2, label=site))+
  scale_shape_manual(values=c(21,22,23))+
  nmdstheme+
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        legend.title=element_blank(),
        legend.background = element_blank())+
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5, 
           label = paste("Stress =", stresses_grouped%>%filter(method=="bray")%>%pull(stress), "k = 10")) +
  labs(shape="",fill="")+
  coord_equal()+
  scale_fill_manual(values=fill_pal)

p1_grouped_pa <- ggplot()+
  geom_polygon(data=hull_group_data%>%filter(method=="pa"),aes(x=NMDS1,y=NMDS2,fill=group),alpha=0.30,col="grey30",lwd=0.2)+
  geom_point(data=data.scores.grouped%>%filter(method=="pa"),aes(x=NMDS1,y=NMDS2,shape=group,fill=group),col="black",size=3)+
  geom_text_repel(data=data.scores.grouped%>%filter(method=="pa"),aes(x=NMDS1,y=NMDS2, label=site))+
  scale_shape_manual(values=c(21,22,23))+
  nmdstheme+
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.1),
        legend.title=element_blank(),
        legend.background = element_blank())+
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5, 
           label = paste("Stress =", stresses_grouped%>%filter(method=="pa")%>%pull(stress), "k = 10")) +
  labs(shape="",fill="")+
  coord_equal()+
  scale_fill_manual(values=fill_pal)

#save plots
ggsave("output/nmds_grouped_jaccard_updated.png",p1_grouped_pa,height=8,width=8,units="in",dpi=300)
ggsave("output/nmds_grouped_braycurtis_updated.png",p1_grouped_bray,height=8,width=8,units="in",dpi=300)
