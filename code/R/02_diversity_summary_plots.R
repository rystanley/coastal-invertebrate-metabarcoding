#make some summary plots to showcase patterns in diversity amongst methods and geography

#load libraries ----
library(tidyverse)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)

#custom colour palette for plots
palette_cust <- c("#154360", "#FF5733", "#1ABC9C")

#load the species list data from each data and merge with the taxonomic information (cleaned in 01_data_processing.R)


#phyla plot ---------
comp_df <- read.csv("data/comparison_data_april2025.csv")

bar_data <- comp_df%>%
              mutate(Status = case_when(
                method == "both" ~ "Shared",
                method == "edna" ~ "eDNA Only",
                method == "traditional" ~ "Traditional Only",
                TRUE ~ NA))%>%
              group_by(Phylum,Status)%>%
              summarise(Count=n())%>%
              ungroup()%>%
              mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
              data.frame()



# bar_data <- comp_df%>%
#             group_by(aphiaID)%>%
#             mutate(Status = case_when(
#                   sum(method == "traditional") == 1 & sum(method == "edna") == 0 ~ "Traditional Only",
#                   sum(method == "traditional") == 0 & sum(method == "edna") == 1 ~ "eDNA Only",
#                   TRUE ~ "Shared"
#                 ))%>%
#             ungroup()%>%
#             distinct(aphiaID,.keep_all=TRUE)%>%
#             group_by(Phylum,Status)%>%
#             summarise(Count=n())%>%
#             ungroup()%>%
#             mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
#             data.frame()

#set up a eDNA centric plot order
phylum_order <- bar_data %>%
                group_by(Phylum) %>%
                summarise(eDNA_Only_Proportion = sum(Status == "eDNA Only") / n(),
                          eDNA_Only_count = sum(Status == "eDNA Only")) %>%
                arrange(desc(eDNA_Only_count)) %>%
                pull(Phylum)  # Extract the ordered phyla

#set up plotting data

count_levels <- bar_data%>%
                filter(Status=="eDNA Only")%>%
                arrange(-Count)%>%
                pull(Phylum)%>%
                c(.,setdiff(bar_data$Phylum,bar_data%>%
                              filter(Status=="eDNA Only")%>%
                              arrange(-Count)%>%
                              pull(Phylum)))

    df1 <- bar_data%>%
          mutate(Phylum = factor(Phylum,levels=count_levels),
                 Status = gsub("Only","only",Status),
                 Status = factor(Status,levels=c("Traditional only","Shared","eDNA only")))%>%
        group_by(Phylum)%>%
        mutate(Proportion = Count / sum(Count),
               TotalCount = sum(Count))%>%
        ungroup()
    
    #order by eDna only proprotions 
      prop_phyla_ord <- df1%>%
                        filter(Status == "eDNA only")%>%
                        arrange(-Proportion)%>%
                        pull(Phylum)%>%
                        as.character()%>%
                        c(.,setdiff(unique(bar_data$Phylum),.))
      
      #plot with the Phylum ordered according to the proportional fraction of eDNA only within each phyla
      df2 <- df1%>%
              mutate(Phylum = factor(Phylum,levels=prop_phyla_ord))%>%
              group_by(Phylum)%>%
              mutate(cumulative_proportion = cumsum(Proportion) - Proportion / 2)%>%
              ungroup()
      
      df2_labs <- df2%>%
        distinct(Phylum,.keep_all = TRUE)
      
      #Assemble data that is for all taxa to join to this plot. 
      df3 <- bar_data%>%
             group_by(Status)%>%
             summarise(Count=sum(Count))%>%
             ungroup()%>%
             mutate(Phylum = "All",
                    Proportion = Count/sum(Count),
                    Status = gsub("Only","only",Status),
                    Status = factor(Status,levels=c("Traditional only","Shared","eDNA only")))
      
      df4 <- df2%>%
             dplyr::select(names(df3))%>%
             rbind(.,df3)%>%
             mutate(Phylum = factor(Phylum,levels=c(levels(df2$Phylum),"All")),
                    Status = gsub(" only","",Status),
                    Status = factor(Status,levels=c("Traditional","Shared","eDNA")),
                    group = ifelse(Phylum == "All","All","Phyla"),
                    group = factor(group,levels=c("Phyla","All")))
        
    
#Stacked bar plot of counts
count_plot <- ggplot(df1, aes(x = Phylum, y = Count, fill = Status)) +
              geom_bar(stat = "identity", position = "stack",col="black") + 
              labs(x = "", y = "Total Count",
                 fill = "") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position = "inside",
                    legend.position.inside = c(0.9, 0.9),
                    legend.background = element_blank())+
              scale_fill_manual(values=palette_cust)+
              scale_y_continuous(expand = c(0.02,0));count_plot


#Proportional stacked bar plot 
prop_plot <- ggplot(df2, aes(x = Phylum, y = Proportion, fill = Status)) +
  # Stacked bar plot with proportions
  geom_bar(stat = "identity", position = "fill", col = "black") +  
  # Add text for counts within each segment of the stacked bar using `position_stack(vjust)`
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5), size = 4, fontface = 2) +
  # Add total count on top of each bar
  # geom_text(data = df2_labs, aes(x = Phylum, y = 1.05, label = TotalCount), 
  #           size = 4, fontface = 2) + 
  # Customize labels and theme
  labs(x = "", y = "Proportion", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.justification = "right",
        legend.direction = "horizontal") +
  scale_y_continuous(expand = c(0.02, 0),labels = scales::percent) +
  scale_fill_manual(values = palette_cust);prop_plot


#save plots
ggsave("output/diversity_bar_count.png",count_plot,width=6,height=4.5,units="in",dpi=300)
ggsave("output/diversity_bar_prop.png",prop_plot,width=6,height=4.5,units="in",dpi=300)

#do a proportional plot with the 'all' added in 

# Filter out the "all" category
df_phyla <- df4 %>% filter(Phylum != "All")
df_all <- df4 %>% filter(Phylum == "All")

# Count the number of phyla excluding 'all'
num_phyla <- n_distinct(df_phyla$Phylum)

# Plot for the phylum categories excluding 'all'
plot_phyla <- ggplot(df_phyla, aes(x = Phylum, y = Proportion, fill = Status)) +
  geom_bar(stat = "identity", position = "fill", col = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4, fontface = 2) +
  labs(x = "", y = "Proportion", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.justification = "right",
        legend.direction = "horizontal") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::percent) +
  scale_fill_manual(values = palette_cust)

# Plot for the 'all' category
plot_all <- ggplot(df_all, aes(x = Phylum, y = Proportion, fill = Status)) +
  geom_bar(stat = "identity", position = "fill", col = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4, fontface = 2) +
  labs(x = "", y = "Proportion", fill = "") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::percent) +
  scale_fill_manual(values = palette_cust)

# Combine the two plots with proportional widths
# num_phyla is used to ensure the widths are proportional
combined_plot <- (plot_phyla + plot_all) + 
                  plot_layout(ncol = 2, widths = c(num_phyla, 1), guides = "collect")&
                  theme(legend.position = "top", legend.justification = "right")

# Display the plot
combined_plot

ggsave("output/proportional_comparison_all.png",combined_plot,width=6,height=4.5,units="in",dpi=300)


#Lifestyle plot  -----

lifestyle_df <- read.csv("data/lifestyle_formatted.csv")

bar_data_ls <- lifestyle_df%>%
  mutate(Status = case_when(
    method == "both" ~ "Shared",
    method == "edna" ~ "eDNA Only",
    method == "traditional" ~ "Traditional Only",
    TRUE ~ NA))%>%
  group_by(lifestyle,Status)%>%
  summarise(Count=n())%>%
  ungroup()%>%
  mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
  data.frame()


#set up a eDNA centric plot order
plot_order_ls <- bar_data_ls %>%
  group_by(lifestyle) %>%
  summarise(eDNA_Only_Proportion = sum(Status == "eDNA Only") / n(),
            eDNA_Only_count = sum(Status == "eDNA Only")) %>%
  arrange(desc(eDNA_Only_count)) %>%
  pull(lifestyle)  # Extract the ordered phyla

#set up plotting data

count_levels_ls <- bar_data_ls%>%
  filter(Status=="eDNA Only")%>%
  arrange(-Count)%>%
  pull(lifestyle)%>%
  c(.,setdiff(bar_data_ls$lifestyle,bar_data_ls%>%
                filter(Status=="eDNA Only")%>%
                arrange(-Count)%>%
                pull(lifestyle)))

df1_ls <- bar_data_ls%>%
  mutate(lifestyle = factor(lifestyle,levels=count_levels_ls),
         Status = gsub("Only","only",Status),
         Status = factor(Status,levels=c("Traditional only","Shared","eDNA only")))%>%
  group_by(lifestyle)%>%
  mutate(Proportion = Count / sum(Count),
         TotalCount = sum(Count))%>%
  ungroup()

#order by eDna only proprotions 
prop_ls_ord <- df1_ls%>%
  filter(Status == "eDNA only")%>%
  arrange(-Proportion)%>%
  pull(lifestyle)%>%
  as.character()%>%
  c(.,setdiff(unique(bar_data_ls$lifestyle),.))

#plot with the Phylum ordered according to the proportional fraction of eDNA only within each phyla
df2_ls <- df1_ls%>%
  mutate(lifestyle = factor(lifestyle,levels=prop_ls_ord))%>%
  group_by(lifestyle)%>%
  mutate(cumulative_proportion = cumsum(Proportion) - Proportion / 2)%>%
  ungroup()

df2_labs_ls <- df2_ls%>%
  distinct(lifestyle,.keep_all = TRUE)

#Assemble data that is for all taxa to join to this plot. 
df3_ls <- bar_data_ls%>%
  group_by(Status)%>%
  summarise(Count=sum(Count))%>%
  ungroup()%>%
  mutate(lifestyle = "All",
         Proportion = Count/sum(Count),
         Status = gsub("Only","only",Status),
         Status = factor(Status,levels=c("Traditional only","Shared","eDNA only")))

df4_ls <- df2_ls%>%
  dplyr::select(names(df3_ls))%>%
  rbind(.,df3_ls)%>%
  mutate(lifestyle = factor(lifestyle,levels=c(levels(df2_ls$lifestyle),"All")),
         Status = gsub(" only","",Status),
         Status = factor(Status,levels=c("Traditional","Shared","eDNA")),
         group = ifelse(lifestyle == "All","All","lifestyle"),
         group = factor(group,levels=c("lifestyle","All")))



# Filter out the "all" category
df_lifestyle <- df4_ls %>% filter(lifestyle != "All")
df_all_ls <- df4_ls %>% filter(lifestyle == "All")

# Count the number of phyla excluding 'all'
num_ls <- n_distinct(df_lifestyle$lifestyle)

# Plot for the phylum categories excluding 'all'
plot_lifestyle <- ggplot(df_lifestyle, aes(x = lifestyle, y = Proportion, fill = Status)) +
  geom_bar(stat = "identity", position = "fill", col = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4, fontface = 2) +
  labs(x = "", y = "Proportion", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.justification = "right",
        legend.direction = "horizontal") +
  scale_y_continuous(expand = c(0.02, 0), labels = scales::percent) +
  scale_fill_manual(values = palette_cust)

ggsave("output/diversity_lifestyle_proportion.png",plot_lifestyle,width=6,height=4.5,units="in",dpi=300)


