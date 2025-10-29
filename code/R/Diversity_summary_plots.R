#make some summary plots to showcase patterns in diversity amongst methods and geography

#load libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)

#custome colour palette for plots
palette_cust <- c("#154360", "#FF5733", "#1ABC9C")

# #load the csv with jus the 'combined-reduced data'. 
# comp_data <- read.csv("data/comp_data_raw.csv")%>%
#              mutate(Phylum = ifelse(Phylum == "Nermertea","Nemertea",Phylum)) #spelling error
# 
# ## stacked bar plot ------
# 
# bar_data <- comp_data%>%
#             mutate(sp = case_when(
#               grepl("\\(", Species) ~ sub(" \\(.*\\)", "", Species),  # Remove text in parentheses
#               grepl(",", Species) ~ sub(",.*", "", Species),          # Remove text after a comma
#               TRUE ~ Species),
#               sp = trimws(sp),
#               sp =  sub("\\.$", "", sp),
#               Traditional = ifelse(Traditional =="X",TRUE,FALSE),
#               eDNA = ifelse(eDNA == "X",TRUE,FALSE))%>%
#             pivot_longer(cols = c("Traditional", "eDNA"), 
#                          names_to = "Method", 
#                          values_to = "Detected") %>%
#             filter(Detected == TRUE) %>%
#   group_by(sp) %>% # Create a status column to indicate shared or exclusive detections
#   mutate(Status = case_when(
#     sum(Method == "Traditional") == 1 & sum(Method == "eDNA") == 0 ~ "Traditional Only",
#     sum(Method == "Traditional") == 0 & sum(Method == "eDNA") == 1 ~ "eDNA Only",
#     TRUE ~ "Shared"
#   ))%>%
#   ungroup()%>%
#   group_by(Phylum, Method, Status) %>%# Create a count column for alluvial plotting
#   summarise(Count = n()) %>%
#   ungroup()%>%
#   group_by(Phylum,Status)%>%
#   distinct(Status,.keep_all = TRUE)%>%
#   ungroup()%>%
#   mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
#   arrange(Phylum,Status,Method)

bar_data <- read.csv("data/comparison_bar_df_formatted.csv")%>%
            mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))

comp_df <- read.csv("data/comparison_df_formatted.csv")

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
      df3 <- comp_df%>%
             gather(key="Method",value="value",Traditional:eDNA)%>%
             #mutate(value = ifelse(value=="X",TRUE,FALSE))%>%
             filter(value)%>%
             group_by(Species)%>%
             summarise(Status=paste0(unique(Method),collapse="-"))%>%
             ungroup()%>%
             mutate(Status = gsub("Traditional-eDNA","Shared",Status))%>%
             group_by(Status)%>%
             summarise(Count=n())%>%
             ungroup()%>%
             mutate(Phylum = "All",
                    Proportion = Count/sum(Count))
      
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

### now do the regional alluvial plots

alluv_df_raw <- read.csv("data/Combined_FurtherFilteredJuly2024.csv")

site_grouping <- data.frame(site=colnames(alluv_df_raw[5:length(colnames(alluv_df_raw))]))%>%
  mutate(site_simple = gsub("\\..*", "", site))%>%
  mutate(group=case_when(grepl("CHB",site) ~ "South",
                         site_simple %in% c("LIN","NOR","ESB","MIR",
                                            "CHT","ASB") ~ "CB",
                         TRUE ~ "East"))

alluv_df <- alluv_df_raw%>%
            gather(key="site",value,5:length(alluv_df_raw))%>%
            mutate(site = gsub("\\..*", "", site))%>%
            group_by(site,Species)%>%
            summarise(count=sum(value))%>%
            ungroup()%>%
            mutate(pa = ifelse(count != 0,1,0))%>%
            left_join(.,alluv_df_raw%>%dplyr::select(1:4))%>%
            left_join(.,site_grouping%>%
                        dplyr::select(-site)%>%
                        rename(site=site_simple)%>%
                        distinct(site,.keep_all=TRUE))

ggplot(alluv_df%>%filter(pa==1),
       aes(axis1 = Phylum, axis2 = Class, axis3 = group, y = pa)) +
  geom_alluvium(aes(fill = group)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phylum", "Class", "Group"), expand = c(0.1, 0.1)) +
  labs(x = "", y = "Count of Species Presence",
       title = "Alluvial Plot: Phylum to Class to Geographic Grouping") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(alluv_df%>%filter(pa==1),
       aes(axis1 = Phylum, axis2 = group, y = pa)) +
  geom_alluvium(aes(fill = group)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phylum", "Group"), expand = c(0.1, 0.1)) +
  labs(x = "", y = "Count of Species Presence",
       title = "Alluvial Plot: Phylum to Class to Geographic Grouping") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# #Alluvial plot ----
# #clean up the data so that it is just species names. 
# # Convert data to long format to show detections by method
comp_data_df <- comp_data%>%
  mutate(sp = case_when(
    grepl("\\(", Species) ~ sub(" \\(.*\\)", "", Species),  # Remove text in parentheses
    grepl(",", Species) ~ sub(",.*", "", Species),          # Remove text after a comma
    TRUE ~ Species),
    sp = trimws(sp),
    sp =  sub("\\.$", "", sp),
    Traditional = ifelse(Traditional =="X",TRUE,FALSE),
    eDNA = ifelse(eDNA == "X",TRUE,FALSE))%>%
  pivot_longer(cols = c("Traditional", "eDNA"),
               names_to = "Method",
               values_to = "Detected") %>%
  filter(Detected == TRUE) %>%
  group_by(sp) %>% # Create a status column to indicate shared or exclusive detections
  mutate(Status = case_when(
    sum(Method == "Traditional") == 1 & sum(Method == "eDNA") == 0 ~ "Traditional Only",
    sum(Method == "Traditional") == 0 & sum(Method == "eDNA") == 1 ~ "eDNA Only",
    TRUE ~ "Shared"
  )) %>%
  ungroup()%>%
  group_by(Phylum, Method, Status) %>%# Create a count column for alluvial plotting
  summarise(Count = n()) %>%
  ungroup()%>%
  mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")),
         Phylum = factor(Phylum,levels = phylum_order))%>%
  arrange(Phylum,Status,Method)

# Create the alluvial plot with separate legends
ggplot(comp_data_df, aes(axis1 = Phylum, axis2 = Status, y = Count)) +
  # Alluvium lines by Status
  geom_stratum(aes(fill = Status), width = 0.2) +
  geom_alluvium(aes(fill = Status), width = 0.2) +
  scale_fill_manual(name = "Status", values = palette_cust) + # Custom colors for Status
  ggnewscale::new_scale_fill() +  # Allows for a new fill scale
  # Strata fill by Phylum
  #geom_stratum(aes(fill = Phylum), width = 0.2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(name = "Phylum", values = c("Phylum1" = "#d62728",
                                                "Phylum2" = "#9467bd",
                                                "Phylum3" = "#8c564b")) + # Custom colors for Phylum
  scale_x_discrete(limits = c("Phylum", "Method")) +
  labs(title = "Species Detection Across Methods",
       x = "Categories", y = "Species Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

