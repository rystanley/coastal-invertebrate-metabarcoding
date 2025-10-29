## In the combined dataset that was compiled some species where removed from the eDNA account. These
#were captured in a re-assessment of the eDNA data in the 'Combined_FurtherFilteredJuly2024.csv'. To reconcile we add
# the missing species as Traditional = FALSE. That data is provided below

#libraries
library(dplyr)

edna_df <- read.csv("data/Combined_FurtherFilteredJuly2024.csv")%>%
            mutate(sp=gsub("_"," ",Species),
                            sp=gsub(" BOLD:.*", "", sp),
                            sp=gsub(" cf.","",sp),
                            sp=gsub(" CMC02","",sp),
                            sp = case_when(sp == "Amphitrite figulus" ~ "Neoamphitrite figulus", #neo is the accepted version https://www.marinespecies.org/aphia.php?p=taxdetails&id=335363
                                           TRUE ~ sp))

edna_sp <- edna_df%>%pull(sp)%>%unique()

comp_df <- read.csv("data/comp_data_raw.csv")%>%
  mutate(Phylum = ifelse(Phylum == "Nermertea","Nemertea",Phylum))%>%#spelling error
  mutate(sp = case_when(
    grepl("\\(", Species) ~ sub(" \\(.*\\)", "", Species),  # Remove text in parentheses
    grepl(",", Species) ~ sub(",.*", "", Species),          # Remove text after a comma
    TRUE ~ Species),
    sp = trimws(sp),
    sp =  sub("\\.$", "", sp),
    sp = gsub(",","",sp),
    Traditional = ifelse(Traditional =="X",TRUE,FALSE),
    eDNA = ifelse(eDNA == "X",TRUE,FALSE))%>%
  filter(eDNA,
         sp != "Polycirrus dubius") # this is not in the eDNA dataset generated after the last clean based on a thershold of eDNA detections so it is removed.  
  
comp_sp <- comp_df%>%pull(sp)%>%unique()

#identify the eDNA species that were trimmed from the original comparison list and add them in. 
missing_species <- data.frame(sp = setdiff(edna_sp,comp_sp))%>%
                   left_join(.,edna_df%>%select(Phylum,Class,sp))%>%
                   rename(Species=sp)%>%
                   select(Phylum,Class,Species)%>%
                   filter(Species != "Clitellio sp.")%>% #in this case there is a species "Clitellio sp." and that family is not identified by traditional methods
                   mutate(Traditional=FALSE,eDNA=TRUE)


#this is formatted using the code for bar_data in Diversity_summary_plots.R
combined_comparison <- read.csv("data/comp_data_raw.csv")%>%
                       mutate(Phylum = ifelse(Phylum == "Nermertea","Nemertea",Phylum))%>%#spelling error
                       mutate(sp = case_when(
                       grepl("\\(", Species) ~ sub(" \\(.*\\)", "", Species),  # Remove text in parentheses
                       grepl(",", Species) ~ sub(",.*", "", Species),          # Remove text after a comma
                       TRUE ~ Species),
                       sp = trimws(sp),
                       sp =  sub("\\.$", "", sp),
                       sp = gsub(",","",sp),
                       Traditional = ifelse(Traditional =="X",TRUE,FALSE),
                       eDNA = ifelse(eDNA == "X",TRUE,FALSE),
                       eDNA = ifelse(sp == "Polycirrus dubius",FALSE,eDNA))%>% #note in the comparison dataset this is 'TRUE' for eDNA but this species is absent from the curated eDNA list, likely due to low yeilds
                       select(-Species)%>%
                       rename(Class=Subphylum.Class,Species=sp)%>%
                       select(names(missing_species))%>%
                       rbind(.,missing_species)

combined_comparison_bar <- combined_comparison%>%
                      pivot_longer(cols = c("Traditional", "eDNA"), 
                                   names_to = "Method", 
                                   values_to = "Detected") %>%
                      filter(Detected == TRUE) %>%
                      group_by(Species) %>% # Create a status column to indicate shared or exclusive detections
                      mutate(Status = case_when(
                        sum(Method == "Traditional") == 1 & sum(Method == "eDNA") == 0 ~ "Traditional Only",
                        sum(Method == "Traditional") == 0 & sum(Method == "eDNA") == 1 ~ "eDNA Only",
                        TRUE ~ "Shared"
                      ))%>%
                      ungroup()%>%
                      group_by(Phylum, Method, Status) %>%# Create a count column for alluvial plotting
                      summarise(Count = n()) %>%
                      ungroup()%>%
                      group_by(Phylum,Status)%>%
                      distinct(Status,.keep_all = TRUE)%>%
                      ungroup()%>%
                      mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
                      arrange(Phylum,Status,Method)

write.csv(combined_comparison_bar,"data/comparison_bar_df_formatted.csv",row.names=FALSE)
write.csv(combined_comparison,"data/comparison_df_formatted.csv",row.names=FALSE)

