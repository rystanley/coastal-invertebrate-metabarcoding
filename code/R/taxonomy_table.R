## Code to create a summary table for the paper

#load libraries
library(dplyr)
library(tidyr)
library(taxize)


#load the data 
comp_df <- read.csv("data/comparison_df_formatted.csv")%>%
  mutate(sp=gsub("sp1","",Species),
         sp=gsub("\\s+sp$","",sp),
         sp=trimws(sp))

#Run a classification on that data - output from this step is saved, so unless the dataset is changed then there is no need to re-run this analysis

  #Note also that some manual interventions are required when multiple matches are returned from WORMS using the classification function
  #In this case, the 'accepted' closest matching search result was retained (in all cases either option 1 or 2)

  # edna_taxonomy <- classification(comp_df$sp, db = 'worms',return_id = TRUE) #this can take a while note 'decisions' below
  
  # #save intermediate output 
  # save(edna_taxonomy,file="data/edna_taxonomy.RData")

load("edna_taxonomy.RData")

edna_flat <- do.call("rbind", edna_taxonomy)%>%
  mutate(call = gsub("\\..*","",rownames(.))) #this creates a column with what itis recieved as input

PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")      

edna_format <- edna_flat%>%
               dplyr::select(name,rank,call)%>%
               group_by(call)%>%
               spread(rank,name)%>%
               ungroup()%>%
               dplyr::select(call,all_of(PhyloNames))%>%
               left_join(.,edna_flat%>%
                           group_by(call)%>%
                           filter(row_number()==n())%>% #only want the AphiaID for the lowest classification
                           select(call,id))%>%
               rename(AphiaID = id)

#create table for publication 
edna_table <- edna_format%>%
              rename(sp=call)%>%
              left_join(.,comp_df%>%select(sp,Traditional,eDNA))%>%
              rowwise()%>%
              mutate(Detection =case_when(Traditional & !eDNA ~ "Traditional",
                                          !Traditional & eDNA ~ "eDNA",
                                          Traditional & eDNA ~ "Shared",
                                          TRUE ~ "something went wrong"))%>%
              select(c(AphiaID,all_of(PhyloNames),Detection))%>%
              mutate(Detection = factor(Detection,levels=c("eDNA","Shared","Traditional")))%>%
              arrange(Phylum,Detection)
  
#save the final table
write.csv(edna_table,file="output/taxonomy_comparison_formatted.csv",row.names=FALSE)
