### this is code for doing data prep for latter plotting

#load libraries
library(tidyverse)
library(taxize)
library(RCurl)
library(worrms)
library(bold)

#taxonomy levels for the final table
PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#Clean the traditional sampling data -----------
#load the final species list from the traditional sampling compilation
trad_df <- read.csv("data/traditional_sampling_species_list.csv")%>%
           mutate(Phylum = ifelse(Phylum == "Nermertea","Nemertea",Phylum),
                  sp = case_when(
                    grepl("\\(", Species) ~ sub(" \\(.*\\)", "", Species),  # Remove text in parentheses
                    grepl(",", Species) ~ sub(",.*", "", Species),          # Remove text after a comma
                    TRUE ~ Species),
                  sp = trimws(sp),
                  sp =  sub("\\.$", "", sp),
                  sp = tolower(sp),
                  sp = gsub(" sp","",sp),
                  sp = ifelse(sp == "circeisirillum","circeis spirillum",sp), #these are errors in the formatting that need to be fixed manually
                  sp = ifelse(sp == "dipolydora1","dipolydora",sp),
                  sp = ifelse(sp == "spirorbisirorbis","spirorbis spirorbis",sp)) 

trad_species <- trad_df$sp

tax_list_trad <- data.frame()

for(i in 1:length(trad_species)){
  
  sp <- trad_species[i]
  
  message(paste0("Working on ",sp," ",i," of ",length(trad_species)))
  
  temp <- classification(sp,db="worms") #run classification
  
  temp2 <- temp[[1]]%>% #unpack classification
    data.frame()
  
  if(nrow(temp2[1])>1){
    
    temp2 <- temp2%>%
      select(rank,name)%>%
      spread(rank,name)%>%
      mutate(aphiaID = temp[[1]]%>%data.frame()%>%slice(n())%>%pull(id))
    
    temp3 <- temp2%>% #trim classification
      select(all_of(c(names(temp2)[names(temp2)%in%PhyloNames],"aphiaID")))
    
    #if data is missing or a certain taxonomic level isn't identified.
    missing_cols <- setdiff(c(PhyloNames,"aphiaID"),names(temp3))
    
    if(length(missing_cols)>0){
      
      temp3[,missing_cols] <- NA
      
      temp3 <- temp3[,c(PhyloNames,"aphiaID")]
      
    }
  }else{temp3 <- data.frame(matrix(NA, nrow = 1, ncol = length(PhyloNames)))
  names(temp3) = PhyloNames
  temp3$aphiaID = NA}
  
  temp3$species_filter <- sp # this is for linking back in the original code
  
  tax_list_trad <- rbind(tax_list_trad,temp3)
  
}
          
trad_formatted <- tax_list_trad%>%
                  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                  mutate(method="traditional")

## clean the eDNA data  ----------------
#load the final (merged) species list from the eDNA sampling
edna_df <- read.csv("data/Merged_Species_Detections.csv")

edna_species <- edna_df%>%
                mutate(species = tolower(gsub("_"," ",Species)),
                       species = case_when(species == "chamaedrilus sphagnetorum" ~ "cognettia sphagnetorum", # chamaedrilus is a invalid
                                           species == "lepidonotus squamatus cmc02" ~ "lepidonotus squamatus",
                                           #species == "stylodrilus heringianus" ~ 1040258, #note that this is a non-marine species
                                           species == "clitellio sp. bold:aag8543" ~ "clitellio",
                                           species == "boonea cf. bisuturalis bold:abw1340" ~ "boonea bisuturalis",
                                           TRUE ~ species))%>%
                distinct(species)%>%
                filter(species != "stylodrilus heringianus")%>% #pull this out using aphiaID later
                pull(species)

tax_list <- data.frame()

for(i in 1:length(edna_species)){
  
  sp <- edna_species[i]
  
  message(paste0("Working on ",sp," ",i," of ",length(edna_species)))
  
  temp <- classification(sp,db="worms") #run classification
  
  temp2 <- temp[[1]]%>% #unpack classification
    data.frame()
  
  if(nrow(temp2[1])>1){
    
    temp2 <- temp2%>%
      select(rank,name)%>%
      spread(rank,name)%>%
      mutate(aphiaID = temp[[1]]%>%data.frame()%>%slice(n())%>%pull(id))
    
    temp3 <- temp2%>% #trim classification
      select(all_of(c(names(temp2)[names(temp2)%in%PhyloNames],"aphiaID")))
    
    #if data is missing or a certain taxonomic level isn't identified.
    missing_cols <- setdiff(c(PhyloNames,"aphiaID"),names(temp3))
    
    if(length(missing_cols)>0){
      
      temp3[,missing_cols] <- NA
      
      temp3 <- temp3[,c(PhyloNames,"aphiaID")]
      
    }
  }else{temp3 <- data.frame(matrix(NA, nrow = 1, ncol = length(PhyloNames)))
  names(temp3) = PhyloNames
  temp3$aphiaID = NA}
  
  temp3$species_filter <- sp # this is for linking back in the original code
  
  tax_list <- rbind(tax_list,temp3)
  
}

#final formatted species list from eDNA
edna_formatted <- tax_list%>%
                  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                  mutate(method="edna")


#now output the cleaned taxonomy for later analysis
comp_df <- rbind(trad_formatted,edna_formatted)

#final cleaned output
save(comp_df,file="data/taxonomy_cleaned.RData")


