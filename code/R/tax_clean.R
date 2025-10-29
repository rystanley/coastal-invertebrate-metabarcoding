#Code to assess which species have refernce databases

https://www.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=2

#load libraries
library(dplyr)
library(bold)
library(tidyr)
library(tibble)
library(taxize)
library(RCurl)
library(worrms)

#load function for extracting taxonomic information
script <- RCurl::getURL("https://raw.githubusercontent.com/rystanley/offshore_edna/main/code/ClassifyFunction.R", ssl.verifypeer = FALSE)
eval(parse(text = script),envir=.GlobalEnv)
rm(script)  

capitalize_first <- function(text) { #from ChatGPT!
  words <- strsplit(text, " ")[[1]] 
  words[[1]] <- paste0(toupper(substr(words[[1]], 1, 1)), substr(words[[1]], 2, nchar(words[[1]])))
  capitalized_text <- paste(words, collapse = " ")  # Recombine the words into a single string
  return(capitalized_text)
}


#load taxon detected
load("data/taxtable.RData")

species <- taxtable%>%
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
  
PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

tax_list <- data.frame()

for(i in 1:length(species)){
  
  sp=species[i]
  
  message(paste0("Working on ",sp," ",i," of ",length(species)))
  
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

#the species name for "stylodrilus heringianus" doesn't work but the aphiaID does. THis will append this to the end. 
temp4 <- classification(1040258,db="worms")

temp4 <- temp4[[1]]%>%
         data.frame()%>%
         select(rank,name)%>%
         spread(rank,name)%>%
         mutate(aphiaID = temp[[1]]%>%data.frame()%>%slice(n())%>%pull(id))

temp4 <- temp4%>%
         select(all_of(c(names(temp2)[names(temp2)%in%PhyloNames],"aphiaID")))%>%
         mutate(species_filter = "stylodrilus heringianus")

tax_list <- rbind(tax_list,temp4)%>%
            select(all_of(c(PhyloNames,"aphiaID","species_filter")))

save(tax_list,file="data/tax_list.Rdata")

#now do bold extraction

tax_bold <- tax_list%>%
            rowwise()%>%
            mutate(sp = capitalize_first(species_filter))%>% #because bold has a insanely narrow search where capitalization matters. 
            data.frame()%>%
            pull(sp)

bold_out <- NULL
for(i in 1:length(tax_bold)){
  
  sp <- tax_bold[i]
  
  message(paste0("Working on '",sp,"' ",i," of ",length(tax_bold)))
  
  #these are the bits of info you want from BOLD
  info_cols <- c("input","taxid","taxon","tax_rank","tax_division","parentid","parentname","specimenrecords")
  
  #check to see if there is records on bold for that species
  temp <- bold_tax_name(sp)
  
  #fill in NULL if no recoreds
  if(is.na(temp$taxid)){
    
    temp <- data.frame(matrix(NA, nrow = 1, ncol = length(info_cols)))
    names(temp) <- info_cols
    temp$input  <- sp
    
  }
  
  #select only target info -- some species have extra image metadata, which is not needed. 
  temp2 <- temp%>%select(input,taxid,taxon,tax_rank,tax_division,parentid,parentname,specimenrecords)
  
  #bind it together
  bold_out <- rbind(bold_out,temp2)
  
}

save(bold_out,file="data/bold_id.RData")

#combine the data together
tax_out <- tax_list%>%
           rowwise()%>%
           mutate(input=capitalize_first(species_filter))%>%
           data.frame()%>%
           left_join(.,bold_out)%>%
           mutate(is_bold = ifelse(is.na(taxid),FALSE,TRUE))

save(tax_out,file="data/tax_process.RData")

