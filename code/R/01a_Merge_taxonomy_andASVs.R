#Merging taxonomy and read count tables, and filtering taxonomy
library(dplyr)

#load in our two data tables - these paths will only work if you open the Courtney-Trask github R project, or set the working directory to your local github folder
tax <- read.table(file = "data/2023_PostBioinformaticsData/Trask_RDPv5.tsv", header = F, sep = "\t") #here the \t is short for tab separation
tax<- rename(tax, OTU.ID= V1) #this is just so there is a matching column between our taxon and asvs files so we can merge them easier (OTU.ID)
                
asvs <- read.table(file="data/2023_PostBioinformaticsData/Trask_feature_table_export.tsv", header = F, sep="\t")
              
taxatable <- full_join(x= asvs, y= tax, by = "V1")  #this will reorganize the new table by OTU alphabetically

colnames(taxatable)
#Clean up this table and remove columns we don't need - here we are 'selecting' the columns from "taxatable" to keep, AND filtering by 0.90 probability of species being kept, AND removing rows with very few ASV counts

taxtable.final <- taxatable %>% filter(V14.y %in% c("Actinopteri", "Anthozoa", "Ascidiacea", "Asteroidea", "Bivalvia", "Branchiopoda", "Calcarea", "Cephalopoda", "Cestoda", "Chondrichthyes", "Echnioidea", "Gastropoda","Gymnolaemata", "Hexactinellida", "Holothuroidea", "Hydrozoa" ,"Malacostraca", "Mammalia", "Ophiuroidea" ,"Ostracoda", "Palaeonemertea", "Pilidiophora", "Pinopsida", "Polychaeta","Polyplacophora", "Pycnogonida", "Scyphozoa", "Staurozoa", "Tentaculata","Thecostraca") & V28.y > 0.90) %>% as.data.frame()

#now we will write this new table to the Github repo
write.table(x = taxtable.final, file = "data/Full_COI_TaxonomyPerSite_2023.txt", quote = F)




# Now load and filter the 2019 data ---------------------------------------


#Read in the csv file - csv file is just a "comma separated" file which is one of the types of data R likes
asvs <- read.csv(file = "data/COI_TaxonomyPerSite_Coastal2019_modified.csv",header = TRUE)
#look at the data
head(asvs)
#look at the dimensions of our data
dim(asvs) #9810 rows and 52 columns

#let's filter out human and insect data first
asvs2 <- asvs %>% filter(Class!="Insecta" & Probability >0.90) 
dim(asvs2) #now by filtering by 90% probability we're down to 1045 rows and 52 columns

#Now let's filter by rows with all zeroes - this may be a bit complicated looking but basically in our 
#asvs2 dataset columns 1 through 5 are non-numeric data, so we start at column 6 in the code below and look at all
#rows in columns 6 through 52 and see if the sum of all values is >0
asvs3 <- asvs2[rowSums(asvs2[,6:length(colnames(asvs2))])>0,]
dim(asvs3) #now we're down to 676 rows and 52 columns

#Let's look at the unique Class names now and see if we want to filter out any other weird species
unique(asvs3$Class)
#We have 43 classes including some diatoms, algage, and vertebrates (cats, dogs, Canada geese) so we could filter those
#we also have class Arachnida in the data, however these are all marine water mites so we should keep them!
asvs4 <- asvs3 %>% filter(!Class %in% c("Aves", "Mammalia", "undef_Discosea","undef_Evosea", "undef_Oomycota"))
dim(asvs4) #now down to 668 by 52 dataframe!

write.csv(x = asvs, file = "data/COI_TaxonomyPerSite_Coastal2019_SUBSET.csv", quote = FALSE)
