## Make a quick map of sampling sites 

#1. Load libraries - make sure these are installed locally 
library(sf)
library(dplyr) #
library(ggplot2)
library(rnaturalearth) #if this doesn't install on its own use devtools::install_github("ropensci/rnaturalearthhires")
library(patchwork)

#create a custom fill palette for the regions
fill_pal <- c("#d1495b","#edae49","#00798c")

#2. Read in the Nova Scotia shapefile and plot it 
## Projections you can use for various mapping elements
latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#3. Make a basemap to use for a figure using the rnaturalearth package

basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
  dplyr::select(name_en,geometry)%>%
  st_as_sf()%>%
  st_union()%>%
  st_transform(latlong)%>%
  st_as_sf()%>%
  mutate(country="Canada")%>%
  rbind(.,ne_states(country = "United States of America",returnclass = "sf")%>%
          dplyr::select(name_en,geometry)%>%
          st_as_sf()%>%
          st_union()%>%
          st_transform(latlong)%>%
          st_as_sf()%>%
          mutate(country="US"),
        ne_states(country = "Greenland",returnclass = "sf")%>%
          dplyr::select(name_en,geometry)%>%
          st_as_sf()%>%
          st_union()%>%
          st_transform(latlong)%>%
          st_as_sf()%>%
          mutate(country="Greenland"))%>%
  st_transform(CanProj)

#outer boundary for the 'zoom out' map - eastern canada
can_bound <- ne_states(country = "Canada",returnclass = "sf")%>%
                filter(name_en %in% c("Ontario","Quebec","New Brunswick","Newfoundland and Labrador","Prince Edward Island"))%>%
                st_transform(CanProj)%>%
                st_bbox()

##basemap of Nova Scotia - this is where we use the 'shapefile' in the Github data folder (the .shp one)
ns_coast <- read_sf("data/NS_coastline_project_Erase1.shp")%>%
            st_transform(latlong)%>%
            mutate(name="Nova Scotia")%>%
            dplyr::select(name,geometry)

#extent of focus area for mapping
focal_bound <- ns_coast%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_transform(utm)%>%
                st_buffer(20)%>% #  km buffer around the focal zone
                st_transform(CanProj)%>%
                st_bbox()

#Read in the sample site coordinates and convert to a spatial 'sf' object - need to add the 2023 Cape Breton site coordinates to this
sample_sites <-read.csv("data/Map.csv")%>%
              st_as_sf(coords=c("Longitude","Latitude"),crs=latlong)%>% #had to change lon and lat to Longitude/Latitude since that's what the columns are called in the Map.csv file
              st_transform(CanProj)%>%
              mutate(year=ifelse(Code %in% c("CHB","CHT","ASB","NOR","LIN","MIR","ESB"),"2023","2019"),
                     group=case_when(Code %in% c("CHB","FAI","RSE") ~ "South",
                                     Code %in% c("LIN","NOR","ESB","MIR",
                                                "CHT","ASB") ~ "CB",
                                     TRUE ~ "East"))

p1 <- ggplot()+
  geom_sf(data=basemap%>%filter(country=="Canada"),fill="grey")+
  geom_sf(data=basemap%>%filter(country!="Canada"),fill="grey90")+
  geom_sf(data=ns_coast,fill="grey")+
  geom_sf(data=sample_sites,aes(shape=year,fill=group),size=3)+
  coord_sf(expand=0,xlim=focal_bound[c(1,3)],ylim=focal_bound[c(2,4)])+
  theme_bw()+
  scale_fill_manual(values=fill_pal)+
  scale_shape_manual(values=c(21,23))+
  theme(legend.position="none")

p2 <- ggplot()+
  geom_sf(data=basemap%>%filter(country=="Canada"),fill="grey")+
  geom_sf(data=basemap%>%filter(country!="Canada"),fill="grey90")+
  geom_sf(data=ns_coast,fill="grey")+
  geom_sf(data=focal_bound%>%st_as_sfc(),fill=NA)+
  theme_bw()+
  theme(axis.text = element_blank())+
  coord_sf(expand=0,xlim=can_bound[c(1,3)],ylim=can_bound[c(2,4)])

p3 <- p1 + p2 + plot_layout(ncol=2)

ggsave("output/readme_plot_1.png",p3,height=5,width=10,units="in",dpi=300)
