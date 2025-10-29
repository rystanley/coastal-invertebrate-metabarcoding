# Load necessary libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(marmap)
library(ggspatial)

#map projections
latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
           dplyr::select(name_en,geometry)%>%
           st_as_sf()%>%
           st_union()%>%
           st_transform(CanProj)%>%
           st_as_sf()

#sample site coordinates
sample_sites <- data.frame(code=c("ASB","CHB"),
                           lon=c(-60.450273,-66.105158), #these are immedately adjacent to the samples sites. Based on the resolution of the marmap depth raster, these have negative depths and are needed for the analysis
                           lat=c(46.910745,43.725419))%>%
                st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)

all_sites <-read.csv("data/Map.csv")%>%
            st_as_sf(coords=c("Longitude","Latitude"),crs=latlong)%>% #had to change lon and lat to Longitude/Latitude since that's what the columns are called in the Map.csv file
            st_transform(CanProj)
  
range_lims <- sample_sites%>%
             st_transform(utm)%>%
             st_buffer(100)%>%
             st_transform(latlong)%>%
             st_bbox()

#get bathymetric data for the study region
bathydata<-marmap::getNOAA.bathy(lon1 = range_lims[1], lon2 = range_lims[3], 
                                 lat1 = range_lims[2], lat2 = range_lims[4],
                                 resolution = 1,keep=FALSE)


#confirm negative depths
marmap::get.depth(bathydata,
                  x=sample_sites$lon,
                  y=sample_sites$lat,
                  locator=F)

#create the transition object for the distance calculation
trans <- marmap::trans.mat(bathydata,max.depth = -70) #hug the coastline 

#format sites
sites<-sample_sites%>%data.frame%>%.[,c("lon","lat")]
rownames(sites)<-sample_sites$code

#calcualte the distance
ASB_CHB <- marmap::lc.dist(trans, 
                sites, 
                res="dist")%>%
           as.numeric();ASB_CHB

#vizualize the path
lc.dists <- marmap::lc.dist(trans, 
                            sites, 
                            res="path")

lc.line <- as.data.frame(lc.dists)%>%
            rename(lon=1,lat=2)%>%
            st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
            mutate(id=1:n())%>%
            st_union()%>%
            st_cast("LINESTRING")


lc_sf_line <- st_sfc(st_linestring(lc.dists[[1]]))
st_crs(lc_sf_line) <- latlong
lc_sf <- lc_sf_line%>%
         st_as_sf()%>%
         st_transform(CanProj)

#make the final plot

plot_lims <- sample_sites%>%
              st_transform(utm)%>%
              st_buffer(100)%>%
              st_transform(CanProj)%>%
              st_bbox()
             
p1 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=lc_sf,lwd=0.75,lty=4)+
  geom_sf(data=all_sites,size=0.75)+
  geom_sf(data=sample_sites%>%st_transform(CanProj),fill="white",shape=21,size=2)+
  theme_bw()+
  annotation_scale()+
  coord_sf(xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])+
  labs(title=paste0("Distance from Aspy Bay to Chebogue - ",ASB_CHB,"km"));p1

ggsave("output/coastal_distance_plot.png",p1,height=8,width=5,units="in",dpi=300)













