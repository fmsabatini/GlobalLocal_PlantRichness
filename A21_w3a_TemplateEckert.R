library(tidyverse)
library(sp)
library(rgdal)
library(viridis)
library(sf)
library(rnaturalearth)
library(dggridR)
library(cowplot)

#Template of the world map
#### Get geographic data
## Downloaded from NaturalEarth
countries <- readOGR("/data/sPlot/users/Francesco/Ancillary_Data/naturalearth/ne_110m_admin_0_countries.shp") %>% 
  st_as_sf() %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()
graticules <- readOGR("/data/sPlot/users/Francesco/Ancillary_Data/naturalearth/ne_110m_graticules_15.shp") %>% 
  st_as_sf() %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()

bb <- readOGR("/data/sPlot/users/Francesco/Ancillary_Data/naturalearth/ne_110m_wgs84_bounding_box.shp") %>% 
  st_as_sf() %>% 
  st_transform(crs = "+proj=eck4") %>% 
  st_geometry()


label_parse <- function(breaks) {
  parse(text = breaks)
}

## basic graph of the world in Eckert projection
w3a <- ggplot() +
  #geom_sf(data = graticules, col = "grey20", lwd = 0.1) +
  geom_sf(data = countries, fill = "grey90", col = NA, lwd = 0.3) +
  geom_sf(data = bb, col = "grey20", fill = NA) +
  coord_sf(crs = "+proj=eck4") +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
        legend.key.height = unit(1.1, "cm"), 
        legend.key.width = unit(1.1, "cm")) +
  scale_fill_viridis(label=label_parse) + 
  scale_color_viridis(label=label_parse)
