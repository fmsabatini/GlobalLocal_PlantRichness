library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(openxlsx)
#install.packages("remotes")
#remotes::install_github("SEEG-Oxford/seegSDM")
library(seegSDM)
library(cleangeo)



rasterOptions(tmpdir="/data/sPlot/users/Francesco/_tmp")
write("TMPDIR = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('TMPDIR'), '.Renviron'))
write("R_USER = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# Import mydata as template
load("../sPlot/_input/Mydata_global.RData")


#### Import New data from literature ####
newdata <- read.xlsx("../sPlot/_Ancillary/NatComm_R1_AdditionalPlots.xlsx", sheet = 1)
newdata <- newdata %>% 
  dplyr::select(-X19) %>% 
  filter(Use) %>% 
  mutate(Reference=coalesce(Reference_complete, Reference)) %>% 
  rename(POINT_X="Long", POINT_Y="Lat") %>% 
  rename(PlotObservationID="PseudoPlotObservationID", 
         rich.complete="Richness_all.(no.epi)",
         rich.woody_all="Richness_woody_all", 
         rich.woody_large="Richness_woody_large") %>% 
  mutate_at(.vars=vars(starts_with("rich")), 
            .funs=~as.numeric(ifelse(.=="-", NA, .))) %>% 
  mutate(plants_recorded=factor(Sampled_vegetation)) %>% 
  mutate(plants_recorded=fct_collapse(plants_recorded, 
                                     complete=c("complete_vegetation", "complete"))) %>% 
  mutate(isforest=1) %>% 
  rename(Rel.area="Plot_size") %>% 
  mutate(rich=ifelse(plants_recorded=="complete", rich.complete, 
                     ifelse(plants_recorded=="woody_all", rich.woody_all, 
                            ifelse(plants_recorded=="woody_large", rich.woody_large, NA)))) %>% 
  # temp - filter out newdata too large or too small
  filter(Rel.area <= 25000) %>% 
  filter(Rel.area >= 10) 




#### Load Predictors ####
##### Vector Predictors ####
# 1. Ecoregions
ecoreg <- readOGR("../sPlot/_Ancillary/official", layer="wwf_terr_ecos")
EcoRichness <- xlsx::read.xlsx("../sPlot/_Ancillary/Kier_et_al_SI_jbi_1272_sm_sa2.xls", 
                               sheetIndex = 1, startRow = 4, header=T)
ecoreg@data <- ecoreg@data %>% 
  dplyr::select(OBJECTID, ECO_NAME, REALM, BIOME, ECO_NUM, ECO_ID, eco_code)

# 2. Continent
sPDF <- rworldmap::getMap(resolution="coarse")
continent <- sPDF[,"continent"]
continent@data[243,"continent"] <- "South America" ## Manually correct missing data
continent = clgeo_Clean(continent)

##### Raster Predictors ####
# 1. Biomes
sBiomes.shp <- readOGR(dsn="../../Ancillary_Data/Biomes_sPlot/", layer="sBiomes")
sBiomes.labels <- sBiomes.shp@data
## checked: BiomeIDs match with sBiomes
sBiomes <- raster("../sPlot/_Ancillary/sBiome_raster01/sBiomes_raster01.tif")

# 2. Climate
CHELSA.stack <- stack(paste("/data/sPlot/users/Francesco/Ancillary_Data/CHELSA/CHELSA_bio10_", 
                            stringr::str_pad(1:19, width=2, side="left", pad="0"), ".tif", sep=""))

# 3. Soil 
ISRIC.paths <- list.files("/data/sPlot/users/Francesco/Ancillary_Data/ISRIC/", 
                          pattern = ".tif$", full.names = T)
ISRIC.paths <- ISRIC.paths[!str_detect(ISRIC.paths, pattern = "Generalized")]
bulk <- which(str_detect(ISRIC.paths, pattern="BLDFIE"))
ISRIC.stack <- stack(ISRIC.paths[-bulk])
ISRIC.bulk <- raster(ISRIC.paths[bulk])

# 4. Topography
landform1km.maj <- raster("../sPlot/_predictors/landform1km.maj.masked.tif")
tri50km <- raster("/data/sPlot/users/Francesco/Ancillary_Data/EarthEnv/tri_50KMmd_GMTEDmd.tif")
landform50km.count <- raster("/data/sPlot/users/Francesco/Ancillary_Data/EarthEnv/geom_50KMcount_GMTEDmd.tif")

# 5. Climate Change velocity
CCV.Prec <- raster("../sPlot/_Ancillary/ClimateChangeVelocity/LGM to Present Precipitation Velocity.tif")
CCV.Temp <- raster("../sPlot/_Ancillary/ClimateChangeVelocity/LGM to Present Temperature Velocity.tif")
CCV.stack <- raster::stack(CCV.Prec, CCV.Temp)
names(CCV.stack) <- c("CCVPre", "CCVTem")

##### Recalculate World PCAs for Soil & Climate #####
load("../sPlot/_predictors/predictors_world.RData")

#Climate
chelsa.world.avg <- apply(chelsa.world %>% 
                            filter(complete.cases(.)), MARGIN=2, "mean")
chelsa.world.sd <- apply(chelsa.world %>% 
                           filter(complete.cases(.)), MARGIN=2, "sd")
chelsa.pca <- prcomp(chelsa.world %>% 
                       filter(complete.cases(.)), center = T, scale = T)

#Soil
isric.world.avg <- apply(isric.world %>% 
                           filter(complete.cases(.)), MARGIN=2, "mean")
isric.world.sd <- apply(isric.world %>% 
                          filter(complete.cases(.)), MARGIN=2, "sd")
isric.pca <- prcomp(isric.world %>% 
                      filter(complete.cases(.)), center = T, scale = T)





#### Extract Predictors for additional plots ####
##### Prepare spatial data #####
shp.newdata <- newdata %>% 
  filter(!is.na(POINT_X) | !is.na(POINT_Y))
newdata.shp <- SpatialPointsDataFrame(coords=shp.newdata %>%
                                       dplyr::select(POINT_X, POINT_Y),
                                     proj4string = crs(ecoreg), 
                                     data=data.frame(PlotID=shp.newdata$PlotID))
rm(shp.newdata)


##### Vectors #####
# 1. Ecoregions
ecoreg.out <- sp::over(x=newdata.shp, y=ecoreg) # ~5 mins
# 103 unmatched plots - find nearest neighbour
toassign <- newdata.shp[which(is.na(ecoreg.out$ECO_NAME)),]

library(parallel)
library(doParallel)
cl <- makeCluster(4, outfile="" )
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(rgdal)
  library(raster)
  library(sp)
  library(dplyr)
  library(rgeos)
})
nearestEcoregion <- foreach(i=1:length(toassign), .packages=c('raster'), .combine=rbind) %dopar% {  
  #  print(i)
  ## create a subset of geoentities based on a 5° buffer radius around each target plot.
  tmp.buff <- gBuffer(toassign[i,], width=5) 
  tryCatch(
    tmp.mypredictor <- spatialEco::spatial.select(
      x = tmp.buff,
      y = ecoreg,
      distance = 0.1,
      predicate = "intersect"
    ),
    error = function(e) {
      print(paste("Nothing close enough for plot", toassign@data$PlotObservationID[i]))
    }
  )
  # find nearest neighbour  
  nearest.tmp <- tryCatch(tmp.mypredictor@data[geosphere::dist2Line(toassign[i,],
                                                                    tmp.mypredictor)[,"ID"],],
                          error = function(e){
                            ee <- ecoreg@data[1,, drop=F]
                            ee[1,] <- rep(NA, ncol(ecoreg))
                          }
  ) %>% 
    mutate_all(~as.character(.))
  return(nearest.tmp)
}


ecoreg.out[is.na(ecoreg.out$ECO_NAME),] <- nearestEcoregion
rm(toassign)

# 2. Continent
continent.out <- sp::over(x=newdata.shp, y=continent)
#correct unassigned points to closest continent
toassign <- newdata.shp[which(is.na(continent.out$continent)),] #146 unassigned
crs(toassign) <- crs(continent)

nearestContinent <- foreach(i=1:length(toassign), .packages=c('raster'), .combine=rbind) %dopar% {  
  #  print(i)
  ## create a subset of geoentities based on a 5° buffer radius around each target plot.
  tmp.buff <- gBuffer(toassign[i,], width=5) 
  tryCatch(
    tmp.mypredictor <- spatialEco::spatial.select(
      x = tmp.buff,
      y = continent,
      distance = 0.1,
      predicate = "intersect"
    ),
    error = function(e) {
      print(paste("Nothing close enough for plot", toassign@data$PlotObservationID[i]))
    }
  )
  # find nearest neighbour  
  nearest.tmp <- tryCatch(tmp.mypredictor@data[geosphere::dist2Line(toassign[i,],
                                                                    tmp.mypredictor)[,"ID"],],
                          error = function(e){
                            ee <- continent@data[1,, drop=F]
                            ee[1,] <- rep(NA, ncol(continent))
                          }
  ) %>% 
    as.character()
  return(nearest.tmp)
}
continent.out[which(is.na(continent.out$continent)),] <- nearestContinent
rm(toassign)






##### Rasters #####
# 1. Biomes
sbiomes.out <- raster::extract(sBiomes, newdata.shp)
toassign <- newdata.shp[which(sbiomes.out==0),] #28 still to assign

## find nearest biome for unmatched plots
nearestBiome <- foreach(i=1:length(toassign), .packages=c('raster'), .combine=rbind) %dopar% {  
  ## create a subset of geoentities based on a 5° buffer radius around each target plot.
  tmp.buff <- gBuffer(toassign[i,], width=5) 
  tryCatch(
    tmp.mypredictor <- spatialEco::spatial.select(
      x = tmp.buff,
      y = sBiomes.shp,
      distance = 0.1,
      predicate = "intersect"
    ),
    error = function(e) {
      print(paste("Nothing close enough for plot", toassign@data$PlotObservationID[i]))
    }
  )
  # find nearest neighbour  
  nearest.tmp <- tryCatch(tmp.mypredictor@data[geosphere::dist2Line(toassign[i,],
                                                                    tmp.mypredictor)[,"ID"],],
                          error = function(e){
                            ee <- sBiomes.shp@data[1,, drop=F]
                            ee[1,] <- rep(NA, ncol(sBiomes.shp))
                          }
  ) %>% 
    mutate_all(~as.character(.))
  return(nearest.tmp)
}
sbiomes.out[which(sbiomes.out==0)] <- as.numeric(nearestBiome$BiomeID)


# 2. Climate 
climate.out <- raster::extract(CHELSA.stack, newdata.shp)
toassign <- newdata.shp[which(!complete.cases(climate.out)),]  # 56 unmatched plots
## find nearest neighbour
near.land <- nearestLand(toassign@coords, CHELSA.stack[[1]], max_distance=10000)
assigned <- raster::extract(CHELSA.stack, near.land)
# two plots remain unassigned
climate.out[which(!complete.cases(climate.out)),] <- assigned
sum(!complete.cases(climate.out)) # 0 still unassigned

# 3. Soil
ISRIC.out <- raster::extract(ISRIC.stack, newdata.shp)
toassign <- newdata.shp[which(!complete.cases(ISRIC.out)),]  # 177 unmatched plots
## find nearest neighbour
near.land <- nearestLand(toassign@coords, ISRIC.stack[[1]], max_distance=10000)
assigned <- raster::extract(ISRIC.stack, near.land)
ISRIC.out[which(!complete.cases(ISRIC.out)),] <- assigned
sum(!complete.cases(ISRIC.out)) # 0 still unassigned  

ISRIC.bulk.out <- raster::extract(ISRIC.bulk, newdata.shp)
toassign <- newdata.shp[which(!complete.cases(ISRIC.bulk.out)),]  # 177 unmatched plots
## find nearest neighbour
near.land <- nearestLand(toassign@coords, ISRIC.bulk, max_distance=10000)
assigned <- raster::extract(ISRIC.bulk, near.land)
ISRIC.bulk.out[which(!complete.cases(ISRIC.bulk.out))] <- assigned
sum(!complete.cases(ISRIC.bulk.out)) # 0 still unassigned  

  
# 4. Topography
tri50km.out <- raster::extract(tri50km, newdata.shp)
#landform1km.maj.out <- foreach(i=1:length(newdata.shp), .packages=c('raster'), .combine=rbind) %dopar% { 
landform1km.maj.out <- raster::extract(landform1km.maj, newdata.shp) #}
landform1km.maj.out <- factor(landform1km.maj.out, levels=1:10, labels=c("flat", "peak", "ridge", 
                                                                         "shoulder", "spur", "slope",
                                                                         "hollow", "footslope", "valley",
                                                                         "pit"))
landform50km.count.out <- raster::extract(landform50km.count, newdata.shp)

# 5. Climate change velocity
CCV.out <- raster::extract(CCV.stack, newdata.shp)
toassign <- newdata.shp[which(!complete.cases(CCV.out)),]  # 31 unmatched plots
## find nearest neighbour
near.land <- nearestLand(toassign@coords, CCV.stack[[1]], max_distance=30000)
assigned <- raster::extract(CCV.stack, near.land)
CCV.out[which(!complete.cases(CCV.out)),] <- assigned
sum(!complete.cases(CCV.out)) # 0 still unassigned  


#### Recombine predictors in newdata  #####
newdata2 <- newdata %>% 
  rename(RELEVE_NR=PlotObservationID) %>% 
  bind_cols(data.frame(CONTINENT=continent.out$continent,
                       climate.out,
                       data.frame(BLDFIE=ISRIC.bulk.out),
                       ISRIC.out, 
                       CCV.out, 
                       tri50km.out, 
                       landform1km.maj.out, 
                       landform50km.count.out, 
                       (data.frame(sBiomeID=sbiomes.out) %>%
                         left_join(sBiomes.labels %>% 
                                     dplyr::select(BiomeID, Name) %>%
                                     rename(sBiomeID=BiomeID, sBiomeName=Name),
                                   by="sBiomeID")),  
                       ecoreg.out %>% 
                         dplyr::select(ECO_NAME, eco_code, REALM) %>% 
                         left_join(EcoRichness %>%
                                     dplyr::select(eco_id2, sp_wfig) %>%
                                     dplyr::rename(eco_code=eco_id2), by="eco_code")
                       )) %>% 
  rename_at(.vars=vars(starts_with("CHELSA")), 
            .funs=~stringr::str_sub(., -8,-1)) %>% 
  rename_at(.vars=vars(starts_with("bio10_")),
            .funs=~gsub(pattern="bio10_", replacement="bio_", x=.)) %>% 
  rename_at(.vars=vars(ends_with("_M_sl2_250m_ll")),
            .funs=~gsub(pattern="_M_sl2_250m_ll", replacement="", x=.)) %>% 
  rename_at(.vars=vars(ends_with(".out")),
            .funs=~gsub(pattern=".out", replacement="", x=.)) %>% 
  mutate_at(.vars=vars(starts_with("rich")), 
            .funs=~round(.))





#### Project points in PCA space ####
#Climate
standardized.chelsa <- t(as.matrix( (t(newdata2 %>% 
                                         dplyr::select(bio_01:bio_19)) - chelsa.world.avg) / chelsa.world.sd))
projected.chelsa <- as.data.frame(standardized.chelsa %*% as.matrix(chelsa.pca$rotation))


#Soil
standardized.isric <- t(as.matrix( (t(newdata2 %>% 
                                        dplyr::select(BLDFIE:SNDPPT)) - isric.world.avg) / isric.world.sd))
projected.isric <- as.data.frame(standardized.isric %*% as.matrix(isric.pca$rotation))

#Merge to newdata2 object
newdata2 <- newdata2 %>% 
  bind_cols(as.data.frame(projected.chelsa) %>% 
              dplyr::select(PC1:PC5) %>% 
              rename_all(.funs = list(~ paste0(.,"_chelsa")))) %>% 
  bind_cols(as.data.frame(projected.isric) %>% 
              dplyr::select(PC1:PC4) %>% 
              rename_all(.funs = list(~ paste0(.,"_isric"))))

## final corrections
newdata2 <- newdata2 %>% 
  mutate(CONTINENT=fct_recode(CONTINENT,
                              AF ="Africa", 
                              AN = "Antarctica", 
                              AU = "Australia",
                              EU = "Eurasia", 
                              "NA"="North America", 
                              SA ="South America"))


save(newdata2, file = "../sPlot/_Ancillary/NatComm_newdata2.RData")

## double check PCAs
#aa <- mydata %>% 
#  filter(sBiomeName %in% c("Subtrop. with year-round rain",
#                        "Temperate midlatitudes", 
#                        "Tropics with summer rain", 
#                        "Tropics with year-round rain")) %>% 
#  sample_n(3000) 
#
#ggplot(data=aa) + 
#  geom_point(data=newdata2, 
#             aes(x=BLDFIE, y=PC1_isric), col="black", alpha=0.5) +
#  geom_point(aes(x=BLDFIE, y=PC1_isric), col="red", alpha=0.5)
#



#### Import data from sPlot 3.0 and reshape them ####

## Explore data from sPlot 3.0, not includes in sPlot 2.0 from the tropics.
## Evaluate whether it makes sesnse to add the to the analysis
load("/data/sPlot/releases/sPlot3.0/header_sPlot3.0.RData")


## select only NEW datasets from outside Europe

new.datasets <- c("AF-00-010", "AS-ID-XXX", "SA-CO-003", "US-NA-016", 
                  "SA-EC-002", "SA-UY-001", "AF-00-011", "AF-NA-001", 
                  "AF-CM-001", "SA-PE-001", "AS-ID-002",
                  #"AF-EG-XXX", ### need to double check for duplicate records
                  "AS-CN-008", "SA-AR-003")
                  # "NA-CU-XXX", can't use for lack of Rel.ara


dbs <- read_csv("/data/sPlot/users/Francesco/_sPlot_Management/Consortium/Databases.out.csv")
# make summary
header %>% 
  filter(`GIVD ID` %in% new.datasets) %>% 
  count(`GIVD ID`, Dataset) %>% 
  left_join(dbs %>% dplyr::select(`GIVD ID`, Custodian))

"   `GIVD ID` Dataset          n Custodian             
   <chr>     <chr>        <int> <chr>                 
 1 AF-00-010 Afro-Alpine    250 Petr Sklenar          
 2 AF-00-011 Ivory_Coast    105 Bruno Hérault         
 3 AF-CM-001 Camaroon       172 Jiri Dolezal          
 4 AF-NA-001 Namibia       1308 Ben Strohbach         
 5 AS-CN-008 China_hainan   472 Hua-Feng Wang         
 6 AS-ID-002 Sumatra        160 Holger Kreft          
 7 AS-ID-XXX Ladakh        4623 Jiri Dolezal          
 8 NA-CU-XXX Cuba           613 Ute Jandt             # can't use
 9 SA-AR-003 Patagonia      147 Karina Speziale       
10 SA-CO-003 Colombia       207 Esteban Alvarez-Davila
11 SA-EC-002 Galapagos      105 Gonzalo Rivas-Torres  
12 SA-PE-001 Peru           152 Antonio Galán-de-Mera 
13 SA-UY-001 Uruguay        308 Felipe Lezama         "

## Filter & reshape header

header3 <- header %>% 
  filter(`GIVD ID` %in% new.datasets) %>% 
  rename(GIVD_ID=`GIVD ID`, 
         RELEVE_NR=PlotObservationID,
         plants_recorded=`Plants recorded`, 
         POINT_X=Longitude, 
         POINT_Y=Latitude, 
         Rel.area=`Relevé area (m²)`, 
         ECO_NAME=Ecoregion,
         eco_code=EcoregionID, 
         sBiomeName=sBiome) %>% 
  mutate(isforest=coalesce(Forest, is.forest)) %>% 
  left_join(EcoRichness %>%
              dplyr::select(eco_code=eco_id, eco_id2, sp_wfig),
            by="eco_code") %>%  
  dplyr::select(-eco_code, eco_code=eco_id2) %>% 
  dplyr::select(any_of(colnames(mydata))) %>% 
  left_join(ecoreg@data %>% 
              dplyr::select(ECO_NAME, REALM) %>% 
              distinct(ECO_NAME, .keep_all = T), 
            by="ECO_NAME") %>% 
  mutate(plants_recorded=fct_recode(plants_recorded, 
                                    woody_all = "Woody plants >= 2 cm dbh",
                                    complete = "All vascular plants")) %>% 
  mutate(plants_recorded=as.character(plants_recorded)) %>% 
  mutate(plants_recorded=ifelse(GIVD_ID %in% c("AF-CM-001"), 
                                "woody_large", 
                                plants_recorded)) %>% 
  #mutate(plants_recorded=ifelse(GIVD_ID %in% c("AF-00-010"), 
  #                              "woody_all", 
  #                              plants_recorded)) %>% 
  mutate(plants_recorded=ifelse(GIVD_ID %in% c("AF-00-010", "AF-NA-001", 
                                               "AF-00-011",
                                               "AS-ID-XXX", "NA-CU-XXX", 
                                               "SA-CO-003", "SA-EC-002", 
                                               "SA-UY-001"), 
                                "complete", 
                                plants_recorded)) %>% 
  mutate(plants_recorded=ifelse(GIVD_ID =="AS-CN-008" & is.na(plants_recorded), 
                                "complete", 
                                plants_recorded)) %>% 
  mutate(plants_recorded=factor(plants_recorded, 
                                levels=c("complete","woody_all", "woody_large"))) %>% 
  mutate(Rel.area=replace(Rel.area, 
                          list = (GIVD_ID=="AS-CN-008" & is.na(Rel.area)), 
                          values = 400)) %>% 
  mutate(Rel.area=replace(Rel.area, 
                          list = (GIVD_ID=="AF-00-010" & is.na(Rel.area)), 
                          values = 200)) %>% 
  mutate(Rel.area=replace(Rel.area, 
                          list = (GIVD_ID=="SA-CO-003" & is.na(Rel.area)), 
                          values = 200)) %>% 
  mutate(Rel.area=ifelse( (GIVD_ID=="SA-AR-003" & Rel.area<0) ,
                          -Rel.area, 
                          Rel.area)) %>% 
  mutate(Rel.area=replace(Rel.area, 
                          list = (GIVD_ID=="AS-ID-XXX" & is.na(Rel.area)), 
                          values = 10000)) %>% 
  mutate(isforest=ifelse(GIVD_ID %in% c("AF-00-011", "AF-CM-001", 
                                        "AS-CN-008", "AS-ID-002", 
                                        "SA-CO-003"), 
                         T, 
                         isforest)) %>% 
  mutate(isforest=ifelse(GIVD_ID %in% c("AF-00-010", "AS-ID-XXX"), 
                         F, 
                         isforest)) %>% 
  mutate(isforest=ifelse(GIVD=="SA-EC-002" & !is.na(is.forest), 
                         is.forest, 
                         F))
         

# Add species richness 
load("/data/sPlot/releases/sPlot3.0/DT_sPlot3.0.RData")
load("/data/sPlot/releases/sPlot3.0/Traits_CWMs_sPlot3.RData")

sel.plots <- header3 %>% 
  pull(RELEVE_NR)

(Sp.richness <- DT2 %>% 
  filter(PlotObservationID %in% sel.plots) %>% 
  mutate(Species=coalesce(Species, Species_original)) %>% 
  filter(!Taxon_group %in% c("Alga_Stonewort", "Lichen", "Moss")) %>% 
  distinct(PlotObservationID, Species) %>% 
  # add growth forms
  left_join(sPlot.traits %>% 
              dplyr::select(Species, GrowthForm), 
            by="Species") %>% 
  # calculate species richness by GF
  left_join({.} %>% 
              filter(str_detect(GrowthForm, "shrub|tree")) %>% 
              count(PlotObservationID) %>% 
              rename(rich.woody_all=n), 
            by="PlotObservationID") %>% 
  left_join({.} %>% 
              filter(GrowthForm %in% c("tree", "shrub/tree", "herb/shrub/tree")) %>% 
              count(PlotObservationID) %>% 
              rename(rich.woody_large=n), 
            by="PlotObservationID") %>% 
  count(PlotObservationID, rich.woody_all, rich.woody_large) %>% 
  rename(rich=n) %>% 
  left_join(header3 %>% 
              dplyr::select(RELEVE_NR, plants_recorded), 
            by=c("PlotObservationID"="RELEVE_NR")) %>% 
  mutate(rich.woody_large=ifelse(plants_recorded=="woody_large", rich, rich.woody_large)) %>% 
  replace_na(list(rich.woody_large=0, rich.woody_all=0)) %>% 
  mutate(rich.woody_all=ifelse(plants_recorded=="woody_large", NA, rich.woody_all)) %>% 
  mutate(rich.woody_all=ifelse(plants_recorded=="woody_all", rich, rich.woody_all)) %>% 
  mutate(rich.complete=ifelse(plants_recorded=="complete", rich, NA)) %>% 
  dplyr::select(-rich)
)


## Attach environmental predictors
load("/data/sPlot/releases/sPlot3.0/SoilClim_sPlot3.RData")
soilclim3 <- soilclim %>% 
  filter(PlotObservationID %in% sel.plots) %>% 
  rename_at(.vars=vars(starts_with("bio")), 
            .fun=~str_replace(., "bio", "bio_")) %>% 
  dplyr::select(any_of(colnames(mydata)))

## extract remaining predictors
header3.shp <- SpatialPointsDataFrame(coords=header3 %>%
                                        dplyr::select(POINT_X, POINT_Y),
                                      proj4string = crs(ecoreg), 
                                      data=data.frame(RELEVE_NR=header3$RELEVE_NR))

# Topography
tri50km.out <- raster::extract(tri50km, header3.shp)
#landform1km.maj.out <- foreach(i=1:length(newdata.shp), .packages=c('raster'), .combine=rbind) %dopar% { 
landform1km.maj.out <- raster::extract(landform1km.maj, header3.shp) #}
landform1km.maj.out <- factor(landform1km.maj.out, levels=1:10, labels=c("flat", "peak", "ridge", 
                                                                         "shoulder", "spur", "slope",
                                                                         "hollow", "footslope", "valley",
                                                                         "pit"))
landform50km.count.out <- raster::extract(landform50km.count, header3.shp)

# Climate change velocity
CCV.out <- raster::extract(CCV.stack, header3.shp)
toassign <- header3.shp[which(!complete.cases(CCV.out)),]  # 61 unmatched plots
## find nearest neighbour
near.land <- nearestLand(toassign@coords, CCV.stack[[1]], max_distance=80000)
assigned <- raster::extract(CCV.stack, near.land)
CCV.out[which(!complete.cases(CCV.out)),] <- assigned
sum(!complete.cases(CCV.out)) # 0 still unassigned  

#### Project points in PCA space ####
#Climate
standardized.chelsa <- t(as.matrix( (t(soilclim3 %>% 
                                         dplyr::select(bio_01:bio_19)) - chelsa.world.avg) / chelsa.world.sd))
projected.chelsa <- as.data.frame(standardized.chelsa %*% as.matrix(chelsa.pca$rotation))


#Soil
standardized.isric <- t(as.matrix( (t(soilclim3 %>% 
                                        dplyr::select(BLDFIE:SNDPPT)) - isric.world.avg) / isric.world.sd))
projected.isric <- as.data.frame(standardized.isric %*% as.matrix(isric.pca$rotation))

## Reassemble
header3 <- header3 %>% 
  left_join(Sp.richness %>% 
              dplyr::select(-plants_recorded), by=c("RELEVE_NR"="PlotObservationID")) %>% 
  bind_cols(soilclim3) %>% 
  bind_cols(data.frame(tri50km=tri50km.out, 
                       landform1km.maj=landform1km.maj.out, 
                       landform50km.count=landform50km.count.out, 
                       CCV.out)) %>% 
  bind_cols(as.data.frame(projected.chelsa) %>% 
              dplyr::select(PC1:PC5) %>% 
              rename_all(.funs = list(~ paste0(.,"_chelsa")))) %>% 
  bind_cols(as.data.frame(projected.isric) %>% 
              dplyr::select(PC1:PC4) %>% 
              rename_all(.funs = list(~ paste0(.,"_isric")))) %>% 
  #change PlotObservationID
  mutate(RELEVE_NR=RELEVE_NR+10000000) %>% 
  mutate(rich=ifelse(plants_recorded=="complete", rich.complete, 
                     ifelse(plants_recorded=="woody_all", rich.woody_all, 
                            ifelse(plants_recorded=="woody_large", rich.woody_large, NA)))) %>% 
  mutate(CONTINENT=fct_recode(CONTINENT,
                              "NA"="N-A", 
                              SA="S-A"))


colnames(mydata)[which(!colnames(mydata) %in% colnames(header3))]

save(header3, file = "../sPlot/_Ancillary/NatComm_header3.RData")


### Create mydata2 ####
load("../sPlot/_input/Mydata_global.RData")
load("/data/sPlot/releases/sPlot2.1/sPlot_header_20161124.RData")
header.v21 <- header
load("../sPlot/_Ancillary/NatComm_header3.RData")
load("/data/sPlot/releases/sPlot3.0/header_sPlot3.0.RData")
header.v3 <- header
load("../sPlot/_Ancillary/NatComm_newdata2.RData")


mydata2 <- mydata %>% 
  ## exclude Naturalness==3 plots
  filter(RELEVE_NR %in% (header.v21 %>% 
                           filter(Naturalness!=3 | is.na(Naturalness)) %>% 
                           pull(PlotObservationID))) %>% 
  bind_rows(newdata2) %>% 
  ## exclude Naturalness==3 plots
  bind_rows(header3 %>% 
              filter(RELEVE_NR %in% (header.v3 %>% 
                                       filter(Naturalness!=3 | is.na(Naturalness)) %>%
                                       mutate(RELEVE_NR=10000000 + PlotObservationID) %>% 
                                       pull(RELEVE_NR)))) %>% 
  dplyr::select(all_of(colnames(mydata))) %>% 
  mutate(sBiomeName=factor(sBiomeName, levels=levels(mydata$sBiomeName))) %>% 
  mutate(ECO_NAME=factor(ECO_NAME, levels=levels(mydata$ECO_NAME))) %>% 
  mutate(isforest=as.logical(isforest)) 
  #dplyr::select(all_of(selected.predictors))

mydata <- mydata2 %>% 
  filter(Rel.area >= 10) 

save(mydata, file = "../sPlot/_input/Mydata_global_NatCommR1.RData")










## Load BRT model and check how predictions match the new data
library(dismo)
library(gbm)

iter <- 10
metric <- "sr10" # index.list$Var1[index]
size <- 10
which.model <- "all"
#  for(metric in c("sr10", "sr100", "sr400", "sr1ha", "Asymp")){
brt.models.path <- "../sPlot/_versions/NatCommR1/BRTglobal"
listf <- list.files(brt.models.path, 
                    pattern=paste0("^",which.model,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T)

order.i <- as.numeric(gsub(pattern=paste0("_",size,"m\\.RData"), replacement="", 
                           x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))

listf <- listf[order(order.i)]

load(listf[iter])
#- load regression coefficient to correct BRT prediction bias
load("../sPlot/_versions/NatCommR1/_world_predictions/BiasCorrect_regrCoefs_all_.RData")



## Predict & Add BiasCorrection 
mydata.out <- mydata
mydata.out$pred.rich <- predict.gbm(object = modello, 
                                  newdata = mydata, type = "response")    

mydata.out2 <- mydata.out %>% 
  dplyr::rename(pred.rich.raw=pred.rich) %>% 
  mutate(pred.rich.bc = (pred.rich.raw-mycoefs[[1]][iter])/mycoefs[[2]][iter]) %>% 
  mutate(pred.rich.bc = ifelse(pred.rich.bc<0, 0, pred.rich.bc)) %>% 
  mutate(rich=ifelse(plants_recorded=="complete", rich.complete, 
                     ifelse(plants_recorded=="woody_all", rich.woody_all, 
                            ifelse(plants_recorded=="woody_large", rich.woody_large, NA))))

mygg <- ggplot(mydata.out2 %>% 
         #filter(!is.na(CONTINENT)) %>%  ### to delete
         #filter(rich<700) %>% 
         mutate(Rel.area.cut=cut(Rel.area, c(0,150, 600, 1200, Inf)))) +
  geom_point(aes(x=rich, y=pred.rich.bc, col=CONTINENT)) + 
  geom_abline(intercept=0, slope=1, col=2) +
  theme_bw() + 
  coord_equal() + 
  xlim(0,400) + 
  ylim(0,400) + 
  facet_grid(.~plants_recorded)

cor(mydata.out2$rich, mydata.out2$pred.rich.bc, use="complete.obs" )



#### alternative plotting
library(sf)
library(viridis)
source("A21_w3a_TemplateEckert.R")

mydata2.sf <- mydata2 
coordinates(mydata2.sf) <- ~POINT_X+POINT_Y
crs(mydata2.sf) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
mydata2.sf <- mydata2.sf %>% 
  st_as_sf() %>% 
  st_transform("+proj=eck4")

w3a + 
  geom_sf(data=mydata2.sf, aes(col=plants_recorded), alpha=0.7, cex=0.7)
#HAHAHAHHAHA we are messing up with your coooode XXOXOXOXOXO 







