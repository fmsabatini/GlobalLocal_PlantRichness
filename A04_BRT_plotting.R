write("TMPDIR = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('TMPDIR'), '.Renviron'))
write("R_USER = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

library(gbm)
library(foreign)
library(data.table)
library(tidyverse)
library(dismo)
library(sp)
library(rgdal)
library(matrixStats)
library(viridis)
library(ggpubr)
library(pals)
library(geosphere)
library(sf)
library(dggridR)
library(rnaturalearth)
library(ggpattern)
library(cowplot)
library(furrr)
library(patchwork)
source("A20_AncillaryFunctions.R")
source("A21_w3a_TemplateEckert.R") ## import map template

### Set paths
path.to.pics <- "../sPlot/_versions/NatCommR2/_pics"
path.to.BRTglobal <- "../sPlot/_versions/NatCommR2/BRTglobal"
path.to.world.predictions <- "../sPlot/_versions/NatCommR2/_world_predictions"
path.to.output <- "../sPlot/_versions/NatCommR2"
path.to.input <- "../sPlot/_input"

##### Ancillary functions ####
{
  ShapeData <- function(out.for=NULL, out.nonfor=NULL, world0, grass.threshold=30){
  if(!is.null(out.for)){
    tmp.for <- world0 %>%
      filter(potforest==T) %>%
      dplyr::select(RAST_ID) %>%
      bind_cols(data.frame(predicted.for=out.for))
  } else tmp.for <- NULL
  
  if(!is.null(out.nonfor)){
    tmp.nonfor <- world0 %>%
      filter(isgrassland==T) %>%
      dplyr::select(RAST_ID, totgrassland) %>%
      bind_cols(data.frame(predicted.nonfor=out.nonfor))# %>% 
    #filter(totgrassland>= grass.threshold) %>% ### !!! ### 
    #dplyr::select(-totgrassland)
  } else tmp.nonfor <- NULL
  return(bind_rows(tmp.for, tmp.nonfor)) # Avoids dealing with pixels close to date line
}


## Ancillary function for plotting hotspots\coldspots
  tominmax <- function(x){
    minmax <- ifelse(x<quantile(world.out$value.out, 0.05, na.rm=T), "min", 
                     ifelse(x>quantile(world.out$value.out, 0.95, na.rm=T), "max", NA))
    minmax <- factor(minmax, levels=c( "max", "min", NA), labels=c("> 95%", "< 5%"))
    return(minmax)
  }
  
  ## Ancillary function to calculate robust mode
  robust.mode.factor <- function(x){if(any(!is.na(x))) {
    a <- x[which(!is.na(x))] #exclude 
    return((names(sort(table(a), decreasing=T))[1]))} else
      return(NA)        
  }
  
  ### Ancillary function for stretching color scale when mapping species richness
  ### and create breaks
  logbreaks <- function(max.value){
    log2series <- 2^(2:12)
    log2series <- log2series[1:which.min(abs(max.value - log2series))]
    if(length(log2series)>5){
      log2series <- log2series[seq(1, length(log2series), by=2)]
    }
    out <- data.frame(breaks=log2series, labels=log2series) %>%
      mutate(labels=replace(labels, list = breaks==4, values="<4"))
    return(out)
  }
  ## calculate minimum distance between each world point and splot points
  minDistGeo <- function(x, y){return(min(distm(x, y, fun = distGeo)))}
}



## Create labels and indices for subsequent work
all.metrics <- c("sr10", "sr100",  "sr400", "sr1000","sr1ha")
vegtypes <- c( "all")#, "for", "nonfor")
grain <- gsub(all.metrics, pattern = "sr", replacement="")
grain <- gsub(grain, pattern = "ha", replacement="~ha")
grain[1:(length(grain)-1)] <- paste(grain[1:(length(grain)-1)], "~m2", sep="")
grain <- gsub(pattern="m2", replacement = "m^2", x=grain)
selected.predictors <- c("PC1_chelsa", "PC2_chelsa", "PC3_chelsa", "PC4_chelsa", "PC5_chelsa",
                         "PC1_isric", "PC2_isric", "PC3_isric", "PC4_isric", 
                         "CCVPre", "CCVTem", 
                         "tri50km", "landform1km.maj", "landform50km.count", 
                         "sBiomeName", 
                         "sp_wfig",  "REALM", "isforest", "plants_recorded", 
                         "Rel.area")#, "interpl.dist")
var.labs0=c("ClimPC1 (Annual T)","ClimPC2 (Prec)","ClimPC3 (P Season)",
           "ClimPC4 (T warm & wet Q)","ClimPC5 (P coldest Q)",
           "SoilPC1 (Bulk density)", "SoilPC2 (Sand)", "SoilPC3 (Coarse frag)", "SoilPC4 (pH)",
           "CCVelocity (Prec)", "CCVelocity (Temp)", 
           "TRI (50km)", "Landform (1km)", "No. landforms (50km)",
           "Biome", "Ecoregion sp. pool (n)", "Realm", "Forest", "Plants recorded", "Plot size (m²)")#, "Interplot dist.")#,"# plots used", "elevation"))

var.labs.long0 <- c("Climate PC1 - Annual mean temperature",
                   "Climate PC2 - Precipitation",
                   "Climate PC3 - Precipitation seasonality",
                   "Climate PC4 - Temp. warmest/wettest q.",
                   "Climate PC5 - Precip. coldest quarter",
                   "Soil PC1 - Bulk density", 
                   "Soil PC2 - Sand %", 
                   "Soil PC3 - Coarse fragments %", 
                   "Soil PC4 - pH",
                   "Climate Change Velocity - Precipitation", 
                   "Climate Change Velocity - Temperature", 
                   "Terrain Ruggedeness Index", 
                   "Dominant landform", 
                   "Number of landforms (50km)",
                   "Biome", 
                   "Ecoregion species pool", 
                   "Realm", "Forest", "Plants recorded", "Plot size")


# Labelling
biome.labs <- data.frame(x= c("Alpine", "Boreal zone", "Dry midlatitudes",
                              "Dry tropics and subtropics", "Subtropics with winter rain", 
                              "Subtrop. with year-round rain", "Temperate midlatitudes", 
                              "Tropics with summer rain","Tropics with year-round rain", "Polar and subpolar zone"),
                         labels=c("ALP", "BOR", "DML", "DTR","STW","STY", "TEM","TRS","TYR", "POL"))
isforest.labs <- data.frame(x=c("for", "nonfor"), labels=c("FOR", "N-F"))
landform.labs <- data.frame(x=c("flat","peak","ridge","shoulder","spur",
                                "slope","hollow","footslope","valley","pit"), 
                            labels=c("flt", "pk", "rdg", "shl", "spr", "slp", "hlw", "fsl", "vll", "pit"))
plant_recorded.labs <- data.frame(x=c("complete", "woody_all", "woody_large"), 
                                  labels=c("all", "w_a", "w_l"))


# Import world raster 
load(file.path(path.to.input, "world.over.RData"))
world.data <- world.over
rm(world.over)


index.list <- expand.grid(all.metrics, vegtypes)
index.list$Var3 <- factor(index.list$Var1, labels=c(10, 100, 400, 1000, 10000))

require(parallel)
require(doParallel)
cl <- makeForkCluster(5, outfile="mylog.log")
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(matrixStats)
  library(dplyr)
  library(purrr)
  library(gbm)
})

#### Step 1 - summarize output #####
#### get bias-correction regression coefficients
index.list.short <- index.list %>% 
  filter(Var3==10) ## only size 10 has the BRT model stored
foreach(index = 1:nrow(index.list.short)) %dopar% {
  size <- index.list$Var3[index]
  metric <- index.list$Var1[index]
  fornonf <- index.list$Var2[index]
  (listf <- list.files(path.to.BRTglobal, 
                       pattern=paste0("^",fornonf,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T))
  mycoefs <- vector(mode = "list", length = 2)
  for(i in 1:length(listf)){
    load(listf[i])
    ##bias correction - regression as in Zhou et al 2016
    mypredicted <- predict.gbm(modello, type="response")
    mylm <- lm(mypredicted ~ modello$data$y)
    mycoefs[[1]] <- c(mycoefs[[1]], coef(mylm)[[1]])
    mycoefs[[2]] <- c(mycoefs[[2]], coef(mylm)[[2]])
    #Sbcfit=max[(Sfit−a)∕b,0
  }
  save(mycoefs, file=paste(file.path(path.to.world.predictions, "BiasCorrect_regrCoefs"),fornonf, ".RData", sep="_"))
}




## imported predictions for world.raster and summarize across resamplings
foreach(index=1:nrow(index.list)) %dopar% {
  index.list$Var3 <- factor(index.list$Var1, labels=c(10, 100, 400, 1000, 10000))
  size <- index.list$Var3[index]
  metric <- index.list$Var1[index]
  fornonf <- index.list$Var2[index]
  #foreach(fornonf=vegtypes) %dopar% { 
  ### Import data
  (listf <- list.files(path.to.BRTglobal, 
                       pattern=paste0("^",fornonf,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T))
  load(file=paste(file.path(path.to.world.predictions, "BiasCorrect_regrCoefs"),fornonf, ".RData", sep="_"))
  
  for(i in 1:length(listf)){
    load(listf[i])
    a <- mycoefs[[1]][i]
    b <- mycoefs[[2]][i]
    
    if(i==1) {
      out <- matrix(NA, nrow=length(p[[1]]), ncol=length(listf))
    }
    #out[,i] <- p[[1]]  #without bias-corr
    out[,i] <- map_dbl(p[[1]], function(x){max((x-a)/b, 0)}) #with bias corr
    if(fornonf=="all"){
      if(i==1) {
        out2 <- matrix(NA, nrow=length(p[[2]]), ncol=length(listf))
      }
      #out2[,i] <- p[[2]]
      out2[,i] <- map_dbl(p[[2]], function(x){max((x-a)/b, 0)}) #with bias corr
    }
    print(i)
  }
  out.summary <- matrix(NA, nrow=length(p[[1]]), ncol=8, dimnames = list(NULL, c("min", "perc025", "perc25","mean", "median","perc75", "perc975", "max")))
  out.summary[,1] <- rowMins(out)
  out.summary[,c(2,3,5,6,7)] <- rowQuantiles(out, probs = c(0.025,0.25, 0.5, 0.75, 0.975))
  out.summary[,4] <- rowMeans(out, na.rm=T)
  out.summary[,8] <- rowMaxs(out)
  if(fornonf=="all"){
    out.summary2 <- matrix(NA, nrow=length(p[[2]]), ncol=8, dimnames = list(NULL, c("min", "perc025", "perc25","mean", "median","perc75", "perc975", "max")))
    out.summary2[,1] <- rowMins(out2)
    out.summary2[,c(2,3,5,6,7)] <- rowQuantiles(out2, probs = c(0.025,0.25, 0.5, 0.75, 0.975))
    out.summary2[,4] <- rowMeans(out2, na.rm=T)
    out.summary2[,8] <- rowMaxs(out2)
    #not saving invididual results from out and out2
    save(out.summary, out.summary2, file=paste(file.path(path.to.world.predictions, "BRT_out"),fornonf,metric, ".RData", sep="_"))
  } else  save(out.summary, file=paste(file.path(path.to.world.predictions, "BRT_out"),fornonf,metric, ".RData", sep="_"))
  #}
}
stopCluster(cl)





#### Step 2 - BRT predicting #### 
#### Maps of sp richness at different scales & species pool size - Forest vs nonforest  ####

vegtypes <- c("for", "nonfor")
index.list.all <- expand.grid(all.metrics, vegtypes, "all" )
index.list.for <- expand.grid(all.metrics, "for", "for")
index.list.nonfor <- expand.grid(all.metrics, "nonfor", "nonfor")
index.list <-rbind(index.list.all, index.list.for,index.list.nonfor)
save <- F
which.summary.metric <- "median"
which.summary.metric.lab <- ifelse(which.summary.metric!="median", which.summary.metric, "")


## crashes if run as a loop. works one by one
## (Probably for lack of memory)
save.raster <- T
for(index in 1:20)
  {
  metric <- index.list$Var1[index]
  fornonf <- index.list$Var2[index]
  which.model <- index.list$Var3[index]
  mygrain <- grain[which(levels(index.list$Var1)==metric)]
  mygrain0 <- as.numeric(gsub("([0-9]+).*$", "\\1", mygrain))
  
  legpos <- c(0.160, .24) # legend position

  load(paste0(file.path(path.to.world.predictions, "BRT_out_"), which.model, "_", metric, "_.RData"))
  #grass.threshold <- 0
  
  if(which.model=="all") tmp.metric <- ShapeData(out.summary, out.summary2, world.data)
  if(which.model=="for") tmp.metric <- ShapeData(out.for=out.summary, out.nonfor = NULL, world.data)
  if(which.model=="nonfor") tmp.metric <- ShapeData(out.for=NULL, out.nonfor = out.summary, world.data)
  world.data2 <- world.data %>%
    left_join(tmp.metric, by="RAST_ID") %>% 
    filter(!(abs(POINT_X) >172.5 & abs(POINT_Y>60))) # Avoid dealing with pixels close to date line
  
  ##Plot separately for for & nonfor
  world.data2 <- world.data2 %>%
    mutate(value.out = !!rlang::sym(paste("predicted.", fornonf, ".", which.summary.metric, sep="")))
  
  #### build hexagon grid
  dggs          <- dgconstruct(spacing=150, metric=T, resround='down')
  #Get the corresponding grid cells for each earthquake epicenter (lat-long pair)
  world.data2$cell <- dgGEO_to_SEQNUM(dggs, world.data2$POINT_X, world.data2$POINT_Y)$seqnum
#  
#  ### W3 - Map of Species richness 
  {
  #Calculate mean metric for each cell
  world.out   <- world.data2 %>% 
    dplyr::select(cell, value.out) %>% 
    filter(!is.na(value.out)) %>% 
    group_by(cell) %>% 
    summarise(value.out=mean(value.out, na.rm=T), n=n()) %>% 
    filter(n> ifelse(fornonf=="for", 80, 400)) #80 for forest
  #Get the grid cell boundaries for cells 
  grid   <- dgcellstogrid(dggs, world.out$cell, frame=F) %>%
    st_as_sf() %>% 
    mutate(cell = world.out$cell) %>% 
    mutate(value.out=world.out$value.out) %>% 
    st_transform("+proj=eck4") %>% 
    st_wrap_dateline(options = c("WRAPDATELINE=YES"))
  
  ## plotting
  #define new env
  {env <- new.env(parent = globalenv())
    env$w3a <- w3a
    env$countries <- countries
    env$metric <- metric
    env$legpos <- legpos
    env$grid <- grid
    env$mygrain0 <- mygrain0
    env$tominmax <- tominmax
    } 
  
  # create w3 in new environment to reduce file size when saving
    w3 <- with(env, {w3a + 
         geom_sf(data=grid, aes(fill=value.out),lwd=0, alpha=0.9)    +
         #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
         labs(fill= if(!metric%in% c("sr1ha", "Asym.gomp")) { 
           bquote(atop('Sp. Richness', 
                       '('~.(mygrain0) ~ m^{2} ~")"))} else {
                         if(metric=="sr1ha") {
                           bquote(atop('Sp. Richness', '(1 ha)'))} else  {
                             bquote(atop("Sp. pool", "size"))}
                       }
         ) +
         theme(legend.position = legpos +c(-0.06, 0.25)) 
    })
    #w3
    assign(paste("w3", metric, fornonf, sep="."), w3)
    save(list=paste("w3", metric, fornonf, sep="."), 
         file=paste(paste(path.to.pics,which.model,"w3", sep="/"), metric, fornonf, which.summary.metric.lab, "RData", sep="."))
    
  ### W3 - Sp. Richness - same as above but with geom_tile instead of hexagons
    world.data3 <- world.data2 %>% 
      filter(!is.na(value.out))
    
    # save raster
    if(save.raster){
    pred.raster <- rasterFromXYZ(world.data2 %>% 
                    dplyr::select(x=POINT_X, y=POINT_Y, z=value.out), 
                    res=c(NA,NA),
                  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", digits=2)
    writeRaster(x = pred.raster, overwrite=T,
                filename = paste(file.path(path.to.pics,which.model,"w3_tile"), 
                                 metric, fornonf, which.summary.metric.lab, "tif", sep="."))
    }
    #
    world.data3 <- st_as_sf(x = world.data3 %>%
                              dplyr::select(POINT_X, POINT_Y, value.out),                         
                     coords = c("POINT_X", "POINT_Y"),
                     crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
      st_transform("+proj=eck4")
    world.data3 <- world.data3 %>% 
      st_coordinates()%>% 
      as.data.frame() %>% 
      bind_cols(world.data3 %>%
                  dplyr::select(value.out)) %>% 
      rename(x=1, y=2, value.out=3) %>% 
      mutate(value.out=ceiling(value.out+0.01))
    ## log color scale 
    mybreaks <- logbreaks(max(world.data3$value.out))
    w3_tile <- #with(env, {
      w3a + 
        geom_tile(data=world.data3 %>%
                    mutate(value.out=ifelse(value.out<4,4,value.out)), 
                  aes(x=x, y=y, fill=value.out, col=value.out)) +
        xlab(NULL) + ylab(NULL) + 
        #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
        labs(fill= if(!metric%in% c("sr1ha", "Asym.gomp")) { 
          bquote(atop('Sp. Richness', 
                      '('~.(mygrain0) ~ m^{2} ~")"))} else {
                        if(metric=="sr1ha") {
                          bquote(atop('Sp. Richness', '(1 ha)'))} else  {
                            bquote(atop("Sp. pool", "size"))}
                      }
        ) +
        theme(legend.position = legpos +c(-0.06, 0.25)) + 
        scale_color_viridis(guide= "none", trans="log2") + 
        scale_fill_viridis( trans="log2", 
                            breaks=mybreaks$breaks,
                            labels=mybreaks$labels,
                            limits=c(4,
                                     min(max(world.data3$value.out))))
    #})
    
    assign(paste("w3_tile", metric, fornonf, sep="."), w3_tile)
    save(list=paste("w3_tile", metric, fornonf, sep="."), 
         file=paste(file.path(path.to.pics,which.model,"w3_tile"), metric, fornonf, which.summary.metric.lab, "RData", sep="."))
  }


    
  ###W3 minmax - Map of hotspots\coldspots
  {
    ### Map of hotspots
    # Recalculate grid values, not as the max of the mean, but as the max of the original values
    #Calculate mean metric for each cell
    world.data2.minmax   <- world.data2 %>% 
      dplyr::select(cell, POINT_X, POINT_Y, value.out) %>% 
      filter(!is.na(value.out)) %>% 
      mutate(q05=quantile(value.out, 0.05)) %>% 
      mutate(q95=quantile(value.out, 0.95)) %>% 
      mutate(minmax = ifelse(value.out>q95, "> 95%", 
                             ifelse(value.out<q05, "< 5%", 
                             NA))) %>% 
      filter(!is.na(minmax))
    
    
    #Calculate mean metric for each cell
    world.out.minmax   <- world.data2.minmax %>% 
      dplyr::select(cell, minmax) %>% 
      group_by(cell) %>% 
      summarise(minmax=robust.mode.factor(minmax), n=n()) %>% 
      filter(n> ifelse(fornonf=="for", 80, 400)) # reduce salt and pepper
    
    #Get the grid cell boundaries for cells 
    grid.minmax   <- dgcellstogrid(dggs, world.out.minmax$cell, frame=F) %>%
      st_as_sf() %>% 
      mutate(cell = world.out.minmax$cell) %>% 
      mutate(minmax=factor(world.out.minmax$minmax, levels=c("< 5%", "> 95%"))) %>% 
      st_transform("+proj=eck4") %>% 
      st_wrap_dateline(options = c("WRAPDATELINE=YES"))
    
    
    
  w3minmax <- with(env, {w3a +  
      geom_sf(data=grid.minmax,
              aes(fill=minmax),lwd=0)    +
      #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
      scale_fill_brewer(palette="Set1", direction = -1) + 
      labs(fill= if(!metric %in% c("sr1ha", "Asym.gomp")) { 
        bquote(atop('Sp. Richness', 
                    '('~.(mygrain0) ~ m^{2} ~")"))} else {
                      if(metric=="sr1ha") {
                        bquote(atop('Sp. Richness', '(1 ha)'))} else  {
                          bquote(atop("Sp. pool", "size"))}
                    }
      ) +
      theme(legend.position = legpos +c(-0.06, 0.25))
  })
  assign(paste("w3minmax", metric, fornonf, sep="."), w3minmax)
  save(list=paste("w3minmax", metric, fornonf, sep="."), 
       file=paste(file.path(path.to.pics,which.model,"w3minmax"), metric, fornonf, which.summary.metric.lab,"RData", sep="."))
  
  
  
  #### MinMax - Tiles
  world.data3.minmax   <-SpatialPointsDataFrame(coords = world.data2.minmax %>% 
                                                  dplyr::select(x=POINT_X, y=POINT_Y), 
                                                data = world.data2.minmax %>% dplyr::select(value.out, minmax),
                                                proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
    st_as_sf() %>% 
    st_transform("+proj=eck4") 
  world.data3.minmax <- world.data3.minmax %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    bind_cols(world.data3.minmax %>% 
                st_drop_geometry()) %>% 
    rename(x=1, y=2, value.out=3)
  
  w3minmax_tile <- with(env, {w3a +  
      geom_tile(data=world.data3.minmax, aes(x=x, y=y, fill=minmax, col=minmax)) +
      xlab(NULL) + ylab(NULL) + 
      labs(fill= if(!metric %in% c("sr1ha", "Asym.gomp")) { 
        bquote(atop('Sp. Richness', 
                    '('~.(mygrain0) ~ m^{2} ~")"))} else {
                      if(metric=="sr1ha") {
                        bquote(atop('Sp. Richness', '(1 ha)'))} else  {
                          bquote(atop("Sp. pool", "size"))}
                    }
      ) +
      theme(legend.position = legpos +c(-0.06, 0.25)) + 
      scale_fill_brewer(palette="Set1", direction=-1) + 
      scale_color_brewer(palette="Set1", guide = FALSE, direction=-1)
  })
  assign(paste("w3minmax_tile", metric, fornonf, sep="."), w3minmax_tile)
  save(list=paste("w3minmax_tile", metric, fornonf, sep="."), 
       file=paste(file.path(path.to.pics,which.model,"w3minmax_tile"), metric, fornonf, which.summary.metric.lab, "RData", sep="."))
  }
  
  
  ####W4 - Map of Interquartile ranges 
  {
  legpos <- c(0.160, .24)
  if(which.summary.metric=="median"){
    ##### IQR 
    pred25.name <- paste("predicted.", fornonf, ".perc25", sep="")
    predmedian.name <- paste("predicted.", fornonf, ".median", sep="")
    pred75.name <- paste("predicted.", fornonf, ".perc75", sep="")
    world.data2$value.out <- pull(abs(world.data2[, pred25.name] - world.data2[, pred75.name])/world.data2[,predmedian.name] * 100)
    ### GRID!!!!
    #Calculate mean metric for each cell
    world.out   <- world.data2 %>% 
      dplyr::select(cell, value.out) %>% 
      group_by(cell) %>% 
      summarise(value.out=mean(value.out, na.rm=T), n=n()) %>% 
      filter(!is.na(value.out) & n > ifelse(fornonf=="for", 80, 400))
    grid2   <- dgcellstogrid(dggs, world.out$cell, frame=F) %>%
      st_as_sf() %>% 
      mutate(cell = world.out$cell) %>% 
      mutate(value.out=world.out$value.out) %>% 
      st_transform("+proj=eck4") %>% 
      st_wrap_dateline(options = c("WRAPDATELINE=YES"))
    
    {env <- new.env(parent = globalenv())
    env$w3a <- w3a
    env$countries <- countries
    env$metric <- metric
    env$legpos <- legpos
    env$grid2 <- grid2
    env$mygrain0 <- mygrain0
    env$tominmax <- tominmax
      }
    w4 <- with(env, {w3a + 
        geom_sf(data=grid2, aes(fill=value.out),lwd=0, alpha=1)    +
        #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
        scale_fill_viridis(option="plasma") +
        labs(fill= bquote("IQR \n (% of \nmedian)")) +
        theme(legend.position = legpos+c(-0.1, 0.2), 
              legend.key.height = unit(.6, "cm"), 
              legend.key.width = unit(1.0, "cm")) 
    })
    assign(paste("w4", metric, fornonf, sep="."), w4)
    save(list=paste("w4", metric, fornonf, sep="."), 
         file=paste(file.path(path.to.pics,which.model,"w4"), metric, fornonf, "RData", sep="."))
  
    w4minmax <- with(env, {w3a + 
        geom_sf(data=grid2 %>% 
                  mutate(minmax=tominmax(value.out)) %>% 
                  filter(!is.na(minmax)),
                aes(fill=minmax),lwd=0)    +
        #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
        scale_fill_brewer(palette="Set1", direction=-1) + 
        labs(fill= bquote("IQR \n (% of \nmedian)")) +
        theme(legend.position = legpos +c(-0.06, 0.25))
    })
    assign(paste("w4minmax", metric, fornonf, sep="."), w4minmax)
    save(list=paste("w4minmax", metric, fornonf, sep="."), 
         file=paste(file.path(path.to.pics,which.model,"w4minmax"), metric, fornonf, "RData", sep="."))
  }
    #### W4 IQR - same as above but based on tiles
    # save raster
    pred.raster <- rasterFromXYZ(world.data2 %>% 
                                   dplyr::select(x=POINT_X, y=POINT_Y, z=value.out), 
                                 res=c(NA,NA),
                                 crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", digits=2)
    writeRaster(x = pred.raster, overwrite=T,
                filename = paste(file.path(path.to.pics,which.model,"w4_tile"), 
                                 metric, fornonf, which.summary.metric.lab, "tif", sep="."))
    #
    world.data3 <- world.data2 %>% 
      filter(!is.na(value.out))
    world.data3 <- SpatialPointsDataFrame(coords = world.data3 %>% 
                                            dplyr::select(x=POINT_X, y=POINT_Y), 
                                          data = world.data3 %>% dplyr::select(value.out),
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
      st_as_sf() %>% 
      st_transform("+proj=eck4") 
    world.data3 <- world.data3 %>% 
      st_coordinates()%>% 
      as.data.frame() %>% 
      bind_cols(world.data3 %>% 
                  st_drop_geometry()) %>% 
      rename(x=1, y=2, value.out=3)
    w4_tile <- with(env, {w3a + 
        geom_tile(data=world.data3, aes(x=x, y=y, fill=value.out, col=value.out)) +
        #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
        scale_fill_viridis(option="plasma") +
        scale_color_viridis(option="plasma", guide=F) +
        labs(fill= bquote("IQR \n (% of \nmedian)")) +
        theme(legend.position = legpos+c(-0.1, 0.2), 
              legend.key.height = unit(.6, "cm"), 
              legend.key.width = unit(1.0, "cm")) 
    })
    
    
    assign(paste("w4_tile", metric, fornonf, sep="."), w4_tile)
    save(list=paste("w4_tile", metric, fornonf, sep="."), 
         file=paste(file.path(path.to.pics,which.model,"w4_tile"), metric, fornonf, "RData", sep="."))
  }
    
  
  ### W5 - Map of ignorance 
    ## calculate only once per fornonf, as they are identical across metrics
    if(metric==all.metrics[1]){
      load(file.path(path.to.input, "/Mydata_global_NatCommR2.RData"))
    
    if(fornonf=="for"){
      mydata.ign <- bind_rows(world.data2 %>% 
                                filter(isforest==T) %>%
                                dplyr::select(POINT_X, POINT_Y) %>% 
                                mutate(isplot=0),
                              mydata %>% 
                                #filter(isforest==ifelse(fornonf=="for", 1, 0)) %>%
                                dplyr::select(POINT_X, POINT_Y) %>%
                                distinct() %>%
                                mutate(isplot=1)) 
    }
    if(fornonf=="nonfor"){
      mydata.ign <- bind_rows(world.data2 %>% 
                                filter(isgrassland==T) %>%
                                #filter(totgrassland >= grass.threshold) %>%
                                dplyr::select(POINT_X, POINT_Y) %>% 
                                mutate(isplot=0),
                              mydata %>% 
                                #  filter(isforest!=T) %>%   ## changed here from (isgrassland==T) 14.10.2019
                                dplyr::select(POINT_X, POINT_Y) %>%
                                distinct() %>%
                                mutate(isplot=1)) 
    }
    
    ## round to 0.5 degree
    mydata.ign <- mydata.ign %>% 
      mutate_at(.vars=vars(POINT_X:POINT_Y), 
                .funs=list(round05=~(round(.*2,0)/2))) %>% 
      #.funs=~round(.,1)) %>% 
      group_by(POINT_X_round05, POINT_Y_round05) %>% 
      arrange(desc(isplot)) %>% 
      slice(1) %>% 
      ungroup() %>% 
      arrange(desc(isplot))
    ## split splot vs world points
    mydata.ign.all <- mydata.ign %>% 
      filter(isplot==0) %>% 
      dplyr::select(-isplot)
    
    mydata.ign.splot <- mydata.ign %>% 
      filter(isplot==1) %>% 
      dplyr::select(-isplot)

    #transform to list of vectors for future_map
    datL2 <- mydata.ign.all %>% 
      dplyr::select(POINT_X_round05, POINT_Y_round05) %>% 
      mutate(xy=map2(POINT_X_round05, POINT_Y_round05, c)) %>% 
      dplyr::select(xy)
    plan(multisession, workers=5)
    minDist.ign <- future_map_dbl(datL2$xy, 
                         #MARGIN=1, 
                         minDistGeo, 
                         mydata.ign.splot %>% 
                           dplyr::select(POINT_X_round05, POINT_Y_round05)
                         )
    rm(datL2)
    plan(sequential)
    
    
    mydata.ign <- world.data2 %>% 
      dplyr::select(POINT_X, POINT_Y, value.out) %>%
      #filter(abs(POINT_X) <176)  %>% # Avoid dealing with pixels close to date line
      mutate_at(.vars=vars(POINT_X:POINT_Y), 
                .funs=list(round05=~(round(.*2,0)/2))) %>% 
      left_join(mydata.ign.all %>% 
                  bind_cols(data.frame(mindist=minDist.ign)) %>% 
                  bind_rows(mydata.ign.splot %>% 
                              mutate(mindist=0)) %>% 
                  dplyr::select(POINT_X_round05, POINT_Y_round05, mindist),
                by=c("POINT_X_round05", "POINT_Y_round05")) %>% 
      filter(!(abs(POINT_X) >175 & POINT_Y>45)) %>% 
      #get a rid of polygons too close to changing date line
      filter(!(abs(POINT_X) > 160 & POINT_Y < -15 & POINT_Y > -30))
    
    
    # Split in hexagons and calculate distance for each hexagon
    mydata.ign$cell <- dgGEO_to_SEQNUM(dggs, mydata.ign$POINT_X, mydata.ign$POINT_Y)$seqnum
    mydata.ign.out   <- mydata.ign %>% 
      dplyr::select(cell, mindist, value.out) %>% 
      #filter(!is.na(value.out)) %>% 
      group_by(cell) %>% 
      summarise(mindist=median(mindist/1000, na.rm=T), n=n()) %>% 
      filter(!is.na(mindist)) %>% 
      filter(n> ifelse(fornonf=="for", 80, 400)) # reduce salt and pepper
    #Get the grid cell boundaries for cells 
    grid   <- dgcellstogrid(dggs, mydata.ign.out$cell, frame=F) %>%
      st_as_sf() %>% 
      mutate(cell = mydata.ign.out$cell) %>% 
      mutate(mindist=mydata.ign.out$mindist) %>% 
      st_transform("+proj=eck4") %>%  
      st_wrap_dateline(options = c("WRAPDATELINE=YES")) 
    
    env <- new.env(parent = globalenv())
    env$w3a <- w3a
    env$countries <- countries
    env$metric <- metric
    env$legpos <- legpos
    env$grid <- grid
    env$mygrain0 <- mygrain0
    
    w5 <-  with(env, {w3a + 
        geom_sf(data=grid %>% 
                  mutate(mindist=ifelse(mindist>2000, 2000, mindist)), aes(fill=mindist),lwd=0, alpha=1)    +
        #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
        labs(fill= "Dist. from  \n nearest \n plot (km)") + 
        xlab("") + ylab(NULL) +
        theme(legend.position = legpos+c(-0.1, 0.2), 
              legend.background = element_rect(size=0.1, linetype="solid", colour = 1), 
              legend.key.height = unit(.6, "cm"), 
              legend.key.width = unit(1.0, "cm")) +
        scale_fill_gradient(low="#0571b0", high="#ca0020", 
                            limits = c(-5, 2010), 
                            breaks=seq(0,2000, by=500), labels=c(0, 500, 1000, 1500,">2000" ))
    })
    
    ## create raster of ignorance
    e <- extent(mydata.ign %>% 
                  dplyr::select(x=POINT_X, y=POINT_Y))
    r <- raster(e, ncol=720, nrow=360)
    ign.raster <- rasterize(mydata.ign %>% dplyr::select(x=POINT_X, y=POINT_Y), 
                            r, 
                            mydata.ign %>% dplyr::select(mindist), 
                            fun=mean, na.rm=T)
    crs(ign.raster) <- crs("+proj=longlat +datum=WGS84 +no_defs ")
    threshold.ign <- 500000 #500 km
    ign.raster.t <- ign.raster < threshold.ign
    ign.raster.t[ign.raster.t] <- NA
    ### Create polygon where we are ignorant
    pp <- rasterToPolygons(ign.raster.t, dissolve=T)
    pp <- spatialEco::remove.holes(pp)
    
    pp.sf <- pp %>% 
      st_as_sf() %>% 
      st_transform("+proj=eck4") %>% 
      st_cast("POLYGON")  
    pp.sf$area <- as.numeric(st_area(pp.sf))/1000000 # to km2
    pp.sf <- pp.sf %>% 
      filter(area > 100000) #select only polygons larger than 100000 km2
    #### redundant back and forth transformation to fix a bug in one of the polygons
    pp.sf <- as(pp.sf, 'Spatial')
    pp.sf <- spatialEco::remove.holes(pp.sf)
    pp.sf <- pp.sf %>% 
      st_as_sf() %>% 
      st_cast("MULTIPOLYGON")
      
    save(mydata.ign, pp.sf, 
         file=paste0(file.path(path.to.pics, "ignorance_"), which.model, "_", fornonf, ".RData"))
    
    assign(paste("w5", metric, fornonf, sep="."), w5)
    save(list=paste("w5", metric, fornonf, sep="."), 
         file=paste(file.path(path.to.pics,which.model,"w5"), metric, fornonf, "RData", sep="."))
    }
    
  if(save) {
    gglist <- list()
    gglist[[1]] <- w3
    gglist[[2]] <- w3minmax
    gglist[[3]] <- w4
    gglist[[4]] <- w4minmax
    names(gglist) <- paste(metric, fornonf, c("data","data.minmax", "IQR", "IQR.minmax"))
    if(metric==all.metrics[1]){
      gglist[[5]] <- w5
      names(gglist)[5] <- paste(metric, fornonf,"ignorance")}
    save(gglist, file = paste(file.path(path.to.pics, "gglist"), metric, fornonf, "RData", sep="."))
  }
    #return(gglist)
} ; index <- 3

#stopCluster(cl)


#####  RELOAD GROBS for FIGURES #####
which.model <- "all"
index.list0 <- index.list %>% 
  filter(Var3==which.model) %>% 
  # to redue amount of memory needed
  filter(  (Var2 == "for" & Var1 %in% c("sr10", "sr400", "sr1000", "sr1ha")) | 
             (Var2=="nonfor"  & Var1 %in% c("sr10", "sr100" , "sr1000")))
for(i in 1:nrow(index.list0)){
  metric <- index.list0$Var1[i]
  fornonf <- index.list0$Var2[i]
  listw <- list.files(file.path(path.to.pics,which.model), # pattern="^w[0-9].[A-Za-z]*")
                     pattern=paste(paste(c("^w[0-9]+[A-Za-z]*", "^w[0-9]+[A-Za-z]*_tile"), metric, fornonf, paste0(which.summary.metric.lab, "RData$"), sep="\\."),collapse="|"),
                     full.names =T)
  #inelegant patch because of the .. double dots in w3 objects
  listw2 <- list.files(file.path(path.to.pics,which.model), # pattern="^w[0-9].[A-Za-z]*")
                      pattern=paste(paste(c("^w[0-9]+[A-Za-z]*", "^w[0-9]+[A-Za-z]*_tile*"),metric, fornonf, which.summary.metric.lab, "RData$", sep="\\."),collapse="|"),
                      full.names =T)
  lapply(listw, load,.GlobalEnv)
  lapply(listw2, load,.GlobalEnv)
}
load(file=paste0(file.path(path.to.pics, "ignorance_"),which.model,"_for.RData"))
stripes.for <- geom_sf_pattern(data=pp.sf, pattern="stripe", pattern_fill="black", pattern_colour=NA, 
                               fill=NA,colour=NA,pattern_density=0.15,pattern_spacing=0.01)
gray.for <- geom_sf(data=pp.sf, fill="grey90",colour=NA) ## NatComm R2

load(file=paste0(file.path(path.to.pics, "ignorance_"),which.model,"_nonfor.RData"))
stripes.nonfor <- geom_sf_pattern(data=pp.sf, pattern="stripe", pattern_fill="black", pattern_colour=NA, 
                                  fill=NA,colour=NA,pattern_density=0.15,pattern_spacing=0.01)
gray.nonfor <- geom_sf(data=pp.sf, fill="grey90",colour=NA) ## NatComm R2


## New version of plots - 08.05.2020
### Figure 1 Compare scales - Forests ####
tile <- T
if(which.model %in% c("all", "for")){
  leftside <- plot_grid(w3.sr400.for ,
                        w3.sr1000.for,
                        w3.sr1ha.for , ncol=1, labels=c("A","B", "C"))#, "D", "E")) #w3.Asym.gomp.for, 
  central <- plot_grid(w3minmax.sr400.for+ stripes.for, 
                       w3minmax.sr1000.for+ stripes.for, 
                       w3minmax.sr1ha.for+ stripes.for,  ncol=1) #w3minmax.Asym.gomp.for,
  
  Fig1 <- plot_grid(leftside,  central,ncol=2) #rightside, 
  ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig1_for", which.summary.metric.lab, ".png"), 
         Fig1, height=10, width=12, unit="in", dpi=600, bg="white")
  ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig1_for", which.summary.metric.lab, ".pdf"), 
         Fig1, height=10, width=12, unit="in", bg="white")
  if(tile==T){
    leftside <- plot_grid(w3_tile.sr400.for + stripes.for, #+  scale_color_viridis(guide = FALSE) + scale_fill_viridis(), 
                          w3_tile.sr1000.for+ stripes.for, #+ scale_color_viridis(guide = FALSE) + scale_fill_viridis(),  
                          w3_tile.sr1ha.for + stripes.for, #+ scale_color_viridis(guide = FALSE) + scale_fill_viridis(),  
                          ncol=1, labels=c("A","B", "C"))#, "D", "E")) #w3.Asym.gomp.for, 
    central <- plot_grid(w3minmax_tile.sr400.for+ stripes.for, 
                         w3minmax_tile.sr1000.for+ stripes.for, 
                         w3minmax_tile.sr1ha.for+ stripes.for,  ncol=1) #w3minmax.Asym.gomp.for,
    
    Fig1 <- plot_grid(leftside,  central,ncol=2) #rightside, 
    ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig1_for_tile", which.summary.metric.lab, ".png"), Fig1, 
           height=10, width=12, unit="in", dpi=600, bg="white")
    ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig1_for_tile", which.summary.metric.lab, ".pdf"), Fig1, 
           height=10, width=12, unit="in", bg="white")
    
  }
}

### Figure 2 Compare scales - NonForests ####
if(which.model %in% c("all", "nonfor")){
  leftside2 <- plot_grid(w3.sr10.nonfor + stripes.nonfor, 
                         w3.sr100.nonfor+ stripes.nonfor, 
                         w3.sr1000.nonfor+ stripes.nonfor, ncol=1, labels=c("A","B", "C"))#, "d", "e")) #w3.Asym.gomp.for, 
  central2 <- plot_grid(w3minmax.sr10.nonfor+ stripes.nonfor, 
                        w3minmax.sr100.nonfor+ stripes.nonfor, 
                        w3minmax.sr1000.nonfor+ stripes.nonfor,  ncol=1) #w3minmax.Asym.gomp.for,
  Fig2 <- plot_grid(leftside2, central2,  ncol=2) #rightside2,
  ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig2_nonfor", which.summary.metric.lab, ".png"), 
         Fig2, height=10, width=12, unit="in", dpi=600, bg="white")
  ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig2_nonfor", which.summary.metric.lab, ".pdf"), 
         Fig2, height=10, width=12, unit="in", bg="white")
  
  if(tile){
    leftside2 <- plot_grid(w3_tile.sr10.nonfor + stripes.nonfor, 
                           w3_tile.sr100.nonfor+ stripes.nonfor, 
                           w3_tile.sr1000.nonfor+ stripes.nonfor, ncol=1, labels=c("A","B", "C"))#, "d", "e")) #w3.Asym.gomp.for, 
    central2 <- plot_grid(w3minmax_tile.sr10.nonfor+ stripes.nonfor, 
                          w3minmax_tile.sr100.nonfor+ stripes.nonfor, 
                          w3minmax_tile.sr1000.nonfor+ stripes.nonfor,  ncol=1) #w3minmax.Asym.gomp.for,

    Fig2 <- plot_grid(leftside2, central2,  ncol=2) 
    ggsave(filename=paste0(path.to.pics,"/",  which.model, "/Fig2_nonfor_tile", which.summary.metric.lab, ".png"), 
           Fig2, height=10, width=12, unit="in", dpi=600, bg="white")
    ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig2_nonfor_tile", which.summary.metric.lab, ".pdf"), 
           Fig2, height=10, width=12, unit="in", bg="white")
    
  }
}

### NatComms R2 - Replace patterns with gray ####
## Figure 1R2 Compare scales - Forests
if(which.model %in% c("all", "for")){
    leftside <- plot_grid(w3_tile.sr400.for + gray.for + stripes.for, #+  scale_color_viridis(guide = FALSE) + scale_fill_viridis(), 
                          w3_tile.sr1000.for+ gray.for + stripes.for, #+ scale_color_viridis(guide = FALSE) + scale_fill_viridis(),  
                          w3_tile.sr1ha.for + gray.for + stripes.for, #+ scale_color_viridis(guide = FALSE) + scale_fill_viridis(),  
                          ncol=1, labels=c("A","B", "C"))#, "D", "E")) #w3.Asym.gomp.for, 
    central <- plot_grid(w3minmax_tile.sr400.for+ gray.for + stripes.for, 
                         w3minmax_tile.sr1000.for+gray.for + stripes.for, 
                         w3minmax_tile.sr1ha.for+ gray.for + stripes.for,  ncol=1) #w3minmax.Asym.gomp.for,
    
    Fig1R2 <- plot_grid(leftside,  central,ncol=2) #rightside, 
    ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig1R2_for_tile", which.summary.metric.lab, ".png"), Fig1R2, 
           height=10, width=12, unit="in", dpi=600, bg="white")
    ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig1R2_for_tile", which.summary.metric.lab, ".pdf"), Fig1R2, 
           height=10, width=12, unit="in", bg="white")
}
rm(leftside, central, Fig1R2)

### Figure 2 Compare scales - NonForests 
if(which.model %in% c("all", "nonfor")){
    leftside <- plot_grid(w3_tile.sr10.nonfor + gray.nonfor + stripes.nonfor, 
                           w3_tile.sr100.nonfor+ gray.nonfor + stripes.nonfor, 
                           w3_tile.sr1000.nonfor+gray.nonfor + stripes.nonfor, ncol=1, labels=c("A","B", "C"))#, "d", "e")) #w3.Asym.gomp.for, 
    central <- plot_grid(w3minmax_tile.sr10.nonfor+  gray.nonfor+ stripes.nonfor, 
                          w3minmax_tile.sr100.nonfor+ gray.nonfor+ stripes.nonfor, 
                          w3minmax_tile.sr1000.nonfor+gray.nonfor+ stripes.nonfor,  ncol=1) #w3minmax.Asym.gomp.for,
    
    Fig2R2 <- plot_grid(leftside, central,  ncol=2) 
    ggsave(filename=paste0(path.to.pics,"/",  which.model, "/Fig2R2_nonfor_tile", which.summary.metric.lab, ".png"), 
           Fig2R2, height=10, width=12, unit="in", dpi=600, bg="white")
    ggsave(filename=paste0(path.to.pics, "/", which.model, "/Fig2R2_nonfor_tile", which.summary.metric.lab, ".pdf"), 
           Fig2R2, height=10, width=12, unit="in", bg="white")
}



#### Figure S2 IQR - forest vs non Forest ####
### only intermediate grain
{
leftside.s2b <- w4.sr1000.for+ 
  stripes.for +
  ggtitle("Forest") + 
  theme(plot.title=element_text(hjust=0.5),
        axis.title = element_blank())

rightside.s2b <- w4.sr100.nonfor + 
  stripes.nonfor + 
  ggtitle("Non Forest") + 
  theme(plot.title=element_text(hjust=0.5),
        axis.title = element_blank())

if(which.model=="all"){
  FigS2 <- plot_grid(leftside.s2b, 
                     rightside.s2b + 
                       ## add color limits to avoid effect of outlier (R2)
                       scale_fill_viridis(option="plasma", limits=c(12,98)) ,  
                     nrow=2, labels = c("A","B")) #rightside2,
  ggsave(filename=paste0(path.to.pics, "/", which.model, "/FigS2b_IQRs", which.summary.metric.lab,".png"), 
         FigS2, height=9, width=8, unit="in", dpi=600, bg="white")
  }

}


### Figure S11 Maps of ignorance ####

if(which.model=="all"){
rightside.s8 <- plot_grid(w5.sr10.for, # + 
                          w5.sr10.nonfor, #+ 
                       ncol=1, labels=c("A", "B"))
FigS8 <- plot_grid(rightside.s8, ncol=1)
ggsave(filename=file.path(path.to.pics, which.model, "FigS8_ignorance.png"),  
      height=6.5, width=6.5, unit="in", dpi=600, FigS8, bg="white")
}


### Figure 5 bivariate map of congruence across grains #####
### input ideas:
## https://timogrossenbacher.ch/2019/04/bivariate-maps-with-ggplot2-and-sf/
## http://www.joshuastevens.net/cartography/make-a-bivariate-choropleth-map/
## http://lenkiefer.com/2017/04/24/bivariate-map/ # continuous chloropleth
#### Bivariate maps # V2
### ONLY FOR ALL MODEL!!!

#### Shape data into a useful form 
tmp.metric <- list()
for(m in 1:length(all.metrics)){
  metric <- all.metrics[m]
  load(paste0(file.path(path.to.world.predictions, "BRT_out_"), which.model, "_", metric, "_.RData"))
  if(which.model=="for") out.summary2 <- NULL
  if(which.model=="nonfor") {out.summary2 <- out.summary; out.summary2 <- NULL}
  tmp.metric[[m]]  <- ShapeData(out.summary, out.summary2, world.data) %>% 
    dplyr::select(RAST_ID, predicted.for.median, predicted.nonfor.median) %>% 
    rename_at(.vars=vars(starts_with("predicted")), 
              .funs=list(~paste0(., ".", metric, sep="")))
  if(m>1) tmp.metric[[m]] <- tmp.metric[[m]] %>% dplyr::select(-RAST_ID)
}
names(tmp.metric) <- all.metrics

world.data3 <- world.data %>%
  left_join(do.call("cbind", tmp.metric) %>% 
              rename(setNames(names(.), sub(paste(paste0(all.metrics, "\\."), 
                                                         collapse="|"), "", names(.)))),
            by="RAST_ID") %>% 
  filter(!(abs(POINT_X) >172 & abs(POINT_Y>60)))


### Create legend first
## continuos color scheme 4x4  
d <- expand.grid(x=1:4,y=1:4)
w6b.legend <- ggplotGrob(
  ggplot(d, aes(x,y,fill=atan(y/x),alpha=(x+y)))+
    geom_tile()+
    scale_fill_viridis()+ 
    labs(y = "Fine SR ⟶️",
         #expression("Sp. rich." ~ (100*m^2), "⟶️"),#"Higher local sp. richness (100m2) ⟶️",
         x = "Coarse SR ⟶️") +
    theme_minimal() +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.background = element_rect(colour = "black"), 
          legend.position="none") + 
    theme(
      axis.title = element_text(size = 12)
    ) +
    coord_fixed()
)



#### build hexagon grid
dggs          <- dgconstruct(spacing=150, metric=T, resround='down')
#Get the corresponding grid cells for each earthquake epicenter (lat-long pair)
world.data3$cell <- dgGEO_to_SEQNUM(dggs, world.data3$POINT_X, world.data3$POINT_Y)$seqnum

#Calculate mean metric for each cell
world.out   <- world.data3 %>% 
  dplyr::select(predicted.for.median.sr10:cell) %>% 
  group_by(cell) %>% 
  summarize_at(.vars = vars(starts_with("predicted")), 
               .funs = list(~median(., na.rm=T), 
                            n=~sum(!is.na(.)))) %>% 
  rename_all(.funs=list(~gsub(pattern="_median$",  replacement="", x =.))) %>% 
  rename_all(.funs=list(~gsub(pattern="predicted",  replacement="p", x =.))) %>% 
  dplyr::select(-p.for.median.sr10_n, -p.for.median.sr100_n, -p.for.median.sr1000, #-p.for.median.sr1ha_n, 
                -p.nonfor.median.sr10_n, -p.nonfor.median.sr100_n, -p.for.median.sr1000, #-p.nonfor.median.sr1ha_n
                ) %>% # redundant columns
  rename_at(.vars=vars(ends_with("_n")),
            .funs=list(~gsub(pattern=".median.sr1ha", replace="", x=.))) %>% 
  left_join(world.data3 %>% 
              dplyr::select(cell) %>% 
              dplyr::group_by(cell) %>% 
              dplyr::summarize(n=n()),
            by="cell")

#Get the grid cell boundaries for cells 
grid   <- dgcellstogrid(dggs, world.out$cell, frame=F) %>%
  st_as_sf() %>% 
  #mutate(cell = world.out$cell) %>% 
  bind_cols(world.out) %>% 
  st_transform("+proj=eck4") 


##Plot separately for for & nonfor
gglist.bivariate <- list()
tick <- 1
breakss <- seq(0,1, length.out=5)
for(fornonf in c("for", "nonfor")){
  grid2 <- grid %>%
    mutate(local = !!rlang::sym(paste("p.", fornonf, ".median.",ifelse(fornonf=="for", 
                                                                       "sr400", "sr10"), sep=""))) %>% 
    mutate(coarse = !!rlang::sym(paste("p.", fornonf, ".median.",ifelse(fornonf=="for", 
                                                                       "sr1ha", "sr1000"), sep=""))) %>% 
    mutate(coarse=as.numeric((cut(coarse, 
                                     breaks=quantile(coarse, breakss, na.rm=T), 
                                     labels=breakss[-1])))) %>%
    mutate(local=as.numeric((cut(local, 
                                 breaks=quantile(local, breakss, na.rm=T), 
                                 labels=breakss[-1])))) %>% 
    filter(complete.cases(local, coarse) & n>ifelse(fornonf=="for", 80,400)) 
  
  
  #load data on ignorance areas
  load(file=paste0(path.to.pics, "/ignorance_",which.model, "_",fornonf, ".RData"))
                
  
  (w6b <-  w3a + 
      geom_sf(data = countries, fill = "white", col = "grey10", lwd = 0.7) +
      geom_sf(data = countries, fill = "white", col = NA) +
      geom_sf(data=grid2, aes(fill=atan(local/coarse), alpha=(local+coarse)), lwd=0)    +
      geom_sf_pattern(data=pp.sf, pattern="stripe", pattern_fill="black", pattern_colour=NA, 
                      fill=NA,colour=NA,pattern_density=0.15,pattern_spacing=0.01) + 
      theme(legend.position="none", 
            plot.title = element_text(hjust=.5)) +
      scale_fill_viridis() #option="magma"
  )
  
  if(fornonf=="nonfor"){
    (w6b <- w6b + 
       annotation_custom(
         grob = w6b.legend,
         xmin =  -45e6,
         ymax =  -1.2e6,
         ymin =  -8.5e6
       ))
  }
  gglist.bivariate[[tick]] <- w6b
  tick <- tick + 1
}

w6.grid <- plot_grid(gglist.bivariate[[1]] + 
                       ggtitle("Forest"), 
                     gglist.bivariate[[2]] + 
                       ggtitle("Non-Forest"), 
                     ncol=1, labels=c("A", "B"))
ggsave(file.path(path.to.pics,  which.model, "Fig5_World_SpPool_hex_bivariate_discrete.png"), 
       height=8, width=8, units="in", dpi=600, w6.grid, bg="white")
save(gglist.bivariate, file = file.path(path.to.pics, which.model, "gglist.bivariate.Rdata"))








#### Table S2 - Calculate summaries ####
## ONLY for ALL MODEL

summary.global <- world.data3 %>%
  filter(!is.na(sBiomeName)) %>% 
  #group_by(sBiomeName) %>% 
  summarize_at(.vars=vars(predicted.for.median.sr10:predicted.nonfor.median.sr1ha), 
               .funs=list(min=~min(., na.rm=T),
                          #                          q1=~quantile(., 0.25, na.rm=T), 
                          med=~median(., na.rm=T),
                          #                          q3=~quantile(., 0.75, na.rm=T), 
                          IQR=~quantile(., 0.75, na.rm=T)-quantile(., 0.25, na.rm=T),
                          max=~max(., na.rm=T))) %>% 
  gather(key="metrics", value=median.richness) %>% 
  mutate(metrics=gsub(pattern="predicted.", replacement="", x=metrics)) %>% 
  mutate(metrics=gsub(pattern="median.", replacement="", x=metrics)) %>% 
  #mutate(metrics=gsub(pattern="gomp.", replacement="", x=metrics)) %>% 
  mutate(metrics=gsub(pattern="\\.", replacement="_", x=metrics)) %>% 
  separate(metrics, into=c("fornonf", "metrics", "stats"), sep="_") %>% 
  mutate(metrics=fct_relevel(metrics, c("sr10", "sr100","sr400", "sr1000", "sr1ha"))) %>% 
  spread(key=stats, value="median.richness") %>% 
  dplyr::select(fornonf, metrics, min, med, max, IQR)

summary.biome <- world.data3 %>%
  filter(!is.na(sBiomeName)) %>% 
  group_by(sBiomeName) %>% 
  summarize_at(.vars=vars(predicted.for.median.sr10:predicted.nonfor.median.sr1ha), 
               .funs=list(min=~min(., na.rm=T),
                          #                          q1=~quantile(., 0.25, na.rm=T), 
                          med=~median(., na.rm=T),
                          #                          q3=~quantile(., 0.75, na.rm=T), 
                          IQR=~quantile(., 0.75, na.rm=T)-quantile(., 0.25, na.rm=T),
                          max=~max(., na.rm=T))) %>% 
  gather(key="metrics", value=median.richness, -sBiomeName) %>% 
  mutate(metrics=gsub(pattern="predicted.", replacement="", x=metrics)) %>% 
  mutate(metrics=gsub(pattern="median.", replacement="", x=metrics)) %>% 
  #mutate(metrics=gsub(pattern="gomp.", replacement="", x=metrics)) %>% 
  mutate(metrics=gsub(pattern="\\.", replacement="_", x=metrics)) %>% 
  separate(metrics, into=c("fornonf", "metrics", "stats"), sep="_") %>% 
  mutate(metrics=fct_relevel(metrics, c("sr10", "sr100","sr400", "sr1000", "sr1ha"))) %>% 
  spread(key=stats, value="median.richness") %>% 
  dplyr::select(sBiomeName, fornonf, metrics, min, med, max, IQR) %>% 
  arrange(fornonf, metrics)

## format as output table
summary.for <- summary.biome %>% 
  bind_rows(summary.global %>% 
          mutate(sBiomeName="Global")) %>% 
  rename(grain=metrics) %>% 
  filter( (fornonf=="for" & grain %in% c("sr400", "sr1000", "sr1ha"))) %>% 
            #(fornonf=="nonfor" & grain %in% c("sr10", "sr100", "sr1000"))) %>% 
  pivot_longer(names_to = "stat", cols = min:IQR) %>% 
  mutate(value=round(value)) %>% 
  pivot_wider(names_from = grain:stat)#, 
    #id_cols = sBiomeName:stat)

summary.nonfor <- summary.biome %>% 
  bind_rows(summary.global %>% 
              mutate(sBiomeName="Global")) %>% 
  rename(grain=metrics) %>% 
  filter(fornonf=="nonfor" & grain %in% c("sr10", "sr100", "sr1000")) %>% 
  pivot_longer(names_to = "stat", cols = min:IQR) %>% 
  mutate(value=round(value,1)) %>% 
  pivot_wider(names_from = grain:stat)#, 
              #id_cols = sBiomeName:stat)

write_csv(x = summary.for, path=file.path(path.to.output, "Summary.for.csv"))
write_csv(x = summary.nonfor, path=file.path(path.to.output, "Summary.nonfor.csv"))



### Figure S3 - Joint Map Forest NonForests ####
# only for 1000 m2
which.model <- "all"
metric <- "sr1000"
mygrain0 <- 1000
load(paste0(file.path(path.to.world.predictions, "BRT_out_"), which.model, "_", metric, "_.RData"))
grass.threshold <- 50

if(which.model=="all") tmp.metric <- ShapeData(out.summary, out.summary2, world.data)
world.data2 <- world.data %>%
  left_join(tmp.metric %>% dplyr::select(-totgrassland), by="RAST_ID") %>% 
  filter(!(abs(POINT_X) >172.5 & abs(POINT_Y>60)))  %>% # Avoid dealing with pixels close to date line
  #mutate(value.out = ifelse(totgrassland>50, predicted.nonfor.median, predicted.for.median)) 
  mutate(value.out = coalesce(predicted.for.median, predicted.nonfor.median))

{  ### W3 - Sp. Richness
  legpos <- c(0.160, .24)
  world.data3 <- world.data2 %>% 
    dplyr::select(POINT_X, POINT_Y, totgrassland, value.out) %>% 
    filter(!is.na(value.out))
  
  world.data3 <- st_as_sf(x = world.data3 %>%
                            dplyr::select(POINT_X, POINT_Y, value.out),                         
                          coords = c("POINT_X", "POINT_Y"),
                          crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
    st_transform("+proj=eck4") 
  world.data3 <- world.data3 %>% 
    st_coordinates()%>% 
    as.data.frame() %>% 
    bind_cols(world.data3 %>%
                dplyr::select(value.out)) %>% 
    rename(x=1, y=2, value.out=3) %>% 
    mutate(value.out=ceiling(value.out+0.01))
  
  
  ## test color scale by quantiles
  w3_tile_joint <- w3a + 
    geom_tile(data=world.data3 %>% 
                mutate(value.out=ifelse(value.out<4, 4, value.out)), 
              aes(x=x, y=y, fill=value.out, col=value.out)) +
    xlab(NULL) + ylab(NULL) + 
    #geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
    labs(fill= bquote(atop('Sp. Richness', 
                           '('~.(mygrain0) ~ m^{2} ~")"))
    ) +
    theme(legend.position = legpos +c(-0.06, 0.25)) + 
    scale_color_viridis(guide = "none", trans="log2") + 
    scale_fill_viridis( trans="log2", 
                        breaks=logbreaks(max(world.data3$value.out))$breaks,
                        labels=logbreaks(max(world.data3$value.out))$labels,
                        limits=c(max(logbreaks(max(world.data3$value.out))$breaks),
                                 min( max(world.data3$value.out))))
  }

###W3 minmax
{
  ### Map of hotspots
  world.data3.minmax   <- world.data3 %>% 
    mutate(q05=quantile(value.out, 0.05)) %>% 
    mutate(q95=quantile(value.out, 0.95)) %>% 
    mutate(minmax = ifelse(value.out>q95, "> 95%", 
                           ifelse(value.out<q05, "< 5%", 
                                  NA))) %>% 
    filter(!is.na(minmax))
  
  #### MinMax - Tiles
  w3minmax_tile_joint <- w3a +  
    geom_tile(data=world.data3.minmax, aes(x=x, y=y, fill=minmax, col=minmax)) +
    xlab(NULL) + ylab(NULL) + 
    labs(fill= bquote(atop('Sp. Richness', 
                           '('~.(mygrain0) ~ m^{2} ~")"))
    ) +
    theme(legend.position = legpos +c(-0.06, 0.25)) + 
    scale_fill_brewer(palette="Set1", direction=-1) + 
    scale_color_brewer(palette="Set1", guide = FALSE, direction=-1)
  }
  
w3joint <- w3_tile_joint/w3minmax_tile_joint
ggsave(file.path(path.to.pics, 
              which.model, "FigS3_World_SpPool_1000_Joint.png"), 
       height=7, width=6, units="in", dpi=600, w3joint, bg="white")




#### PDPs - Calculate data ######
load(file.path(path.to.input, "Mydata_global_NatCommR1.RData"))
impp.list <- list()
pdp.list <- list()
impp.summary.list <- list()
tick <- 1
for(which.model in c( "all")){#,"nonfor", "for")){
  for(metric in "sr10"){ ## only first needed
    size <- factor(metric, 
                   levels=c("sr10","sr100","sr400","sr1000","sr1ha" ), 
                   labels=c(10,100,400,1000, 10000))
    ### Import data
    (listf <- list.files(path.to.BRTglobal, 
                         #pattern =paste0(which.model, "_", metric, "_[0-9]*.\\.Rdata$"), full.names =T))
                         pattern=paste0("^",which.model,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T))
    out <- list()
    order.i <- as.numeric(gsub(pattern=paste0("_",size,"m\\.RData"), replacement="", 
                               x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))
    ## create list of variable importances across all iterations
    impp.out <- NULL
    for(i in 1:length(order.i)){
      load(listf[i])
      if(i==1){
        var.labs <- data.frame(var.names=modello$gbm.call$predictor.names, 
                               var.labs=var.labs0)#,"# plots used", "elevation"))
        var.lab.df <-data.frame(variable=modello$gbm.call$predictor.names, 
                                label=var.labs$var.labs[match(modello$gbm.call$predictor.names, var.labs$var.names)])
        impp <- modello$contributions
        varr <- modello$gbm.call$predictor.names
        var.order <- order(impp[varr,"rel.inf"], decreasing=T)
        impp.out <- data.frame(impp, iter=order.i[i], n=1:nrow(impp))
      } else {
        impp.out <- rbind(impp.out, 
                          data.frame(modello$contributions, iter=order.i[i],n=1:nrow(impp)))
      }
    }
    ##reorder and summarize impp.out
    impp.summary0 <- impp.out %>%
      group_by(var) %>%
      dplyr::summarize(rel.inf.mean=mean(rel.inf)) %>%
      arrange(desc(rel.inf.mean))
    impp.summary <- as.character((impp.summary0)$var)
    
    if(which.model %in% c("for", "nonfor")){
      var.order <- (left_join(data.frame(var=impp.summary), data.frame(impp, n=1:nrow(impp)), by="var") %>%
                      filter(var!="isforest"))$n
      impp.summary <- impp.summary[-which(impp.summary=="isforest")]
    }
    impp.out$var.labs <- var.labs$var.labs[match(impp.out$var, var.labs$var.names)]
    impp.list[[tick]] <- impp.out
    names(impp.list)[tick] <- paste(metric, which.model)
    impp.summary.list[[tick]] <- impp.summary
    names(impp.summary.list)[tick] <- paste(metric, which.model)
    
    ### extract list of predicted values for each iteration and each variable
    list.final <- list()
    
    for(i in 1:length(order.i)){  
      print(i)
      load(listf[i])
      ordered.iter.labs <- as.numeric(gsub(pattern=paste0("_",size,"m\\.RData"), replacement="", 
                                           x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))
      numvar <- 20
      ## calculate means and ranges for each predictor to create newdata for predictions
      #original.mean.response <- mean(modello$data$y)
      original.data <- reconstructGBMdata(modello) 
      original.ranges <- apply(original.data[,-1], MARGIN=2, "range", na.rm=T)
      original.means <- apply(original.data[,-1], MARGIN=2, "mean", na.rm=T)
      original.means <- data.frame(data.frame(original.means) %>% 
                                     rownames_to_column("var") %>% 
                                     pivot_wider(names_from=var, values_from=original.means))
      original.data <- original.data %>%
        dplyr::select(-y.data) %>% 
        mutate(REALM=factor(REALM, labels=modello$var.levels[[which(modello$var.names=="REALM")]])) %>% 
        mutate(plants_recorded=factor(plants_recorded, labels=modello$var.levels[[which(modello$var.names=="plants_recorded")]])) %>% 
        mutate(sBiomeName=factor(sBiomeName, labels=modello$var.levels[[which(modello$var.names=="sBiomeName")]])) %>% 
        mutate(isforest=factor(isforest, labels=modello$var.levels[[which(modello$var.names=="isforest")]])) %>% 
        mutate(landform1km.maj=factor(landform1km.maj, labels=modello$var.levels[[which(modello$var.names=="landform1km.maj")]])) 

      for(j in 1:numvar){
        var.to.predict.name <- colnames(original.means)[j]
        original.means.j <- original.means %>% 
          dplyr::select(-all_of(j), -Rel.area, -isforest, -REALM, -sBiomeName, -landform1km.maj)
        
        if(var.to.predict.name=="Rel.area") {
          mylength <- 1000 
        } else {
            mylength <- 100
          }
        
        temp <- plot.gbm(modello, i.var=c(j,18, 20), 
                          return.grid = T, 
                          continuous.resolution = mylength)
        if(var.to.predict.name %in% c("Rel.area", "isforest")){
          temp <- temp[,-which(colnames(temp)==var.to.predict.name)[2]] %>% 
            distinct()
        }
        temp <- temp %>% 
          as_tibble() %>% 
          mutate(y=y -mean(y)) %>% 
          #filter(Rel.area %in% c(10,100,400,1000, 10000)) %>% 
          mutate(x= !!rlang::sym(var.to.predict.name)) %>% 
          mutate(var=var.to.predict.name) %>% 
          mutate(which.model=which.model) %>% 
          mutate(iter=ordered.iter.labs[i]) %>% 
          mutate(type=class(x)) %>% 
          mutate(x=as.character(x)) %>% 
          dplyr::select(which.model, iter, 
                        Rel.area, isforest, var, x, y, type)    
        if(var.to.predict.name %in% c("Rel.area", "isforest")){
          temp <- temp %>% 
            mutate( !!var.to.predict.name := NA)
          }
        
        if(j==1){list.final[[i]] <- list()}
        list.final[[i]][[j]] <- temp
        names(list.final[[i]])[j] <- colnames(original.data)[j]
        #jj <- jj+1
        }
      }
#    }
    names(list.final) <-  ordered.iter.labs
    pdp.list[[tick]] <- list.final
    names(pdp.list)[tick] <- paste(metric, which.model)
    print(tick)
    tick <- tick +1
  }
}

save(impp.list, impp.summary.list, var.labs, 
     file = file.path(path.to.pics, which.model, "impp.list.RData"))
save(pdp.list, 
     file = file.path(path.to.pics, which.model, "pdp.list.RData"))




### Figure 4 - PDPs - plotting - new version #####
load(file.path(path.to.pics, which.model, "impp.list.RData"))
load(file.path(path.to.pics, which.model, "pdp.list.RData"))
pdp.df0 <- bind_rows(lapply(pdp.list[[1]], function(x){bind_rows(x)}))

## correct some labeling  
## not neededed when rerunning everything
var.labs$var.labs <- var.labs0


pdp.df0 <- pdp.df0 %>% 
  left_join(biome.labs, by="x") %>% 
  mutate(x=ifelse(is.na(labels), x, as.character(labels))) %>% 
  dplyr::select(-labels) %>% 
  left_join(isforest.labs, by="x") %>% 
  mutate(x=ifelse(is.na(labels), x, as.character(labels))) %>% 
  dplyr::select(-labels) %>%
  left_join(landform.labs, by="x") %>% 
  mutate(x=ifelse(is.na(labels), x, as.character(labels))) %>% 
  dplyr::select(-labels) %>% 
  left_join(plant_recorded.labs, by="x") %>% 
  mutate(x=ifelse(is.na(labels), x, as.character(labels))) %>% 
  dplyr::select(-labels) %>% 
  mutate(isforest=factor(isforest, levels=c("nonfor", "for"), labels=c("Non-For", "For")))


mymetric <- "sr10"
impp.index <- which(names(impp.summary.list)==paste(mymetric, which.model))
impp.m <- impp.list[[impp.index]] %>% 
  group_by(var) %>% 
  dplyr::summarize(imp=median(rel.inf), 
                   q025=quantile(rel.inf, 0.025)) %>% 
  filter(q025>5) %>% 
  arrange(desc(imp)) %>% 
  #filter(!var %in% c("isforest", "Rel.area")) %>% 
  left_join(var.labs, by=c("var"="var.names"))

mypalette <- c("#000000", #black
               "#e7298a", #magenta
               "#e6ab02", #dark yellow
               "#91bfdb")


#plotranges <- mydata %>% 
#  dplyr::select(isforest, Rel.area, all_of(impp.m$var)) %>% 
#  pivot_longer(Rel.area:sp_wfig, values_to="x", names_to="var") %>% 
#  group_by(isforest, var) %>% 
#  summarize(q05=quantile(x, 0.0001, na.rm=T),
#            q95=quantile(x, 0.9999, na.rm=T)) %>% 
#  mutate(which.model=which.model) %>% 
#  mutate(isforest=as.factor(ifelse(isforest, "For", "Non-For"))) %>% 
#  mutate(var=factor(var, levels=impp.m$var, 
#                    labels=impp.m$var.labs))
  

pdp.df <- pdp.df0 %>% 
  replace_na(list(Rel.area=-999)) %>% 
  mutate(var=factor(var, levels=impp.m$var, 
                    labels=impp.m$var.labs)) %>% 
  #filter(!is.na(var)) %>% 
  arrange(which.model, iter, Rel.area, isforest, var, x) %>% 
  left_join({.} %>% 
              filter(type=="numeric") %>% 
              filter(var!="Plot size (m²)") %>% 
              group_by(which.model, iter, Rel.area, isforest, var) %>% 
              summarize(q05=quantile(as.numeric(x), 0.05),
                        q95=quantile(as.numeric(x), 0.95)), 
            by=c("which.model", "iter", "Rel.area", "isforest", "var")) %>%
  filter( var != "Plot size (m²)" | 
             (isforest=="For" & var=="Plot size (m²)" & as.numeric(x)>95 & as.numeric(x)<25100) |
            (isforest=="Non-For" & var=="Plot size (m²)" & as.numeric(x)>7 & as.numeric(x) <1010)) %>% 
  mutate(grain.lab="none") %>% 
  mutate(grain.lab=ifelse( (isforest=="For" & Rel.area==400) | 
                             (isforest=="Non-For" & Rel.area==10), "Fine", grain.lab)) %>% 
  mutate(grain.lab=ifelse( (isforest=="For" & Rel.area==1000) | 
                           (isforest=="Non-For" & Rel.area==100), "Intermediate", grain.lab)) %>% 
  mutate(grain.lab=ifelse( (isforest=="For" & Rel.area==10000) | 
                             (isforest=="Non-For" & Rel.area==1000), "Coarse", grain.lab)) %>%
  mutate(grain.lab=factor(grain.lab, levels=c("none", "Fine", "Intermediate", "Coarse"))) %>% 
  mutate(grain.col=factor(grain.lab, labels=mypalette)) %>% 
  mutate(grain.col=as.character(grain.col))



nvars <- length(impp.m$var)
pdp.gg <- list()
tick <- 1
for(ff in 1:2){
  f <- c("For", "Non-For")[ff]
  pdp.gg[[ff]] <- list()
  for(v in 1:nvars){
    if(f=="For"){
      area.filter <- c(400, 1000, 10000, -999)
    } else {
      area.filter <- c(10, 100, 1000, -999)
    }
    
    tmp.data <- pdp.df %>% 
      mutate(Rel.area=as.factor(Rel.area)) %>% 
      filter(isforest==f & Rel.area %in% area.filter)  %>% 
      filter(var %in% impp.m$var.labs[v]) %>% 
      mutate(x=as.numeric(x)) %>% 
      filter(!is.na(x)) 
    
    if(impp.m$var.labs[v] != "Plot size (m²)") {
      tmp.data <- tmp.data %>% 
        filter((x>q05 & x<q95))
    }
    
    rug.data <- mydata %>% 
      filter(isforest== (f=="For")) %>% 
      dplyr::select(x=all_of(impp.m$var[v])) %>%
      mutate(x=jitter(x)) %>% 
      mutate(y=2.51) #%>% 
      #sample_frac(0.25)
    
    mylim <- range(tmp.data$x)*c(1,1.03)
   if(impp.m$var.labs[v]=="Plot size (m²)" & f=="Non-For"){
      mylim <- c(0, 1010)
   }
    if(impp.m$var.labs[v]=="Plot size (m²)" & f=="For"){
      mylim <- c(0, 25010)
    }

    
    pdp.gg[[ff]][[v]] <- ggplot(data=tmp.data) +
        geom_line(aes(x=x, y=y, 
                      group=interaction(iter, Rel.area), 
                      col=grain.col), 
                  alpha=1/10) + 
        geom_rug(data=rug.data, aes(x=x, y=y), alpha=1/5, sides="b") +
        scale_color_identity() +
        theme_bw() + 
        scale_y_continuous(limits= c(-1, 1.6), breaks = NULL, labels = NULL, name=NULL)  +
        scale_x_continuous(limits = mylim, breaks = scales::pretty_breaks(n = 3), name=NULL) + 
        theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
              axis.title = element_text(size=9),
              axis.text = element_text(size=9, hjust=0.65), 
              panel.grid = element_blank(), 
              plot.margin=margin(l=-0.8, unit="cm"))
        theme(legend.position = "none") 
      if(v==1) {
        pdp.gg[[ff]][[v]] <-pdp.gg[[ff]][[v]] + 
          scale_y_continuous(limits= c(-1,1.5),
                             breaks = scales::pretty_breaks(n = 4), 
                             name="Fitted function")
      }
    if(f=="Non-For") {
      myxlab <- ifelse(impp.m$var.labs[v]=="ClimPC1 - Mean yr T",
                     "ClimPC1 - Mean T",
                     as.character(impp.m$var.labs[v])
      )
      pdp.gg[[ff]][[v]] <-pdp.gg[[ff]][[v]] + 
        scale_x_continuous(limits = mylim, 
                           breaks = scales::pretty_breaks(n = 3), 
                           name=myxlab) + 
        theme(axis.text.x = element_text(hjust=0.65))

    }
  }
}

mylegend <- get_legend(ggplot(data=pdp.df %>% 
         dplyr::filter(iter==1) %>% 
         filter(!is.na(grain.lab)) %>% 
         filter(grain.lab != "none") %>% 
         mutate(grain.lab=factor(grain.lab)))+
  geom_line(aes(x=x, y=y, 
                group=interaction(iter, Rel.area), 
                col=grain.lab)) + 
  scale_color_manual(values=mypalette[-1], drop=F, name="Grain:") +
  theme(legend.position="bottom"))

# now add the title
title.for <- ggdraw() +
  draw_label("Forest",fontface = 'bold',x = 0,hjust = 0.5) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 310)
  )

title.nonfor <- ggdraw() +
  draw_label("Non-Forest",fontface = 'bold',x = 0,hjust = 0.5) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 310)
  )



mywidths <- c(1.02, 1,1,1,1) #c(1, 1,1,1)
myheights <- c(0.05, .38, 0.05, 0.44, 0.08)

pdp.panel <- plot_grid(
  plot_grid(plotlist=list(rep(zeroGrob(),5)), 
                       rel_heights = myheights, ncol=1),
  plot_grid(title.for,
    plot_grid(plotlist=pdp.gg[[1]], nrow=1, align = "v", axis = "b", rel_widths=mywidths),
    title.nonfor,
    plot_grid(plotlist=pdp.gg[[2]], nrow=1, align = "v", axis = "b", rel_widths=mywidths), 
    mylegend, nrow=5, rel_heights = myheights, labels = c("A", "", "B", "", "")),
  plot_grid(plotlist=list(rep(zeroGrob(),5)), 
            rel_heights = myheights, ncol=1),
  ncol=3, rel_widths = c(0.04, 0.92, 0.04))


ggsave(paste0(file.path(path.to.pics, which.model), "/Fig4b_pdp_new_", nvars, "vars_horiz2.png"), 
       width = 9, height = 5, dpi=600, pdp.panel,bg = "white")
ggsave(paste0(file.path(path.to.pics, which.model), "/Fig4b_pdp_new_", nvars, "vars_horiz.pdf"), 
       width = 9, height = 5, pdp.panel, bg = "white")



#### Figure 3 compare metrics with forest-plots ######
var.labs.long <- data.frame(var.names=selected.predictors, 
                            var.labs=var.labs.long0)


impp.gg <- rbindlist(impp.list[1:6], use.names=T, idcol=T) %>% 
  mutate(.id=sapply(strsplit(.id, " "), function(x) x[[1]][1])) %>% 
  dplyr::rename(metric=.id) %>% 
  mutate(metric=factor(metric,levels=all.metrics, labels=grain)) %>% 
  mutate(var.labs=fct_rev(var.labs)) %>% 
  as_tibble() %>% 
  dplyr::select(-var.labs) %>% 
  left_join(var.labs.long, by=c(var="var.names"))

## calculate non-parametric 95% c.i. 
impp.gg.SE <- impp.gg %>% 
  group_by(metric, var, var.labs) %>% 
  dplyr::summarize(rel.inf_mean=mean(rel.inf),
            sd=sd(rel.inf),
            q025=quantile(rel.inf, 0.025),
            rel.inf_median=median(rel.inf),
            q095=quantile(rel.inf, 0.975)) %>% 
#impp.gg.SE <- summarySE(impp.gg, measurevar = "rel.inf", 
#                        groupvars=c( "var.labs"), na.rm=T) %>% #"metric",
  as_tibble() %>% 
  ungroup() %>% 
  mutate_at(.vars=vars(var.labs),
            .funs=list(~factor(.)))
impp.gg.SE <- impp.gg.SE %>% 
  mutate(var.labs=factor(var.labs, levels=(impp.gg.SE %>% 
                                             #filter(metric=="Sp.~pool") %>% 
                                             arrange(rel.inf_median) %>% 
                                             pull(var.labs)))
  )

#order of variables in plotting 
impp.gg <- impp.gg %>% 
  mutate(var.labs=factor(var.labs, levels=levels(impp.gg.SE$var.labs)))



(Fig3.metrics <- ggplot(impp.gg, 
                        aes(x = var.labs, y = rel.inf)) + #, fill = metric
    geom_hline(yintercept = 5, col=gray(0.2), lty=2) + 
    geom_flat_violin(aes(fill = NULL), position = position_nudge(x = 0.0, y = 0), 
                     scale="area", adjust = 2, width=1.5, trim = F, alpha = .5, lwd=0.3, col=NA) +
    geom_point(data = impp.gg.SE, 
               aes(x = as.numeric(var.labs), y = rel.inf_median), shape = 18) +
    geom_errorbar(data = impp.gg.SE, 
                  aes(x = as.numeric(var.labs), y=rel.inf_mean, 
                      ymin=q025, ymax=q095), width = .05) +

    coord_flip() +
    scale_color_brewer(palette="Dark2", label=label_parse) + 
    ylim(c(0,24)) + 
    xlab(NULL) + 
    ylab("Relative Influence (%)") +
    theme_bw() + 
    theme(#legend.title=element_blank(), 
      legend.position=c(0.76, 0.5),
      legend.text.align = 0,
      axis.text.x=element_text(size=8),
      legend.box.background = element_rect(color="black", size=1),
      legend.spacing.y = unit(.01, 'cm'), 
      panel.grid.minor = element_blank()
    ) +  
    guides(colour=guide_legend(ncol=1,byrow=TRUE, title = "Spatial grain"))
)

ggsave(filename=paste0(file.path(path.to.pics, which.model), "/Fig3_compareMetrics_median.png"), 
       height=4, width=5,dpi = 600, unit="in", Fig3.metrics)
ggsave(filename=paste0(file.path(path.to.pics, which.model), "/Fig3_compareMetrics_median.pdf"), 
       height=4, width=5,dpi = 600, unit="in", Fig3.metrics)







