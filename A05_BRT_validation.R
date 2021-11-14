library(tidyverse)
library(gbm)
library(dismo)
library(sp)
library(rgdal)
#library(pals)
library(blockCV)
library(sf)

# set folders for storing temporary files
#rasterOptions(tmpdir="/data/sPlot/users/Francesco/_tmp")
#write("TMPDIR = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('TMPDIR'), '.Renviron'))
#write("R_USER = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))


validate.BRT <- function(mydata.path, 
                         world.data.path, 
                         brt.models.path, 
                         output,
                         iteration, 
                         verbose=T){
  biome.labs <- data.frame(name= c("Alpine", "Boreal zone", "Dry midlatitudes",
                                "Dry tropics and subtropics", "Subtropics with winter rain", 
                                "Subtrop. with year-round rain", "Temperate midlatitudes", 
                                "Tropics with summer rain","Tropics with year-round rain", "Polar and subpolar zone"),
                           labels=c("ALP", "BOR", "DML", "DTR","STW","STY", "TEM","TRS","TYR", "POL"))
  selected.predictors <- c("PC1_chelsa", "PC2_chelsa", "PC3_chelsa", "PC4_chelsa", "PC5_chelsa",
                           "PC1_isric", "PC2_isric", "PC3_isric", "PC4_isric", 
                           "CCVPre", "CCVTem", 
                           "tri50km", "landform1km.maj", "landform50km.count", 
                           "sBiomeName", 
                           "sp_wfig",  "REALM", "isforest", "plants_recorded", "Rel.area")
  
  
  if(verbose){print("load data")}
  # world raster 
  load(world.data.path)

  # load mydata
  load(mydata.path)
  mydata <- mydata %>%
      #mutate(isforest=as.numeric(isforest)) %>% ##classification based on land cover
      mutate(isforest=ifelse(isforest, "for", "nonfor")) %>%
      mutate(isforest=factor(isforest)) %>% 
      mutate(CONTINENT=factor(CONTINENT))
  
  all.metrics <- c("sr10", "sr100", "sr400", "sr1000", "sr1ha")
  grain <- gsub(all.metrics, pattern = "sr", replacement="")
  grain <- gsub(grain, pattern = "ha", replacement="~ha")
  grain[1:(length(grain)-2)] <- paste(grain[1:(length(grain)-2)], "~m2", sep="")
  grain <- gsub(pattern="m2", replacement = "m^2", x=grain)
  
  
  vegtypes <- c( "all")#, "for", "nonfor")
  index.list <- expand.grid("sr10", vegtypes) %>% 
    mutate_all(~as.character(.))
  
  #cor.out1 <- data.frame(NULL)
  #foreach(index=1:nrow(index.list), .combine=rbind) %do% {
  index <- 1
  metric <- "sr10" # index.list$Var1[index]
  size <- 10
  which.model <- index.list$Var2[index]
  #  for(metric in c("sr10", "sr100", "sr400", "sr1ha", "Asymp")){
  listf <- list.files(brt.models.path, 
                      pattern=paste0("^",which.model,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T)
  
  order.i <- as.numeric(gsub(pattern=paste0("_",size,"m\\.RData"), replacement="", 
                             x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))
  
  
  
  
  
  
  if(verbose){print("Re-calculate list of plots to be used in each iteration")}
  ### Relevées used in each iteration for BRTs
  nrows <- 99
  set.seed(999)
  ## stratified by biome, forest, realm, releve area - each group is capped to 100 relevées
  rel.list <- lapply(1:nrows, function(x){mydata %>% 
      group_by(sBiomeName, isforest, REALM, cut(Rel.area, c(0,150, 600, 1200, Inf))) %>% 
      sample_frac(1) %>% 
      slice(1:100) %>% 
      #sample_n(100, replace=T) %>% 
      ungroup() %>% 
      dplyr::select(RELEVE_NR) %>% 
      distinct() %>%
      pull(RELEVE_NR)})
  
  all.used.plots <- unique(unlist(rel.list))
  
  ### Select biomes with >10000 complete plots for forests
  ### 20% of thse plots in these biomes will be transformed to either 
  ### woody_all or woody_large
  sel.biomes <- mydata %>% 
    filter(plants_recorded =="complete") %>% 
    filter(isforest=="for") %>% 
    group_by(sBiomeName) %>% 
    summarize(n=n()) %>% 
    filter(n>10000) %>% 
    mutate(sBiomeName=as.character(sBiomeName)) %>% 
    pull(sBiomeName)
  
  ## For each plot in forest, check whether it is in a well represented biomes, 
  ## if yes, reduce to woody_all or woody_large. Loop over all resamples
  p.rec.list <- lapply(rel.list, function(x){
    mydata %>% 
      filter(RELEVE_NR %in% x) %>% 
      rowwise() %>% 
      mutate(p.rec=ifelse(sBiomeName %in% sel.biomes &
                            isforest=="for" & plants_recorded=="complete", 
                          as.character(base::cut(x = runif(1), 
                                                 breaks=c(0, 0.1, 0.2, 1), 
                                                 labels=c("woody_large", "woody_all", "complete"))),
                          as.character(plants_recorded))) %>% 
      ungroup() %>% 
      dplyr::select(RELEVE_NR, sBiomeName, plants_recorded, p.rec) %>% 
      pull(p.rec)
  })
  
  
  
  
  ### Prepare spatial data for spatial CV
  ### Trasnform world data to raster with 0.1° res
  if(verbose){print("Prepare spatial data for spatial block CV")}
  selected.variables.quantitative <- c("PC1_chelsa", "PC2_chelsa", "PC3_chelsa", "PC4_chelsa", "PC5_chelsa",
                                       "PC1_isric", "PC2_isric", "PC3_isric", "PC4_isric", 
                                       "CCVPre", "CCVTem", 
                                       "tri50km", "landform1km.maj", "landform50km.count")
  
  # Create world rasters of PCA values and extract plot values by geographic intersection
  # raster at half a degree resolution (cf. 30 arc minute resolution)
  rgeo <- brick(nrows=360, ncols=720, xmn=-180, xmx=180, ymn=-90, ymx=90, nl=14) 
  rgeo <- disaggregate(rgeo, fact=12) # raster at 2.5 arc minute resolution
  world.over$cellID <- cellFromXY(rgeo, cbind(world.over$POINT_X, world.over$POINT_Y))
  
  ### create rasters from PCA
  posit <- world.over$cellID
  temp <- getValues(rgeo)
  temp <- as.data.frame(temp)
  
  temp[posit,] <- world.over[, selected.variables.quantitative]
  colnames(temp) <- selected.variables.quantitative
  world.raster <- setValues(rgeo, as.matrix(temp))
  world.raster <- raster::projectRaster(world.raster, crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs") 
    
  rm(world.over)
  
  
  
  #### Calculate cv correlations across models, and biomes ####
  if(verbose){print("Start main loop")}
  cor.out <- NULL
  #cor.out0 <- data.frame(NULL)
  #cor.out <- data.frame(NULL)
  #cor.out0 <- foreach(iteration=1:length(listf), .combine=rbind) %do% {    
#for(iteration in 1:length(listf)){
  i <- order.i[iteration]
  #i <- iteration
  print(paste(which.model, metric, i))
  load(listf[iteration])
  cor.out <- rbind(cor.out, 
                   data.frame(model=which.model, 
                              metric=metric, 
                              cor=modello$cv.statistics$correlation.mean, 
                              var="cv.model", 
                              type="all", 
                              iter=i, 
                              n=nrow(modello$gbm.call$dataframe)))
  
  ## Spatial CrossValidation using BlockCV package
  # Re-create dataset i and create spatial sf object
  print(paste("calculate BlockCV for iteration", iteration))
  mydata.i <- mydata %>%
    filter(RELEVE_NR %in% rel.list[[i]]) %>%
    mutate(p.rec = p.rec.list[[i]]) %>%
    mutate(richness = ifelse(
      p.rec == "complete", rich.complete, ifelse(
        p.rec == "woody_all", rich.woody_all, ifelse(
          p.rec == "woody_large", rich.woody_large, NA)))) %>%
    as.data.frame()
  

  pa_data <- st_as_sf(mydata.i, coords = c("POINT_X", "POINT_Y"), crs = crs(rgeo)) %>% 
    st_transform(., sf::st_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
  

  
  ## Calculate spatial autocorrelation of predictors based on plot locations only
  # Need to exclude plots with NA in world.raster, first
  # mydata.world <- raster::extract(world.raster, pa_data)
  sac <- spatialAutoRange(rasterLayer = world.raster,
                          sampleNumber = 5000,
                          #speciesData=pa_data %>% 
                          #  filter(complete.cases(mydata.world)),
                          doParallel = F,
                          showPlots = TRUE)
  autocorrelation.i <- sac$range
  print(paste("Autocorrelation in run", i, "is", round(autocorrelation.i/1000), "km"))
  
  #plotting
  #pdf(file = "../sPlot/_derived_data/Resample1/_pics/Figv6_autocorrelation.pdf", width = 8, height = 6, bg = "white")
  #plot(sac)
  #dev.off()
  
  # spatial blocking with randomly assigned grid cells, 
  # having the of autocorrelation.i
  sb2 <- spatialBlock(speciesData = pa_data, # presence-background data
                      species = "richness",
                      rasterLayer = world.raster,
                      theRange=autocorrelation.i, ## median of the range of the variograms of the predictors
                      k = 5,
                      iteration=99,
                      selection = "random", 
                      seed=1,
                      progress = T)
  # plotting
  #sb2$plots + 
  #  geom_sf(data = pa_data, alpha = 0.5, cex=0.6)
  #ggsave("../sPlot/_derived_data/Resample1/_pics/Figv7_SpatialBlocks_iteration4.png", width = 8, height=6,dpi = 300, bg = "white")
  #ggsave("../sPlot/_derived_data/Resample1/_pics/Figv7_SpatialBlocks_iteration4.pdf", width = 8, height=6,bg = "white")
  
  rm(world.raster)
  
  # Run BRT on spatialBlocks
  # get model specifications
  
  
  gbm.step.wrapper <- function(lr){gbm.step(
    data = mydata.i,
    gbm.x = modello$gbm.call$predictor.names,
    gbm.y = "richness",
    family = "poisson",
    fold.vector = sb2$foldID,
    n.folds=5,
    step.size = 100,
    tree.complexity = 5,
    learning.rate = lr, 
    bag.fraction = 0.5,
    verbose = F )}
  
  rm(mod.blockCV)
  mylr0 <- modello$gbm.call$learning.rate
  set.seed(15)
  mod.blockCV <- gbm.step.wrapper(mylr0)
  
  mylr <- mylr0
  while(is.null(mod.blockCV) & mylr>0.001){
    set.seed(15)
    mylr <- mylr/2
    print(paste("Retry with smaller lr = ", mylr))
    mod.blockCV <- gbm.step.wrapper(mylr)}
  
  cor.out <- rbind(cor.out, 
                   data.frame(
                     model = which.model,
                     metric = metric,
                     cor = mod.blockCV$cv.statistics$correlation.mean,
                     var = "cv.block",
                     type = "all",
                     iter = i,
                     n = nrow(modello$gbm.call$dataframe)
                   ))
  
  
  mylr0 <- mod.blockCV$gbm.call$learning.rate
  gbm.step.wrapper2 <- function(data, lr, verbose=F)
    {gbm.step(
      data = data,
      gbm.x = modello$gbm.call$predictor.names,
      gbm.y = "richness",
      family = "poisson",
      n.trees = 100,
      step.size = 100,
      n.folds = 2,  # keep number of folds low to decrease computing time
      tree.complexity = 5,
      learning.rate = lr, 
      bag.fraction = 0.5, 
      verbose=verbose)
    }
  ## calculate correlation by continent
  for(cc in levels(mydata.i$CONTINENT)){
    for(fold.i in 1:5){
      print(paste("iteration = ", iteration,"continent =", cc, "foldID = ", fold.i))
      traindata.i.cc <- mydata.i %>% 
        bind_cols(foldID=sb2$foldID) %>% 
        filter(!(foldID==fold.i & CONTINENT == cc)) # use all plots EXCEPT those in continent cc and fold ID i
      
      testdata.i.cc <- mydata.i %>% 
        bind_cols(foldID=sb2$foldID) %>% 
        filter(foldID==fold.i & CONTINENT == cc)
      
     if(nrow(testdata.i.cc)>100){
       mod.cc <- gbm.step.wrapper2(data=traindata.i.cc, lr=mylr0, verbose=F)
       # in case it fails
       mylr <- mylr0
       while(is.null(mod.cc) & mylr>0.001){
         set.seed(15)
         mylr <- mylr/2
         print(paste("Retry with smaller lr = ", mylr))
         mod.cc <- gbm.step.wrapper2(data=traindata.i.cc, mylr)}
        
        if(!is.null(mod.cc)){
          p.valid.cont <- predict(mod.cc, newdata=testdata.i.cc  ,
                                n.trees=mod.cc$gbm.call$best.trees, type="response")
          mycor <- (cor.test(testdata.i.cc$richness, p.valid.cont))$estimate
        } else {mycor <- NA}

      } else {
        print(paste("iteration = ", iteration,"continent =", cc, "foldID = ", fold.i, "Not enough test plots, skip"))
        mycor <- NA}
      cor.out <- rbind(cor.out, 
                           data.frame(model=which.model, 
                                      metric=metric, 
                                      cor=mycor,
                                      var="cv.continent", 
                                      type=cc, 
                                      iter=i, 
                                      n=nrow(testdata.i.cc)))
    }
  }
  
  #calculate cor across biomes
  for(bb in levels(mydata.i$sBiomeName)){
    for(fold.i in 1:5){
      print(paste("iteration = ", iteration, "Biome =", bb, "foldID = ", fold.i))
      traindata.i.bb <- mydata.i %>% 
        bind_cols(foldID=sb2$foldID) %>% 
        filter(!(foldID==fold.i & sBiomeName == bb)) # use all plots EXCEPT those in biome bb and fold ID i
      
      testdata.i.bb <- mydata.i %>% 
        bind_cols(foldID=sb2$foldID) %>% 
        filter(foldID==fold.i & sBiomeName == bb)
      
      if(nrow(testdata.i.bb)>100){
        mod.bb <- gbm.step.wrapper2(data=traindata.i.bb, lr=mylr0)
        # in case it fails
        mylr <- mylr0
        while(is.null(mod.bb) & mylr>0.001){
          set.seed(15)
          mylr <- mylr/2
          print(paste("Retry with smaller lr = ", mylr))
          mod.bb <- gbm.step.wrapper2(data=traindata.i.bb, mylr)}
        
        if(!is.null(mod.bb)){
          p.valid.cont <- predict(mod.bb, newdata=testdata.i.bb  ,
                                  n.trees=mod.bb$gbm.call$best.trees, type="response")
          mycor <- (cor.test(testdata.i.bb$richness, p.valid.cont))$estimate
        } else {mycor <- NA}
        
      } else {
        print(paste("iteration = ", iteration,"biome =", bb, "foldID = ", fold.i, "Not enough test plots, skip"))
        mycor <- NA}
      cor.out <- rbind(cor.out, 
                       data.frame(model=which.model, 
                                  metric=metric, 
                                  cor=mycor,
                                  var="cv.biome", 
                                  type=bb, 
                                  iter=i, 
                                  n=nrow(testdata.i.bb)))
    }
  }
  
  
  
  
  ## subset data for cross validation with data not used in model calibration (neither test nor training)
  valid.data <- mydata %>% 
    filter(RELEVE_NR %in% all.used.plots) %>% 
    filter(!RELEVE_NR %in% rel.list[[i]]) %>%
    sample_n(length(rel.list[[i]])) 
  
  if(which.model %in% c("for", "nonfor")){
    valid.data <- valid.data %>%
      filter(isforest==which.model)
  }
  
  valid.data <- valid.data %>%
######mutate(value.out = !!rlang::sym(as.character(metric)))
  mutate(value.out = ifelse(plants_recorded=="complete", rich.complete, 
                            ifelse(plants_recorded=="woody_all", rich.woody_all,
                                   rich.woody_large)))
    
  if(which.model %in% c("for", "nonfor")){
    valid.data <- valid.data %>% 
      filter(isforest=which.model)
  }
  
  ### calculate overall correlation 
  print("Calculate overall correlation with data from other iterations")
  p.valid <- predict(modello, newdata=valid.data,
                       n.trees=modello$gbm.call$best.trees, type="response")
  cor.out <- rbind(cor.out, 
                   data.frame(model=which.model,
                              metric=metric, 
                              cor=( cor.test(valid.data$value.out, p.valid))$estimate, 
                              var="cv.other", 
                              type=which.model, 
                              iter=i, 
                              n=length(p.valid)))

  #return(cor.out)
    
  #  return(cor.out0)
  #}
  
  #cor.out <- cor.out0
  save(cor.out, file= output)
}

