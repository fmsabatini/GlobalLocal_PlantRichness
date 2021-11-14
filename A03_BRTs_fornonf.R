library(gbm)
library(foreign)
library(data.table)
library(tidyverse)
library(dismo)
library(sp)
library(rgdal)


BRTs_sPlot <- function(mydata, world.data, output, verbose, nrows=99, fornonf, index) {
  selected.predictors <- c("PC1_chelsa", "PC2_chelsa", "PC3_chelsa", "PC4_chelsa", "PC5_chelsa",
                           "PC1_isric", "PC2_isric", "PC3_isric", "PC4_isric", 
                           "CCVPre", "CCVTem", 
                           "tri50km", "landform1km.maj", "landform50km.count", 
                           "sBiomeName", 
                           "sp_wfig",  "REALM", "isforest", "plants_recorded", 
                           #"interpl.dist", 
                           "Rel.area"
                           #"CONTINENT",
                           #"n.plots.area", "elevation"
                           #BIOGEOGRAPHICAL REGION; ...what else? 
                           #"n.plots.area",  "GIVD_ID"
                           )
  
  if (verbose) {
    print("starting to load data ...")
  }
  ##Load data
  load(mydata)
  mydata <- mydata %>%
	  #mutate(isforest=as.numeric(isforest)) %>% ##classification based on land cover
    mutate(isforest=ifelse(isforest, "for", "nonfor")) %>%
    mutate(isforest=factor(isforest)) # %>% 
    ### log transform releve area
    #mutate(Rel.area=log10(Rel.area))
  
  
  if (verbose) {
    print(paste("split into", nrows,"resamples ..."))
  }
  set.seed(999)
  ## create resamples # 2k plots per biome per round
  #rel.list <- lapply(1:nrows, function(x){mydata %>% 
  #    group_by(sBiomeName) %>% 
  #    sample_n(2000, replace=F) %>% 
  #    pull(RELEVE_NR)})
  
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
  
  
  if (verbose) {
    print("Balance plots with different recorded plants")
  }
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
                            isforest=="for" & plants_recorded=="complete" & !is.na(rich.woody_all) & !is.na(rich.woody_large), 
                          as.character(base::cut(x = runif(1), 
                                    breaks=c(0, 0.1, 0.2, 1), 
                                    labels=c("woody_large", "woody_all", "complete"))),
                          as.character(plants_recorded))) %>% 
      ungroup() %>% 
      dplyr::select(RELEVE_NR, sBiomeName, plants_recorded, p.rec) %>% 
      pull(p.rec)
  })

#  require(parallel)
#  require(doParallel)
#  cl <- makeForkCluster(ncores, outfile="" )
#  registerDoParallel(cl)
  
  
#  if (verbose) {
#    print("starting main foreach loop ...")
#  }
  
  if(nrows==0) nrows <- 99
  if(verbose){
  print(paste("run BRT on iteration", index))
  }
  
#  foreach(i = 1:nrows) %dopar% {
  i <- index
    if (verbose)
      print(paste("fitting",  "BRT", i))
    
    ## Create response variable
    ## assign appropriate richness estimation in mydata based on p.rec
    mydata.i <- mydata %>%
      filter(RELEVE_NR %in% rel.list[[i]]) %>%
      mutate(p.rec = p.rec.list[[i]]) %>%
      mutate(richness = ifelse(
        p.rec == "complete", rich.complete, ifelse(
          p.rec == "woody_all", rich.woody_all, ifelse(
            p.rec == "woody_large", rich.woody_large, NA)))) %>%
      as.data.frame()
    
    if (fornonf != "all") {
      mydata.i <- mydata.i %>%
        filter(isforest == fornonf)
    }
    print(head(mydata$RELEVE_NR)) #double check
    
    set.seed(15)
    lrs <-
      c(0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025)
    models.m.all <- list()
    cv.m <- rep(NA, length(lrs))
    ntrees.m <- cv.m
    j <- 1
    while (j <= length(lrs) &
           ifelse(j > 1, ifelse(is.na(ntrees.m[j - 1]), T, ntrees.m[j - 1] < 10000), T)) {
      #print(paste("lr", lrs[j]))
      mod.brt <- gbm.step(
        data = mydata.i,
        gbm.x = selected.predictors,
        gbm.y = "richness",
        family = "poisson",
        tree.complexity = 5,
        learning.rate = lrs[j],
        bag.fraction = 0.5,
        silent = F
      )
      
      cv.m[j] <-
        ifelse(
          length(mod.brt$cv.statistics$correlation.mean) > 0 ,
          mod.brt$cv.statistics$correlation.mean,
          NA
        )
      ntrees.m[j] <-
        ifelse(length(mod.brt$n.trees) > 0 , mod.brt$n.trees, NA)
      models.m.all[[j]] <- mod.brt
      j <- j + 1
    }
    modello <-
      models.m.all[[which(cv.m == max(cv.m[which(ntrees.m != 10000)]))]]
    model.name <- paste("mod_brt_", Sys.Date(), i, sep = "")
    rm(models.m.all)
    
    if(verbose) {print("Load world data")}
    ##make predictions
    load(world.data)
    world.data <- world.over %>%
      # use potential forest from WRI, instead of actual forest from FRA
      #mutate(isforest = ifelse(isforest, "for", "nonfor")) %>%
      #mutate(isforest = factor(isforest))
      mutate(isforest = ifelse(potforest, "for", "nonfor")) %>%
      mutate(isforest = factor(isforest))
    
    rm(world.over)
    
    for (metric in c(10, 100, 400, 1000, 10000)) {
      if (verbose) {
        print(paste("predicting at grain = ", metric))
      }
      
      world.data <- world.data %>%
        mutate(Rel.area = metric) 
      p <- list()
      if (fornonf == "for") {
        p[[1]] <- predict(
          modello,
          newdata = world.data %>%
            filter(isforest == "for"),
          n.trees = modello$gbm.call$best.trees,
          type = "response"
        )
      }
      if (fornonf == "nonfor") {
        p[[1]] <- predict(
          modello,
          newdata = world.data %>%
            filter(isgrassland == T) %>%
            mutate(isforest != "for"),
          n.trees = modello$gbm.call$best.trees,
          type = "response"
        )
      }
      if (fornonf == "all") {
        p[[1]] <- predict(
          modello,
          newdata = world.data %>%
            filter(isforest == "for"),
          n.trees = modello$gbm.call$best.trees,
          type = "response"
        )
        p[[2]] <- predict(
          modello,
          newdata = world.data %>%
            filter(isgrassland == T) %>%
            mutate(isforest != "for"),
          n.trees = modello$gbm.call$best.trees,
          type = "response"
        )
        names(p) <- c("for", "nonfor")
      }
      if (metric == 10) {
        save(modello, model.name, p, cv.m,#models.m.all,
               file = paste(output, "_", metric, "m.RData", sep = ""))
      } else {
        save(model.name, p,
             file = paste(output, "_", metric, "m.RData", sep = ""))
      }
    }
  }
#  stopCluster(cl)
#}

