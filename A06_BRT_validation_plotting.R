# set folders for storing temporary files
rasterOptions(tmpdir="/data/sPlot/users/Francesco/_tmp")
write("TMPDIR = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('TMPDIR'), '.Renviron'))
write("R_USER = /data/sPlot/users/Francesco/_tmp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

library(ggpubr)
library(tidyverse)
library(gbm)
library(dismo)
library(sp)
library(rgdal)
library(matrixStats)
library(viridis)
library(cowplot)
library(pals)
library(blockCV)
library(sf)
library(patchwork)


source("A20_AncillaryFunctions.R")
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


### Set paths
path.to.pics <- "../sPlot/_versions/NatCommR1/_pics"
path.to.BRTglobal <- "../sPlot/_versions/NatCommR1/BRTglobal"
path.to.world.predictions <- "../sPlot/_versions/NatCommR1/_world_predictions"
path.to.output <- "../sPlot/_versions/NatCommR1"
path.to.input <- "../sPlot/_input"
path.to.CValidation <- "../sPlot/_versions/NatCommR1/CValidation/"

### Ancillary functions ####
## extract and bind all BRT predictions into the same object
{
  ShapeData <- function(out.for, out.nonfor, world0, grass.threshold=30){
    tmp.for <- world0 %>%
      filter(isforest==T) %>%
      dplyr::select(RAST_ID) %>%
      bind_cols(data.frame(predicted.for=out.for))
    
    tmp.nonfor <- world0 %>%
      filter(isgrassland==T) %>%
      dplyr::select(RAST_ID, totgrassland) %>%
      bind_cols(data.frame(predicted.nonfor=out.nonfor)) %>% 
      filter(totgrassland>= grass.threshold) %>% ### !!! ### 
      dplyr::select(-totgrassland)
    return(bind_rows(tmp.for, tmp.nonfor)) # Avoids dealing with pixels close to date line
  }
  
  ShapeData2 <- function(x){
    load(x)
    pattern="BRT_out_all_(.*?)_.RData"
    metric <- regmatches(x,regexec(pattern,x))[[1]][2]
    tmp.metric <- ShapeData(out.summary, out.summary2, world.over) %>% 
      dplyr::select(RAST_ID, ends_with(".mean")) %>% 
      rename_at(.vars=vars(starts_with("predicted.")), 
                .funs=list(~gsub(., pattern="mean", replacement=as.character(metric)))) %>%
      rename_at(.vars=vars(starts_with("predicted.")), 
                .funs=list(~gsub(., pattern="predicted.", replacement="p.")))
    return(tmp.metric)}
  
  ## define function to match spatially mydata plots with the NN in world data
  ## to decrease computing time, a previous version of min dist was calculated only WITHIN continent (12 hr with 8 cores!)
  ## now it crops the world.file in 1x1 tiles. (~30 min with 8 cores)
  ## the function returns the RAST_ID of the matching world.over 
  which.minDistGeo <- function(x0, y0){
    require(geosphere)
    lon.rang <- range(x0$POINT_X)  #c(floor(x0[1]), ceiling(x0[1]+0.000001)) ##add a little jitter to avoid problems with round numbers
    lat.rang <- range(x0$POINT_Y)  #c(floor(x0[2]), ceiling(x0[2]+0.000001))
    x <- data.frame(POINT_X=x0[1], POINT_Y=x0[2])
    y <- y0 %>% 
      filter(POINT_X>=lon.rang[1] & POINT_X<lon.rang[2] & 
               POINT_Y>=lat.rang[1] & POINT_Y<lat.rang[2]) %>% 
      dplyr::select(POINT_X, POINT_Y, RAST_ID)
    rastid <- rep(NA, nrow(x))
    for(k in 1:nrow(x)){
      md <- tryCatch(which.min(distm(x[k,], y %>% dplyr::select(-RAST_ID), fun = distGeo)),
                     error = function(e){NA}
      )
      if(!is.na(md)) {
        rastid[k] <- y %>% 
          slice(md) %>% 
          pull(RAST_ID)
      }
    }
    return(data.frame(x, RAST_ID=rastid))
  }
  
  
} # end of ancillary functions


#### Data Preparation ####
## Set labels
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
var.labs0=c("ClimPC1 - Annual T","ClimPC2 - Prec","ClimPC3 - P Season",
            "ClimPC4 - T warm/wet Q","ClimPC5 - P coldest Q",
            "SoilPC1 - Bulk density", "SoilPC2 - Sand", "SoilPC3 - Coarse frag", "SoilPC4 - pH",
            "CCVelocity - Prec", "CCVelocity- Temp", 
            "TRI (50km)", "Landform (1km)", "No. landforms (50km)",
            "Biome", "Ecoregion sp. pool", "Realm", "Forest", "Plants recorded", "Plot size")#, "Interplot dist.")#,"# plots used", "elevation"))

# load world raster 
load(file.path(path.to.input, "world.over.RData"))

# load and adjust mydata
load(file.path(path.to.input, "Mydata_global_NatCommR1.RData"))
mydata <- mydata %>%
  mutate(isforest=ifelse(isforest, "for", "nonfor")) %>%
  mutate(isforest=factor(isforest)) %>% 
  mutate(CONTINENT=factor(CONTINENT))

## create list of grains with labels
all.metrics <- c("sr10", "sr100", "sr400", "sr1000", "sr1ha")
grain <- gsub(all.metrics, pattern = "sr", replacement="")
grain <- gsub(grain, pattern = "ha", replacement="~ha")
grain[1:(length(grain)-2)] <- paste(grain[1:(length(grain)-2)], "~m2", sep="")
grain <- gsub(pattern="m2", replacement = "m^2", x=grain)


vegtypes <- c( "all")#, "for", "nonfor")
index.list <- expand.grid("sr10", vegtypes) %>% 
  mutate_all(~as.character(.))



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


## only needed for block CV
#do.blocks <- F
#if(do.blocks)
{
  ## select biomes with enough samples
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
}



## Import CV correlations and summarize 
allcvs <- list.files(path=path.to.CValidation, pattern = ".RData", full.names = T)
cor.out0 <- list()
for(i in 1:length(allcvs)){
  load(allcvs[[i]])
  cor.out0[[i]] <- cor.out
}

cor.out <- bind_rows(cor.out0) %>% 
  group_by(model, metric, var, type, iter) %>% 
  dplyr::summarize(cor=mean(cor, na.rm=T), 
                   n=median(n, na.rm=T), 
                   n.cors=sum(!is.na(cor)))
rm(cor.out0)


### Fig S4 - Model validation ####
##### V1 - Top - spider plot correlations ###### 
## change labelling
fornonf.lab <- data.frame(label=c("all", "for", "nonfor"),
                          type=c("All", "Forest", "Non Forest"))
cordata2 <- cor.out %>%
  #filter(n>20) %>%
  filter(var != "cv.continent") %>% 
  mutate_at(.vars=vars(var, type),
            .funs=list(~as.character(.))) %>%
  mutate(type=replace(type, list=(var=="cv.model"), values="modelCV")) %>%
  mutate(type=replace(type, list=(var=="cv.block"), values="blockCV")) %>%
  mutate(type=replace(type, list=(var=="cv.other"), values="other")) %>%
  filter(model=="all") %>%
  left_join(biome.labs %>% 
              dplyr::rename(`label`=`labels`, type=name) %>% 
              bind_rows(fornonf.lab) %>%
              mutate_all(~as.character(.)),
            by="type") %>%
  mutate(label=ifelse(is.na(label), type, label)) %>%
  mutate_at(.vars=vars(type, var, label),
            .funs=list(~as.factor(.))) %>%
  mutate(var=factor(var, levels=c("cv.model", "cv.block", "cv.other", "cv.biome"),
                    labels=c("ModelCV", "BlockCV", "Other", "Biome"))) %>% 
  mutate(type=label) %>% 
  mutate(type=fct_rev(type)) %>% 
  mutate(type=fct_relevel(type,   "blockCV","other",  "modelCV"))

cordata2 <- as.data.frame(cordata2)
cordata.SE <- summarySE(cordata2, measurevar = "cor", 
                        groupvars=c("metric", "type"), na.rm=T) 


FigV1.template <- ggplot(cordata2, 
                         aes(x = type, y = cor)) + #, fill = metric
  scale_color_brewer(palette="Set1", label=label_parse) + 
  ylim(c(0,1)) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=8),
        legend.spacing.y = unit(.01, 'cm')
  ) +  
  coord_flip() 


(FigV1.metrics.top <- FigV1.template %+% 
    (cordata2 %>% 
       filter(type %in% c("modelCV","blockCV", "other"))) + 
    geom_point(data = cordata.SE %>% filter(type %in% c("modelCV","blockCV", "other")), 
               aes(x = as.numeric(type), y = cor_mean, 
                   group = type), shape = 18) +
    geom_linerange(data = cordata.SE %>% filter(type %in% c("modelCV","blockCV", "other")) %>% mutate(type=factor(type)), 
                   aes(x = as.numeric(type), y=cor_mean, 
                       ymin=cor_mean-sd, ymax=cor_mean+sd,
                       group = type), width = .05) +
    scale_x_continuous(breaks=1:3, labels = levels(cordata.SE%>% 
                                                     filter(type %in% c("modelCV","blockCV", "other")) %>%  
                                                     mutate(type=factor(type, 
                                                                        levels=c("blockCV", "other", "modelCV"), 
                                                                        labels=c("bCV", "other", "mCV"))) %>% 
                                                     pull(type)),
                       limits = c(0.1, 3.9)) + 
    theme(legend.position="none",
          legend.direction = "horizontal", 
          panel.grid.minor = element_blank()
    )
)

cordata2.tmp <- cordata2 %>% 
  filter(!type %in% c("modelCV","blockCV", "other")) %>% 
  mutate(type=factor(type))
cordata.SE.tmp <- cordata.SE %>% 
  filter(!type %in% c("modelCV","blockCV", "other")) %>% 
  mutate(type=droplevels(type))

(FigV1.metrics.bottom <- FigV1.template %+% (cordata2.tmp) + 
    geom_flat_violin(aes(fill = NULL), position = position_nudge(x = 0.0, y = 0), 
                     scale="area", adjust = 2, width=1.5, trim = F, alpha = .5, lwd=0.3, col=NA) +
    geom_point(data = cordata.SE.tmp, 
               aes(x = as.numeric(type), y = cor_mean, 
                   group = type), shape = 18) +
    geom_errorbar(data = cordata.SE.tmp, 
                  aes(x = as.numeric(type), y=cor_mean, 
                      ymin=cor_mean-sd, ymax=cor_mean+sd,
                      group = type), width = .05) +
    ylab("Correlation") +
    theme(legend.title=element_blank(), 
          legend.position="bottom",
          legend.direction = "horizontal",
          legend.text.align = 0,
          axis.text.x=element_text(size=8),
          #legend.spacing.y = unit(.01, 'cm')
    ) +  
    guides(colour=guide_legend(nrow=2, byrow=TRUE, parse=T))
)

FigV1.cor_panel <- plot_grid(
  plot_grid(NULL, NULL, labels = c("a", "b"), nrow=2,  rel_heights = c(0.28, 0.72)),
  plot_grid(FigV1.metrics.top, FigV1.metrics.bottom, nrow=2,  
            rel_heights = c(0.28, 0.72), align="v", axis = "l"), 
  rel_widths = c(0.05, 0.95))


#### Calculate data for graphs Predicted vs observed 
## comparing predictions of SR, or SAR predictions of individual plots to predictions of respective pixel in world grid
mydata.predobs0 <- mydata %>% 
  as_tibble() %>% 
  dplyr::select(RELEVE_NR, POINT_X, POINT_Y, CONTINENT, Rel.area, plants_recorded, rich.complete:rich.woody_large, sBiomeName, isforest) %>% 
  dplyr::rename(isforest.my=isforest) %>% 
  mutate(RELEVE_NR=factor(RELEVE_NR)) %>% 
  mutate(rich=ifelse(plants_recorded=="complete", rich.complete, 
                     ifelse(plants_recorded=="woody_all", rich.woody_all, 
                            ifelse(plants_recorded=="woody_large", rich.woody_large, NA))))




#### find nearest neighbour between mydata and world.over 


## go parallel
### skip here and reload
require(parallel)
require(doParallel)
ncores <- 12
cl <- makeForkCluster(ncores, outfile="" )
registerDoParallel(cl)

### to reduce computing time, I round coordinates to 1/50 of degree, and only find NN for each distinct combination of coordinates
mydata.predobs.distinct <- mydata.predobs0 %>% 
  mutate(POINT_X=round(POINT_X*50)/50) %>% 
  mutate(POINT_Y=round(POINT_Y*50)/50) %>% 
  distinct(POINT_X, POINT_Y) %>% 
  ## divide in chunks 1x1°
  rowwise() %>% 
  mutate(Class5degree=factor(paste(round(POINT_X/5)*5, round(POINT_Y/5)*5, sep="_"))) %>% 
  ungroup()


system.time(
  MyMinDist.out <- foreach(i=1:nlevels(mydata.predobs.distinct$Class5degree), .combine = rbind) %dopar%{
    tile <- levels(mydata.predobs.distinct$Class5degree)[i]
    MyMinDist.tmp <- which.minDistGeo(mydata.predobs.distinct %>% 
                                        filter(Class5degree==tile), 
                                      world.over)
    return(MyMinDist.tmp)
  }
)
stopCluster(cl)
#user  system elapsed 
#23.231  33.681 655.806 
save(MyMinDist.out, file = file.path(path.to.output, "MyMinDist.out.RData"))

## Reload here
load(file.path(path.to.output, "MyMinDist.out.RData"))
## predict 
mydata.predobs <- mydata.predobs0 %>% 
  mutate(POINT_X.r=round(POINT_X*50)/50) %>% 
  mutate(POINT_Y.r=round(POINT_Y*50)/50) %>% 
  left_join(MyMinDist.out %>% 
              dplyr::select(POINT_X.r=POINT_X, POINT_Y.r=POINT_Y, RAST_ID), 
            by=c("POINT_X.r", "POINT_Y.r")) %>% 
  dplyr::select(-POINT_X.r, -POINT_Y.r) 

### create new dataset to make predictions, based on world.over environmental features, and 
### plot level intrinsic features (e.g., rel.Area, isforest, plants_recorded)

### predict from models themselves
require(parallel)
require(doParallel)
cl <- makeForkCluster(9, outfile="" )
registerDoParallel(cl)
# only takes a couple of minutes

mydata.predobs3 <- NULL
index <- 1
metric <- index.list$Var1[index]
which.model <- index.list$Var2[index]
listf <- list.files(path.to.BRTglobal,
                    pattern=paste0("^",which.model,"BRTs_direct99-[0-9]*-[0-9]*_10m\\.RData$"), full.names =T)
mydata.predobs.tmp <- list()
mydata.predobs3 <- foreach(i=1:length(listf), .combine = rbind) %dopar% {
  mydata.predobs.i <- mydata.predobs %>% 
    filter(RELEVE_NR %in% rel.list[[i]]) %>% 
    #sample_n(10000) %>% 
    dplyr::select(RAST_ID, Rel.area, plants_recorded, isforest=isforest.my, obs.rich=rich) %>% 
    left_join(world.over %>% dplyr::select(-isforest, -plants_recorded), by="RAST_ID") %>% 
    dplyr::select(RAST_ID, POINT_X, POINT_Y, obs.rich, all_of(selected.predictors)) %>% 
    mutate(which.model=which.model)
  
  load(listf[i])
  world.over.i <- mydata.predobs.i
  world.over.i$iter <- as.numeric(gsub(pattern=paste0("_10m\\.RData"), replacement="", 
                                       x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))[i]
  world.over.i$pred.rich <- predict.gbm(object = modello, 
                                        newdata = world.over.i, type = "response")              
  mydata.predobs.tmp[[i]] <- world.over.i
  return(bind_rows(mydata.predobs.tmp))
}
stopCluster(cl)

mydata.predobs3 <- mydata.predobs3 %>% 
  mutate(sBiomeName=factor(sBiomeName, levels=biome.labs$name,
                           labels=biome.labs$labels)) %>% 
  mutate(isforest=factor(isforest, levels=c("for", "nonfor"), labels=c("Forest", "Non Forest"))) %>% 
  mutate(metric=factor(metric, levels=all.metrics, labels=grain)) %>% 
  distinct() ##### Some plots are saved multiple times - Temporary patch to solve

save(mydata.predobs3, 
     file = file.path(path.to.output, "PredObs_Validation.RData"))




#### Plotting Residuals 
load(file.path(path.to.output, "PredObs_Validation.RData"))
## Add BiasCorrection - load regression coefficient to correct BRT prediction bias
load(file.path(path.to.world.predictions, "BiasCorrect_regrCoefs_all_.RData"))

mydata.predobs3 <- mydata.predobs3 %>% 
  dplyr::rename(pred.rich.raw=pred.rich) %>% 
  mutate(pred.rich.bc = (pred.rich.raw-mycoefs[[1]][iter])/mycoefs[[2]][iter]) %>% 
  mutate(pred.rich.bc = ifelse(pred.rich.bc<0, 0, pred.rich.bc))



mysample <- 1 #sample(1:99, 10)
## Which predicted values should I plot? Raw or bias.corrected?
which.pred <- "pred.rich.bc" # "pred.rich.raw" # 

mydata.predobs3 <- mydata.predobs3 %>% 
  dplyr::mutate(pred.rich = !!rlang::sym(which.pred))

mylimits <- c(0,250)
mycuts <- c(0, 20, 150, 500, 1200, Inf)
cutlabels <- c(10, 100, 400, 1000, 10000)



#### V3 - Center - BRTs vs OBS (Global & Biomes) ####
ggpredobs <- list()
(ggpredobs[[1]] <- ggplot(data=mydata.predobs3 %>% 
                            filter(which.model=="all") %>% 
                            filter(iter %in% mysample)) + 
    geom_point(aes(x=obs.rich, y=pred.rich, color=isforest), alpha=0.3, cex=0.6) + 
    geom_abline(slope=1, intercept=0, lty=2, lwd=0.5)+
    scale_color_brewer(palette="Set1", direction=-1) + 
    ylab("Predicted species richness") +
    xlab("Observed species richness") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits) +
    theme_bw() + 
    theme(legend.position = c(0.8,0.85), 
          legend.title = element_blank(), 
          legend.box.background = element_rect(colour = "black")) + 
    guides(color = guide_legend(override.aes= list(alpha = 1, cex=1.5))) + 
    coord_fixed()
)

for(br in 1:nlevels(mydata.predobs3$sBiomeName)){
  bi <- levels(mydata.predobs3$sBiomeName)[br]
  ggpredobs[[br+1]] <- ggpredobs[[1]] %+% (mydata.predobs3%>% 
                                             filter(which.model=="all") %>% 
                                             filter(sBiomeName==bi) %>% 
                                             filter(iter %in% mysample)) + 
    geom_text(data=data.frame(x=max(mylimits)*0.45, y=max(mylimits)*0.85, label=bi), 
              aes(x, y, label=label)) + 
    theme(axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  if(br %in% 1:4){
    ggpredobs[[br+1]] <- ggpredobs[[br+1]] + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits, position="right") + 
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
  }
  if(br == 5){
    ggpredobs[[br+1]] <- ggpredobs[[br+1]] + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits, position="right") + 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits) 
  }
  if(br %in% 6:9){
    ggpredobs[[br+1]] <- ggpredobs[[br+1]] + 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits, position="bottom") + 
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
  }
}


rightside.V3b <- plot_grid( 
  ggpredobs[[2]] + theme(axis.text = element_blank()), 
  ggpredobs[[3]] + theme(axis.text = element_blank()),  
  ggpredobs[[4]] + theme(axis.text = element_blank()),  
  ggpredobs[[5]] + theme(axis.text = element_blank()),  
  ggpredobs[[6]] + theme(axis.text = element_blank()),
  NULL,
  ncol=1, 
  rel_heights = c(rep(0.17,5), 0.025))#, labels=c("b", "c", "d", "e"))
bottomside.V3b <- plot_grid(NULL, 
                            ggpredobs[[7]] + theme(axis.text = element_blank()), 
                            ggpredobs[[8]]+  theme(axis.text = element_blank()), 
                            ggpredobs[[9]]+ theme(axis.text = element_blank()), 
                            ggpredobs[[10]]+ theme(axis.text = element_blank()), 
                            ncol=5, rel_widths = c(0.1, rep(0.225,4)))#,
(main2 <- plot_grid(plot_grid(ggpredobs[[1]] +
                                geom_text(aes(x=max(mylimits)*0.5, y=max(mylimits)*0.95, label="Global")) +
                                theme(legend.position = "none", 
                                      axis.title=element_text(size=8)) + 
                                ylab("Predicted species richness") + 
                                xlab("Observed species richness"),
                              bottomside.V3b,
                              nrow = 2, rel_heights=c(0.75, 0.25)),
                    rightside.V3b, ncol=2, 
                    rel_widths =c(.77,.23)))


######## V9 - Bottom - error (i.e. pred-observed) ########
#### graph predicted-observed vs plot size
ggerror <- list()
(ggerror[[1]] <- ggplot(data=mydata.predobs3 %>% 
                          filter(which.model=="all") %>% 
                          filter(iter %in% mysample) %>% 
                          filter(!is.na(sBiomeName)) %>% 
                          mutate(metric=cut(Rel.area, 
                                            breaks=mycuts, 
                                            labels=c("10~m^2", "100~m^2", "400~m^2", "1000~m^2",  "1~ha")))
) +
    geom_point(aes(x=Rel.area, y=pred.rich-obs.rich, col=isforest), alpha=0.3, cex=0.6) + 
    geom_abline(slope=0, intercept=0, lty=2) +
    scale_color_brewer(palette="Set1", direction=-1) + 
    scale_y_continuous(limits=c(-200,200)) +
    scale_x_log10(limits = c(10, 25010)) +
    xlab(parse(text="Plot~size~(m^2)")) +
    ylab("Pred. - obs. richness") + 
    theme_bw() + 
    theme(legend.position = c(0.17,0.13), 
          legend.title = element_blank(), 
          legend.box.background = element_rect(colour = "black")) + 
    guides(color = guide_legend(override.aes= list(alpha = 1, cex=1.5)))
)

## create separate legend

predobs.legend <- cowplot::get_legend(
  ggerror[[1]] + 
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.box.background = element_rect(color = NA)))
  
  
#### Panel V1-V5-V9  #########
V5V9 <- cowplot::plot_grid(plot_grid(FigV1.metrics.top + 
                                       theme(plot.margin = margin(r=0.5, t=0.5, unit="cm"), 
                                             axis.text.x = element_blank()), 
                                     FigV1.metrics.bottom + 
                                       guides(colour=guide_legend(parse=T, byrow=T, nrow=3))+
                                       theme(axis.title=element_text(size = 8), 
                                             plot.margin = margin(r=0.5, t=0.4, unit="cm")),
                                     nrow=2, rel_heights = c(0.25, 0.75), 
                                     align="v", axis = "l",labels = c("A", "B"), vjust = 1),
                           main2, 
                           ggerror[[1]] + 
                             theme(legend.position = "none", 
                                   text = element_text(size=9), 
                                   plot.margin = margin(l=0.2, r=0.5, b=0.4, unit="cm")
                             ),
                           predobs.legend,  
                           nrow=4, align="h", axis="l", 
                           rel_heights = c(.52, .7, .4,.03),
                           labels=c("", "C", "D"), vjust=1)
ggsave(filename=paste0(path.to.pics, "/FigS4_Error_", which.pred, ".png"), 
       width=4.6, height=10, dpi=600, units="in", V5V9, bg="white")
ggsave(filename=paste0(path.to.pics, "/FigS4_Error_", which.pred, ".pdf"), 
       width=4.6, height=10,  units="in", V5V9, bg="white")

### Additional Plots Not in Paper ####
# #### Fig V2 - BRT vs Obs - Grains 
# ##### alternative for pred obs to circumvent the locked aspect ratio in facet_grid
# ##### Horizontal version - Not in paper
# 
# ggpredobs3b <- list()
# index2.list <- data.frame(isforest=rep(c("Forest", "Non Forest"), each=3), 
#                           grain=(c(400, 1000, 10000, 10, 100, 1000)), 
#                           graintoparse=c("400~m^2", "1000~m^2", "1~ha", "10~m^2", "100~m^2",  "1000~m^2"))
# 
# for(m in 1:nrow(index2.list)){
#   mygrain <- index2.list$grain[m]
#   mygraintoparse <- as.character(index2.list$graintoparse[m])
#   fornonf <- index2.list$isforest[m]
#   tmp.data <- mydata.predobs3 %>% 
#     filter(which.model=="all") %>%
#     mutate(metric=cut(Rel.area, breaks=mycuts, labels=cutlabels)) %>% 
#     filter(metric==mygrain) %>% 
#     filter(isforest==fornonf) %>% 
#     filter(iter %in% mysample)
#   lim0 <- c(0,
#             ceiling((max(tmp.data %>% 
#                            dplyr::select(obs.rich, pred.rich) %>% pull()))/10)*10) 
#   ggpredobs3b[[m]] <- ggplot(data=tmp.data) + 
#     geom_point(aes(x=obs.rich, y=pred.rich, color=plants_recorded), cex=0.6, alpha=1/5) + 
#     geom_abline(slope=1, intercept=0, alpha=.5, cex=0.3, linetype=2)+
#     ggtitle(parse(text = mygraintoparse)) + 
#     scale_color_brewer(palette="Pastel1", drop=FALSE) + 
#     scale_x_continuous(breaks = scales::pretty_breaks(n = 2), limits = lim0, name=NULL) +
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 2), limits = lim0, name=NULL) +
#     theme_bw() + 
#     theme(legend.position = "none", 
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), 
#           strip.text.y = element_blank() , 
#           strip.background = element_blank()
#     ) + 
#     guides(color = guide_legend(override.aes= list(alpha = 1, cex=1.5))) + 
#     coord_equal()
# }
# 
# ggpredobs3b[[m+1]] <- get_legend( ggpredobs3b[[m]] + 
#                                     theme(legend.position = "bottom",
#                                           legend.direction = "horizontal",
#                                           legend.title = element_blank()))
# y.grob <- text_grob(label="Predicted SR", rot = 90)
# x.grob <- text_grob(label="Observed SR")
# 
# 
# FigV2_Predobs_h <- plot_grid(y.grob, plot_grid( ggpredobs3b[[m+1]],
#                                                 plot_grid(plotlist=ggpredobs3b[1:m], nrow=2, align="v", axis="t"),
#                                                 x.grob,
#                                                 nrow=3, rel_heights = c(0.1, 0.8, 0.1)), rel_widths = c(0.05,.95))
# ggsave(filename=paste0("../sPlot/_versions/NatCommR1/_pics/FigV2_ObsBRT_grains_horizontal_", which.pred, ".png"), 
#        width=4.5, height=4, dpi=300, units="in", FigV2_Predobs_h, bg="white")
# 
# 
# #### Fig V2b - BRT vs OBS - grains VERTICAL 
# ## A better horizontal version of this graph is below. This vertical is meant for the panel of all validation graphs
# # not in paper
# ggpredobs3 <- list()
# ggpredobs3[[1]] <- ggplot(data=mydata.predobs3 %>% 
#                             filter(which.model=="all") %>% 
#                             filter(iter %in% mysample) %>% 
#                             mutate(metric=cut(Rel.area, 
#                                               breaks=mycuts, 
#                                               labels=c("10~m^2", "100~m^2", "400~m^2", "1000~m^2",  "1~ha")))
# ) + 
#   geom_point(aes(x=obs.rich, y=pred.rich, color=isforest), cex=0.6, alpha=1/10) + 
#   geom_abline(slope=1, intercept=0, alpha=.5, cex=0.3, linetype=2)+
#   #stat_smooth(geom='line', aes(x=predicted.sar, y=predicted.brt, group=iter), 
#   #            method="lm", se=F, na.rm=T, linetype=1, alpha=1/5,
#   #            lwd=0.5, col=1, fullrange=T) + 
#   scale_color_brewer(palette="Set1", direction = -1) + 
#   xlab("Observed species richness") +
#   ylab("Predicted species richness") + 
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 2), limits = mylimits) +
#   theme_bw() + 
#   theme(legend.position = "top",
#         legend.direction = "horizontal",
#         legend.title = element_blank(),
#         strip.background = element_blank(),
#         strip.text.y.right = element_text(angle=0, hjust = 0),
#         panel.border = element_rect(colour = "black"), 
#   ) + 
#   facet_grid(metric~., labeller = labeller(metric = label_parsed)) + 
#   guides(color = guide_legend(override.aes= list(alpha = 1, cex=1.5))) + 
#   coord_equal()
# 
# 
# ggpredobs3.legend <- cowplot::get_legend(ggpredobs3[[1]])
# ggsave(filename=paste0("../sPlot/_versions/NatCommR1/_pics/FigV2b_SarBRT_grains_vert_", which.pred, ".png"), 
#        width=2.5, height=6, dpi=300, units="in", ggpredobs3[[1]], bg="white")
# 


### Standalone version of V3 with axis labels
# rightside.V3 <- plot_grid(NULL, 
#                           ggpredobs[[2]], 
#                           ggpredobs[[3]], 
#                           ggpredobs[[4]], 
#                           ggpredobs[[5]], 
#                           ggpredobs[[6]], ncol=1, rel_heights = c(0.0015, rep(0.17,4), .2))#, labels=c("b", "c", "d", "e"))
# bottomside.V3 <- plot_grid(NULL, 
#                            ggpredobs[[7]], 
#                            ggpredobs[[8]], 
#                            ggpredobs[[9]],
#                            ggpredobs[[10]], ncol=5, rel_widths = c(0.1, rep(0.225,4)))#,
# 
# (main <- plot_grid(plot_grid(ggpredobs[[1]] +
#                                geom_text(aes(x=max(mylimits)*0.5, y=max(mylimits)*0.95, label="Global")),
#                              bottomside.V3,
#                              NULL, nrow = 3, rel_heights=c(0.79, 0.21, 0.0095)),
#                    rightside.V3, ncol=2, 
#                    rel_widths =c(.78,.22)))

#ggsave(filename=paste0(path.to.pics, "/FigV3_ObsPred_biomes_", which.pred, ".png"), 
#       width=5, height=5, dpi=300, units="in", main, bg="white")
#ggsave(filename=paste0(path.to.pics, "/FigV3_ObsPred_biomes_", which.pred, ".pdf"), 
#       width=5, height=5, units="in", main, bg="white")


### Fig V3c - BRTs vs OBS - compared raw vs Bias-corrected for global
# ggpredobs.raw <- ggpredobs[[1]] %+% (mydata.predobs3 %>% 
#                                        filter(iter==mysample) %>% 
#                                        mutate(pred.rich = pred.rich.raw)) + 
#   theme(legend.position="none")
# ggpredobs.bc <-  ggpredobs[[1]] %+% (mydata.predobs3 %>% 
#                                        filter(iter==mysample) %>% 
#                                        mutate(pred.rich = pred.rich.bc)) #+ 
# ggpredobs.compare <- ggpredobs.raw + ggpredobs.bc + 
#   plot_annotation(tag_levels = 'A')
# #ggsave(filename=file.path(path.to.pics, "FigV3c_ObsBRT_compare_nolm.png"), 
# #       width=7, height=4, dpi=600, units="in", ggpredobs.compare)
# #ggsave(filename=file.path(path.to.pics, "FigV3c_ObsBRT_compare_nolm.pdf"), 
# #       width=7, height=4, units="in", ggpredobs.compare)
# 
# 

#### Fig V4 - BRT vs OBS - Grains x Biomes
# # Not in paper
# ggpredobs3[[2]] <- ggplot(data=mydata.predobs3 %>% 
#                             filter(which.model=="all") %>% 
#                             filter(iter %in% mysample) %>% 
#                             filter(!is.na(sBiomeName)) %>% 
#                             mutate(metric=cut(Rel.area, 
#                                               breaks=mycuts, 
#                                               labels=c("10~m^2", "100~m^2", "400~m^2", "1000~m^2",  "1~ha")))
# ) + 
#   geom_point(aes(x=obs.rich, y=pred.rich, color=isforest), cex=0.6, alpha=1/10) + 
#   geom_abline(slope=1, intercept=0, alpha=.5, cex=0.3, linetype=2)+
#   #    stat_smooth(geom='line', aes(x=predicted.sar, y=predicted.brt, group=iter), 
#   #                method="lm", se=F, na.rm=T, linetype=1, alpha=1/5,
#   #                lwd=0.5, col=1, fullrange=T) + 
#   scale_color_brewer(palette="Set1", direction=-1) + 
#   coord_equal() + 
#   xlab("Observed species richness") +
#   ylab("Predicted species richness") + 
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 1), limits = mylimits) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 1), limits = mylimits) +
#   theme_bw() + 
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.title = element_blank(),
#         strip.background = element_blank(),
#         panel.border = element_rect(colour = "black")
#   ) + 
#   facet_grid(metric~sBiomeName, labeller = labeller(metric = label_parsed)) + 
#   guides(color = guide_legend(override.aes= list(alpha = 1, cex=1.5)))
# 
# 
# ggsave(filename=paste0("../sPlot/_versions/NatCommR1/_pics/FigV4_ObsBRT_Grains_x_Biomes_",which.pred, ".png"), 
#        width=7, height=4, dpi=600, units="in", ggpredobs3[[2]])
# ggsave(filename=paste0("../sPlot/_versions/NatCommR1/_pics/FigV4_ObsBRT_Grains_x_Biomes_",which.pred, ".pdf"), 
#        width=7, height=4, units="in", ggpredobs3[[2]])
# 
# 


#### V5 -  Panel of all validation plots 

#fig4.legend <- get_legend(Fig4.metrics.bottom + 
#             theme(legend.position = "right", 
#                   legend.direction = "vertical") + 
#             guides(colour=guide_legend(parse=T, byrow=F)))
#
# rightside.up.V5 <- plot_grid(NULL, plot_grid(FigV1.metrics.top, FigV1.metrics.bottom + 
#                                                guides(colour=guide_legend(parse=T, byrow=T, nrow=3))+
#                                                theme(axis.title=element_text(size = 8)),
#                                              nrow=2, rel_heights = c(0.28, 0.72), 
#                                              align="v", axis = "l",labels = c("B", "C"), hjust=c(2.3,2.4)),
#                              NULL, ncol=3, 
#                              rel_widths=c(0.1, .8, 0.1))
# rightside.V5 <- plot_grid(rightside.up.V5,  
#                           main2, 
#                           ggpredobs3.legend,  
#                           nrow=3, align="h", axis="l", 
#                           rel_heights = c(.45, .5, .1),
#                           labels=c("", "D"), hjust=1)
# 
# panel.V5 <- plot_grid(ggpredobs3[[1]] + 
#                         theme(axis.title=element_text(size = 8), 
#                               legend.position="none"), 
#                       rightside.V5,
#                       nrow=1, ncol=2, rel_widths = c(0.4,0.6) ,labels=c("A", ""))





#### _______________  ####
### Fig S14 - Map of residuals ####
library(sf)
library(dggridR)
source("A21_w3a_TemplateEckert.R")

mydata.predobs3.i <- mydata.predobs3 %>% 
  filter(iter==mysample) %>% 
  filter(!is.na(POINT_X))

mydata.predobs3.i <- SpatialPointsDataFrame(coords = mydata.predobs3.i %>% 
                                              filter(!is.na(POINT_X))%>% 
                                              dplyr::select(x=POINT_X, y=POINT_Y) , 
                                            data = mydata.predobs3.i %>% 
                                              filter(!is.na(POINT_X)) %>% 
                                              dplyr::select(obs.rich,  isforest:pred.rich),
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
  st_as_sf() %>% 
  st_transform("+proj=eck4") 
mydata.predobs3.i <- mydata.predobs3.i %>% 
  st_coordinates()%>% 
  as.data.frame() %>% 
  bind_cols(mydata.predobs3.i %>% 
              st_drop_geometry()) %>% 
  ## calculate Obs-pred 
  mutate(resid=pred.rich - obs.rich)

mydata.predobs.i.toplot <- mydata.predobs3.i %>% 
  mutate(resid=ifelse(resid> 20,  20, resid)) %>% 
  mutate(resid=ifelse(resid< -20, -20, resid)) %>% 
  slice(sample(1:n())) ## random order

(ggresid.for <- w3a + 
    geom_point(data=mydata.predobs.i.toplot %>% 
                 filter(isforest=="Forest"), 
               aes(x=X, y=Y, col=resid), pch=16, cex=0.2, alpha=1, 
               position = position_jitter()) +
    ggtitle("Forest") + 
    theme(legend.position="none")
)
(ggresid.nonfor <- w3a + 
    geom_point(data=mydata.predobs.i.toplot %>% 
                 filter(isforest=="Non Forest"), 
               aes(x=X, y=Y, col=resid), pch=16, cex=0.2, alpha=1, 
               position = position_jitter()) +
    ggtitle("Non-Forest") 
)
ggresid <- ggresid.for / ggresid.nonfor

ggsave(filename=paste0(path.to.pics, "/FigS14_Residual_space_", which.pred, ".png"), 
       width=8, height=10, dpi=300, units="in", ggresid)
ggsave(filename=paste0(path.to.pics, "/FigS14_Residual_space_", which.pred, ".pdf"), 
       width=5, height=5, units="in", ggresid)


### Fig S15 - Check autocorrelation of residuals ####
library(gstat)
library(automap)
which.model <- "all"
size <- 10
## import BRT models
(listf <- list.files(path.to.BRTglobal,
                     pattern=paste0("^",which.model,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T))
order.i <- as.numeric(gsub(pattern=paste0("_",size,"m\\.RData"), replacement="", 
                           x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))
listf <- listf[order(order.i)]

myboundaries <- c(30, 50, 100, 250, 500, 750, 1500, 2000, 3000)
vg <- gstat::variogram(res~1, mydata.i.sp, boundaries=myboundaries)
ggplot(data=vg, aes(x=dist, y=gamma)) + 
  geom_point() + 
  ylim(0,3.5) + 
  theme_bw() + 
  ylab("Semivariance") + 
  xlab("Distance (km)")
ggsave(filename = file.path(path.to.pics, "FigS15_AutocorrelationResiduals_iter1.pdf"), 
       width = 6, height = 4, bg="white")
ggsave(filename = file.path(path.to.pics,"FigS15_AutocorrelationResiduals_iter1.png"), 
       width = 6, height = 4, bg="white", dpi=300, units="in")



### Fig S5 - Compare distributions of observed and predicted species richness ####
## only for iteration 1 ###

i <- mysample
load(listf[i])

(gbmdata <- reconstructGBMdata(modello) %>% 
    as_tibble() %>% 
    mutate(RELEVE_NR=mydata.i.sp$RELEVE_NR) %>% 
    mutate(pred = modello$fitted) %>% 
    #function(x){max((x-a)/b, 0)} #bias correction
    mutate(pred=(pred-mycoefs[[1]][i])/mycoefs[[2]][i]) %>% 
    mutate(pred=ifelse(pred<0,0,pred)) %>% 
    # end of bias correction
    dplyr::select(RELEVE_NR, 
                  Rel.area,
                  obs=y.data, 
                  pred, isforest, sBiomeName,
                  plants_recorded) %>% 
    mutate(plants_recorded=factor(plants_recorded, labels=levels(mydata.i.sp$plants_recorded)),
           sBiomeName=factor(sBiomeName, labels=levels(mydata.i.sp$sBiomeName)),
           isforest=factor(isforest, labels=levels(mydata.i.sp$isforest))) %>%
    mutate(isforest = fct_recode(isforest, For="for", `Nonfor`="nonfor")) %>% 
    mutate(Rel.area=cut(Rel.area, breaks=c(0,150, 600, 1200, Inf))) %>%
    mutate(Rel.area = fct_recode(Rel.area, `(600,1200]`="(600,1.2e+03]", `(1200-Inf]` = "(1.2e+03,Inf]")) %>% 
    mutate(plants_recorded = factor(plants_recorded,
                                    levels=c("complete", "woody_all", "woody_large"), 
                                    labels=c("compl", 
                                             "woody", 
                                             "trees"))) %>% 
    pivot_longer(cols = c(obs, pred), names_to = "value", values_to = "Richness") %>% 
    mutate(value = fct_recode(value, Observed="obs", Predicted="pred")) %>% 
    left_join(biome.labs, by=c("sBiomeName"="name")) %>% 
    dplyr::rename("Biome"=labels)
)

xlims <- quantile(gbmdata$Richness, c(0.01, 0.99))
ggtemplate <- ggplot(data=gbmdata, aes(x=Richness, group=value, fill=value)) + 
  geom_density( alpha=0.55, col=NA)+
  scale_fill_brewer(palette = "Set1", name=NULL)+
  xlim(xlims) + 
  theme_classic() + 
  theme(legend.position = "bottom")
#by for nonfor
gg1 <- ggtemplate + 
  facet_grid(isforest ~ Rel.area) + 
  theme(legend.position = "none", 
        axis.title.x=element_blank())
#by biome
gg2 <- ggtemplate + 
  facet_grid(Biome ~ Rel.area) + 
  theme(legend.position = "none", 
        axis.title.x=element_blank())
#by completeness
gg3 <- ggtemplate + 
  facet_grid(plants_recorded ~ Rel.area) + 
  ylim(c(0, 0.2))

ggtot <- gg1/gg2/gg3 + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = unit(c(2, 10, 3), c('null')))

ggsave(filename=file.path(path.to.pics, "FigS5_ComparePredObsDistr.png"), 
       width=8, height=10, dpi=300, plot = ggtot)
ggsave(filename=file.path(path.to.pics, "FigS5_ComparePredObsDistr.pdf"), 
       width=8, height=10, plot = ggtot)



### Table S5 - Create Summary showing observer Species richness at different grains #### 
#### Calculate summaries
## ONLY for complete plots effectively used in the analysis

mydata.summ <- mydata %>%
  filter(RELEVE_NR %in% all.used.plots) %>% 
  filter(plants_recorded=="complete") %>% 
  dplyr::rename(rich=rich.complete) %>% 
  mutate(Rel.area.cut=cut(Rel.area, breaks=c(0,20, 150, 600, 1200, Inf)))  %>% 
  filter( (isforest=="nonfor" & Rel.area.cut %in% c("(0,20]", "(20,150]", "(600,1.2e+03]")) | 
            (isforest=="for" & Rel.area.cut %in%    c("(150,600]", "(600,1.2e+03]", "(1.2e+03,Inf]"))) %>% 
  unite(isforest, Rel.area.cut, col="Rel.area.class", remove=F) %>% 
  mutate(Rel.area.class=factor(Rel.area.class)) %>% 
  mutate(Rel.area.class = fct_collapse(Rel.area.class, 
                                       small = c( "for_(150,600]", "nonfor_(0,20]"),
                                       medium = c("for_(600,1.2e+03]", "nonfor_(20,150]"), 
                                       large = c("for_(1.2e+03,Inf]", "nonfor_(600,1.2e+03]"))) %>% 
  mutate(Rel.area.class=factor(Rel.area.class, levels=c("small","medium", "large")))

summary.global <- mydata.summ %>%
  group_by(isforest, Rel.area.class) %>% 
  dplyr::summarize_at(.vars=vars(rich), 
                      .funs=list(n = ~n(), 
                                 min=~min(., na.rm=T),
                                 med=~median(., na.rm=T),
                                 max=~max(., na.rm=T), 
                                 IQR=~quantile(., 0.75, na.rm=T)-quantile(., 0.25, na.rm=T))) %>% 
  pivot_longer(cols = n:IQR, names_to="metrics", values_to = "Richness") %>% 
  unite(metrics, Rel.area.class, sep = "_", col = "metrics", remove=T) %>% 
  mutate(sBiomeName="zGlobal") %>% 
  pivot_wider(id_cols = c(sBiomeName, isforest), names_from = "metrics", values_from = "Richness")

summary.biome <- mydata.summ %>%
  group_by(isforest, sBiomeName, Rel.area.class) %>% 
  dplyr::summarize_at(.vars=vars(rich), 
                      .funs=list(n = ~n(), 
                                 min=~min(., na.rm=T),
                                 med=~median(., na.rm=T),
                                 max=~max(., na.rm=T), 
                                 IQR=~quantile(., 0.75, na.rm=T)-quantile(., 0.25, na.rm=T))) %>% 
  pivot_longer(cols = n:IQR, names_to="metrics", values_to = "Richness") %>% 
  unite(metrics, Rel.area.class, sep = "_", col = "metrics", remove=T) %>% 
  pivot_wider(id_cols = c("sBiomeName", isforest), names_from = "metrics", values_from = "Richness")

## format as output table
(summary.out <- summary.biome %>% 
    bind_rows(summary.global) %>% 
    arrange(isforest, sBiomeName))
write_csv(x = summary.out, file = file.path(path.to.output,"TableS3_SummaryObserved.csv"))




### SpatialCV graphs ####
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



index <- 1
metric <- "sr10" # index.list$Var1[index]
size <- 10
which.model <- index.list$Var2[index]
listf <- list.files(path.to.BRTglobal, 
                    pattern=paste0("^",which.model,"BRTs_direct99-[0-9]*-[0-9]*_", size, "m\\.RData$"), full.names =T)
order.i <- as.numeric(gsub(pattern=paste0("_",size,"m\\.RData"), replacement="", 
                           x=str_extract(listf, pattern="[0-9]*_[0-9]*m\\.RData$")))

iteration <- 1
#### Calculate cv correlations across models, and biomes
cor.out <- NULL
i <- order.i[iteration]
print(paste(which.model, metric, i))
load(listf[iteration])
## Spatial CrossValidation using BlockCV package
# Re-create dataset i and create spatial sf object
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
##rename world.raster first
renaming.string  <- data.frame(selected.predictors, var.labs0) %>% 
  filter(selected.predictors %in% selected.variables.quantitative) %>% 
  #mutate(selected.predictors=paste0("'", selected.predictors, "'")) %>% 
  #mutate(var.labs0=paste0("'", var.labs0, "'")) %>% 
  #unite(var.labs0, selected.predictors, col = "string", sep = " = ") %>% 
  pull(var.labs0) %>% 
  #renaming.string <- renaming.string %>% 
  str_replace(pattern=" - ", "_") %>% 
  str_replace(pattern="-", "_") %>% 
  str_replace(pattern=" \\(", "_") %>%
  str_replace(pattern="\\)", "") %>%
  str_replace(pattern="\\/", "_") %>% 
  str_replace(pattern=" \\.|\\. ", ".") %>% 
  str_replace_all(pattern=" ", "_") %>% 
  str_squish()

names(world.raster) <- renaming.string

### Fig S12 - Spatial autocorrelation of predictors ####
sac <- spatialAutoRange(rasterLayer = world.raster,
                        sampleNumber = 5000,
                        #speciesData=pa_data %>% 
                        #  filter(complete.cases(mydata.world)),
                        doParallel = F,
                        showPlots = TRUE)
autocorrelation.i <- sac$range
print(paste("Autocorrelation in run", i, "is", round(autocorrelation.i/1000), "km"))

#plotting
pdf(file = file.path(path.to.pics, "FigS12_autocorrelation.pdf"), width = 8, height = 6, bg = "white")
plot(sac)
dev.off()

png(file = file.path(path.to.pics, "FigS12_autocorrelation.png"), width = 8, height = 6, bg = "white", res =300, units = "in")
plot(sac)
dev.off()


### Fig S13 - SpatialCV blocks ####
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
sb2$plots + 
  geom_sf(data = pa_data, alpha = 0.5, cex=0.6)
ggsave(file.path(path.to.pics, "FigS13_SpatialBlocks_iteration4.png"), width = 8, height=5,dpi = 300, bg = "white")
ggsave(file.path(path.to.pics, "FigS13_SpatialBlocks_iteration4.pdf"), width = 8, height=5,bg = "white")




