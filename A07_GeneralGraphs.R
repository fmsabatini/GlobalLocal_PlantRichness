## Additional Figures and Tables 
library(tidyverse)
library(sp)
library(rgdal)
library(viridis)
library(sf)
library(rnaturalearth)
library(dggridR)
library(cowplot)

source("A21_w3a_TemplateEckert.R")

### Set paths
path.to.pics <- "../sPlot/_versions/NatCommR1/_pics"
path.to.output <- "../sPlot/_versions/NatCommR1"
path.to.input <- "../sPlot/_input"



# load mydata
load(file.path(path.to.input, "Mydata_global_NatCommR1.RData"))
mydata <- mydata %>%
  #mutate(isforest=as.numeric(isforest)) %>% ##classification based on land cover
  mutate(isforest=ifelse(isforest, "for", "nonfor")) %>%
  mutate(isforest=factor(isforest)) %>% 
  mutate(CONTINENT=factor(CONTINENT))

# 417180 pre-seltected plots 

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

length(unique(unlist(rel.list))) #
'[1] 170700' # total number of plots effectively used

## calculate leverage of plots
res <- data.frame(RELEVE_NR=unlist(rel.list)) %>% 
  group_by(RELEVE_NR) %>% 
  summarize(n=n())


### Fig S6 - No. of plots per grain across biomes ######
biome.order <- c('Polar and subpolar zone' ,'Alpine' ,'Boreal zone' ,'Temperate midlatitudes' ,
                 'Dry midlatitudes' ,'Dry tropics and subtropics' ,'Subtrop. with year-round rain' ,
                 'Subtropics with winter rain' ,'Tropics with summer rain' ,'Tropics with year-round rain')
biome.labs <- c('Polar &\n subpolar' ,'Alpine' ,'Boreal zone' ,'Temperate\n midlatitudes' ,
                'Dry midlatitudes' ,'Dry tropics\n & subtropics' ,'Subtropics -\n year-round rain' ,
                'Subtropics -\n winter rain' ,'Tropics -\n summer rain' ,'Tropics -\n year-round rain')

plotdistr <- mydata %>% 
  mutate(sBiomeName=factor(sBiomeName, levels=biome.order, labels=biome.labs)) %>% 
  mutate(plot_size=cut(Rel.area, c(0,150, 600, 1200, Inf), labels=c("(0, 150]", "(150, 600]", "(600, 1200]", ">1200"))) %>% 
  group_by(isforest, sBiomeName, plot_size) %>% 
  mutate(isforest=factor(isforest, levels=c("for", "nonfor"), labels=c("Forest", "Non forest"))) %>% 
  summarize(n=n()) %>% 
  complete(isforest, sBiomeName, plot_size, fill=list(n=NA))

ggplotsize <- ggplot(data=plotdistr) + 
  geom_bar(aes(x=sBiomeName, y=n, group=plot_size, fill=plot_size), stat="identity", position ="dodge") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1), 
        legend.position = "bottom") +  
  facet_grid(isforest~.) +
  scale_y_log10(labels = function(x) format(x, scientific = F)) +
  xlab("Biome") + 
  ylab("Number of plots") + 
  scale_fill_discrete(name= bquote(Plot~size~(m^{2})))

ggsave(filename=file.path(path.to.pics, "FigS6_plotsize.png"), device = "png", dpi=400, 
       width = 7, height=5, ggplotsize)

### Fig S7 - No. of plots per completeness across biomes ######
completeness_distr <- mydata %>% 
  mutate(plants_recorded=factor(plants_recorded, 
                                levels=c("woody_large", "woody_all", "complete"),
                                labels=c("only trees", "trees & shrubs", "complete vegetation"))) %>% 
  mutate(sBiomeName=factor(sBiomeName, levels=biome.order, labels=biome.labs)) %>% 
  group_by(sBiomeName, plants_recorded, isforest) %>% 
  mutate(isforest=factor(isforest, levels=c("for", "nonfor"), labels=c("Forest", "Non forest"))) %>% 
  summarize(n=n()) %>% 
  complete(isforest, sBiomeName, plants_recorded, fill=list(n=NA))

ggcompleteness <- ggplot(data=completeness_distr) + 
  geom_bar(aes(x=sBiomeName, y=n, group=plants_recorded, fill=plants_recorded), stat="identity", position ="dodge") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1), 
        legend.position = "bottom") + 
  facet_grid(isforest~.) +
  scale_y_log10(labels = function(x) format(x, scientific = F)) +
  xlab("Biome") + 
  ylab("Number of plots") + 
  scale_fill_discrete(name= "Plants recorded")

ggsave(filename=file.path(path.to.pics, "FigS7_completeness.png"), device = "png", dpi=400, 
       width = 7, height=5, ggcompleteness)




### Table S1 - Plot no. by database ##### 
### create table of databases with plot number and biblio
#Import databases and create reference tags
library(bib2df)
bib.db <- bib2df("/data/sPlot/users/Francesco/_sPlot_Management/Consortium/sPlot_References.bib")
#Import database-level information
databases <- read_csv("/data/sPlot/users/Francesco/_sPlot_Management/Consortium/Databases.out.csv")

# create citation tags that can be picked up by Manubot
databases <- databases %>% 
  left_join(bib.db %>% 
              dplyr::select(BIBTEXKEY, DOI, URL), 
            by="BIBTEXKEY") %>%
  dplyr::select(-DOI, -URL)


# Create Table S1
databases21 <- read_csv("/data/sPlot/users/Francesco/_sPlot_Management/Consortium/Databases.out2.1.csv")
load("/data/sPlot/releases/sPlot2.1/sPlot_header_20161124.RData")
source("A96_fixheaderPaper11.R")
header.fix <- fix.header11(header)
header.fix <- header.fix %>% 
  dplyr::select(PlotObservationID, Dataset) %>% 
  mutate(Dataset=replace(Dataset, list=Dataset %in% c("Argentina_Chaco_Espinal","Argentina_Cordoba"), 
                         "Argentina_Central")) %>%
  mutate(Dataset=replace(Dataset, list=Dataset=="SIVIM wetlands", 
                         "Spain_sivim_wetlands")) %>%
  left_join(databases21 %>% 
              dplyr::select(`GIVD ID`, label), 
            by=c("Dataset"="label"))


table1 <- mydata %>%
  left_join(header.fix %>% 
              dplyr::select(RELEVE_NR=PlotObservationID, `GIVD ID`)) %>% 
  mutate(GIVD_ID=ifelse( (Country=="Peru" & plants_recorded=="complete" & RELEVE_NR>9000000), 
                         "00-00-001", 
                         GIVD_ID)) %>% 
  mutate(GIVD_ID=coalesce(GIVD_ID, `GIVD ID`)) %>% 
  dplyr::select(-`GIVD ID`) %>% 
  left_join(databases, by=c("GIVD_ID"="GIVD ID")) %>% 
  group_by(GIVD_ID) %>% 
  summarize(contributed_plots=n(), .groups = 'drop') %>% 
  left_join(databases %>% 
              filter(`Still in sPlot`==T, 
                     Via!="Aggregator") %>% 
              dplyr::select(-Via, -`Still in sPlot`) %>% 
              distinct(),
            by=c("GIVD_ID"="GIVD ID")) %>% 
  filter(!is.na(contributed_plots)) %>% 
  dplyr::select(GIVD_ID, 
                `Dataset name`=`DB_name GIVD`,  
                `Nr. of unique plots used` = contributed_plots, Ref=BIBTEXKEY) %>% 
  distinct() %>% 
  arrange(GIVD_ID) %>% 
  replace_na(list(Ref=""))


write_csv(table1, file.path(path.to.output, "TableS1_Databases.csv"))




### Fig S1 - PLOT distribution of resampled plots ###############
#### Transform to spatial data.frame
mycrs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
header.filtered.sf <- SpatialPointsDataFrame(coords=mydata %>% 
                                               dplyr::select(POINT_X, POINT_Y),
                                             proj4string = mycrs, 
                                             data=mydata %>% 
                                               dplyr::select(-POINT_X, -POINT_Y)) %>% 
  st_as_sf() %>% 
  st_transform("+proj=eck4")

header.resampled.sf <- header.filtered.sf[which(header.filtered.sf$RELEVE_NR %in% res$RELEVE_NR),]

## resampled plots
header.resampled.sf <- header.resampled.sf %>% 
  mutate(isforest=fct_recode(isforest, Forest="for", `Non Forest`="nonfor")) 
ggRes.for <- w3a + 
  geom_sf(data=header.resampled.sf %>% 
            filter(isforest=="Forest"), 
          aes(col=isforest), pch="+", cex=0.9, alpha=1/6) + 
  scale_color_brewer(palette="Set1", direction=-1, name = "Vegetation type") +
  #guides(color = guide_legend(override.aes= list(alpha = 1, cex=3))) + 
  theme(legend.position="none", 
        legend.background = element_rect(color = NA), 
        plot.title = element_text(hjust=0.5)) + 
  ggtitle(paste("Forest (n = ", header.resampled.sf %>% 
                  filter(isforest=="Forest") %>% nrow(), ")", sep=""))

ggRes.nonfor <- w3a + 
  geom_sf(data=header.resampled.sf %>% 
            filter(isforest=="Non Forest"), 
          aes(col=isforest), pch="+", cex=0.9, alpha=1/6) + 
  scale_color_brewer(palette="Set1", direction=-1, name = "Vegetation type") +
  #guides(color = guide_legend(override.aes= list(alpha = 1, cex=3))) + 
  theme(legend.position="none", 
        legend.background = element_rect(color = NA), 
        plot.title = element_text(hjust=0.5))+ 
  ggtitle(paste("Non-Forest (n = ", header.resampled.sf %>% 
                  filter(isforest=="Non Forest") %>% nrow(), ")", sep=""))


header.resampled.sf.res <- header.resampled.sf %>% 
  left_join(res)
ggRes.leverage <- w3a + 
  geom_sf(data=header.resampled.sf.res, aes(col=n), pch="+", cex=0.9) + 
  scale_color_viridis(name="Plot leverage (%)", limits=c(0,100), breaks=scales::pretty_breaks(n=5)) +
  theme(legend.position="bottom", 
        legend.background = element_rect(color = NA), 
        plot.title = element_text(hjust=0.5))+ 
  ggtitle("Plot Leverage")

FigS1 <- plot_grid(ggRes.for, ggRes.nonfor, ggRes.leverage, 
                   nrow=3, labels = c("A", "B", "C"), rel_heights = c(1,1,1.4))
ggsave(filename=file.path(path.to.pics, "FigS1_AllResampledPlots_and_leverage.png"), device = "png", dpi=600, 
       width = 5, height=9, FigS1, bg="white")



### Additional plot, not in paper ####
## PLOT density of all used plots (for vs nonfor) in hexagons 
# #### build hexagon grid
# dggs          <- dgconstruct(spacing=200, metric=T, resround='down')
# #Get the corresponding grid cells for each earthquake epicenter (lat-long pair)
# mydata.resampled <- mydata %>% 
#   filter(RELEVE_NR %in% res$RELEVE_NR)
# mydata.resampled$cell <- dgGEO_to_SEQNUM(dggs, mydata.resampled$POINT_X, mydata.resampled$POINT_Y)$seqnum
# 
# #Calculate mean metric for each cell
# mydata.resampled.out   <- mydata.resampled %>% 
#   dplyr::select(cell, isforest) %>% 
#   group_by(cell) %>% 
#   summarise(count=log(n(),10),
#             count.for=log(sum(isforest=="for"),10),
#             count.nonfor=log(sum(isforest=="nonfor"),10))
# #Get the grid cell boundaries for cells 
# grid   <- dgcellstogrid(dggs, mydata.resampled.out$cell, frame=F) %>%
#   st_as_sf() %>% 
#   mutate(cell = mydata.resampled.out$cell) %>% 
#   mutate(value.out=mydata.resampled.out$count) %>% 
#   st_transform("+proj=eck4") %>% 
#   st_wrap_dateline(options = c("WRAPDATELINE=YES"))
# 
# ## plotting
# ggCount <- w3a + 
#   geom_sf(data=grid, aes(fill=value.out),lwd=0, alpha=0.9)    +
#   geom_sf(data = countries, col = "grey10", fill=NA, lwd = 0.3) + 
#   labs(fill = "# of Plots\n(Log 10)")
# ggsave(filename=file.path(path.to.pics, "FigSXXV_AllResampledPlots.png"), 
#        device = "png", dpi=400, width = 6, height=3.5, ggCount, bg="white")



### Calculate share of complete\incomplete plots after resampling ####
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


get.incomplete <- function(x){mydata %>%
  filter(RELEVE_NR %in% rel.list[[x]]) %>% 
  mutate(plants_recorded=fct_collapse(plants_recorded, 
                                      woody=c("woody_all", "woody_large"))) %>% 
  dplyr::count(sBiomeName, plants_recorded) %>% 
  mutate(well.sampled=ifelse(sBiomeName %in% sel.biomes, T, F)) %>% 
  pivot_wider(values_from=n, names_from=plants_recorded) %>% 
  group_by(well.sampled) %>% 
  summarize(n=mean(woody))
}
prop.incomplete <- lapply(1:99, get.incomplete) %>% 
  bind_rows() %>% 
  group_by(well.sampled) %>% 
  summarize(mean(n))

