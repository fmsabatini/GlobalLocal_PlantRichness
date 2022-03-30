library(tidyverse)

load("../sPlot/_input/Mydata_global_NatCommR1.RData")

load("/data/sPlot/releases/sPlot2.1/DT2_20161025.RData")
DT2 -> DT2.1
load("/data/sPlot/releases/sPlot3.0/DT_sPlot3.0.RData")
DT2 -> DT3
rm(DT2)


### Check VegBank Data
#between 40-50° N and 103-125 ° W
vb <- read_csv("../sPlot/_input/USA_VegBank_check/VegBank_to_check.csv")
inw <- vb %>% 
  filter(str_detect(authorplotcode_vb, "^INW")) %>% 
  pull(RELEVE_NR)

inw3 <- mydata %>% 
  filter(POINT_Y>40 & POINT_Y<50 & POINT_X< -103 & POINT_X> -125) %>% 
  filter(RELEVE_NR %in% inw) %>% 
  filter(rich.complete<5) %>% 
  pull(RELEVE_NR)

DT2.1 %>% filter(PlotObservationID == sample(inw3, 1))

rejects <- openxlsx::read.xlsx("../sPlot/_input/USA_VegBank_check/Rejects2F.xlsx", sheet = 1) %>% 
  filter(!is.na(authorplotcode_vb))

rejects.vegbank <- vb %>% 
  filter(authorplotcode_vb %in% rejects$authorplotcode_vb) %>% 
  pull(RELEVE_NR)



## Check species poor plots
aa <- mydata %>% 
  #filter(RELEVE_NR %in% res$RELEVE_NR) %>% 
  filter(rich.complete<3) %>% 
  filter(plants_recorded=="complete") %>% 
  pull(RELEVE_NR)

# check species in species poor plots
DT2.1 %>% 
  as_tibble() %>% 
  dplyr::filter(PlotObservationID %in% aa) %>% 
  count(species) %>% 
  arrange(desc(n))

### Exclude wetland plots, based on hydrophytes
wetland.species <- c('Lemna minor', 
                     'Lemna trisulca', 
                     'Spirodela polyrrhiza', 
                     'Glyceria fluitans',
                     "Glyceria notata",
                     "Myriophyllum implicatum",
                     "Myriophyllum verticillatum",
                     "Myriophyllum sibiricum", 
                     "Myriophyllum spicatum", 
                     "Myriophyllum alterniflorum", 
                     "Stuckenia pectinata", 
                     "Lemna", 
                     "Lemna gibba", 
                     "Alisma plantago-aquatica", 
                     "Callitriche hermaphroditica", 
                     "Ceratophyllum demersum", 
                     "Potamogeton crispus",
                     "Butomus umbellatus",
                     "Lemna turionifera",
                     "Hottonia palustris", 
                     "Hippuris vulgaris",
                     "Callitriche cophocarpa",
                     "Nymphaea candida",
                     "Pistia stratiotes",
                     "Potamogeton pusillus", 
                     "Azolla filiculoides", 
                     "Sagittaria latifolia", 
                     "Sagittaria sagittifolia", 
                     "Nymphoides peltata", 
                     "Potamogeton berchtoldii", 
                     "Potamogeton natans", 
                     "Utricularia vulgaris", 
                     "Utricularia intermedia", 
                     "Callitriche obtusangula",
                     "Elodea canadensis", 
                     "Zannichellia palustris", 
                     "Wolffia arrhiza", 
                     "Eleocharis acicularis", 
                     "Menyanthes trifoliata", 
                     "Potamogeton gramineus")
      
bb <- DT2.1 %>% 
  as_tibble() %>% 
  filter(PlotObservationID %in% mydata$RELEVE_NR) %>% 
  mutate(wetsp=species %in% wetland.species) %>% 
  group_by(PlotObservationID) %>% 
  summarize(wet.cover=sum(Relative.cover*wetsp), rich=n()) %>% 
  filter(wet.cover>=.50) %>% 
  pull(PlotObservationID)

mydata.R2 <- mydata %>% 
  filter(!RELEVE_NR %in% bb) %>% 
  filter(!RELEVE_NR %in% rejects.vegbank)

mydata <- mydata.R2
save(mydata, file = "../sPlot/_input/Mydata_global_NatCommR2.RData")


# check species in species poor plots
aa2 <- mydata.R2 %>% 
  filter(rich.complete<3) %>% 
  filter(plants_recorded=="complete") %>% 
  pull(RELEVE_NR)

DT2.1 %>% 
  as_tibble() %>% 
  dplyr::filter(PlotObservationID %in% aa2) %>% 
  count(species) %>% 
  arrange(desc(n))


