library(tidyverse)

load("/data/sPlot/releases/sPlot2.1/DT2_20161025.RData")
DT2.0 <- DT2
##### Filter out all species entry of non-vascular plants #####
#use taxon groups from sPlot 3.0
load("/data/sPlot/releases/sPlot3.0/DT_sPlot3.0.RData")
DT3.0 <- DT2
rm(DT2)
dt3.species <- DT3.0 %>% 
  dplyr::rename(Taxon.group=taxon_group) %>% 
  distinct(species, species_original, Taxon.group) %>% 
  filter(Taxon.group != "Unknown") %>% 
  separate(species, sep=" ", into = c("genus", "species")) %>% 
  dplyr::select(genus, Taxon.group) %>% 
  mutate(Taxon.group=replace(Taxon.group, 
                             list=genus=="Friesodielsia", 
                             values="Vascular plant")) %>% 
  mutate(Taxon.group=replace(Taxon.group, 
                             list=genus=="Lamprothamnus", 
                             values="Alga_Stonewort")) %>% 
  filter(!is.na(genus)) %>% 
  filter(!genus %in% c("Hepatica", "×")) %>% 
  distinct()
  

dt2.species <- DT2.0 %>% 
  distinct(species) %>% 
  mutate(species0 = species) %>% 
#  mutate(Taxon.group=replace(Taxon.group, 
#                             list=Taxon.group=="Unknown", 
#                             values=NA)) %>% 
#  mutate(Taxon.group=replace(Taxon.group, 
#                             list=Taxon.group %in% c("Alga", "Stonewort"), 
#                             values="Alga_Stonewort")) %>% 
  separate(species, sep=" ", into = c("genus", "species")) %>% 
  left_join(dt3.species, by="genus") %>% 
#  mutate(Taxon.group=coalesce(Taxon.group.x, as.character(Taxon.group.y))) %>% 
  dplyr::select(species0, Taxon.group) %>% 
  rename(species=species0)

DT2 <- DT2.0 %>% 
  mutate(Taxon.group=replace(Taxon.group, 
                               list=Taxon.group=="Unknown", 
                             values=NA)) %>% 
  mutate(Taxon.group=replace(Taxon.group, 
                               list=Taxon.group %in% c("Alga", "Stonewort"), 
                               values="Alga_Stonewort")) %>% 
  left_join(dt2.species, by="species") %>% 
  #coalesce prioritizing sPlot 3.0
  mutate(Taxon.group=coalesce(as.character(Taxon.group.y), as.character(Taxon.group.x))) %>% 
  dplyr::select(-Taxon.group.x, -Taxon.group.y) %>% 
  filter(!Taxon.group %in% c("Alga_Stonewort", "Lichen", "Moss"))


### exclude all taxa matched at family or higher level #ending "...aceae or ...psida"
DT2.out <- DT2 %>% 
  filter(!grepl(pattern="[A-Za-z]*aceae$|[A-Za-z]*opsida$|^[A-Za-z]*ales$|^[A-Za-z]*phyta$", species))
  

DT2 <- DT2.out
save(DT2, file = "../sPlot/_derived_data/DT2_20161025_filtered.RData")





#### Add column tree/shrub  ####
load("../sPlot/_derived_data/DT2_20161025_filtered.RData")
load("/data/sPlot2.0/TRY.all.mean.sd.3.by.genus.species.tree.Rdata")

gf <- TRY.all.mean.sd.3.by.genus.species.tree %>% 
  dplyr::select(species=StandSpeciesName,is.tree.or.tall, is.shrub) %>% 
  filter(complete.cases(.))
## only cover ~half of the species

load("/data/sPlot/releases/sPlot3.0/Traits_CWMs_sPlot3.RData")
gf2 <- sPlot.traits %>% 
  mutate(is.shrub=NA) %>% 
  mutate(is.shrub=replace(is.shrub, 
                          list=str_detect(GrowthForm, "shrub"), 
                          values=T)) %>% 
  dplyr::select(species, is.shrub, is.tree.or.tall.shrub) %>% 
  filter(complete.cases(.))
  
  
DT2 <- DT2 %>% 
  left_join(gf, by="species") %>%
  left_join(gf2, by="species") %>% 
  mutate(is.tree.or.tall=coalesce(is.tree.or.tall, is.tree.or.tall.shrub)) %>% 
  mutate(is.shrub=coalesce(is.shrub.x, is.shrub.y)) %>% 
  dplyr::select(-is.tree.or.tall.shrub, -is.shrub.x, -is.shrub.y)

tmp <- DT2 %>% distinct(species, is.tree.or.tall, is.shrub)
table(tmp$is.tree.or.tall, tmp$is.shrub, exclude=NULL)
"        FALSE  TRUE  <NA>
  FALSE 15222  3535     0
  TRUE   6571  1498     0
  <NA>      0     0 31062"  
  
save(DT2, file = "../sPlot/_derived_data/DT2_20161025_filtered.RData")

### use info from header. Assign to is.tree.or.tall all entries from plots
### where only woody species were sampled
load("/data/sPlot/releases/sPlot2.1/sPlot_header_20161124.RData")
source("A96_fixheaderPaper11.R")
header <- fix.header11(header)
plants.recorded <- header %>% 
  dplyr::select(PlotObservationID, plants_recorded2) %>% 
  filter(plants_recorded2 != "complete")

gf.complement <- DT2 %>% 
  filter(PlotObservationID %in% (plants.recorded  %>% 
                                   pull(PlotObservationID))) %>% 
  left_join(plants.recorded, by="PlotObservationID") %>% 
  dplyr::select(species, plants_recorded2) %>% 
  distinct(species, .keep_all = T) %>% 
  mutate(is.tree.or.tall2 = ifelse(plants_recorded2=="woody_large", T, NA)) %>% 
  mutate(is.shrub2 = ifelse(plants_recorded2=="woody_large", F, NA)) %>% 
  mutate(is.shrub2 = ifelse(plants_recorded2=="woody_all", T, is.shrub2))
  
  
DT2 <- DT2 %>% 
  left_join(gf.complement, by="species") %>% 
  mutate(is.tree.or.tall=coalesce(is.tree.or.tall, is.tree.or.tall2)) %>% 
  mutate(is.shrub=coalesce(is.shrub, is.shrub2)) %>% 
  dplyr::select(-is.tree.or.tall2, -is.shrub2)

tmp <- DT2 %>% distinct(species, is.tree.or.tall, is.shrub)
table(tmp$is.tree.or.tall, tmp$is.shrub, exclude=NULL)

"        FALSE  TRUE  <NA>
  FALSE 15222  3535     0
  TRUE   9684  1498     0
  <NA>      0  3852 24097"

### check the most species rich genera 
tmp2 <- tmp %>% 
  filter(is.na(is.tree.or.tall) & is.na(is.shrub)) %>% 
  separate(species, into = c("genus", "species")) %>% 
  group_by(genus) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  slice(1:200)
sum(tmp2$n) ## 10528 species in the 200 most speciose genera

genera.gf <- c('Carex' = 'h','Astragalus' = 'h','Acacia' = 't','Euphorbia' = NA,'Silene' = 'h',
               'Eucalyptus' = 't','Ranunculus' = 'h','Rubus' = 's','Senecio' = 'h','Viola' = 'h',
               'Galium' = 'h','Erica' = 's','Festuca' = 'h','Ficus' = 't','Hypericum' = 'h',
               'Hieracium' = 'h','Cyperus' = 'h','Centaurea' = 'h','Miconia' = 't','Salix' = 't',
               'Solanum' = 'h','Allium' = 'h','Quercus' = 't','Potentilla' = 'h','Juncus' = 'h',
               'Dianthus' = 'h','Veronica' = 'h','Indigofera' = 's','Syzygium' = 't','Artemisia' = 'h',
               'Campanula' = 'h','Inga' = 't','Saxifraga' = 'h','Poa' = 'h','Taraxacum' = 'h',
               'Eugenia' = 't','Trifolium' = 'h','Cirsium' = 'h','Eragrostis' = 'h','Diospyros' = 't',
               'Alyssum' = 'h','Alchemilla' = 'h','Stipa' = 'h','Geranium' = 'h','Polygala' = 's',
               'Pedicularis' = 'h','Psychotria' = 't','Ocotea' = 't','Piper' = 's','Panicum' = 'h',
               'Aconitum' = 'h','Draba' = 'h','Ilex' = 't','Oxytropis' = 'h','Aristida' = 'h',
               'Polygonum' = 'h','Cerastium' = 'h','Minuartia' = 'h','Erigeron' = 'h',
               'Acantholimon' = 'h','Prunus' = 't','Asplenium' = 'h','Croton' = NA,'Limonium' = 'h',
               'Plantago' = 'h','Thymus' = 'h','Rosa' = 's','Rumex' = 'h','Verbascum' = 'h',
               'Astracantha' = 'h','Atriplex' = NA,'Restio' = NA,'Salvia' = 'h','Aspalathus' = 's',
               'Crotalaria' = NA,'Helichrysum' = 'h','Pouteria' = 't','Sedum' = 'h','Crassula' = 'h',
               'Vicia' = 'h','Rhynchospora' = 'h','Gentiana' = 'h','Elymus' = 'h','Oxalis' = 'h',
               'Linum' = 'h','Erysimum' = 'h','Dryopteris' = 'h','Eleocharis' = 'h',
               'Delphinium' = 'h','Ipomoea' = 'h','Melaleuca' = 't','Clematis' = 's',
               'Epilobium' = 'h','Saussurea' = 'h','Crepis' = 'h','Lathyrus' = 'h','Ribes' = 's',
               'Agrostis' = 'h','Pelargonium' = 'h','Tephrosia' = NA)

genera.gf2 <- tmp2$genus
genera.gf2 <- genera.gf2[which(!genera.gf2 %in% names(genera.gf))]
#paste(paste("'", genera.gf2, "' = ''", sep=""), collapse=",")
genera.gf2 <- c('Gypsophila' = 'h','Asperula' = 'h','Paronychia' = 'h','Arenaria' = 'h',
                'Cliffortia' = 's','Thesium' = 's','Achillea' = 'h','Scutellaria' = 'h',
                'Calamagrostis' = 'h','Elaphoglossum' = 'h','Stachys' = 'h','Huperzia' = 'h',
                'Linaria' = 'h','Gentianella' = 'h','Valeriana' = 'h','Phylica' = 's',
                'Stellaria' = 'h','Cardamine' = 'h','Ficinia' = 'h','Fimbristylis' = 'h',
                'Lobelia' = 's','Ornithogalum' = 'h','Heliotropium' = 'h','Iris' = 'h',
                'Diplostephium' = 't','Scorzonera' = 'h','Baccharis' = 's','Asparagus' = 's',
                'Cheilanthes' = 'h','Solidago' = 'h','Thlaspi' = 'h','Lupinus' = 'h',
                'Anemone' = 'h','Aster' = 'h','Genista' = 's','Peucedanum' = 'h',
                'Selaginella' = 'h','Isatis' = 'h','Aethionema' = 'h','Agathosma' = 's',
                'Erodium' = 'h','Parrya' = 'h','Polystichum' = 'h','Tragopogon' = 'h',
                'Wahlenbergia' = 'h','Acanthophyllum' = 's','Anthemis' = 'h','Athyrium' = 'h',
                'Corydalis' = 'h','Eryngium' = 'h','Hermannia' = 's','Hesperis' = 's',
                'Justicia' = 's','Pentacalia' = NA,'Scleria' = 'h','Scrophularia' = 'h',
                'Berberis' = 's','Bupleurum' = NA,'Castilleja' = 'h','Eremogone' = NA,
                'Gnidia' = NA,'Goodenia' = 'h','Myosotis' = 'h','Papaver' = 'h',
                'Pentaschistis' = 'h','Salsola' = 's','Tetraria' = 's','Teucrium' = 's',
                'Crataegus' = 't','Dioscorea' = 'h','Lachemilla' = NA,'Othonna' = NA,
                'Zygophyllum' = NA,'Carduus' = 'h','Colchicum' = 'h','Convolvulus' = 'h',
                'Desmodium' = 'h','Moraea' = 'h','Muraltia' = 's','Nepeta' = 'h',
                'Onosma' = 'h','Xyris' = 'h','Clinopodium' = 'h','Hibiscus' = NA,
                'Hydrocotyle' = 'h','Knautia' = 'h','Leandra' = NA,'Luzula' = 'h',
                'Packera' = 'h','Pteris' = 'h','Ruschia' = 'h','Androsace' = 'h',
                'Centella' = 'h','Consolida' = 'h','Grevillea' = 't','Lampranthus' = 'h',
                'Lilium' = 'h','Phyllanthus' = 'h','Puccinellia' = 'h','Rhamnus' = 't',
                'Sporobolus' = 'h','Thalictrum' = 'h','Angelica' = 'h','Armeria' = 'h',
                'Calceolaria' = 'h','Lactuca' = 'h','Lotus' = 'h','Nassella' = 'h',
                'Platanthera' = 'h','Scabiosa' = 'h','Symphyotrichum' = 'h','Arabis' = 'h',
                'Asclepias' = 'h','Barleria' = 'h','Begonia' = 'h','Crocus' = 'h',
                'Eriocaulon' = 'h','Eupatorium' = 'h','Gynoxys' = NA,'Jurinea' = 'h',
                'Lepidium' = 'h','Orobanche' = 'h','Paspalum' = 'h','Potamogeton' = 'h',
                'Pterostylis' = 'h','Reseda' = 'h','Seseli' = 'h','Chusquea' = 's',
                'Eriogonum' = 'h','Herniaria' = 'h','Melica' = 'h','Myrcia' = NA,
                'Penstemon' = 'h','Plectranthus' = 'h','Selago' = 'h','Tanacetum' = 'h',
                'Tillandsia' = 'h','Bomarea' = NA,'Boronia' = 's','Bromus' = 'h',
                'Echinops' = 'h','Euphrasia' = 'h','Gagea' = 'h','Hymenophyllum' = 'h',
                'Metalasia' = NA,'Phlox' = 'h','Primula' = 'h','Pultenaea' = 's',
                'Scirpus' = 'h','Sida' = 's','Sideritis' = 'h','Swainsona' = NA,
                'Brachyotum' = 's','Cleome' = NA,'Commelina' = 'h','Cuscuta' = 'h',
                'Cynanchum' = 'h','Espeletia' = 's','Halenia' = 'h','Lonicera' = 'h',
                'Muhlenbergia' = 'h','Puya' = NA,'Sesleria' = 'h','Trigonella' = NA, 
                'Ageratina' = 's','Bartsia' = 'h','Blechnum' = 'h','Chlorophytum' = 'h',
                'Conyza' = 'h','Cousinia' = 'h','Dendrobium' = 'h','Ehrharta' = 'h',
                'Epidendrum' = 'h','Epipactis' = 'h','Helianthemum' = 'h','Hibbertia' = 't',
                'Inula' = 'h','Isoetes' = 'h','Lotononis' = NA,'Oenothera' = 'h',
                'Rhododendron' = 's','Rhynchosia' = 'h','Rorippa' = 'h','Sorbus' = 't',
                'Stylidium' = 'h','Symplocos' = 't','Valerianella' = 'h','Asarum' = 'h',
                'Boechera' = 'h','Bulbophyllum' = 'h','Drimia' = 'h','Drosanthemum' = NA,
                'Gnaphalium' = 'h','Liatris' = 'h')

genera.gf <- c(genera.gf, genera.gf2)

tmp3 <- tmp %>% 
  mutate(species2=species) %>% 
  separate(species2, into = c("genus", "sp")) %>% 
  left_join(data.frame(genus=names(genera.gf), gf=genera.gf), 
            by="genus") %>% 
  mutate(is.tree.or.tall=ifelse(gf=="t", T, F)) %>% 
  mutate(is.shrub=ifelse(gf=="s", T, F))

DT2 <- DT2 %>% 
  left_join(tmp3 %>% 
              dplyr::select(species, is.tree.or.tall, is.shrub), 
            by="species") %>% 
  mutate(is.tree.or.tall=coalesce(is.tree.or.tall.x, is.tree.or.tall.y)) %>% 
  mutate(is.shrub=coalesce(is.shrub.x, is.shrub.y)) %>% 
  dplyr::select(-is.tree.or.tall.x, -is.tree.or.tall.y, -is.shrub.x, -is.shrub.y, -plants_recorded2) 

tmp <- DT2 %>% distinct(species, is.tree.or.tall, is.shrub)
table(tmp$is.tree.or.tall, tmp$is.shrub, exclude=NULL)

"        FALSE  TRUE  <NA>
  FALSE 24841  5155     0
  TRUE  10242  1784     0
  <NA>      0  3256 12610"

save(DT2, file = "../sPlot/_derived_data/DT2_20161025_filtered2.RData")


