### A03 - BRT
source("A03_BRTs_fornonf.R")
# ------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------

fornonf     <- "all"
mydata      <- paste("../Input//Mydata_global_NatCommR1.RData", sep="")
world.data  <- "../Input/world.over.RData"
output      <- paste("../sPlot/_derived_data/Resample1/", "BRTglobal_", fornonf, "_", sep="") #intermediate output
verbose     <- T
nrows       <- 99
ncores      <- 1#detectCores()

BRTs_sPlot(mydata, world.data, output, verbose, nrows, fornonf, ncores) 




### A05 - Validation
source("A05_BRT_validation.R")
# ------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------
mydata.path <- "../sPlot/_derived_data/Resample1/Mydata_global.RData"
world.data.path <- "../sPlot/_derived_data/world.over.RData"
brt.models.path <- "../sPlot/_derived_data/Resample1/BRTglobal"
iteration <- 5
output <- paste0("../sPlot/_derived_data/Resample1/cv.validation_", iteration, ".RData")
verbose <- T


validate.BRT(mydata.path, world.data.path, brt.models.path, output, verbose)



### A98 
source("A98_PredictorsExtract.R")
# ------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------

x.shp     <- "../sPlot/_derived_data/header.shp.RData"
myfunction<- "robust.mean"
output    <- "../sPlot/_derived_data/ecoreg.out.csv"
toextract <- "../Spatial/TEOW/official/wwf_terr_ecos.shp"
typp      <- "shp"
ncores    <- 2

PredExtr(x.shp, myfunction, output, toextract, typp, ncores)
  

