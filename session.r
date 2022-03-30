source("A01_SpeciesPool_sPlot.R")

# ------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------
vegtype <- "for"

data     <- "/data/splot/_data/DT2_20161025.RData"
#dt_beals <- "DT_beals_sPlot.rData"
Mij      <- "/data/splot/_data/Mij_tot_AllPlots_sparse_2019.rData"
header   <- "/data/splot/_data/sPlot_header_20161124.RData"
output   <- paste("/data/splot/_output/test_result6.splot", Sys.Date(), vegtype,"csv", sep=".")
#subset   <- "/data/splot/_data/BufferedRandomSubset20km.RData"
subset <- ifelse(vegtype=="for", 
                 "/work/sabatini/subset20km/Forests99SpeciesPool_sPlot-4926365.RData", 
                 "/work/sabatini/subset20km/NonForests99SpeciesPool_sPlot-4926378.RData"
                 )

ncores   <- detectCores()
verbose  <- TRUE
nrows    <- 0
index    <- 1


# ------------------------------------------------------------------------------
# actual program
# ------------------------------------------------------------------------------

SpeciesPool_sPlot(data, Mij, header, subset, output, ncores, verbose, nrows, index, vegtype)




############### same but for iDiv rStudio SERVER
source("A01_SpeciesPool_sPlot.R")

# ------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------

vegtype <- "nonfor"

data     <- "../sPlot/_derived_data/DT2_20161025_filtered2.RData"
header   <- "/data/sPlot/releases/sPlot2.1/sPlot_header_20161124.RData"
Mij      <- "../sPlot/_derived_data/Mij_tot_AllPlots_sparse_filtered.rData"
output   <- paste("../sPlot/_derived_data/Resample1/result6.splot", Sys.Date(), vegtype, "RData", sep=".")
subset <- ifelse(vegtype=="for", 
                 "../sPlot/_derived_data/Resample1/Resample1.for.RData",
                 "../sPlot/_derived_data/Resample1/Resample1.nonfor.RData")

#nthreads <- detectCores()
ncores <- 6 #detectCores()
verbose  <- TRUE
nrows    <- 0
index    <- 1
# ------------------------------------------------------------------------------
# actual program
# ------------------------------------------------------------------------------

SpeciesPool_sPlot(data, Mij, header, subset, output, ncores, verbose, nrows, index, vegtype)
##############################















### RUN RESAMPLING
source("A99_ResamplingByRadius.R")
# -----"-------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------

header    <- "/data/splot/_data/sPlot_header_20161124.RData"
output    <- "/data/splot/_data/test.RData"

#for idiv rstudio
header   <- "/data/sPlot/releases/sPlot2.1/sPlot_header_20161124.RData"
output   <- "../sPlot/_output/test.RData"



t.radius <- 20000
reps     <- 5
cores    <- detectCores()
verbose  <- T
vegtype  <- "for"
# ------------------------------------------------------------------------------
# actual program
# ------------------------------------------------------------------------------

res <- res.sub(header, output, t.radius, reps, cores, vegtype, verbose)




### A03 - BRT
source("A03_BRTs_fornonf.R")
# ------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------

fornonf     <- "all"
mydata      <- paste("../sPlot/_derived_data/Resample1/Mydata_global_NatCommR1.RData", sep="")
world.data  <- "../sPlot/_derived_data/world.over.RData"
output      <- paste("../sPlot/_derived_data/Resample1/", "BRTglobal_", fornonf, "_", sep="")
verbose     <- T
nrows       <- 17#99
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
  

