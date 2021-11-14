library(optparse)

# ------------------------------------------------------------------------------
# defaults
# ------------------------------------------------------------------------------

default.ncores <- 32
default.verbose <- FALSE

# ------------------------------------------------------------------------------
# parsing arguments
# ------------------------------------------------------------------------------

options <- list (
  make_option(
		opt_str = c("-v", "--verbose"),
		action  = "store_true",
		default = default.verbose,
		help    = "print more output on what's happening"
	),

  make_option(
		opt_str = c("-n", "--nrows"),
		dest    = "nrows",
		type    = "integer",
		default = 0,
		help    = "number of rows from data to use, defaults to all of them",
		metavar = "4"
	),
	
	make_option(
	  opt_str = c("-f", "--fornonf"),
	  dest    = "fornonf",
	  type    = "character",
	  default = NA,
	  help    = "Which vegetation type? all, for or nonfor",
	  metavar = "4"
	), 
	make_option(
	  opt_str = c("-i", "--index"),
	  dest    = "index",
	  type    = "integer",
	  default = 0,
	  help    = "Which realization of resampled dataset",
	  metavar = "4"
	)
)

parser <- OptionParser(
       usage       = "Rscript %prog [options] data dt_beals header output",
       option_list = options,
       description = "\nan awesome R script",
       epilogue    = "use with caution, the awesomeness might slap you in the face!"
)

cli <- parse_args(parser, positional_arguments = 3)

# ------------------------------------------------------------------------------
# assign a few shortcuts
# ------------------------------------------------------------------------------

mydata     <- cli$args[1]
world.data <- cli$args[2]
output     <- cli$args[3]
verbose    <- cli$options$verbose
nrows      <- cli$options$nrows
fornonf    <- cli$options$fornonf
index      <- cli$options$index


# ------------------------------------------------------------------------------
# actual program
# ------------------------------------------------------------------------------

source("A03_BRTs_fornonf.R")


BRTs_sPlot(mydata, world.data, output, verbose, nrows, fornonf, index)
