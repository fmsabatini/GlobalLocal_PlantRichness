
## Function 1 - SummarySE.R   #####

## define function to calculate summary statistics
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and 
## confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
##   adapted from Ryan Hope's function: 
##   https://www.rdocumentation.org/packages/Rmisc/versions/1.5/topics/summarySE




# summarySE function
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = .95, .drop = TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      length(x)
    }
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, median, and sd
  
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N      = length2(xx[[col]], na.rm=na.rm),
                           mean   = mean(xx[[col]], na.rm=na.rm),
                           median = median(xx[[col]], na.rm=na.rm),
                           sd      = sd(xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" and "median" columns    
  datac <- plyr::rename(datac, c("mean" = paste(measurevar, "_mean", sep = "")))
  datac <- plyr::rename(datac, c("median" = paste(measurevar, "_median", sep = "")))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



## Function 2 - Raincloud ####
### This script creates an R function to generate raincloud plots, then simulates
### data for plots. If using for your own data, you only need lines 1-80.
### It relies largely on code previously written by David Robinson
### (https://gist.github.com/dgrtwo/eb7750e74997891d7c20)
### and the package ggplot2 by Hadley Wickham

# Check if required packages are installed 
packages <- c("cowplot", "readr", "ggplot2", "dplyr", "lavaan", "smooth", "Hmisc")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

# Load packages 
library(ggplot2)

# Defining the geom_flat_violin function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )




## Function 3 - cracked plot.gbm ####
## function cracked in a way that the last variable (which should be Rel.area)
## is only predicted for the levels 10,100, 400, 1000, 10000

plot.gbm <- function (x, i.var = 1, n.trees = x$n.trees, continuous.resolution = 100, 
                      return.grid = FALSE, type = c("link", "response"), level.plot = TRUE, 
                      contour = FALSE, number = 4, overlap = 0.1, col.regions = viridis::viridis, 
                      ...) 
{
  type <- match.arg(type)
  if (all(is.character(i.var))) {
    i <- match(i.var, x$var.names)
    if (any(is.na(i))) {
      stop("Requested variables not found in ", deparse(substitute(x)), 
           ": ", i.var[is.na(i)])
    }
    else {
      i.var <- i
    }
  }
  if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
    warning("i.var must be between 1 and ", length(x$var.names))
  }
  if (n.trees > x$n.trees) {
    warning(paste("n.trees exceeds the number of tree(s) in the model: ", 
                  x$n.trees, ". Using ", x$n.trees, " tree(s) instead.", 
                  sep = ""))
    n.trees <- x$n.trees
  }
  if (length(i.var) > 3) {
    warning("plot.gbm() will only create up to (and including) 3-way ", 
            "interaction plots.\nBeyond that, plot.gbm() will only return ", 
            "the plotting data structure.")
    return.grid <- TRUE
  }
  grid.levels <- vector("list", length(i.var))
  for (i in 1:(length(i.var)-1)) {
    if (is.numeric(x$var.levels[[i.var[i]]])) {
      grid.levels[[i]] <- seq(from = min(x$var.levels[[i.var[i]]]), 
                              to = max(x$var.levels[[i.var[i]]]), length = continuous.resolution)
    }
    else {
      grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]], 
                                            levels = x$var.levels[[i.var[i]]])) - 1
    }
  }
  grid.levels[[i+1]] <- c(10,100,400, 1000, 10000)
  X <- expand.grid(grid.levels)
  names(X) <- paste("X", 1:length(i.var), sep = "")
  if (is.null(x$num.classes)) {
    x$num.classes <- 1
  }
  y <- .Call("gbm_plot", X = as.double(data.matrix(X)), cRows = as.integer(nrow(X)), 
             cCols = as.integer(ncol(X)), n.class = as.integer(x$num.classes), 
             i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees), 
             initF = as.double(x$initF), trees = x$trees, c.splits = x$c.splits, 
             var.type = as.integer(x$var.type), PACKAGE = "gbm")
  if (x$distribution$name == "multinomial") {
    X$y <- matrix(y, ncol = x$num.classes)
    colnames(X$y) <- x$classes
    if (type == "response") {
      X$y <- exp(X$y)
      X$y <- X$y/matrix(rowSums(X$y), ncol = ncol(X$y), 
                        nrow = nrow(X$y))
    }
  }
  else if (is.element(x$distribution$name, c("bernoulli", "pairwise")) && 
           type == "response") {
    X$y <- 1/(1 + exp(-y))
  }
  else if ((x$distribution$name == "poisson") && (type == "response")) {
    X$y <- exp(y)
  }
  else if (type == "response") {
    warning("`type = \"response\"` only implemented for \"bernoulli\", ", 
            "\"poisson\", \"multinomial\", and \"pairwise\" distributions. ", 
            "Ignoring.")
  }
  else {
    X$y <- y
  }
  f.factor <- rep(FALSE, length(i.var))
  for (i in 1:length(i.var)) {
    if (!is.numeric(x$var.levels[[i.var[i]]])) {
      X[, i] <- factor(x$var.levels[[i.var[i]]][X[, i] + 
                                                  1], levels = x$var.levels[[i.var[i]]])
      f.factor[i] <- TRUE
    }
  }
  names(X)[1:length(i.var)] <- x$var.names[i.var]
  if (return.grid) {
    return(X)
  }
  nx <- length(i.var)
  if (nx == 1L) {
    plotOnePredictorPDP(X, ...)
  }
  else if (nx == 2) {
    plotTwoPredictorPDP(X, level.plot = level.plot, contour = contour, 
                        col.regions = col.regions, ...)
  }
  else {
    plotThreePredictorPDP(X, nx = nx, level.plot = level.plot, 
                          contour = contour, col.regions = col.regions, number = number, 
                          overlap = overlap, ...)
  }
}


## Function 4 - label_parse ####
## helper for plotting labels with subscripts
label_parse <- function(breaks) {
  parse(text = breaks)
}
