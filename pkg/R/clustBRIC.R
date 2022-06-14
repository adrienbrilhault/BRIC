
suppressPackageStartupMessages(require("CompQuadForm", quietly = TRUE))

#' Bootstrap and Refine Iterative Clustering
#'
#' Robust clustering algorithm based on depth measures and convex body
#' minimizers
#'
#' @param data Matrix or Data-Frame of numerical values containing the
#'   observations (rows correspond to observations, columns to variables)
#' @param maxIterations Maximum number of iterations performed by the algorithm
#'   (i.e. max number of potential clusters encountered). Set to NULL or 0 for
#'   unlimited number (default)
#' @param minUnassigned Numerical value between 0 and 1 (default: 0.1),
#'   providing the proportion of unassigned samples from `data` under which the
#'   algorithm will terminate the search
#' @param nsamp Number of samples randomly selected from `data` for subsampling
#'   calculations, or "best", "exact" or "sample". If "sample" is chosen, the
#'   subset will include up to 2000 observations; with "best" up to 4000
#'   (default); with "exact" (or 0), exhaustive search will be attempted on the
#'   complete dataset (computation in this case might take a long time). When
#'   subsampling, the remaining observations will be assigned to the cluster of
#'   their closest neighbor.
#' @param method Method to use. Valid options are "MCD" and "MVE" for convex
#'   body minimizers, or "L2", "Lui", "Mahalanobis", "Oja", "Projection"
#'   (default), "Spatial" and "Tukey" for depth functions
#' @param alpha Proportion of samples trimmed at each iteration of the recursive
#'   median estimate (numerical value between 0 and 1, default: 0.5), see
#'   [median_rec()]
#' @param testUnimodal Statistical test used for unimodality. Valid options are
#'   "DIP" or a user-defined function (see [filter_outliers()])
#' @param threshUnimodal Threshold of significance for the unimodality test
#'   (numerical value between 0 and 1, default: 0.05)
#' @param distUnimodal Distance metric used for ordering the samples in the
#'   unimodal filtering. Valid options are "Euclidean" (default), or "MCD",
#'   "MVE", and "OGK" for robust distances. "Euclidean" is strongly advised for
#'   unimodality tests.
#' @param testNormal Statistical test used for normality. Valid options are
#'   "Mardia" (default), "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro",
#'   "Lillie", "Chisq", or a user-defined function (see [filter_outliers()])
#' @param threshNormal Threshold of significance for the normality test
#'   (numerical value between 0 and 1, default: 0.05)
#' @param distNormal Distance metric used for ordering the samples in the normal
#'   filtering. Valid options are "Euclidean", or "MCD" (default), "MVE", and
#'   "OGK" for robust distances. Robust distances are strongly advised for
#'   normality tests.
#' @param trimmedPerFilteringIteration Number of samples trimmed at each
#'   iteration of the unimodality and normality filtering (default: 1), see
#'   [filter_outliers()]
#' @param exitWhenUnimodal Logical value. `TRUE` will terminate the execution of
#'   the algorithm as soon as an unimodal subset is encountered on the start of
#'   a global iteration. `FALSE` (default) will let that last iteration proceed
#'   before terminating
#' @param debug Logical value. `TRUE` will compute all p.values in the filtering
#'   steps (even after they exceed the selection threshold, see
#'   [plot.BRIC.Filtering()])
#' @param warnings Logical value, to display the warnings and errors caught
#'
#' @return The function returns an S3 object of type `BRIC` containing the
#'   following values:
#'   \item{`call`}{Parameters of the call (contains `data`,
#'   `maxIterations`, `minUnassigned`, `nsamp`, `method`, `alpha`, `testUnimodal`,
#'   `threshUnimodal`, `distUnimodal`, `testNormal`, `threshNormal`,
#'   `distNormal`, `trimmedPerFilteringIteration`, and `exitWhenUnimodal`)}
#'   \item{`iterations`}{A list with every global iteration of the algorithm,
#'   each containing the two filtering procedures performed: `filteringUnimodal`
#'   and `filteringNormal` (both being S3 object of class `BRIC.Filtering`, see
#'   [filter_outliers()])}
#'   \item{`nbClusters`}{Number of groups encountered}
#'   \item{`labels`}{Labels of the groups encountered (corresponding to the
#'   number of the iteration they were identified in)}
#'   \item{`clustersCenters`}{Matrix containing the coordinates of the centers
#'   of each group (row-wise)}
#'   \item{`clustersSizes`}{Array with the number of
#'   samples in each group}
#'   \item{`mainCluster`}{Index of the group identified as main mode}
#'
#' @references Adrien Brilhault, Sergio Neuenschwander, and Ricardo Rios - A New
#'   Robust Multivariate Mode Estimator for Eye-tracking Calibration - Behavior
#'   Research Methods, 2022 - \href{https://rdcu.be/cI9Pf}{rdcu.be/cI9Pf}
#'
#' @seealso [plot.BRIC()], [print.BRIC()], [filter_outliers()], [median_rec()],
#'   [median_mv()], [depth_values()]
#'
#' @examples
#'
#' # Create a sample distribution and run clustBRIC() function
#' data <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#'
#' res <- clustBRIC(data, debug = TRUE)
#' print(res)
#'
#' # Plot the mode and groups encountered
#' plot(res)
#'
#' # Plot the different iterations (interactive)
#' \dontrun{
#' plot(res, contents = "iterations", asp = 1)
#' }
#'
#' # See ?plot.BRIC() for other plotting examples
#'
#' @export
#'
clustBRIC <- function(data, maxIterations = 0, minUnassigned = 0.1, nsamp = "best",
                 method = "Projection", alpha = 0.5,
                 testUnimodal = "DIP", threshUnimodal = 0.05, distUnimodal = "Euclidean",
                 testNormal = "Mardia", threshNormal = 0.05, distNormal = "MCD",
                 trimmedPerFilteringIteration = 1, exitWhenUnimodal = FALSE,
                 debug = FALSE, warnings = FALSE) {

  DEV_DEBUG <- FALSE

  ######## Check parameters & initialize

  if (missing(data)) {
    stop("Parameter `data` containing the observations is mandatory")
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("Parameter `data` must contain numerical values only")
  }
  if (is.null(ncol(data)) || ncol(data) < 2) {
    stop("Parameter`data` must contain at least 2 columns, and more than one row")
  }
  if (!is.null(maxIterations) && (!is.numeric(maxIterations) || maxIterations < 0)) {
    stop("Parameter `maxIterations` must be either NULL or a positive entire number")
  }
  if (!is.numeric(minUnassigned) || minUnassigned < 0 || minUnassigned > 1) {
    stop("Parameter `minUnassigned` must be a numerical value between 0 and 1 (included)")
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Parameter `alpha` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.numeric(threshUnimodal) || threshUnimodal < 0 || threshUnimodal > 1) {
    stop("Parameter `threshUnimodal` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.numeric(threshNormal) || threshNormal < 0 || threshNormal > 1) {
    stop("Parameter `threshNormal` must be a numerical value between 0 and 1 (excluded)")
  }
  if (!is.numeric(trimmedPerFilteringIteration) || trimmedPerFilteringIteration < 1) {
    stop("Parameter `trimmedPerFilteringIteration` must be an entire number superior or equal to 1")
  }
  trimmedPerFilteringIteration <- round(trimmedPerFilteringIteration)
  if (!is.logical(debug)) {
    stop("Parameter `debug` must be a logical value")
  }
  if (!is.logical(warnings)) {
    stop("Parameter `warnings` must be a logical value")
  }

  opts_test <- c("DIP", "Mardia", "Kurtosis", "Skewness", "KS", "KS-adj", "Shapiro", "Lillie", "Chisq")
  if (is.character(testUnimodal) && (tolower(testUnimodal) %in% tolower(opts_test))) {
    testUnimodal <- opts_test[which(tolower(testUnimodal) == tolower(opts_test))]
  } else if (!is.function(testUnimodal)) {
    stop("Parameter `testUnimodal` should be a function or one of \"", paste(opts_test, collapse = "\", \""), "\"")
  }

  if (is.character(testNormal) && (tolower(testNormal) %in% tolower(opts_test))) {
    testNormal <- opts_test[which(tolower(testNormal) == tolower(opts_test))]
  } else if (!is.function(testNormal)) {
    stop("Parameter `testNormal` should be a function or one of \"", paste(opts_test, collapse = "\", \""), "\"")
  }

  opts_dist <- c("Euclidean", "MCD", "MVE", "OGK")
  if (missing(distUnimodal) || is.null(distUnimodal)) {
    distUnimodal <- "Euclidean"
  } else {
    if (is.character(distUnimodal) && tolower(distUnimodal) %in% tolower(opts_dist)) {
      distUnimodal <- opts_dist[which(tolower(distUnimodal) == tolower(opts_dist))]
    } else {
      stop("Parameter `distUnimodal` should be NULL or one of \"", paste(opts_dist, collapse = "\", \""), "\"")
    }
  }

  if (missing(distNormal) || is.null(distNormal)) {
    distNormal <- "MCD"
  } else {
    if (is.character(distNormal) && tolower(distNormal) %in% tolower(opts_dist)) {
      distNormal <- opts_dist[which(tolower(distNormal) == tolower(opts_dist))]
    } else {
      stop("Parameter `distNormal` should be NULL or one of \"", paste(opts_dist, collapse = "\", \""), "\"")
    }
  }

  if (missing(maxIterations) || is.null(maxIterations)) {
    maxIterations <- 0
  } else {
    maxIterations <- round(maxIterations)
  }

  subsampling <- FALSE

  if (nrow(data) > 5000 && (is.null(nsamp) || nsamp == 0 || nsamp == "exact")) {
    warning("Large number of samples (", nrow(data), "), computations might take a long time,
    consider using nsamp=\"best\" or nsamp=\"sample\"")
  }
  if (is.numeric(nsamp)) {
    if (nsamp < 0) {
      stop("Parameter `nsamp` should be NULL, or a round number superior or equal to 0, or \"sample\", \"best\", or \"exact\".")
    } else if (nsamp > 0 && nsamp < nrow(data) ){
      if (nsamp < 2000 && (nsamp < (ncol(data) * 300)) && (nrow(data) > 1.5 * ncol(data) * 300)) {
        warning("Parameter `nsamp` seems particulary small (nsamp=", nsamp, "), consider using at least ", ncol(data) * 300, ".")
      }
      subsampling <- TRUE
    }
  } else if (is.character(nsamp)) {
    nsamp = tolower(nsamp)
    if (!nsamp %in% c("exact","sample","best")) {
      stop("Parameter `nsamp` should be NULL, or a round number superior or equal to 0, or \"sample\", \"best\", or \"exact\".")
    }
    if (nsamp == "sample" && nrow(data) > 2000) {
      subsampling <- TRUE
      nsamp <- 2000
    } else if (nsamp == "best" && nrow(data) > 5000) {
      subsampling <- TRUE
      nsamp <- 5000
    }
  }

  if (subsampling == TRUE) {
    indexSubset <- sample(1:nrow(data), nsamp, replace = FALSE)
    dataOriginal <- data
    data <- data[indexSubset,]
  }

  ## Internal Parameter
  reorderClustersBySize <- FALSE
  recenterAfterUnimodal <- TRUE
  exitWhenUnimodal <- FALSE
  discardOutliersWithSupsampling <- FALSE

  # Output of the function
  output <- list()
  class(output) <- "BRIC"
  output$call$data <- data
  output$call$maxIterations <- maxIterations
  output$call$minUnassigned <- minUnassigned
  output$call$nsamp <- nsamp
  output$call$subsampling <- subsampling
  output$call$method <- method
  output$call$alpha <- alpha
  output$call$testUnimodal <- testUnimodal
  output$call$threshUnimodal <- threshUnimodal
  output$call$distUnimodal <- distUnimodal
  output$call$testNormal <- testNormal
  output$call$threshNormal <- threshNormal
  output$call$distNormal <- distNormal
  output$call$trimmedPerFilteringIteration <- trimmedPerFilteringIteration
  output$call$exitWhenUnimodal <- exitWhenUnimodal

  # Exit if not enough observations to compute depth values
  minSampleSize <- ncol(data) + 1
  if (nrow(data) < minSampleSize) {
    if (warnings) {
      warning(sprintf(
        "depth_values(method=\"%s\"): sample size of `data``is %d but should be greater than %d, NULL returned",
        method, nrow(data), minSampleSize
      ))
    }
    return(NULL)
  }

  # Initialization
  iterationNB <- 1
  output$iterations <- list()
  output$nbClusters <- 0
  output$clustersCenters <- c()
  output$clustersSizes <- c()
  output$labels <- rep(0, dim(data)[1])

  ######## Main loop, selecting a new group at each iteration

  while ((maxIterations == 0 || iterationNB <= maxIterations) &&
         (sum(output$labels == 0) > max(minSampleSize, minUnassigned * nrow(data)))) {
    unassignedIndices <- which(output$labels == 0)

    if (DEV_DEBUG)
      cat("BRIC - Iteration",iterationNB, "(nb samples:",nrow(output$call$data),"- unprocessed:",length(unassignedIndices),"\n")

    ## Bootstrap : compute first estimate

    bootstrapEstimate <- median_rec(data[unassignedIndices, ], method = method, alpha = alpha, warnings = warnings)$median

    if (DEV_DEBUG)
      cat("* Bootstrap estimate (median_rec): ",bootstrapEstimate,"\n")

    ## Refine estimate

    # Filter furthest outliers until reaching an unimodal subset
    filteringUnimodal <- filter_outliers(
      data[unassignedIndices, ],
      center = bootstrapEstimate, test = testUnimodal,
      threshold = threshUnimodal, distType = distUnimodal, debug = debug, warnings = warnings
    )
    output$iterations[[iterationNB]] <- list("filteringUnimodal" = filteringUnimodal)
    unimodalIndices <- unassignedIndices[filteringUnimodal$selected]

    if (DEV_DEBUG)
      cat("* Filtering unimodal:", length(filteringUnimodal$selected), "selected","\n")

    # If the current distribution, before any filtering, was already unimodal,
    # terminate the overall recursive search now or at the end of this iteration
    if (utils::head(filteringUnimodal$p.values, 1) > threshUnimodal) {
      if (exitWhenUnimodal == TRUE && iterationNB != 1) {
        maxIterations <- iterationNB
        break
      } else {
        maxIterations <- iterationNB
      }
    }

    # Re-estimate central location within this subset
    if (recenterAfterUnimodal == TRUE) {
      bootstrapEstimate <- median_rec(
        data[unimodalIndices, ], method = method, alpha = alpha, warnings = warnings)$median
    }

    if (DEV_DEBUG)
      cat("* Adjusted center:", bootstrapEstimate,"\n")

    # Filter furthest outliers until reaching a normal subset
    filteringNormal <- filter_outliers(
      data[unimodalIndices, ],
      center = bootstrapEstimate, test = testNormal,
      threshold = threshNormal, distType = distNormal, debug = debug, warnings = warnings
    )
    output$iterations[[iterationNB]]$filteringNormal <- filteringNormal
    normalIndices <- unimodalIndices[filteringNormal$selected]

    # Estimate the center of this group (refined estimate)
    center <- colMeans(data[normalIndices, ])

    if (DEV_DEBUG) {
      cat("* Filtering normal:", length(filteringNormal$selected), "selected","\n")
      cat("* Cluster final center:", center,"\n")
    }

    # Update the groups structures
    output$labels[normalIndices] <- iterationNB
    output$nbClusters <- output$nbClusters + 1
    output$clustersCenters <- rbind(output$clustersCenters, center)
    output$clustersSizes <- c(output$clustersSizes, length(normalIndices))
    rownames(output$clustersCenters) <- NULL

    iterationNB <- iterationNB + 1

    # # If the size of the remaining samples is inferior to the largest
    # # group encountered so far, it will not be the main mode, we can exit
    # if (sum(output$labels == 0) < max(output$clustersSizes))
    #   break
  }

  ### Select main cluster/mode

  if (output$nbClusters == 1) {
    output$mainCluster <- 1
    # output$mode <- output$clustersCenters
  } else {
    output$mainCluster <- which.max(utils::head(output$clustersSizes, -1))
    # output$mode <- output$clustersCenters[output$mainCluster, ]
  }

  ### Process samples excluded from the random subset with nsamp

  if (subsampling == TRUE) {

    classes <- sort(unique(output$labels[output$labels>0]))

    ## Compute robust scatter estimate for each cluster
    if (discardOutliersWithSupsampling) {

      covClusters <- list()

      for (i in classes) {

        covCluster <- withCallingHandlers(
          withRestarts({
            robustScatter <- robustbase::covMcd(data, tolSolve = 1e-20)

            # Check values are correct if it restarts from a warning
            if (!exists("robustScatter") || !("cov" %in% names(robustScatter)) ||
                is.null(robustScatter$cov) || all(robustScatter$cov == 0)) {
              covCluster <- NULL
              if (warnings) {
                warning("In clustBRIC(), robust scatter estimate could not be computed",
                        call. = FALSE
                )
              }
            } else {
              covCluster <- robustScatter$cov
            }
            covCluster
          },
          muffleError = function() NULL
          ),
          warning = function(w) {
            if (warnings && !grepl("loss of accuracy", w, fixed = TRUE)) {
              warning("Warning caught in clustBRIC(), when computing clusters' robust scatter estimate:\n*  ", w,
                      call. = FALSE
              )
            }
            invokeRestart("muffleWarning")
          },
          error = function(e) {
            if (warnings) {
              warning("Error caught in clustBRIC(), when computing clusters' robust scatter estimate:\n*  ", e,
                      call. = FALSE
              )
            }
            invokeRestart("muffleError")
          }
        )
        covClusters[[i]] <- covCluster
      }
    }

    output$call$data <- dataOriginal
    output$call$indexSubset <- indexSubset

    labelsComplete <- integer(nrow(dataOriginal))
    labelsComplete[indexSubset] <- output$labels
    indexUnclassified <- setdiff(1:nrow(dataOriginal), indexSubset)

    for (i in indexUnclassified){
      smallestDistance <- NULL
      for (j in indexSubset){
        distance <- sqrt(sum((dataOriginal[i,] - dataOriginal[j,]) ^ 2))
        if (is.null(smallestDistance) || distance < smallestDistance) {
          smallestDistance <- distance
          labelsComplete[i] <- labelsComplete[j]
        }
      }

      # Remove samples assigned by closest neighbor to one of the main clusters
      # if their robust distance to the cluster center exceed a threshold
      if (discardOutliersWithSupsampling) {
        clustNb <- labelsComplete[i]
        if (clustNb != 0 && length(covClusters) >= clustNb
            && !is.null(covClusters[[clustNb]])) {
          if (sqrt(stats::mahalanobis(dataOriginal[i,], output$clustersCenters[clustNb,], covClusters[[clustNb]])) > 3) {
            labelsComplete[i] <- 0
          }
        }
      }
    }

    output$labels <- labelsComplete
    for (i in 1:length(classes)) {
      output$clustersSizes[i] <- length(which(output$labels==classes[i]))
    }
  }

  return(output)
}


#' Print method for `BRIC` objects
#'
#' @param x An object of class `BRIC` (see [clustBRIC()])
#' @param maxDisplayed Number of elements to display in the output. Set to NULL
#'   (or 0) to show all values.
#' @param ... Other arguments passed to or from other methods
#'
#' @seealso [clustBRIC()], [plot.BRIC()]
#'
#' @export
#'
print.BRIC <- function(x, maxDisplayed = NULL, ...) {
  if (class(x) != "BRIC") {
    stop("the object provided is not of class \"BRIC\"")
  }

  if (is.null(maxDisplayed) || maxDisplayed == 0) {
    maxDisplayed <- length(x$labels)
  } else if (!is.numeric(maxDisplayed) || maxDisplayed < 0) {
    stop("Parameter `maxDisplayed` must be a numerical value superior or equal to 1")
  }
  maxDisplayed <- round(maxDisplayed)

  cat("\n=> Results for clustBRIC() using method \"", x$call$method, "\" (alpha=", x$call$alpha,
      "), ", ifelse(is.function(x$call$testUnimodal),"User-defined",x$call$testUnimodal),
      " Unimodality Test (> ", x$call$threshUnimodal,"), and ",
      ifelse(is.function(x$call$testNormal),"User-defined",x$call$testNormal),
      " Normality test (> ", x$call$threshNormal, ")\n", sep = ""
  )

  cat("   ", nrow(x$call$data), " samples: ", x$nbClusters, " clusters identified (sizes ",
      toString(x$clustersSize), ")",
      sep = ""
  )
  if (sum(x$labels == 0) > 0) {
    cat(", and ", sum(x$labels == 0), " samples unassigned", sep = "")
  }

  cat("\n\nClusters Sizes:\n")
  print(x$clustersSizes, ...)

  cat("\n\nClusters Centers:\n")
  print(x$clustersCenters, ...)

  cat("\n\nLabels:\n")
  print(utils::head(x$labels, maxDisplayed), ...)
  if (length(x$labels) > maxDisplayed) {
    cat(paste0(" ... (", length(x$labels) - maxDisplayed, " hidden)"))
  }

  cat("\n\n")

  invisible(x)
}


#' Plot method for `BRIC` objects
#'
#' @param x An object of class `BRIC` (see [clustBRIC()])
#' @param contents Contents to be displayed, options are "scatterplot", or
#'   "iterations" (only one option possible)
#' @param showCenters Logical value used when `contents = "scatterplot"`, to
#'   show or not the clusters' center
#' @param col Default color of samples
#' @param colCenters Color of the center(s) when `contents = "scatterplot"`
#' @param colClusters List (or array) of colors for each of the
#'   clusters/iterations (length must be at least equal to the number of groups
#'   identified by the function [clustBRIC()], i.e. `x$nbClusters`)
#' @param iterationsIndices Numerical value or array of numerical values, used
#'   when `contents = "iterations"`, which provides the indices of the
#'   iterations to be plotted. If more than one iteration is requested, an
#'   interactive menu in the console will be used for the selection. 0 or NULL
#'   (default) will include all the iterations. Values that are negative or
#'   superior to the number of iterations performed by the execution
#'   [clustBRIC()] will be ignored
#' @param iterationsOptions List of additional parameters to be passed to the
#'   [plot.BRIC.Filtering()] function when `contents = "iterations"` is selected
#'   (see [plot.BRIC.Filtering()] for details). Example: `iterationsOptions =
#'   list(xlab = NA, ylab = NA, contents = c("p.values", "scatterplot"), asp =
#'   1)`
#' @param ... Other arguments passed to or from other methods (such as pch for
#'   the symbols, main and sub for title and subtitle, xlab, xmin, ...)
#'
#' @seealso [clustBRIC()], [print.BRIC()], [filter_outliers()],
#'   [print.BRIC.Filtering()]
#'
#' @examples
#'
#' # Create a sample distribution and run clustBRIC() function
#' data <- rbind(
#'   mvtnorm::rmvnorm(300, c(0, 0), diag(2) * 3 - 1),
#'   mvtnorm::rmvnorm(100, c(15, 20), diag(2)),
#'   mvtnorm::rmvnorm(150, c(-10, 15), diag(2) * 2 - 0.5),
#'   mvtnorm::rmvnorm(200, c(5, 5), diag(2) * 200)
#' )
#' res <- clustBRIC(data, debug = TRUE)
#'
#'
#' # Plot the clusters
#' plot(res)
#'
#' # Plot the clusters with extra graphic options
#' plot(res, showCenters = TRUE, main = "Multivariate Mode Estimate",
#'      col = "gray", colCenters = "black", colClusters = c("yellow","cyan",
#'      "purple","red"), asp = 1, pch = 3 )
#'
#' # Plot the second iteration
#' plot(res, contents = "iteration", iterationsIndices = 2)
#'
#' # Plot the second iteration (with arguments to plot.BRIC.filtering())
#' plot(res, contents = "iteration", iterationsIndices = 2,
#'      iterationsOptions = list(
#'        contents = c("scatterplot", "p.values"),
#'        colSelection = "blue", mfrow = c(2,1), asp = 1))
#'
#' \dontrun{
#' # Plot all iterations (interactive mode)
#' plot(res, contents = "iterations")
#'
#' # Plot the 3 first iterations with options (interactive mode)
#' plot(res, contents = "iterations", iterationsIndices = c(1:3),
#'      iterationsOptions = list(
#'        contents = c("scatterplot"),
#'        xlim  = c(-50,50),  ylim  = c(-30,30), asp = 1))
#' }
#'
#' @importFrom graphics plot points
#' @export
#'
plot.BRIC <- function(x, contents = "plot", showCenters = FALSE, col, colCenters,
                      colClusters = NULL, iterationsIndices = NULL,
                      iterationsOptions = NULL, ...) {
  if (class(x) != "BRIC") {
    stop("the object provided is not of class \"BRIC\"")
  }

  # Read extra parameters
  otherArgs <- list(...)

  if (!is.null(iterationsOptions) && !is.list(iterationsOptions)) {
    warning("Parameter `iterationsOptions`, if provided, must be a list")
    iterationsOptions <- NULL
  }
  if (length(contents) > 1) {
    warning("Parameter `contents` can only have one value (among \"scatterplot\" and \"iterations\"")
    contents <- contents[1]
  }
  if (!is.character(contents)) {
    stop("Parameter `contents` must be a string (either \"scatterplot\" or \"iterations\"")
  }
  contents <- tolower(contents)

  if (!is.logical(showCenters)) {
    stop("Parameter `showCenters` must be a logical value")
  }
  if (is.null(iterationsIndices) || anyNA(iterationsIndices) || iterationsIndices == 0) {
    iterationsIndices <- 1:x$nbClusters
  }
  if (!is.numeric(iterationsIndices)) {
    stop("Parameter `iterationsIndices`, if provided, must be numerical value(s)")
  }
  iterationsIndices <- sort(unique(iterationsIndices))
  if (max(iterationsIndices) > x$nbClusters) {
    warning("Only ", x$nbClusters, " iterations, values above discarded")
    iterationsIndices <- iterationsIndices[iterationsIndices <= x$nbClusters]
  }
  if (min(iterationsIndices) < 1) {
    warning("Parameter `iterations` can not take negative values")
    iterationsIndices <- iterationsIndices[iterationsIndices > 0]
  }
  if (missing(col) || is.null(col)) {
    col <- "black"
  }
  if (missing(colCenters) || is.null(colClusters)) {
    colCenters <- "black"
  }

  myColours <- grDevices::rainbow(x$nbClusters)
  if (!missing(colClusters) && !is.null(colClusters)) {
    if (length(colClusters) >= x$nbClusters) {
      myColours <- colClusters
    } else if (length(colClusters) >= length(iterationsIndices)) {
      myColours[iterationsIndices] <- colClusters
    } else {
      warning(
        "Parameter `colClusters` require at least ", length(iterationsIndices),
        " elements, reverting to default colors"
      )
    }
  }

  # SCATTERPLOTS
  if (any(unlist(lapply(c("scater", "scatter", "plot"), grepl, contents)))) {
    do.call(plot, utils::modifyList(
      list(x = x$call$data, xlab = "X", ylab = "Y", col = col),
      otherArgs
    ))

    # plot clusters
    for (i in 1:x$nbClusters) {
      selection <- x$call$data[x$labels == i, ]
      do.call(points, utils::modifyList(
        list(selection[, 1], selection[, 2], col = myColours[i]),
        otherArgs
      ))
    }

    if (showCenters == TRUE) {
      points(x$mode[1], x$mode[2], col = colCenters, pch = 4, lwd = 2
      )
    }

    # ITERATIONS
  } else if (grepl("iter", contents)) {
    xlim <- c(min(x$call$data[, 1]), max(x$call$data[, 1]))
    ylim <- c(min(x$call$data[, 2]), max(x$call$data[, 2]))

    defaultOptions <- list(col = col)

    if (is.null(iterationsOptions) ||
        (is.list(iterationsOptions) && "contents" %in% names(iterationsOptions) &&
         (length(iterationsOptions$contents) > 1 || iterationsOptions$contents == "all"))) {
      defaultOptions <- append(defaultOptions, list(
        xlim = c(min(x$call$data[, 1]), max(x$call$data[, 1])),
        ylim = c(min(x$call$data[, 2]), max(x$call$data[, 2]))
      ))
    }

    if (!is.null(iterationsOptions)) {
      extraOptions <- utils::modifyList(iterationsOptions, otherArgs)
    } else {
      extraOptions <- otherArgs
    }

    if (length(iterationsIndices) == 1) {

      # Plot unimodal filtering
      do.call(plot, utils::modifyList(
        append(defaultOptions, list(
          x = x$iterations[[iterationsIndices]]$filteringUnimodal,
          colSelection = myColours[iterationsIndices]
        )),
        extraOptions
      ))

      # Plot normal filtering
      do.call(plot, utils::modifyList(
        append(defaultOptions, list(
          x = x$iterations[[iterationsIndices]]$filteringNormal,
          colSelection = myColours[iterationsIndices]
        )),
        extraOptions
      ))
    } else {

      # build options list
      options <- c()
      for (i in iterationsIndices) {
        options <- c(options, paste0("Iteration ", i, " - Unimodal Filtering"))
        options <- c(options, paste0("iteration ", i, " - Normal Filtering"))
      }

      # Plot the first figure and present menu menu
      choice <- 1
      repeat{
        if (choice %% 2 == 1) {
          do.call(plot, utils::modifyList(
            append(defaultOptions, list(
              x = x$iterations[[(choice + 1) %/% 2]]$filteringUnimodal,
              colSelection = myColours[((choice + 1) %/% 2)]
            )),
            extraOptions
          ))
        } else {
          do.call(plot, utils::modifyList(
            append(defaultOptions, list(
              x = x$iterations[[(choice + 1) %/% 2]]$filteringNormal,
              colSelection = myColours[((choice + 1) %/% 2)]
            )),
            extraOptions
          ))
        }

        # re-present menu waiting user choice
        cat("\nCURRENTLY DISPLAYED: ", options[choice], "\n\n")
        choice <- utils::menu(options, graphics = FALSE, title = "NEW PLOT SELECTION:")

        if(choice == 0){
          break
        }
      }
    }
  }
  else {
    stop("Value for parameter `contents` is not supported. Valid options are
         \"scatterplot\" or \"iterations\"")
  }
}
