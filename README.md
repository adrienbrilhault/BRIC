<!-- badges: start -->

[![R-CMD-check](https://github.com/adrienbrilhault/BRIL/workflows/R-CMD-check/badge.svg)](https://github.com/adrienbrilhault/BRIC/actions)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](pkg/LICENSE)

<!-- badges: end -->

# BRIC

[BRIC](https://github.com/adrienbrilhault/BRIC) (Bootstrap and Refine
Iterative Clustering) is an [R](https://www.r-project.org) package for
robust data clustering.

The BRIC algorithm searches for the different clusters through three
principal steps:

-   *BOOTSTRAP:* a recursive depth (or convex body minimizers) trimming
    to locate a first central location estimate.
-   *REFINE:* a two-pass outliers filtering, the first relying on
    euclidean distances to the first estimate and unimodality tests, the
    second on robust distances and multinormality tests.
-   *ITERATE:* after removing the samples selected in the REFINE step
    from the global distribution, the same process is reapplied to
    search for additional clusters.

To browse and test this package in an online RStudio environment, follow
the link below:

[![Binder](https://tinyurl.com/badgeRStudio)](https://mybinder.org/v2/gh/adrienbrilhault/BRIC/HEAD?urlpath=rstudio)

To cite this package, please refer to:

*Adrien Brilhault, Sergio Neuenschwander, and Ricardo Rios - A New
Robust Multivariate Mode Estimator for Eye-tracking Calibration -
[Behavior Research Methods, 2022](https://rdcu.be/cI9Pf)* <br><br>

### Installation

To install BRIC from the [GitHub
repository](https://github.com/adrienbrilhault/BRIC), add first the
[remotes](https://github.com/r-lib/remotes) package if missing.

``` r
install.packages("remotes")
```

Then proceed to the installation of the BRIL package and its
dependencies using the `install_github` function from
[remotes](https://github.com/r-lib/remotes).

``` r
remotes::install_github("adrienbrilhault/BRIC", subdir = "pkg")
```

<br>
*Note: As of May 2022, the
[OjaNP](https://cran.r-project.org/web/packages/OjaNP/) package - one of
BRIC's dependencies - has been removed from the main CRAN repository.
If missing, install it manually with the following command, and then launch the BRIC
package installation again:*

``` r
remotes::install_github("cran/OjaNP", subdir = "pkg")
```
<br>

### Documentation

The package documentation is available in
[BRIC-manual.pdf](https://github.com/adrienbrilhault/BRIC/raw/master/BRIC-manual.pdf),
at
[adrienbrilhault.github.io/BRIC/](https://adrienbrilhault.github.io/BRIC/),
and directly in the R environment with the help commands.

``` r
# main documentation page
help(package = "BRIC")

# clustBRIC() function help
?BRIC::clustBRIC()
```

<br>

### Example

If missing, install the
[mvtnorm](https://CRAN.R-project.org/package=mvtnorm) package from CRAN
with the command `install.packages("mvtnorm")`. Then try the following
example, which creates a sample multivariate mixture, and estimates its
main mode.

``` r
library(BRIC)
library(mvtnorm)

XY <- rbind(
  rmvnorm(300, c(0,0), diag(2)*3-1),
  rmvnorm(100, c(15,20), diag(2)),
  rmvnorm(150, c(-10,15), diag(2)*2-0.5),
  rmvnorm(200, c(5,5), diag(2)*200)
)

res <- clustBRIC(XY)

print(res)
plot(res)
```

<br>

### Multivariate Mode

To obtain the main mode of a distribution, as with the package
[adrienbrilhault.github.io/BRIL/](https://adrienbrilhault.github.io/BRIL/),
retrieve the center of the main cluster:

``` r
library(BRIC)
library(mvtnorm)

XY <- rbind(
  rmvnorm(300, c(0,0), diag(2)*3-1),
  rmvnorm(100, c(15,20), diag(2)),
  rmvnorm(150, c(-10,15), diag(2)*2-0.5),
  rmvnorm(200, c(5,5), diag(2)*200)
)

res <- clustBRIC(XY)

mode_coordinates <- res$clustersCenters[res$mainCluster, ]
```

<br>

#### Illustrative notebooks

- For a demonstration of our method with artificial data, please refer to the 
notebook [demo_ArtificialData.ipynb](demo_ArtificialData.ipynb) 
( **[▶️ run binder](https://mybinder.org/v2/gh/adrienbrilhault/BRIC/HEAD?filepath=demo_ArtificialData.ipynb)** )

*(note that the online viewer from Github sometimes needs one or two reload if the notebook was not rendered)*

To run this notebook interactively in JupyterLab, follow this link:

[![Binder](https://tinyurl.com/badgeJupyterLab)](https://mybinder.org/v2/gh/adrienbrilhault/BRIC/HEAD?urlpath=lab)
