
<!-- README.md is generated from README.Rmd. Please edit that file -->

# relMix

relMix makes relationship inference involving DNA mixtures with unknown
profiles and interprets DNA mixtures with related contributors. The main
function is the graphical user interface `relMixGUI`. A tutorial can be
found here: <https://gdorum.github.io/relMix/articles/relMix.html>

## Installation

Install from GitHub as follows:
``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")
#Install relMix from GitHub:
devtools::install_github("gdorum/relMix")
```

To provide pedigree plots, relMix uses the package `tkrplot`. However,
this package has some compatibility issues with MacOS and hence is not
included as a *hard* dependency. Users who wish to see pedigree plots in
the results screen have to install the package `tkrplot` manually with

``` r
install.packages("tkrplot")
```

The `tkrplot` package will be loaded by relMix automatially, so users do
not need to run `library("tkrplot")` in advance.
