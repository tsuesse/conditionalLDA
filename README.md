---
# Introduction

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Here we give an introduction to the package conditionalLDA in R! 
We follow the methods described in: 

Suesse, T., Grupp. V., Brenning, A. "Spatial Linear Discriminant Analysis Approaches for Remote-Sensing Classification" (2023)


We do a sample analysis on the Aconcaguia data set.
The data in the package has been preprocessed. 

```{r setup}
# install.packages(devtools)
# library(devtools)
# devtools::install_github("jcortesr/PerMuTe")
library(conditionalLDA)
```

```{r}
image(temp_gistemp[,,1])
```

Aconcagua data set
d data 
formula formula for LDA
predictors predictors for LDA








# conditionalLDA

<!-- badges: start -->
<!-- badges: end -->

The goal of conditionalLDA is to ...

## Installation

You can install the development version of conditionalLDA from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tsuesse/conditonalLDA")
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(conditionalLDA)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

```{r cars}
summary(cars)
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

```{r pressure, echo = FALSE}
plot(pressure)
```

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
