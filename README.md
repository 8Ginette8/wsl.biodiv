# wsl_biodiv

Developing a toolbox to efficiently handle SDMs in collaborative projects. Package aims at providing varied and flexible model fitting and tools for SDMs. Meta information should be stored (author, date, taxon,...) efficiently. Model evaluation and prediction should flexibly handle various fitted model objects. Evaluation metrics for presence-only and presence-absence models should be mutually implemented. Filtering and testing tools should help the user to find adequate predictors to apply meaningful models. Poisson Point Process models (PPPM) implementation should remain user-friendly. Regularization and variable selection should be implemented in this framework. Efficient methods to correct sampling bias in model fitting should be implements (pseudo-absences sempling strategies, environmental bias correction, bias covariate correction...)

# Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("8Ginette8/wsl.biodiv")
library(wsl.biodiv)
```

# Example

## Ensemble models

``` r
...
```

## Point process models (PPMs) lasso

``` r
...
```

# Mandatory citations

Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O., Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species distribution projections under climate change. Journal of Biogeography, 47(1), 130-142. doi: <a href="https://doi.org/10.1111/jbi.13734">10.1111/jbi.13734</a>

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. doi: <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>
