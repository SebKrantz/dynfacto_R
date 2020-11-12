### Estimator for dynamic factor models

Most of current factor models in macroeonomics are only available in `Matlab` language. This repository contains the equivalent of a small sample of these models rewritten in `R`. Although state space models, of which dynamic factor models are a subset, are already available in `R` in a number of packages, e.g. `MARSS`, it is useful to have this in `R` for speed, understanding and a few other reasons.

Initially, the code was just a line-by-line rewrite in `R` of replication files from Doz, Gianone and Reichlin (2011) available to download from [here](http://homepages.ulb.ac.be/~dgiannon/). From a wider perspective, it will aim to refactor the code as well to make it more readable and maintainable as well as add an estimation option in case of missing data based on Banbura and Modugno (2010). Once it's done, it will be packaged in a simple `R` library.

Original link is not available anymore. However, replication files for `Matlab` have been preserved and are available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=hdl:1902.1/21827).

#### Current status

Currently, the development branch for arbitrary `q`, `r` and `p` as long as `q <= r` and if matrix inversion remains stable. Extensive numerical tests have not been run yet to ensure correctness, but a few tests of the original code on `Octave` on a sample data yield the same results, so it's rather promising. Original code doesn't seem to work when `q != r` and so further investigation is needed.

`em_functions.R` file contains all the functions called for factor estimation, namely Kalman filtering, Kalman smoothing and EM maximimation steps. `em.R` is a wrapper script to estimate the full model and was originally a global `DynFA` function.

Code in the development branch is more than a rewrite. It aims to simplify and improve the notation and readability.

```
devtools::install_github("rbagd/dfm-estimator")
library(dynfactoR)
```
