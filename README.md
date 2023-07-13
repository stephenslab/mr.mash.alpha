# mr.mash.alpha

 [![R-CMD-check](https://github.com/stephenslab/mr.mash.alpha/workflows/R-CMD-check/badge.svg)](https://github.com/stephenslab/mr.mash.alpha/actions)
 [![CircleCI build status](https://circleci.com/gh/stephenslab/mr.mash.alpha.svg?style=svg)](https://circleci.com/gh/stephenslab/mr.mash.alpha)

R package implementing methods for multivariate multiple regression
with adaptive shrinkage priors (mr.mash).

## Quick Start

Install and load the package using [remotes][remotes]:

```R
remotes::install_github("stephenslab/mr.mash.alpha")
library(mr.mash.alpha)
```

Note that installing the package will require a C++ compiler setup
that is appropriate for the version of R installed on your
computer. For details, refer to the documentation on the
[CRAN website][cran].

This command should automatically install all required packages if
they are not installed already.

## Citing this work

If you find the `mr.mash.alpha` package or any of the source code in this
repository useful for your work, please cite:
> Morgante, F., Carbonetto, P., Wang, G., Zou, Y., Sarkar, A. &
> Stephens, M. (2023). A flexible empirical Bayes approach to 
> multivariate multiple regression, and its improved accuracy 
> in predicting multi-tissue gene expression from genotypes.
> *PLoS Genetics* 19(7): e1010539. https://doi.org/10.1371/journal.pgen.1010539

## License

Copyright (c) 2020-2023, Fabio Morgante, Peter Carbonetto and Matthew
Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license].


[remotes]: https://github.com/r-lib/remotes
[cran]: https://cran.r-project.org
[mit-license]: https://opensource.org/licenses/mit-license.html

