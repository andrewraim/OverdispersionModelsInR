# Overview
The OverdispersionModelsInR package was developed for demonstration purposes,
to accompany the technical report
[Modeling Overdispersion in R](/publications/#OverdispersionModelsInR2015).
The objective was to use R to reproduce analyses from the book
[Overdispersion Models in SAS](http://www.sas.com/store/prodBK_62693_en.html)
by Jorge Morel and Nagaraj Neerchal. We were especially interested in analyses involving likelihoods for binomial data with extra variation. For more
information, see the technical report or the ProbStatDay 2014 workshop presentation
_ R Supplement to "Analysis of Overdispersed Data using SAS"_
[[slides](http://www.umbc.edu/~araim1/pub/psday2014-workshop/slides.pdf) |
[handout](http://www.umbc.edu/~araim1/pub/psday2014-workshop/handout.pdf) |
[tech report](http://www.umbc.edu/~araim1/publications#OverdispersionModelsInR2015)].

The original release of the package is [OverdispersionModelsInR_0.1.tar.gz](https://github.com/andrewraim/OverdispersionModelsInR/archive/v0.1.tar.gz). Some improvements and corrections have been made
since then.

Note that the code is currently not documented. This means, for example, that
the R command "?d.rcb" asking for help on the d.rcb function currently returns
nothing. I hope to add this documentation in the future,
and continue to make improvements to the package. If you find the package to be
helpful, or have suggestions, please send me an email and let me know.

# Installation
The package can be installed directly from Github using
``` R
> library(devtools)
> install_github("andrewraim/OverdispersionModelsInR")
```

Or, download the package tarball and install it locally.
``` R
> install.packages("/path/to/OverdispersionModelsInR_0.1.tar.gz", repos = NULL)
```

You may need to install development tools as a prerequisite to build the
package, as some of it is written in C. See this
[link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for instructions specific to Windows, Mac, and Linux users.
