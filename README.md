# Overview
The OverdispersionModelsInR package was developed for demonstration purposes,
to accompany the technical report
[Modeling Overdispersion in R](/publications/#OverdispersionModelsInR2015).
The objective is to develop the capability in R to reproduce analyses from the
book
[Overdispersion Models in SAS](http://www.sas.com/store/prodBK_62693_en.html)
by Jorge Morel and Nagaraj Neerchal,
especially those analyses requiring likelihoods
for binomial data with extra variation. For more information, see the technical
report or the ProbStatDay 2014 workshop presentation
<em>R Supplement to "Analysis of Overdispersed Data using SAS"</em>
[[slides](http://www.umbc.edu/~araim1/pub/psday2014-workshop/slides.pdf) |
[handout](http://www.umbc.edu/~araim1/pub/psday2014-workshop/handout.pdf)].

Note that the code is currently not documented. This means, for example, that
the R command "?d.rcb" asking for help on the d.rcb function currently returns
nothing. I hope to add this documentation in the future,
and continue to make improvements to the package. If you find the package to be
helpful, or have suggestions, please send me an email and let me know. If there
is interest, I may also put the source code into a more collaborative place
so that others can contribute.

# Download
Latest version of the package.
* Source:
[OverdispersionModelsInR_0.1.tar.gz](http://www.umbc.edu/~araim1/pub/OverdispersionModelsInR/OverdispersionModelsInR_0.1.tar.gz)
* Linux x86_64 binary: (TBD)
* Windows binary: (TBD)
* Mac binary: (TBD)

## Installing Binary Packages
Most users can install the binary version specific to their
computing platform. This does not require special development tools, which are
needed to compile the source version of the package.

To install the binary version, save the appropriate file to your computer
(Suppose it is called OverdispersionModelsInR-XYZ) and run the following
command in R.
``` R
> install.packages("/path/to/file/OverdispersionModelsInR-XYZ")
```


## Installing the Source Package
It may be necessary to install the source package if the prepared binary versions
don't work on your computer, or if you make modifications to the one I have
provided.

1. Install [R](http://www.r-project.org) if you haven't already.

2. Install the development tools needed to build R packages. This includes C/C++
compilers and LaTeX. See this [link](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
for instructions specific to Windows, Mac, and Linux users.

3. In R, install the devtools package, if not already installed.
 ``` R
 > install.packages("devtools")
 ```

4. In R, load the devtools package and ensure Rtools can be located
 ``` R
 > library(devtools)
 > find_rtools()
 ```

5. Finally, in R, install the OverdispersionModelsInR package using devtools
 ``` R
 > devtools::install_url("http://www.umbc.edu/~araim1/pub/OverdispersionModelsInR/OverdispersionModelsInR_0.1.tar.gz")
 ```
