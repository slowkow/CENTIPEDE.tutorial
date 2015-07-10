# CENTIPEDE Tutorial

This is a practical tutorial for running [CENTIPEDE] with [DNase-Seq] data. It
explains how to prepare the data and how to run the analysis. The goal is to
predict if a putative transcription factor binding site is actually bound or
not. For details about the statistical models underlying the methods, please
see ([Pique-Regi, et al. 2011][1]).

[Click here][2] to read the tutorial.

This repository has an R code to prepare data for CENTIPEDE. Install and load it:

```r
install.packages("devtools")
library(devtools)

devtools::install_github("slowkow/CENTIPEDE.tutorial")
library(CENTIPEDE.tutorial)
```

[CENTIPEDE]: http://centipede.uchicago.edu/ 
[DNase-Seq]: https://en.wikipedia.org/wiki/DNase-Seq
[1]: http://genome.cshlp.org/content/21/3/447
[2]: https://github.com/slowkow/CENTIPEDE.tutorial/raw/master/vignettes/centipede-tutorial.pdf
