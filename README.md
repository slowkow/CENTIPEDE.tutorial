# CENTIPEDE Tutorial

> CENTIPEDE fits a bayesian hierarchical mixture model to learn TF-specific
> distribution of experimental data on a particular cell-type for a set of
> candidate binding sites described by a motif.

This is a practical tutorial for running [CENTIPEDE] with [DNase-Seq] data. It
explains how to prepare the data and how to run the analysis. The goal is to
predict if a putative transcription factor binding site is actually bound or
not. For details about the statistical models underlying the methods, please
see ([Pique-Regi, et al. 2011][1]).

[CENTIPEDE]: http://centipede.uchicago.edu/ 
[DNase-Seq]: https://en.wikipedia.org/wiki/DNase-Seq
[1]: http://genome.cshlp.org/content/21/3/447

Read the tutorial online or download the PDF:

- [http://slowkow.github.io/CENTIPEDE.tutorial](http://slowkow.github.io/CENTIPEDE.tutorial/)
- [centipede-tutorial.pdf][2]

[2]: https://github.com/slowkow/CENTIPEDE.tutorial/raw/master/vignettes/centipede-tutorial.pdf

This repository has functions to ease the use of CENTIPEDE:

- `centipede_data()` converts data to the format required for CENTIPEDE.
- `parse_region()` parses a string like "chr1:123-456".
- `read_bedGraph()` reads a bedGraph file with 4 columns: chrom, start, end, score.
- `read_fimo()` reads a text file output by FIMO and selects sites that meet a
  significance threshold.

I also provide example data that you can use to follow the tutorial:

- `cen` is a list with two items:
  - `cen$mat` is a matrix of read-start counts for 3,337 genomic regions.
  - `cen$regions` is a dataframe describing those regions.
- `site_cons` is a vector with mean conservation scores for the 3,337 regions,
  computed across 100 vertebrates.

# Installation

Install CENTIPEDE by running this in your shell (not within an R session):

```bash
wget http://download.r-forge.r-project.org/src/contrib/CENTIPEDE_1.2.tar.gz
R CMD INSTALL CENTIPEDE_1.2.tar.gz
```

Next, install the tutorial package:

```r
# This command didn't work for me.
# install.packages("CENTIPEDE", repos="http://R-Forge.R-project.org")

install.packages("devtools")
devtools::install_github("slowkow/CENTIPEDE.tutorial")
```
