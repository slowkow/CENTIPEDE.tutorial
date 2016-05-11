#!/usr/bin/env Rscript
# make_rdata.R
# Kamil Slowikowski
# May 11, 2016
#
# TODO Write a description

library(CENTIPEDE.tutorial)

cen <- centipede_data(
  bam_file = "data-raw/ENCFF000SHS.bam",
  fimo_file = "data-raw/M6496_1.02.fimo.txt.gz",
  pvalue = 1e-4,
  flank_size = 100
)
save(list = c("cen"), file = "data/cen.rda")

