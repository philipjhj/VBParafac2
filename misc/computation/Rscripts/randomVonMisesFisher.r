library("rstiefel", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")

args <- commandArgs(trailingOnly = TRUE)

rmf.matrix(matrix(rnorm(10*5),10,5))