
#' Load Required R Packages

library(readr)
library(parallel)
library(data.table)
library(BTYD)
library(ggplot2)
library(cowplot)
library(beanplot)
library(fst)

release <- '1.1.4'
if (!"BTYDplus" %in% installed.packages()[,"Package"] || packageVersion('BTYDplus')!=release) {
  devtools::install_github("mplatzer/BTYDplus", ref = release, dependencies=TRUE)
}

library(BTYDplus)
