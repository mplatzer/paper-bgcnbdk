
#' Load Required R Packages

library(data.table)
library(BTYD)
library(ggplot2)
library(cowplot)
library(beanplot)

release <- '1.1.4'
if (!"BTYDplus" %in% installed.packages()[,"Package"] || packageVersion('BTYDplus')!=release) {
  devtools::install_github("mplatzer/BTYDplus", ref = release, dependencies=TRUE)
}

library(BTYDplus)
