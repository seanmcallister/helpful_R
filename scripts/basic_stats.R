#! /usr/bin/env Rscript
d<-scan("stdin", quiet=TRUE)

summary(d, digits = 300)
#sd(d[ , 1])
#cat("min", "max", "median", "mean", "\n", sep="\t")
#cat(min(d), max(d), median(d), mean(d), sep="\n")


