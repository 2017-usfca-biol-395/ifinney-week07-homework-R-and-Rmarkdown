
## call rmarkdown on all .Rmd files
f <- list.files(recursive = TRUE)
rmds <- f[grepl(".Rmd$", f)]
lapply(rmds, rmarkdown::render)
lapply(rmds, lintr::lint)
