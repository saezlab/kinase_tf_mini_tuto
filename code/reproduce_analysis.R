if( !require('renv') ) {
  install.packages('renv')
}
renv::restore()
source(knitr::purl("analysis/index.Rmd"), output = tempfile())