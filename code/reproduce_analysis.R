if( !require('renv') ) {
  install.packages('renv')
}
renv::restore()
source(knitr::purl(here("analysis/index.Rmd"), output = tempfile()))
