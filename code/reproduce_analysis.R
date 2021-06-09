if( !require('renv') ) {
  install.packages('renv')
}
renv::restore()
source(rmarkdown::purl(here("analysis/index.Rmd"), output = tempfile()))
