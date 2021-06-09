if( !require('renv') ) {
  install.packages('renv')
}
if( !require('workflowr') ) {
  install.packages('workflowr')
}
renv::restore()
workflowr::wflow_build('analysis/index.Rmd')
