##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/sparseEigen")
# Installation from CRAN
install.packages("sparseEigen")
# Getting help
library(sparseEigen)
help(package="sparseEigen")
?spEigen




##
## Developer commands (http://r-pkgs.had.co.nz/)
##
library(devtools)
#devtools::create("sparseEigen")
devtools::load_all()  #or Ctrl-Shift-L
#devtools::use_package("mvtnorm")
devtools::install()
devtools::build()  # to generate the installation file
#devtools::use_readme_rmd()  # to create the README file
#devtools::use_data_raw()  # to set up the raw-data folder

# Vignettes
#devtools::use_vignette("sparse_eigenvectors")  # to create the folder the first time
#rmarkdown::render("vignettes/sparse_eigenvectors.Rmd", "md_document")  # this is to generate the .md for GitHub
#rmarkdown::render("vignettes/sparse_eigenvectors.Rmd", "pdf_document")
rmarkdown::render("vignettes/sparse_eigenvectors.Rmd", "all")  # this also generates the pdf
#tools::compactPDF("vignettes/sparse_eigenvectors.pdf", gs_quality = "printer")  # this compresses the pdf
devtools::build_vignettes()  # or just install() will do

# Documentation
devtools::document()  #to generate all documentation via roxygen
?spEigen

# README (.md file has to be generated manually)


# code style
lintr::lint_package()


# Code tests
#devtools::use_testthat()  # the first time
devtools::test()
covr::package_coverage()  #coverage of tests
#goodpractice::gp()  # overall checks


# CRAN check
devtools::check()
rcmdcheck::rcmdcheck()
#R CMD check . --as-cran  # this is before submission to CRAN

# to upload to CRAN
devtools::build_win()
#devtools::release()  #for CRAN




