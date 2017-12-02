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
devtools::document()  #to generate all documentation via roxygen
devtools::install()
devtools::build()  # to generate the installation file
#devtools::use_readme_rmd()  # to create the README file
#devtools::use_data_raw()  # to set up the raw-data folder

#for vignettes
#devtools::use_vignette("sparse_eigenvectors")  # to create the folder the first time
#rmarkdown::render("vignettes/sparse_eigenvectors.Rmd", "md_document")  # this is to generate the .md for GitHub
rmarkdown::render("vignettes/sparse_eigenvectors.Rmd", "all")  # this also generates the pdf
tools::compactPDF("vignettes/sparse_eigenvectors.pdf", gs_quality = "printer")  # this compresses the pdf
devtools::build_vignettes() #or just install()


# code checking
lintr::lint_package()
devtools::check()

# code tests
#devtools::use_testthat()
devtools::test()
covr::package_coverage()  #coverage of tests


# overall checks:
#goodpractice::gp()

# CRAN check
rcmdcheck::rcmdcheck()
#R CMD check --as-cran  # this is before submission to CRAN
#R CMD check --as-cran --compact-vignettes

# to upload to CRAN
#devtools::release()  #for CRAN
