#######


### Install packages
install.packages("devtools")
install.packages("roxygen2")

### Load
library(devtools)
library(roxygen2)

###
load_all(".") # Working directory should be in the package SCC_R_package


###
roxygenise()  # Builds the help files

