#install.packages("devtools")
library(devtools)

pkg=file.path("~/git", "aneuploidy_score")
setwd(pkg)


#use_description()

#### Assembling data ####
usethis::use_data_raw()
devtools::load_all()

#### Building ####
#Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/MacOS/pandoc")

devtools::document(pkg)
devtools::check(pkg)
devtools::build_vignettes(pkg)

devtools::build(pkg)
devtools::install(pkg)
#devtools::install_github("quevedor2/aneuploidy_score")


