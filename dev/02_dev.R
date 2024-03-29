# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
usethis::use_package("plotly")
usethis::use_package("shinyWidgets")
usethis::use_package("shinydashboard")
usethis::use_package("Seurat")
usethis::use_package("shinybusy")
usethis::use_package("shinyjs")
usethis::use_package("dashboardthemes")
usethis::use_package("shinyBS")
usethis::use_package("Matrix")
usethis::use_package("scales")
usethis::use_package("RColorBrewer")
usethis::use_package("velocyto.R")
usethis::use_package("ggVennDiagram")
usethis::use_package("pheatmap")
usethis::use_package("DT")
usethis::use_package("ggplot2")

## Add modules ----
## Create a module infrastructure in R/
## Présentation
golem::add_module(name = "Presentation") 

## Patient
golem::add_module(name = "Patient_Expression") 
golem::add_module(name = "Patient_Metadata")
golem::add_module(name = "Patient_Mitochondrie")
golem::add_module(name = "Patient_Reduction_Dimension")
golem::add_module(name = "Patient_Subset")
golem::add_module(name = "Patient_VDJ")
golem::add_module(name = "Patient_Velocity")
golem::add_module(name = "Patient_GOT")




## All
golem::add_module(name = "All_Cell")
golem::add_module(name = "All_DE")
golem::add_module(name = "All_Enrichissement")
golem::add_module(name = "All_Reduction_Dimension")
golem::add_module(name = "All_scRepertoir")
golem::add_module(name = "All_VDJ")
golem::add_module(name = "All_Velocity")

## Other
golem::add_module(name = "Contact")




## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct("helpers", with_test = TRUE)
golem::add_utils("helpers", with_test = TRUE)

## External resources
## Creates .js and .css files at inst/app/www
golem::add_js_file("script")
golem::add_js_handler("handlers")
golem::add_css_file("custom")
golem::add_sass_file("custom")

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw(name = "my_dataset", open = FALSE)

## Tests ----
## Add one line by test you want to create
usethis::use_test("app")

# Documentation

## Vignette ----
usethis::use_vignette("scShiny")
devtools::build_vignettes()

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
##
## (You'll need GitHub there)
usethis::use_github()

# GitHub Actions
usethis::use_github_action()
# Chose one of the three
# See https://usethis.r-lib.org/reference/use_github_action.html
usethis::use_github_action_check_release()
usethis::use_github_action_check_standard()
usethis::use_github_action_check_full()
# Add action for PR
usethis::use_github_action_pr_commands()

# Travis CI
usethis::use_travis()
usethis::use_travis_badge()

# AppVeyor
usethis::use_appveyor()
usethis::use_appveyor_badge()

# Circle CI
usethis::use_circleci()
usethis::use_circleci_badge()

# Jenkins
usethis::use_jenkins()

# GitLab CI
usethis::use_gitlab_ci()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")
