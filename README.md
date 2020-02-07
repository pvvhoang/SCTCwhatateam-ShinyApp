# SCTCwhatateam - Shiny Application

This is the Shiny application for the paper "The winning methods of predicting cellular position in the DREAM single cell transcriptomics challenge".

## Installation

To use this Shiny application, you can use RStudio and you will need to install R and Python packages.

We recommend anaconda3 if you need to install python, pip, and conda. You also need to install magic-impute.

In R, you may need to install packages: shiny, shinythemes, reticulate, Rmagic, SCTCwhatateam, and ggplot2. For the SCTCwhatateam, you can install from github https://github.com/thanhbuu04/SCTCwhatateam

## Quick Start

### Set up the environment

In file server.R, you may need to set up the variables below:

c_python <- ".../anaconda3/python.exe"
c_conda <- ".../anaconda3/Scripts/conda.exe"
c_linmethods <- ".../linmethods.py"
c_run_experiment <- ".../runExperiments.R"
c_shiny_run <- ".../Shiny_run" # Folder for output
c_shiny_run_dash <- ".../Shiny_run/shiny_run_" # Output file

In file runExperiments.R, please set the following variable:

c_linmethods_py <- ".../linmethods.py"

The sample data and user guide can be obtained from here: https://www.dropbox.com/sh/q7fq3y5yw1lmtrq/AAAnVDVAWzdBRS6FJor-Bk1Aa?dl=0
