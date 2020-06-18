# Testing out qPCR package to analyze raw qPCR data using sigmoidal fitting methods
# Author: Prashant Kalvapalle, Stadler lab at Rice University
# Date: 18 Jun 2020

# User inputs ----

flnm <- 'WW4-BRSV_BCoV_Vaccines'  # set the filename

source('./general_functions.R') # Source the general_functions file before running this
library(qpcR)

# Read file ----
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
fl <- readqpcr(flpath) # read excel file exported by Quantstudio

plate_template_raw <- read_sheet('https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118', sheet = 'Plate import setup', range = 'G1:S9')
plate_template <- read_plate_to_column(plate_template_raw, 'Sample Name') %>% separate(`Sample Name`,c(NA, 'Sample Name'),'-')  # ignore target name within the sample name column

# reads named CT data from 
