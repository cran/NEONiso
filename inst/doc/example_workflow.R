## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(NEONiso)

## ---- eval = FALSE------------------------------------------------------------
#  manage_local_EC_archive(file_dir = "~/Desktop",
#                          get = TRUE,
#                          unzip_files = TRUE,
#                          sites = "ONAQ")

## ---- eval = FALSE------------------------------------------------------------
#  for (i in 1:length(fnames.out)) {
#    calibrate_carbon_bymonth(fnames[i],
#                             fnames.out[i],
#                             site=site.code[i],
#                             method = "Bowling_2003")
#  }

## ---- eval = FALSE------------------------------------------------------------
#  data.dir <- "/your/path/here/DP4_00200_001/"
#  
#  fnames <- list.files(path = data.dir,
#                       pattern = ".h5",
#                       recursive = TRUE,
#                       full.names = TRUE)
#  
#  # unselect gz files.
#  fnames <- fnames[!grepl(".gz", fnames)]
#  
#  fname.byfolder <- strsplit(fnames, split=".", fixed = TRUE)
#  site.code  <- sapply(fname.byfolder, "[[", 3)
#  
#  # inspect site.code in the environment: is it a vector with repeated "ONAQ"?
#  fnames.tmp <- gsub(".h5", ".calibrated.h5", fnames)
#  fnames.spt <- strsplit(fnames.tmp, split = "/")
#  fnames.out <- sapply(fnames.spt, "[[", 7)
#  
#  # create new output directory
#  outpaths   <- paste0(your_path, "/ONAQ/output/")
#  # apply function used here to generalize in case you wanted to run all sites
#  sapply(unique(outpaths), dir.create, showWarnings=FALSE)
#  
#  # update fnames.out to include desired output paths.
#  fnames.out <- paste0(outpaths, "/", fnames.out)

