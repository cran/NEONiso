% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_functions.R
\name{write_water_calibration_data}
\alias{write_water_calibration_data}
\title{write_water_calibration_data}
\usage{
write_water_calibration_data(outname, site, calDf)
}
\arguments{
\item{outname}{Output file name.}

\item{site}{NEON 4-letter site code.}

\item{calDf}{Calibration data frame -
this is the output from fit_water_regression}
}
\value{
Nothing to the environment, but writes out the
calibration parameters (e.g.,
regression slopes and intercepts) to the output
hdf5 file.
}
\description{
write_water_calibration_data
}
\author{
Rich Fiorella \email{rfiorella@lanl.gov}
}
