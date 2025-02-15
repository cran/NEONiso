% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isotope_conversions.R
\name{delta_to_R}
\alias{delta_to_R}
\title{Converts delta value to heavy-to-light isotope ratio}
\usage{
delta_to_R(delta_values, element)
}
\arguments{
\item{delta_values}{A vector of isotope ratios in delta notation.}

\item{element}{Which element to return R values -
carbon, oxygen, or hydrogen.}
}
\value{
Vector of isotope ratios (R values).
}
\description{
Converts a delta value (in per mil) to the heavy-to-light
isotope ratio.
}
\examples{
delta_to_R(delta_values = 0, element = 'oxygen') # 2005.2e-6 for VSMOW.

}
\author{
Rich Fiorella \email{rfiorella@lanl.gov}
}
