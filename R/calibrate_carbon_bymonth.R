#' calibrate_carbon_bymonth
#'
#' `r lifecycle::badge("deprecated")`
#' This function drives a workflow that reads in NEON carbon isotope data
#' of atmospheric CO2, calibrates it to the VPDB scale, and (optionally)
#' writes the calibrated data to a new HDF5 file. Two different approaches
#' are possible: a) a calibration on 12CO2 and 13CO2 isotopologues separately,
#' after Bowling et al. 2003 (Agr. For. Met.), or b) a direct calibration
#' of d13C and CO2 values using linear regression. The vast majority of the time
#' the results generated from either method are extremely similar to each other.
#' Wen et al. 2013 compared several different carbon
#' isotope calibration techniques and found this to be the superior method
#' under most circumstances. We also found this to be the case for NEON data
#' (Fiorella et al. 2021; JGR-Biogeosciences).
#'
#' The 'linreg' method simply takes measured and reference d13C and CO2 values
#' and generates a transfer function between them using `lm()`. For the
#' gain-and-offset method, d13C and CO2 values are converted to 12CO2 and 13CO2
#' mole fractions. Gain and offset parameters are calculated for each
#' isotopologue independently, and are analogous to regression slope
#' and intercepts, but jointly correct for CO2 concentration dependence
#' and place d13C values on the VPDB scale.
#' The gain and offset parameters are defined by:
#'
#' \deqn{G = (X_{2,ref}-X_{1,ref})/(X_{2,meas}-X_{1,meas})}
#' \deqn{O = X_{2,ref}- G X_{2,meas}}
#' Calibrated ambient isotopologues are then given as:
#' \deqn{X_{cal} = X_{meas} G + O}
#'
#' Measurements of reference materials were considered "good" if the following
#' conditions were met:
#' \itemize{
#'   \item Measured CO2 concentrations were within 10 ppm
#'         of known "reference" concentrations.
#'   \item Variance of the CO2 concentration in standard peak was < 5 ppm.
#'   \item Measured d13C value must be within 5 per mil
#'         of known "reference" d13C value.
#' }
#' The first two criteria are intended to filter out periods where there is
#' a clear issue with the gas delivery system (i.e., nearly empty gas tank,
#' problem with a valve in the manifold, etc.); the third criterion was adopted
#' after visual inspection of data timeseries revealed that often the first
#' standard measurement following an instrument issue had higher-than-expected
#' error. This criterion clips clearly poor values. Selection of these criteria
#' will become a function argument, and therefore customizable,
#' in a future release.
#'
#' @author Rich Fiorella \email{rfiorella@@lanl.gov}
#'
#' @param inname Name of the input file. (character)
#' @param outname Name of the output file. (character)
#' @param force_cal_to_beginning Extend first calibration to the beginning
#'             of the file? (default true)
#' @param force_cal_to_end Extend last calibration to the end of the file?
#'                         (default true)
#' @param site Four letter NEON site code for site being processed. (character)
#' @param gap_fill_parameters Should function attempt to 'gap-fill' across a
#'            bad calibration by carrying the last known good calibration
#'            forward?
#'            Implementation is fairly primitive currently, as it only carries
#'            the last known good calibration that's available forward rather
#'            than interpolating, etc. Default FALSE.
#' @param filter_ambient Apply the median absolute deviation filter (Brock 86)
#'            to remove impulse spikes in output ambient data?
#'            (logical; default true)
#' @param r2_thres Minimum r2 threshold of an "acceptable" calibration. Acts to
#'            remove calibration periods where a measurement error makes
#'            relationship nonlinear. Default = 0.95
#' @param correct_refData NEON has indicated there are a few instances where
#'            reported d13C or CO2 reference values are wrong. If set to true,
#'            correct known incorrect values. This argument will (hopefully,
#'            eventually) go away after NEON has fixed the reference database.
#'            Users will be warned prior to removal of this argument.
#' @param write_to_file Write calibrated ambient data to file?
#'              (Mostly used for testing)
#' @param method Are we using the Bowling et al. 2003 method
#'              ("Bowling_2003") or direct linear regression of
#'              d13C and CO2 mole fractions ("linreg")?
#' @param calibration_half_width Determines the period (in days)
#'        from which reference data are selected (period
#'        is 2*calibration_half_width).
#'
#' @return Returns nothing to the environment, but creates a new output HDF5
#'         file containing calibrated carbon isotope values.
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{fin <- system.file('extdata',
#' 'NEON.D15.ONAQ.DP4.00200.001.nsae.2019-05.basic.20201020T211037Z.packed.h5',
#'          package = 'NEONiso', mustWork = TRUE)
#' calibrate_carbon_bymonth(inname = fin, outname = 'out.h5',
#'          site = 'ONAQ', write_to_file = FALSE)
#' calibrate_carbon_bymonth(inname = fin, outname = 'out.h5',
#'          site = 'ONAQ', method = 'linreg', write_to_file = FALSE)}
#'
calibrate_carbon_bymonth <- function(inname,
                                     outname,
                                     site,
                                     method = "Bowling_2003",
                                     calibration_half_width = 0.5,
                                     force_cal_to_beginning = TRUE,
                                     force_cal_to_end = TRUE,
                                     gap_fill_parameters = FALSE,
                                     filter_ambient = TRUE,
                                     r2_thres = 0.95,
                                     correct_refData = TRUE,
                                     write_to_file = TRUE) {

  lifecycle::deprecate_warn("0.6.0",
                            "calibrate_carbon_bymonth()",
                            "calibrate_carbon()")

  #-----------------------------------------------------------
  # Extract reference data from input HDF5 file.
  #-----------------------------------------------------------
  # pull all carbon isotope data into a list.
  ciso <- ingest_data(inname, analyte = "Co2", name_fix = FALSE, avg = 9)

  # extract the data we need from ciso list
  refe <- extract_carbon_calibration_data(ciso$refe_stacked)

  # Okay this function now needs some work. *************
  if (correct_refData == TRUE) {

    # do some work to correct the reference data frame
    refe <- correct_carbon_ref_cval(refe,site)

  }

  # get calibration parameters data.frame.
  cal_df <- fit_carbon_regression(ref_data = refe, method = method,
                                  calibration_half_width = calibration_half_width)

#----------------------------------------------------------------------------
#  calibrate ambient data.
#  extract ambient measurements from ciso
  ciso <- rhdf5::h5read(inname, paste0("/", site, "/dp01/data/isoCo2"))
  ciso_logical <- grepl(pattern = "000",
                        x = names(ciso)) & grepl(pattern = "09m",
                                                 x = names(ciso))
  ciso_subset <- ciso[ciso_logical]

  if (method == "Bowling_2003") {
    ciso_subset_cal <- lapply(names(ciso_subset),
                              function(x) {
                                calibrate_ambient_carbon_Bowling2003(
                                  amb_data_list = ciso_subset[[x]],
                                  caldf = cal_df,
                                  site = site,
                                  filter_data = filter_ambient,
                                  force_to_end = force_cal_to_end,
                                  force_to_beginning = force_cal_to_beginning,
                                  r2_thres = r2_thres)
                              })
  } else if (method == "linreg") {
    ciso_subset_cal <- lapply(names(ciso_subset),
                              function(x) {
                                calibrate_ambient_carbon_linreg(
                                  amb_data_list = ciso_subset[[x]],
                                  caldf = cal_df,
                                  site = site,
                                  filter_data = filter_ambient,
                                  force_to_end = force_cal_to_end,
                                  force_to_beginning = force_cal_to_beginning,
                                  r2_thres = r2_thres)
                              })
  }

  names(ciso_subset_cal) <- names(ciso_subset)

  #-----------------------------------------------------------
  # write out these data.frames to a new output file.
  #-----------------------------------------------------------
  if (write_to_file) {
    cal_df$timeBgn <- convert_POSIXct_to_NEONhdf5_time(cal_df$timeBgn)
    cal_df$timeEnd <- convert_POSIXct_to_NEONhdf5_time(cal_df$timeEnd)
    setup_output_file(inname, outname, site, "co2")
    write_carbon_calibration_data(outname, site, cal_df, method = method)
    write_carbon_ambient_data(outname, site, ciso_subset_cal)
    write_carbon_reference_data(inname, outname, site, cal_df)
    write_qfqm(inname, outname, site, "co2")
    write_ucrt(inname, outname, site, "co2")

    validate_output_file(inname, outname, site, "co2")

    # one last invocation of hdf5 close all, for good luck
    rhdf5::h5closeAll()
  }

}
