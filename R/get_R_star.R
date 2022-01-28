#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param day
#' @param amplitude
#' @param wavelength
#' @return
#' @author dhduncan
#' @export
get_R_star <- function(day = .abm_globals$day, 
                       amplitude = 0.1, 
                       wavelength = 0.09) {

  # amplitude and wavelength arbitrarily selected based on viewing curve(exp(amplitude * sin(x * wavelength)), from = 1, to = 365) and looking for R=1 threshold crossing behaviour
  
  exp(amplitude * sin(day * amplitude))

}
