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
get_R_star <- function(day = .abm_globals$day) {

height = log(.abm_parameters$R) # just trying get in the ball park for now
amplitude = 0.2
wavelength = 0.03
# amplitude and wavelength arbitrarily selected based on viewing curve(exp(amplitude * sin(x * wavelength)), from = 1, to = 365) and looking for R=1 threshold crossing behaviour
  
  exp(height + amplitude * sin(day * wavelength))

}
