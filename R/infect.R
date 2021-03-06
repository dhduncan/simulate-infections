#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param infections
#' @param day
#' @return
#' @author Nick Golding
#' @export
infect <- function(infections) {
  
  # add on new infections for this timestep, infecting the vaccinated and
  # unvaccinated populations separately
  rbind(
    infections,
    new_infections(infections, vaccinated = TRUE),
    new_infections(infections, vaccinated = FALSE)
  )
  
}
