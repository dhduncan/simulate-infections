#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param initial_infections
#' @param parameters
#' @return
#' @author Nick Golding
#' @export
sim_abm <- function(
  infections = sim_initial_infections(100),
  parameters = setup_abm(),
  max_days = 365,
  max_infections = 1e4
) {
  
  # put the parameters and some mutable variables in the global environment, to
  # be accessed anywhere
  .abm_parameters <<- parameters
  .abm_globals <<- list(highest_id = 0)
  .abm_tracker <<- tibble(days = integer(0),
                          ind_infectiousness = numeric(0))
  
  days <- seq_len(max_days)
  for (day in days) {
    
    .abm_globals$day <<- day
    .abm_globals$highest_id <<- nrow(infections)
    
    # infect people in the population
    infections <- infect(infections)
    
    
    # do passive detection (before contact tracing, so we can contact trace from them)   
    if (parameters$symptomatic_detections) {
      infections <- symptomatic_presentation(infections)
    }
    
        # do contact tracing for any people put into isolation today
    if (parameters$contact_tracing) {
      infections <- do_contact_tracing(infections)
    }
    
    # do workplace screening for asymptomatic people, so they can be isolated (and their contacts traced??) 
    if (parameters$workplace_screening) {
      infections <- do_screening(infections)
    }
    
    # quit if we hit the maximum total infections
    if (nrow(infections) >= max_infections) {
      break()
    }
    
  }
  
  infections
  
}
