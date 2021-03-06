#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param infections
#' @return
#' @author Nick Golding
#' @export
symptomatic_presentation <- function(infections) {

  # renamed to make it clear that its about the symptomatic presentation in context of passive surveillance
  
  # only detectable if they are not in isolation and are symptomatic
  detectable <- infections$symptomatic &
    !is.finite(infections$isolation_day)
    
  # is the day of test seeking
  days_since_infection <- .abm_globals$day - infections$infection_day 
  
  # a multiplier on detection for the vaccinated
  vaccination_multiplier <- ifelse(
    infections$vaccinated,
    .abm_parameters$vaccination_test_seeking_multiplier,
    1
  )
  
  # they have a 50% probability of detection overall, or possibly lower if
  # vaccinated, scaled by the probability that this is the day of detection
  p_detection <- .abm_parameters$passive_detection_given_symptoms *
    time_to_symptomatic_test_pmf(days_since_infection) *
    vaccination_multiplier *
    detectable
  
  detected <- rbinom(nrow(infections), 1, p_detection)
  
  # put them in isolation on this day if detected
  infections$isolation_day <- ifelse(
    detected,
    .abm_globals$day,
    infections$isolation_day
  )
  
  infections$isolated <- ifelse(
    detected,
    TRUE,
    infections$isolated
  )
  
  
  infections$case_found_by <- ifelse(
    detected,
    "symptomatic_surveillance",
    infections$case_found_by
  )
  
  infections

}
