#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
setup_abm <- function(
  
  # pre-vaccination R to aim for
  R = 3.62,
  
  # vaccination effects on transmission    
  vaccination_coverage = 0.94,
  # incorporate correction for onward transmission to account for reduction in
  # infectiousness due to being symptomatic
  ve_onward = 0.339, # sensitivity analysis values, lower=0.575, upper=0.7029
  ve_susceptibility = 0.435,
  ve_symptoms = 0.375,
  
  # symptomaticity and passive detection
  clinical_fraction = 0.307, #0.307, # was previously - 0.8
  passive_detection_given_symptoms = 0.5,
  asymptomatic_relative_infectiousness = 0.5,
  vaccination_test_seeking_multiplier = 1,
  
  # whether to do routine screening at workplaces
  workplace_screening = TRUE,
  
  screenable_fraction =  0.4 *  107e+05 / 23402e+03, # rounded population stats 
  
  # the first bit is a an arbitrary middle value from: https://www.sgsep.com.au/publications/insights/closing-the-divide-essential-workers-australian-cities-and-covid-19 and the latter is the workforce fraction of the australian population from https://profile.id.com.au/australia/population. (62% of that workforce is in FTE, but that detail not in sim at present).
  
  # switch for passive presentation of symptomatic individuals for testing
  symptomatic_detections = TRUE,
  
   static_R_star = TRUE,

  # probability of an infectee being found by contact tracing from the source
  p_active_detection = 0.5,
  
  # relative probability of active detection for vaccinated individuals
  rel_active_detection_vaccinated_source = 1,
  
  # relative probability of quarantinng/isolating vaccinated contact
  rel_active_detection_vaccinated_contact = 1,
  
  # whether to do downstream contact tracing
  contact_tracing = TRUE,
  
  # placeholder for delays
  isolation_to_interview_samples = get_optimal_isol_interview_samples(),
  isolation_days_vax=14,
  isolation_start_day=c("isolation", "infection")
) {
  
  isolation_start_day <- match.arg(isolation_start_day)
  
  args <- list(
    R = R,
    vaccination_coverage = vaccination_coverage,
    ve_onward = ve_onward,
    ve_susceptibility = ve_susceptibility,
    ve_symptoms = ve_symptoms,
    clinical_fraction = clinical_fraction,
    passive_detection_given_symptoms = passive_detection_given_symptoms,
    asymptomatic_relative_infectiousness = asymptomatic_relative_infectiousness,
    vaccination_test_seeking_multiplier = vaccination_test_seeking_multiplier,
    workplace_screening = workplace_screening,
    screenable_fraction = screenable_fraction,
    symptomatic_detections = symptomatic_detections,
    p_active_detection = p_active_detection,
    rel_active_detection_vaccinated_source = rel_active_detection_vaccinated_source,
    rel_active_detection_vaccinated_contact = rel_active_detection_vaccinated_contact,
    contact_tracing = contact_tracing,
    isolation_to_interview_samples = isolation_to_interview_samples,
    isolation_days_vax=isolation_days_vax,
    isolation_start_day=isolation_start_day,
     static_R_star = static_R_star

    
  )
  
  # correct R for reduced infectiousness of asymptomatics
  args$R_star <- with(
    args,
    R / ((1 - clinical_fraction) + clinical_fraction * asymptomatic_relative_infectiousness)
  )
  
  args
  
}
