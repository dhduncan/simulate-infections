# Freya Shearer's file, embedding the bAByM 

## Loads all packages and defines how to handle NAMESPACE conflicts

# source("./packages.R")
# CTRL + SHFT + ALT + P

set.seed(2)

future::plan(multisession(workers = 8))
#future::plan(sequential, split=TRUE)
sims <- expand_grid(
  vaccination_coverage = 0.9,
  vaccination_test_seeking_multiplier = 1,
  passive_detection_given_symptoms = 0.85,
  rel_active_detection_vaccinated_source = 1,
  rel_active_detection_vaccinated_contact = 1,
  p_active_detection =  0.85, #c(0.5, 0.75, 0.95),
  ve_onward=0.2,
  isolation_days_vax=c(7), 
  isolation_start_day=c('isolation'),
  symptomatic_detections=TRUE, 
  contact_tracing=TRUE,
  workplace_screening=TRUE,
  static_R_star = FALSE
  ) %>%
  mutate(
    # tweak starting R to get optimal reproduction number at about 1
    R = case_when(
      vaccination_coverage <= .79 ~ 5.8, # 4.6 fails < 1.55
      vaccination_coverage >= 0.9 ~ 5, # 4.6 fails <2.58
      TRUE ~ 4.6 # 4.6 fails < 1.9 | assuming this is for vax levels between 8 and 9
    )
  ) %>%
  rowwise() %>%
  mutate(
    parameters = list(
      setup_abm(
        R = R,
        vaccination_coverage = vaccination_coverage,
        vaccination_test_seeking_multiplier = vaccination_test_seeking_multiplier,
        passive_detection_given_symptoms = passive_detection_given_symptoms,
        rel_active_detection_vaccinated_source = rel_active_detection_vaccinated_source,
        rel_active_detection_vaccinated_contact = rel_active_detection_vaccinated_contact,
        ve_onward = ve_onward,
        isolation_days_vax = isolation_days_vax,
        isolation_start_day = isolation_start_day,
        symptomatic_detections = symptomatic_detections,
        contact_tracing = contact_tracing,
        workplace_screening = workplace_screening,
        static_R_star = static_R_star
      )
    )
  ) %>%
  mutate(
    simulations = list(
      get_valid_abm_samples(parameters, n_samples = 1) # we can only use one sample for now
    )
  )  

xlim <- sims %>% unnest(simulations) %>% summarise(max(infection_day)) %>% pull()

# plot simulations to check that they make sense
plot_df <- sims %>%
  unnest(simulations) %>% 
  group_by(simulation, p_active_detection,
           infection_day,
  ) %>% 
  summarise(infections=n()
  ) #%>% 
#filter(simulation %in% paste0('sim_', 1:50))

plot_infections <- ggplot(plot_df,
       aes(
         x=infection_day, 
         y=infections,
         group=simulation)
) +
  scale_x_continuous(limits=c(0,xlim)) +
  # facet_grid(.~p_active_detection) +
  theme_cowplot() +
  geom_line(aes(color=simulation)) +
  theme(legend.position = "none") 

plot_ascertainment <- sims %>% 
  unnest(simulations) %>% 
  filter(
    !is.na(source_id) &  # may be superfluous at present
      infection_day > 1 & # this isn't already done in get_valid_sim_sample?
  infection_day < (max(infection_day) - 14)
  ) %>% 
  group_by(simulation,
           # p_active_detection, 
           infection_day,
           case_found_by) %>% 
  summarise(infections = n(),
            ) %>% 
  mutate(proportion = infections / sum(infections)) %>% 
  filter(is.na(case_found_by)) %>% 
  summarise(ascertainment = 1 - proportion) %>% 
  ggplot(
    aes(
      x=infection_day, 
      y=ascertainment,
      group=simulation)) +
  scale_x_continuous(limits=c(0,xlim)) +
  # facet_grid(.~p_active_detection) +
  theme_cowplot() +
  geom_line() +
  theme(legend.position = "none") 

(sim_plot <- plot_grid(plot_infections, plot_ascertainment, nrow = 2))


# check for superspreading ----
# onward_infections <- sims %>% 
#   unnest(simulations) %>% 
#   filter(!is.na(source_id)) %>% 
#   group_by(simulation, source_id) %>% 
#   summarise(onward_infections = n(),
#             infection_day = first(infection_day))
# 
# onward_infection_hist <- ggplot(
#   onward_infections, 
#   aes(onward_infections)) + 
#   geom_histogram(binwidth = 1)
# 
# onward_thru_time <- sims %>% 
#   unnest(simulations) %>% 
#   left_join(
#     select(
#       ungroup(onward_infections), 
#       onward_infections, source_id), 
#     by = "source_id") %>% 
#   ggplot(
#     aes(x = infection_day,
#         y = onward_infections)
#   ) + 
#   geom_point(aes(alpha = onward_infections), 
#              show.legend = FALSE) +
#   facet_grid(case_found_by~.) +
#   theme_cowplot()
#   
# 
# ggplot(source_check_df,
#        aes(x = infection_day, y = onward_infections)) +
#   geom_point()
  
  
  
# ggsave(
#   'outputs/plots/infections_R5.8_pdetect0.5.png',
#   height=6,
#   width=8,
#   bg='white'
# )

# clean up simulations
trim_sims <- sims %>% 
  unnest(simulations) %>% 
  filter(
  # find sources to consider (exclude those during burn in and last two weeks
  # due to truncation of onward infection)
  !is.na(source_id) &  # may be superfluous at present
    infection_day > 10 & # 70% 20-50, 80% 50-80, 90% 10-300
    infection_day < (max(infection_day) - 14)) %>% 
  group_by(simulation,
           infection_day,
           case_found_by) %>% 
  summarise(infections = n()) %>% 
  group_by(simulation, infection_day) %>%
            mutate(proportion = infections / sum(infections),
                   infections = infections,
                   case_found_by=case_when(
                     is.na(case_found_by)~'undetected',
                     TRUE ~ case_found_by)
                   )


# could try and include the essential params in the file name rather than date which the system takes care of  

# 2)
# check trends in mode of detection of cases over time
# plot the number of cases 
# detected via symptomatic surveillance
# detected via active detection
# neither
trim_sims %>% 
  mutate(
    case_found_by = factor(case_found_by, 
                          levels = c("contact_tracing", 
                          "symptomatic_surveillance",
                          "workplace_screening", 
                          "undetected"),
                          labels = c("contact\ntracing", 
                                     "symptomatic\nsurveillance",
                                     "workplace\nscreening", 
                                     "undetected"),
                          ordered = TRUE)
  ) %>% 
ggplot(
  aes(
    x=infection_day, 
    y=infections,
    group=simulation)) + 
   facet_grid(
     .~case_found_by
   ) +
  theme_cowplot() +
  geom_line(aes(colour = simulation)) +
  theme(legend.position = "none") +
  labs(title="")
    
# ggsave(
#   'outputs/plots/infections_trajectories_by_mode_pdetect0.5.png',
#   height=6,
#   width=12,
#   bg='white'
# )

# save outputs for NG ascertainment model ----
ascertainment_data <- trim_sims %>% 
  #  filter(case_found_by != 'undetected') %>% 
  select(-proportion) %>% 
  mutate(
    case_found_by = case_when(
      case_found_by == "contact_tracing" ~ "contact",
      case_found_by == "workplace_screening" ~ "screening",
      case_found_by == "symptomatic_surveillance" ~ "symptomatic",
      TRUE ~ "undetected"
    )
  ) %>% 
  pivot_wider(names_from = case_found_by,
              values_from = infections,
              names_prefix = "count_tested_") #%>% 
#mutate(across(ends_with("_fraction"), ~if_else(is.na(.), 0, .)))

write_csv(ascertainment_data, "outputs/sim_output.csv")

saveRDS(sims$parameters[[1]], file = "outputs/abm_params.rds") 




sry_sims <- trim_sims %>% 
  group_by(case_found_by,
        infection_day) %>% 
  summarise(mean=mean(infections)) %>% 
  group_by(infection_day)%>% 
  mutate(sum=sum(mean)) %>% 
  mutate(fraction=mean/sum)

ggplot(sry_sims) + 
  geom_col(aes
           (x = infection_day,
             y = fraction, 
             fill = case_found_by
           )
           )+
  scale_fill_manual(values = c('#b3cde3',
                               '#ccebc5',
                               '#decbe4',
                               "grey")) +
  theme_cowplot() +
  ylab('Fraction') +
  xlab("Infection day") 

# ggsave(
#   'outputs/plots/fraction_infections_by_mode_pdetect0.5.png',
#   height=6,
#   width=12,
#   bg='white'
# )
