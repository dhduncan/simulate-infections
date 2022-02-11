# Freya Shearer's file, embedding the bAByM 

## Loads all packages and defines how to handle NAMESPACE conflicts
source("./packages.R")

future::plan(multisession(workers = 8))
#future::plan(sequential, split=TRUE)
sims <- expand_grid(
  vaccination_coverage = 0.74,
  vaccination_test_seeking_multiplier = 1,
  passive_detection_given_symptoms = c(0.5),
  rel_active_detection_vaccinated_source = c(1,0),
  rel_active_detection_vaccinated_contact = c(1,0),
  #ve_onward=0.639*c(0.9, 1, 1.1),
  ve_onward=0.639,
  isolation_days_vax=c(14), 
  isolation_start_day=c('isolation'),
  symptomatic_detections=TRUE, 
  contact_tracing=TRUE,
  workplace_screening=TRUE#,
#  static_R_star = TRUE
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
        workplace_screening = workplace_screening#,
        # static_R_star = static_R_star
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
  group_by(simulation, 
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
  theme_cowplot() +
  geom_line(aes(color=simulation)) +
  theme(legend.position = "none") 

plot_ascertainment <- sims %>% 
  unnest(simulations) %>% 
  filter(
    !is.na(source_id) &  # may be superfluous at present
      infection_day > 1 & 
  infection_day < (max(infection_day) - 14)) %>% 
  group_by(simulation,
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
  theme_cowplot() +
  geom_line() +
  theme(legend.position = "none") 

(sim_plot <- plot_infections + plot_ascertainment)


# check for superspreading ----
onward_infections <- sims %>% 
  unnest(simulations) %>% 
  filter(!is.na(source_id)) %>% 
  group_by(simulation, source_id) %>% 
  summarise(onward_infections = n(),
            infection_day = first(infection_day))

onward_infection_hist <- ggplot(
  onward_infections, 
  aes(onward_infections)) + 
  geom_histogram(binwidth = 1)

onward_thru_time <- sims %>% 
  unnest(simulations) %>% 
  left_join(
    select(
      ungroup(onward_infections), 
      onward_infections, source_id), 
    by = "source_id") %>% 
  ggplot(
    aes(x = infection_day,
        y = onward_infections)
  ) + 
  geom_point(aes(alpha = onward_infections), 
             show.legend = FALSE) +
  facet_grid(case_found_by~.) +
  theme_cowplot()
  

ggplot(source_check_df,
       aes(x = infection_day, y = onward_infections)) +
  geom_point()
  
  
  
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

# save output for NG ascertainment model
ascertainment_output <- trim_sims %>% 
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
  
write_csv(ascertainment_output, "outputs/sim_output.csv")

# could try and include the essential params in the file name rather than date which the system takes care of  

# 2)
# check trends in mode of detection of cases over time
# plot the number of cases 
# detected via symptomatic surveillance
# detected via active detection
# neither
  
ggplot(trim_sims) +
  aes(
    x=infection_day, 
    y=infections,
    group=simulation) + 
  facet_wrap(
    ~case_found_by,
    ncol=4
  ) +
  theme_cowplot() +
  geom_line(aes(color=simulation)) +
  theme(legend.position = "none") +
  labs(title="")
    
# ggsave(
#   'outputs/plots/infections_trajectories_by_mode_pdetect0.5.png',
#   height=6,
#   width=12,
#   bg='white'
# )

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

# 3)

# estimated number of infections = average cases found by contact tracing + average cases found by symptomatic pres * correction factor

# correction factor = 1/(symptomatic fraction * test-seeking fraction)
correction_factor <- 1/(0.307*0.2)

# first check that when contact tracing is switched off - that we can accurately estimate the number of detected infections using the correction factor
case_ascertainment <- trim_sims %>% 
  group_by(case_found_by,
           simulation) %>% # groups by detected/screening and simulation
  summarise(infections_by_sim_mode=sum(infections)) %>% # total infections by detection mode and sim
  group_by(simulation) %>% 
  mutate(all_infections=sum(infections_by_sim_mode)) %>% # total infections by sim
  filter(case_found_by=='symptomatic_surveillance') %>% 
  mutate(estimated_infections=infections_by_sim_mode*correction_factor) %>% 
  rename(detected_infections=infections_by_sim_mode) %>% 
  mutate(case_ascertainment=detected_infections/estimated_infections) %>% 
  mutate(true_case_ascertainment=detected_infections/all_infections) %>% 
  distinct(simulation, case_ascertainment, true_case_ascertainment) %>% 
  pivot_longer(cols=c(case_ascertainment, true_case_ascertainment))


  

# case ascertainment = (average cases found by contact tracing + average cases found by symptomatic pres)/ estimated number of infections

# compare case ascertainment to true case ascertainment

# for each simulation, calculate the estimated number of infections
metrics <- trim_sims %>% 
  group_by(case_found_by,
           simulation) %>% 
  summarise(sim_infections=sum(infections)) %>% 
  group_by(simulation) %>% 
  mutate(all_infections=sum(sim_infections))
  
case_ascertainment <- metrics %>%  
  summarise(detected_infections=all_infections-sim_infections[case_found_by=='undetected']) %>%
  distinct() %>% 
  left_join(metrics, 
            by='simulation') %>% 
  # do some pivot_wider here
  pivot_wider(
    names_from = case_found_by,
    values_from = sim_infections
  ) %>%
  mutate(
    estimated_infections = contact_tracing + workplace_screening + (symptomatic_surveillance * correction_factor)) %>% 
  # mutate(estimated_infections=sim_infections[case_found_by=='contact_tracing']+
  #          (sim_infections[case_found_by=='screening']*correction_factor)) %>% 
  mutate(case_ascertainment=detected_infections/estimated_infections) %>% 
  mutate(true_case_ascertainment=detected_infections/all_infections) %>% 
  distinct(simulation, case_ascertainment, true_case_ascertainment) %>% 
  pivot_longer(cols=c(case_ascertainment, true_case_ascertainment))



# multiple simulations - consider TTIQ response strategies
# generate cloud of true and estimated ascertainments  

ggplot(case_ascertainment) +
  geom_point(aes(x=factor(name),
                 y=value),
             position=position_jitter(width=0.1),
             colour='grey40',
             size=1) +
  # geom_pointrange(data = sim_sry,
  #                 aes(x=factor(1), 
  #                     y=mean, 
  #                     ymin = mean-2*se,
  #                     ymax = mean+2*se
  #                 ),
  #                 size=0.3) +
  #facet_grid(rel_active_detection_vaccinated_contact~rel_active_detection_vaccinated_source) +
  theme_cowplot() + 
  ylim(0,1) +
  theme(panel.grid=element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(size=0.2, color = 'black'),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1,
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background =element_rect(fill='white')) +
  labs(x='',
       y='case ascertainment')  

# ggsave(
#   #paste0('outputs/plots/ca_pdetect_', passive_detection_given_symptoms, '.png'),
#   'outputs/plots/ca_pdetect_0.2_no_tracing.png',
#   height=8,
#   width=5,
#   bg='white'
# )
