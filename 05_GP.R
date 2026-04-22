# Random Walk


# install.packages("EpiNow2",
#                  repos = c("https://epiforecasts.r-universe.dev",
#                            getOption("repos")))

#
library(EpiNow2)
library(tidyverse)
library(data.table)

#
source("00_make_reports_df.R")
setDT(reports_df)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
#
# GP refs:
#  - https://peterroelants.github.io/posts/gaussian-process-kernels/
#  - https://epiforecasts.io/EpiNow2/articles/gaussian_process_implementation_details.html
#  - https://wellcomeopenresearch.org/articles/5-112
#  - https://katbailey.github.io/post/gaussian-processes-for-dummies/
#  - https://epiforecasts.io/EpiNow2/articles/estimate_infections.html
#  - https://thegradient.pub/gaussian-process-not-quite-for-dummies/
#  - https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssa.12919
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

##
gi_pmf         <- Gamma(shape = si_shape, rate = si_rate)
infect_to_test <- NonParametric(pmf = c(1))
sym_report_delay_pmf <- NonParametric(pmf = c(1))

mean(gi_pmf) + mean(infect_to_test) + mean(sym_report_delay_pmf)

setnames(reports_df, 'N', 'confirm')

## -------------------
# many things that you could change
get_EpiNow2output <- function(l = 1.5, b = 0.2, kernel = 'se') {
  
  res_epinow <- epinow(
    data                 = reports_df[, c(1,2)],
    generation_time      = generation_time_opts(gi_pmf),
    delays               = delay_opts(infect_to_test),
    truncation           = trunc_opts(sym_report_delay_pmf),
    rt                   = rt_opts(), # use default
    gp                   = gp_opts(boundary_scale = l,
                                   basis_prop = b,
                                   kernel = kernel),
    backcalc             = backcalc_opts(prior = 'reports'),
    stan                 = stan_opts(chains = 4, cores = 4),
    obs                  = obs_opts(), # ok to use these defaults
    forecast             = forecast_opts(),
    CrIs                 = c(0.2, 0.5, 0.9)
  )
  
  R_df <- subset(summary(res_epinow, type = 'parameters'), variable == 'R')
  
  R_df$model <- paste0("l = ", l)
  R_df$Rt <- R_df$median
  R_df$Rt_lb <- R_df$lower_90
  R_df$Rt_ub <- R_df$upper_90
  return(R_df)
}

R_df <- get_EpiNow2output(l = 3)


##
R_df_all <- R_df

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# PLOT R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

#####
load_abm <- function(model_type, R0_val) {
  path <- file.path(
    'saved_data',
    sprintf("%s_R0_%s_n_1e+05_nsim_100_rt_ci.rds", model_type, R0_val)
  )
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

abm_raw <- load_abm(model_type, R0_value)

abm_summary <- NULL
if (!is.null(abm_raw) && nrow(abm_raw) > 0) {
  abm_summary <- abm_raw %>%
    rename(date = source_exposure_date) %>%
    # filter(date >= date_min, date <= date_max) %>%
    select(date, mean_rt, ci_lower, ci_upper) %>%
    mutate(date = lubridate::make_date(2020, 3, 1) + date)
}

abm_summary
####

first_day <- min(reports_df$date)
last_day <- max(reports_df$date)

library(ggpubr)

ggplot(R_df_all) +
  ##
  theme_classic2() +
  geom_hline(yintercept = 1, linetype = '11') +
  ##
  geom_ribbon(aes(x = date,
                  ymin = ci_lower, ymax = ci_upper),
              alpha = 0.25, fill = 'grey75',
              data = abm_summary) +
  geom_line(aes(x = date, y = mean_rt),
            linewidth = 0.25, 
            data = abm_summary) +
  # ##
  # ##
  geom_ribbon(aes(x = date,
                  ymin = Rt_lb, ymax = Rt_ub,
                  fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date, y = Rt, color = model),
            linewidth = 0.25, show.legend = T) +
  geom_point(aes(x = date, y = Rt, color = model),
             size = 0.5,shape = 1,
             show.legend = T) +
  scale_color_discrete(name = 'Length scale') +
  scale_fill_discrete(name = 'Length scale') +
  ylab(expression(R[t])) +
  xlab(NULL) 

ggsave("img/EPINOW.png", height = 5, width = 5, dpi = 600)
