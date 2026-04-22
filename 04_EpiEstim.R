# Fixed Sliding Windows

# install.packages('EpiEstim',
#                  repos = c('https://mrc-ide.r-universe.dev',
#                            'https://cloud.r-project.org'))
library(EpiEstim)
library(tidyverse)
library(data.table)
library(ggpubr)

#
source("00_make_reports_df.R")
setDT(reports_df)

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

setnames(reports_df, "N", "I")

get_epistim_output <- function(window_size) {
  
  # THEN SET SLIDING WINDOW SIZE
  estimation_window <- window_size
  
  t_start <- seq(2, nrow(reports_df) - estimation_window)
  t_end <- t_start + estimation_window
  
  reports_df$date_int <- 1:nrow(reports_df)
  
  # THEN ESTIMATE R(t) USING SLIDING WINDOWS
  getR <- EpiEstim::estimate_R(
    incid = reports_df,
    # method = "non_parametric_si",
    # config = make_config(list(
    #   si_distr = c(0, serial_interval_pmf),
    #   t_start = t_start,
    #   t_end = t_end
    # )),
    
    method = "parametric_si",
    config = make_config(list(
      mean_si = 6, std_si = 1,
      t_start = t_start,
      t_end = t_end)),
    
    backimputation_window = 15
  )
  
  # **********
  # INCLUDE THE DECONVOLUTION
  # getR$R$date_int <- getR$R$t_end - seeding_time
  getR$R$date_int <- getR$R$t_end
  # **********
  
  # **********
  # MOVE TO WINDOW CENTER
  getR$R$date_int <- getR$R$t_end - floor(estimation_window/2)
  # **********
  
  #
  EpiEstim_R <- getR$R[, c('t_start', 't_end', 'date_int', 'Median(R)',
                           'Quantile.0.05(R)', 'Quantile.0.95(R)')]
  
  names(EpiEstim_R) <- c('t_start', 't_end', "date_int", "Rt", "Rt_lb", "Rt_ub")
  
  EpiEstim_R$model <- paste0(window_size, " days")
  EpiEstim_R$date_int <- as.integer(EpiEstim_R$date_int)
  reports_df <- data.frame(reports_df)
  
  EpiEstim_R <- merge(EpiEstim_R, reports_df[, c('date', 'date_int')])
  
  return(EpiEstim_R)
  
}

EpiEstim_R_full <- get_epistim_output(15)

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

ggplot(EpiEstim_R_full) +
  ##
  theme_classic2() +
  #
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
  geom_ribbon(aes(x = date,
                  ymin = Rt_lb, ymax = Rt_ub,
                  fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date, y = Rt, color = model),
            linewidth = 0.25, show.legend = T) +
  geom_point(aes(x = date, y = Rt, color = model),
             size = 0.5,shape = 1,
             show.legend = T) +
  scale_color_discrete(name = 'Sliding window') +
  scale_fill_discrete(name = 'Sliding window') +
  ylab(expression(R[t])) +
  xlab(NULL) 

ggsave("img/EpiEstim.png", height = 5, width = 5, dpi = 600)

