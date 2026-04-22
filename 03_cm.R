# Filtering

source("00_make_reports_df.R")

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# GET GAMMA SERIAL INTERVAL
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

# RtESTIM assumes that the serial interval is gamma
# here is code that will find a gamma for your discrete distribution

library(nloptr)

# serial interval
x <- serial_interval_pmf

#
x[x < 1e-6] <- 1e-6
x

#
my_eval <- function(params, data) {
  
  si_shape <- params[1]
  si_rate <- params[2]
  
  if (si_shape <= 0 || si_rate <= 0) return(1e10)  # penalize invalid params
  
  w <- sapply(1:(length(data)), function(x){
    pgamma(x, si_shape, si_rate) - pgamma(x - 1, si_shape, si_rate)
  })
  
  cor(w, data) * -1
}


# Initial parameter guesses
init_params <- c(shape = 2, rate = 10)

# Bounds (positive shape and rate)
lower_bounds <- c(1e-6, 1e-6)
upper_bounds <- c(100, 1000)

# Optimize using nloptr
library(nloptr)
fit <- nloptr(
  x0 = init_params,
  eval_f = my_eval,
  lb = lower_bounds,
  ub = upper_bounds,
  opts = list("algorithm" = "NLOPT_LN_SBPLX",
              "xtol_rel" = 1e-8,
              "maxeval" = 1000),
  data = x
)

# View estimated parameters
si_shape = fit$solution[1]
si_rate  = fit$solution[2]

w <- sapply(1:(length(serial_interval_pmf)), function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x - 1, si_shape, si_rate)
})

# looks pretty good
plot(x)
lines(w, col = 'red')

si_shape
si_rate

# ////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------
# MODEL R(t)
# ----------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////////////

get_rtestim_output <- function(k) {
  
  # rt estimation
  rtestim <- cv_estimate_rt(
    dist_gamma      = c(si_shape, si_rate),
    observed_counts = reports_df$N,
    x               = reports_df$date,
    korder          = k,
    nsol            = 1000,
    maxiter         = 1e8
  )
  
  #approximate confidence bands
  rtestim_cb <- confband(rtestim, lambda = "lambda.1se")
  #lambda: the selected lambda. May be a scalar value,
  # or in the case of cv_poisson_rt objects, "lambda.min" or "lambda.max"
  
  # create dataframe
  # if you want to, shift backwards by the seeding_time, aka the sum
  # of the means of the delay distributions
  plot_rtestim <- data.frame(
    model = paste0("k = ",k),
    ## *****
    ## date = reports_df$date - seeding_time,
    date = reports_df$date,
    ###
    Rt = rtestim_cb$fit,
    Rt_lb = rtestim_cb$`2.5%`,
    Rt_ub = rtestim_cb$`97.5%`
  )
  return(plot_rtestim)
}

df1 <- get_rtestim_output(1)
df3 <- get_rtestim_output(3)
rtestim_df <- rbind(df1, df3)

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

date_min    = 12     # first day shown in plots (post-burn-in)
date_max    = 75     # last day shown

abm_raw <- load_abm(model_type, R0_value)

abm_summary <- NULL
if (!is.null(abm_raw) && nrow(abm_raw) > 0) {
  abm_summary <- abm_raw %>%
    rename(date = source_exposure_date) %>%
    filter(date >= date_min, date <= date_max) %>%
    select(date, mean_rt, ci_lower, ci_upper) %>%
    mutate(date = lubridate::make_date(2020, 3, 1) + date)
}

abm_summary

#####

library(ggpubr)

rt_max <- 3
first_day <- min(reports_df$date)
last_day <- max(reports_df$date)
nowcast_start    = last_day - seeding_time
forecast_window  = last_day + 15

ggplot(rtestim_df) +
  ##
  theme_classic2() +
  geom_hline(yintercept = 1, linetype = '11') +
  ##
  geom_ribbon(aes(x = date,
                  ymin = ci_lower, ymax = ci_upper),
              alpha = 0.25, fill = 'grey75',
              data = abm_summary) +
  # ##
  geom_ribbon(aes(x = date,
                  ymin = Rt_lb, ymax = Rt_ub,
                  fill = model),
              alpha = 0.25) +
  geom_line(aes(x = date, y = Rt, color = model),
            linewidth = 0.25, show.legend = T) +
  scale_color_discrete(name = 'Filter degree') +
  scale_fill_discrete(name = 'Filter degree') +
  ylab(expression(R[t])) +
  xlab(NULL) 

plot_rt1


