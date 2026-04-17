# Filtering

# INSTALL:
# remotes::install_github("dajmcdon/rtestim")

library(rtestim)
library(ggplot2)
library(epiworldR)

load_seir_data <- function(file_prefix) {
  cat("Loading data from:", file_prefix, "...\n")
  
  data <- list(
    transitions = readRDS(paste0(file_prefix, "_transitions.rds")),
    reproductive = readRDS(paste0(file_prefix, "_reproductive.rds")),
    rt_ci = readRDS(paste0(file_prefix, "_rt_ci.rds")),
    generation = readRDS(paste0(file_prefix, "_generation.rds")),
    gi_stats = readRDS(paste0(file_prefix, "_gi_stats.rds")),
    transmission = readRDS(paste0(file_prefix, "_transmission.rds")),
    total_hist = readRDS(paste0(file_prefix, "_total_hist.rds")),
    metadata = readRDS(paste0(file_prefix, "_metadata.rds"))
  )
  
  cat("Data loaded successfully!\n")
  return(data)
}

#
R0_values   <- c(1.5, 2.0, 3.0, 5.0)
model_types <- c("full", "partial")
transitions <- c("susceptible_to_exposed", "exposed_to_infected")

## ok start by just picking one
R0_value   <- R0_values[1]
model_type <- model_types[1]
transition <- transitions[2]

## loading the data
label  <- paste(model_type, "| R0 =", R0_value, "|", transition)
prefix <- paste0("saved_data/", model_type, 
                 "_R0_", R0_value, "_n_1e+05_nsim_100")

cat(strrep("-", 60), "\n")
cat("Scenario:", label, "\n")

# Load data
saved_data <- load_seir_data(prefix)

####
data.frame(head(saved_data$transitions))

## I think the column we want to model is exposed_to_infected
reports_df <- data.frame(date = saved_data$transitions$date, 
                         sim_id = saved_data$transitions$id, 
                         N = saved_data$transitions$exposed_to_infected)

## pick one simulation
reports_df <- reports_df %>% filter(sim_id == 1)
reports_df$sim_id <- NULL

# convert date to date
reports_df$date = lubridate::make_date(2020, 3, 1) + reports_df$date

# plot
ggplot(reports_df) + 
  geom_line(aes(x = date, y = N))

## and get the generation stats
gx <- saved_data$generation %>% filter(sim_num == 1)
mean(gx$gentime)
sd(gx$gentime)

hist(gx$gentime)

# so find a gamma that matches this distribution
table(gx$gentime)
dy <- as.numeric(names(table(gx$gentime)))
gxdf <- data.frame(table(gx$gentime))
gxdf$Var1 <- as.numeric(as.character(gxdf$Var1))

serial_interval_pmf <- data.frame(Var1 = 1:max(dy))
serial_interval_pmf <- serial_interval_pmf %>% left_join(gxdf)
serial_interval_pmf[is.na(serial_interval_pmf$Freq), 2] <- 1e-6
serial_interval_pmf$Freq <- serial_interval_pmf$Freq / 
  sum(serial_interval_pmf$Freq)

serial_interval_pmf <- serial_interval_pmf$Freq
serial_interval_pmf

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