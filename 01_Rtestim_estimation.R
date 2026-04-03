# =============================================================================
#  FILE 1 OF 2 — RtEstim Estimation
#
#  Purpose : Estimate the time-varying reproduction number (Rt) using
#            RtEstim across all scenarios, then validate against the
#            Agent-Based Model (ABM) ground truth.
#
#  Run this file FIRST. It writes results to rtestim_results/ which
#  the plotting script (02_rtestim_plots.R) reads.
#
#  Scenarios:
#    Model types : full, partial
#    R0 values   : 1.5, 2.0, 3.0, 5.0
#    Transitions : susceptible_to_exposed (S→E), exposed_to_infected (E→I)
#    → 2 × 4 × 2 = 16 scenario combinations in total
#
#  Outputs (written to rtestim_results/):
#    <model>_R0_<r0>_<transition>.rds   — data frame of per-sim Rt estimates
#    <model>_R0_<r0>_<transition>.csv   — same, in CSV format
#
#  Dependencies:
#    install.packages(c("rtestim", "dplyr", "parallel", "MASS"))
#
#  Author : sima
#  Date   : April 2026
# =============================================================================


# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

library(rtestim)    # Rt estimation via trend filtering
library(dplyr)      # data wrangling
library(parallel)   # mclapply — parallel across simulations
library(MASS)       # fitdistr — gamma fit to generation times

select <- dplyr::select   # prevent clash with MASS::select

source("data001.R")       # provides load_seir_data()


# ── 1. CONFIG ─────────────────────────────────────────────────────────────────
#
#  All scenario switches live here. Change these to subset the run.

CFG <- list(
  R0_values   = c(1.5, 2.0, 3.0, 5.0),
  model_types = c("full", "partial"),
  transitions = c("susceptible_to_exposed", "exposed_to_infected"),
  ncores      = 10,
  out_dir     = "rtestim_results"
)

dir.create(CFG$out_dir, showWarnings = FALSE)


# ── 2. GENERATION INTERVAL PMF ────────────────────────────────────────────────
#
#  Converts raw generation times from the ABM into a discrete PMF over
#  [0, max_day] for use in the RtEstim renewal equation.
#
#  Why exclude seeds and seed offspring?
#    Seed cases are artificially introduced; their generation times are
#    unrepresentative of natural transmission and bias the PMF.
#
#  Strategy:
#   1. Fit a Gamma distribution — smooth and principled.
#   2. Fall back to an empirical PMF if the Gamma fit fails.

build_gi_pmf <- function(gen_times, cap_quantile = 0.95) {
  
  # -- Attempt 1: parametric Gamma ------------------------------------------
  pmf_gamma <- tryCatch({
    
    fit <- MASS::fitdistr(gen_times, "gamma")
    shp <- fit$estimate["shape"]
    rat <- fit$estimate["rate"]
    
    cap     <- ceiling(qgamma(cap_quantile, shape = shp, rate = rat))
    days    <- 0:cap
    raw_pmf <- dgamma(days, shape = shp, rate = rat)
    raw_pmf <- raw_pmf / sum(raw_pmf)
    
    # Sanity checks — reject degenerate fits
    if (any(!is.finite(raw_pmf)) || any(raw_pmf < 0) || sum(raw_pmf) < 0.99) {
      NULL
    } else {
      list(pmf = raw_pmf, source = "gamma", shape = shp, rate = rat)
    }
    
  }, error = function(e) NULL)
  
  if (!is.null(pmf_gamma)) return(pmf_gamma)
  
  # -- Fallback: empirical PMF ----------------------------------------------
  cap       <- ceiling(quantile(gen_times, cap_quantile))
  gt_capped <- gen_times[gen_times <= cap]
  max_gen   <- max(gt_capped)
  
  gen_counts <- table(factor(gt_capped, levels = 0:max_gen))
  raw_pmf    <- as.numeric(gen_counts) / sum(gen_counts)
  raw_pmf    <- raw_pmf / sum(raw_pmf)
  
  list(pmf = raw_pmf, source = "empirical")
}


# ── 3. LAMBDA SELECTION HELPERS ───────────────────────────────────────────────
#
#  RtEstim's cross-validation produces two smoothing parameter candidates:
#    lambda.min — minimises CV error (potentially overfitted)
#    lambda.1se — 1-SE rule (more conservative, more smoothing)
#
#  We evaluate both on held-out Poisson log-likelihood (second half of the
#  incidence series) and pick the winner. When only one candidate produces
#  a valid confidence band, we use that one.

# Poisson log-likelihood: sum_t [ y_t * log(mu_t) - mu_t ]
poisson_loglik <- function(y, mu) {
  keep <- is.finite(y) & is.finite(mu) & mu > 0
  sum(y[keep] * log(mu[keep]) - mu[keep])
}

# Expected incidence under the renewal equation:
#   E[I_t] = Rt * sum_{s=1}^{max_w} w_s * I_{t-s}
compute_expected_counts <- function(rt_fit, incidence, pmf) {
  n     <- length(incidence)
  w     <- pmf
  max_w <- length(w) - 1
  expected <- rep(NA_real_, n)
  for (t in 2:n) {
    lags        <- 1:min(t - 1, max_w)
    weights     <- w[lags + 1]
    past_inc    <- incidence[t - lags]
    expected[t] <- rt_fit[t] * sum(weights * past_inc)
  }
  expected
}

# Safely call confband() — returns NULL on any failure or non-finite output
try_confband <- function(rtestim_fit, lam) {
  if (is.null(lam) || !is.finite(lam) || lam < 0) return(NULL)
  cb <- suppressWarnings(tryCatch(
    confband(rtestim_fit, lambda = lam),
    error = function(e) NULL
  ))
  if (is.null(cb) || any(!is.finite(cb$fit))) return(NULL)
  cb
}

# Pick the better lambda using held-out Poisson log-likelihood
select_best_lambda <- function(rtestim_fit, incidence, pmf) {
  
  cb_min <- try_confband(rtestim_fit, rtestim_fit$lambda.min)
  cb_1se <- try_confband(rtestim_fit, rtestim_fit$lambda.1se)
  
  # Graceful fallback when one candidate fails entirely
  if (is.null(cb_min) && is.null(cb_1se)) return(NULL)
  if (is.null(cb_min)) return(cb_1se)
  if (is.null(cb_1se)) return(cb_min)
  
  # Both valid — compare on held-out second half of the series
  n       <- length(incidence)
  val_idx <- floor(n / 2):n
  y       <- incidence[val_idx]
  
  mu_min <- compute_expected_counts(cb_min$fit, incidence, pmf)[val_idx]
  mu_1se <- compute_expected_counts(cb_1se$fit, incidence, pmf)[val_idx]
  
  ll_min <- poisson_loglik(y, mu_min)
  ll_1se <- poisson_loglik(y, mu_1se)
  
  if (is.finite(ll_min) && is.finite(ll_1se)) {
    if (ll_min >= ll_1se) cb_min else cb_1se
  } else if (is.finite(ll_min)) {
    cb_min
  } else {
    cb_1se
  }
}


# ── 4. SINGLE SIMULATION: ESTIMATE Rt ────────────────────────────────────────
#
#  Processing pipeline for one simulation:
#
#  (a) Extract incidence series; subtract artificially seeded cases.
#  (b) Trim leading zeros and a sparse tail (< 2 cases/day).
#  (c) Drop the first timepoint — at t = 1 there is no convolution
#      history, so weighted_past_counts = 0, which breaks CV silently.
#  (d) Build generation interval PMF (seeds + seed offspring excluded).
#  (e) Fit RtEstim with korder = 2 and 3.
#      Select korder = 3 only if CV error is >5% lower; otherwise use
#      korder = 2 (safer at series edges, especially for low-R0 runs).
#  (f) Select best smoothing lambda via held-out log-likelihood.
#  (g) Apply burn-in of 2 × mean generation interval; cap Rt at 20.

process_single_sim <- function(
    sim_id,
    transitions_data,
    generation_data,
    transmission_data,
    how_infected
) {
  
  tryCatch({
    
    # (a) Filter to this simulation
    sim_transitions  <- transitions_data  %>% filter(id == sim_id)
    sim_generation   <- generation_data   %>% filter(sim_num == sim_id)
    sim_transmission <- transmission_data %>% filter(sim_num == sim_id)
    
    incidence <- if (how_infected == "susceptible_to_exposed") {
      sim_transitions$susceptible_to_exposed
    } else {
      sim_transitions$exposed_to_infected
    }
    
    days <- sim_transitions$date
    
    # Subtract seed introductions from day 0 (S→E only)
    if (how_infected == "susceptible_to_exposed") {
      n_seeds      <- sum(sim_transmission$source == -1)
      incidence[1] <- max(0, incidence[1] - n_seeds)
    }
    
    # (b) Trim
    first_nz    <- which(incidence > 0)[1]
    if (is.na(first_nz)) return(NULL)
    
    last_stable <- max(which(incidence >= 2))
    if (is.na(last_stable) || last_stable < 10) return(NULL)
    
    incidence <- incidence[first_nz:last_stable]
    days      <- days[first_nz:last_stable]
    if (length(incidence) < 10) return(NULL)
    
    # (c) Drop first timepoint (no convolution history)
    incidence <- incidence[-1]
    days      <- days[-1]
    if (length(incidence) < 10) return(NULL)
    
    # (d) Generation interval PMF
    seed_ids       <- sim_transmission %>% filter(source == -1) %>% pull(target)
    seed_offspring <- sim_transmission %>%
      filter(source %in% seed_ids) %>% pull(target)
    exclude_ids    <- union(seed_ids, seed_offspring)
    
    gen_times <- sim_generation %>%
      filter(!(source %in% exclude_ids)) %>%
      pull(gentime) %>%
      { .[!is.na(.) & . > 0] }
    
    if (length(gen_times) < 10) return(NULL)
    
    gi_result   <- build_gi_pmf(gen_times, cap_quantile = 0.95)
    gen_int_pmf <- gi_result$pmf
    mean_gi     <- mean(gen_times)
    
    # (e) Fit korder 2 and 3
    fit_k <- function(k) {
      suppressWarnings(tryCatch(
        cv_estimate_rt(
          observed_counts = incidence,
          delay_distn     = gen_int_pmf,
          korder          = k,
          nsol            = 200,
          maxiter         = 1e6
        ),
        error = function(e) NULL
      ))
    }
    
    fit2 <- fit_k(2)
    fit3 <- fit_k(3)
    
    cvm_valid <- function(fit) {
      !is.null(fit) && length(fit$cvm) > 0 && any(is.finite(fit$cvm))
    }
    
    rtestim_fit <- if (cvm_valid(fit2) && cvm_valid(fit3)) {
      cv2 <- min(fit2$cvm, na.rm = TRUE)
      cv3 <- min(fit3$cvm, na.rm = TRUE)
      if ((cv2 - cv3) / cv2 > 0.05) fit3 else fit2
    } else if (!is.null(fit2)) {
      fit2
    } else if (!is.null(fit3)) {
      fit3
    } else {
      NULL
    }
    
    if (is.null(rtestim_fit)) return(NULL)
    
    # (f) Best lambda
    rtestim_cb <- select_best_lambda(rtestim_fit, incidence, gen_int_pmf)
    if (is.null(rtestim_cb)) return(NULL)
    
    # (g) Burn-in + cap
    burn_in <- ceiling(2 * mean_gi)
    
    data.frame(
      sim_id    = sim_id,
      method    = "RtEstim",
      date      = days,
      median_rt = rtestim_cb$fit,
      q025_rt   = rtestim_cb$`2.5%`,
      q975_rt   = rtestim_cb$`97.5%`
    ) %>%
      filter(date > burn_in) %>%
      mutate(across(c(median_rt, q025_rt, q975_rt), ~ pmin(.x, 20)))
    
  }, error = function(e) NULL)
}


# ── 5. PARALLEL WRAPPER ───────────────────────────────────────────────────────

run_rtestim_parallel <- function(saved_data, how_infected, ncores) {
  
  sim_nums <- unique(saved_data$transitions$id)
  cat("  Simulations:", length(sim_nums), "| Cores:", ncores, "\n")
  
  results <- parallel::mclapply(sim_nums, function(sid) {
    select <- dplyr::select    # re-declare inside forked process
    process_single_sim(
      sim_id            = sid,
      transitions_data  = saved_data$transitions,
      generation_data   = saved_data$generation,
      transmission_data = saved_data$transmission,
      how_infected      = how_infected
    )
  }, mc.cores = ncores)
  
  valid  <- Filter(is.data.frame, results)
  n_fail <- length(sim_nums) - length(valid)
  
  cat("  Succeeded:", length(valid), "/ Failed:", n_fail, "\n")
  if (length(valid) == 0) return(NULL)
  bind_rows(valid)
}


# ── 6. VALIDATION METRICS ─────────────────────────────────────────────────────
#
#  Per-simulation RMSE and MAE, then summarised across simulations.
#  Only days within [date_min, Inf) are compared (post-burn-in).

compute_validation <- function(results, abm_data, date_min = 12) {
  
  abm_clean <- abm_data %>%
    rename(date = source_exposure_date) %>%
    filter(date >= date_min) %>%
    select(date, mean_rt)
  
  per_sim <- results %>%
    filter(date >= date_min) %>%
    left_join(abm_clean, by = "date") %>%
    filter(!is.na(mean_rt)) %>%
    group_by(sim_id) %>%
    summarise(
      RMSE   = sqrt(mean((median_rt - mean_rt)^2)),
      MAE    = mean(abs(median_rt - mean_rt)),
      n_days = n(),
      .groups = "drop"
    )
  
  per_sim
}

print_validation <- function(val, label) {
  cat(sprintf("\n  ┌─ Validation: %s\n", label))
  cat(sprintf("  │  N sims validated : %d\n",       nrow(val)))
  cat(sprintf("  │  Mean RMSE        : %.4f\n",     mean(val$RMSE)))
  cat(sprintf("  │  Median RMSE      : %.4f\n",   median(val$RMSE)))
  cat(sprintf("  │  Mean MAE         : %.4f\n",     mean(val$MAE)))
  cat(sprintf("  └─ Median MAE       : %.4f\n\n", median(val$MAE)))
}


# ── 7. MAIN LOOP ──────────────────────────────────────────────────────────────
#
#  Iterates over all 16 scenario combinations.
#  For each scenario:
#    1. Load the ABM simulation data.
#    2. Run RtEstim in parallel across simulations.
#    3. Save results (.rds + .csv).
#    4. Validate against ABM ground truth and print metrics.

cat("\n", strrep("═", 62), "\n", sep = "")
cat("  RtEstim Estimation — All Scenarios\n")
cat("  Model types : ", paste(CFG$model_types, collapse = ", "), "\n", sep = "")
cat("  R0 values   : ", paste(CFG$R0_values,   collapse = ", "), "\n", sep = "")
cat("  Transitions : ", paste(CFG$transitions,  collapse = ", "), "\n", sep = "")
cat(strrep("═", 62), "\n\n", sep = "")

for (model_type in CFG$model_types) {
  for (R0_val in CFG$R0_values) {
    for (transition in CFG$transitions) {
      
      label  <- sprintf("%s | R0 = %s | %s", model_type, R0_val, transition)
      prefix <- sprintf(
        "saved_data/%s_R0_%s_n_1e+05_nsim_100",
        model_type, R0_val
      )
      
      cat(strrep("─", 62), "\n")
      cat("Scenario:", label, "\n")
      
      if (!file.exists(paste0(prefix, "_metadata.rds"))) {
        cat("  ⚠  Data not found, skipping.\n\n")
        next
      }
      
      # 1. Load
      saved_data <- load_seir_data(prefix)
      
      # 2. Estimate
      rt_results <- run_rtestim_parallel(
        saved_data,
        how_infected = transition,
        ncores       = CFG$ncores
      )
      
      if (is.null(rt_results)) {
        cat("  ⚠  No results produced.\n\n")
        next
      }
      
      # 3. Save
      out_stem <- file.path(
        CFG$out_dir,
        sprintf("%s_R0_%s_%s", model_type, R0_val, transition)
      )
      saveRDS(rt_results, paste0(out_stem, ".rds"))
      write.csv(rt_results, paste0(out_stem, ".csv"), row.names = FALSE)
      cat("  ✓ Results saved →", basename(out_stem), "\n")
      
      # 4. Validate
      abm_file <- paste0(prefix, "_rt_ci.rds")
      if (file.exists(abm_file)) {
        abm_data <- readRDS(abm_file)
        val      <- compute_validation(rt_results, abm_data)
        print_validation(val, label)
      } else {
        cat("  ⚠  ABM ground-truth file not found:", abm_file, "\n\n")
      }
    }
  }
}

cat(strrep("═", 62), "\n", sep = "")
cat("  Estimation complete.\n")
cat("  Results written to: ", CFG$out_dir, "/\n", sep = "")
cat("  → Now run 02_rtestim_plots.R\n")
cat(strrep("═", 62), "\n\n", sep = "")