# =============================================================================
#  RtEstim Estimation — All Scenarios
#
#  Based on Chad's codes:
#    1. Fit a gamma to the generation times  (same method advisor used for SI)
#    2. Run cv_estimate_rt() with dist_gamma = c(shape, rate)
#    3. Extract confidence band with confband()
#    4. Save results and print RMSE / MAE vs ABM ground truth
#
#  Scenarios:
#    Model types : full, partial
#    R0 values   : 1.5, 2.0, 3.0, 5.0
#    Transitions : susceptible_to_exposed (S->E), exposed_to_infected (E->I)
#
#  Run this FIRST, then run 02_rtestim_plots.R
#
#  Dependencies:
#    install.packages(c("rtestim", "dplyr", "parallel", "MASS", "nloptr"))
# =============================================================================


# ── PACKAGES ──────────────────────────────────────────────────────────────────

library(rtestim)
library(dplyr)
library(parallel)
library(MASS)
library(nloptr)

select <- dplyr::select

source("data001.R")   # load_seir_data()


# ── CONFIG ────────────────────────────────────────────────────────────────────

R0_values   <- c(1.5, 2.0, 3.0, 5.0)
model_types <- c("full", "partial")
transitions <- c("susceptible_to_exposed", "exposed_to_infected")
ncores      <- 8
out_dir     <- "rtestim_results"

dir.create(out_dir, showWarnings = FALSE)


# ── STEP 1: FIT GAMMA TO GENERATION TIMES ────────────────────────────────────
#
#  Advisor's approach: fit a gamma distribution to the serial interval
#  using nloptr, then pass c(shape, rate) directly to cv_estimate_rt().
#  We do the same but with generation times from the ABM.
#
#  The objective function maximises correlation between the fitted gamma
#  PMF and the observed discrete distribution.

fit_gamma_to_gentimes <- function(gen_times) {
  
  # Build empirical discrete PMF from generation times
  max_gt  <- max(gen_times)
  emp_pmf <- as.numeric(table(factor(gen_times, levels = 1:max_gt)))
  emp_pmf <- emp_pmf / sum(emp_pmf)
  
  # Remove leading zeros that would confuse the optimizer
  emp_pmf[emp_pmf == 0] <- 1e-6
  
  # Objective: maximise correlation between gamma PMF and empirical PMF
  # (same approach as advisor's code for serial interval)
  my_eval <- function(params, data) {
    si_shape <- params[1]
    si_rate  <- params[2]
    if (si_shape <= 0 || si_rate <= 0) return(1e10)
    w <- sapply(1:length(data), function(x) {
      pgamma(x, si_shape, si_rate) - pgamma(x - 1, si_shape, si_rate)
    })
    cor(w, data) * -1   # negative because nloptr minimises
  }
  
  fit <- tryCatch(
    nloptr(
      x0     = c(shape = 2, rate = 0.5),
      eval_f = my_eval,
      lb     = c(1e-6, 1e-6),
      ub     = c(100,  100),
      opts   = list(
        algorithm  = "NLOPT_LN_SBPLX",
        xtol_rel   = 1e-8,
        maxeval    = 1000
      ),
      data = emp_pmf
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    # Fallback: simple method-of-moments gamma fit via MASS
    fit_mm <- tryCatch(MASS::fitdistr(gen_times, "gamma"), error = function(e) NULL)
    if (is.null(fit_mm)) return(NULL)
    return(c(
      shape = fit_mm$estimate["shape"],
      rate  = fit_mm$estimate["rate"]
    ))
  }
  
  c(shape = fit$solution[1], rate = fit$solution[2])
}


# ── STEP 2: PROCESS ONE SIMULATION ───────────────────────────────────────────
#
#  For each simulation:
#   (a) Get incidence series, subtract seed cases
#   (b) Trim leading zeros and sparse tail
#   (c) Drop first timepoint (no convolution history at t=1)
#   (d) Fit gamma to generation times -> get shape and rate
#   (e) Run cv_estimate_rt() with dist_gamma = c(shape, rate)  <- advisor's way
#   (f) Get confidence band with confband()
#   (g) Apply burn-in, cap Rt at 20

process_single_sim <- function(sim_id,
                               transitions_data,
                               generation_data,
                               transmission_data,
                               how_infected) {
  tryCatch({
    
    # (a) Extract incidence
    sim_tr  <- transitions_data  %>% filter(id == sim_id)
    sim_gen <- generation_data   %>% filter(sim_num == sim_id)
    sim_tx  <- transmission_data %>% filter(sim_num == sim_id)
    
    incidence <- if (how_infected == "susceptible_to_exposed") {
      sim_tr$susceptible_to_exposed
    } else {
      sim_tr$exposed_to_infected
    }
    days <- sim_tr$date
    
    # Subtract artificially seeded cases from day 0
    if (how_infected == "susceptible_to_exposed") {
      n_seeds      <- sum(sim_tx$source == -1)
      incidence[1] <- max(0, incidence[1] - n_seeds)
    }
    
    # (b) Trim leading zeros and sparse tail
    first_nz    <- which(incidence > 0)[1]
    if (is.na(first_nz)) return(NULL)
    
    last_stable <- max(which(incidence >= 2))
    if (is.na(last_stable) || last_stable < 10) return(NULL)
    
    incidence <- incidence[first_nz:last_stable]
    days      <- days[first_nz:last_stable]
    if (length(incidence) < 10) return(NULL)
    
    # (c) Drop first timepoint — weighted_past_counts = 0 at t=1 breaks CV
    incidence <- incidence[-1]
    days      <- days[-1]
    if (length(incidence) < 10) return(NULL)
    
    # (d) Fit gamma to generation times
    #     Exclude seed cases and their direct offspring (unrepresentative GTs)
    seed_ids    <- sim_tx %>% filter(source == -1) %>% pull(target)
    seed_off    <- sim_tx %>% filter(source %in% seed_ids) %>% pull(target)
    exclude_ids <- union(seed_ids, seed_off)
    
    gen_times <- sim_gen %>%
      filter(!(source %in% exclude_ids)) %>%
      pull(gentime) %>%
      { .[!is.na(.) & . > 0] }
    
    if (length(gen_times) < 10) return(NULL)
    
    gamma_params <- fit_gamma_to_gentimes(gen_times)
    if (is.null(gamma_params)) return(NULL)
    
    si_shape <- gamma_params["shape"]
    si_rate  <- gamma_params["rate"]
    mean_gi  <- si_shape / si_rate   # mean of gamma distribution
    
    # (e) Fit RtEstim — using dist_gamma just like advisor's code
    #     Try korder 2 and 3, pick the one with lower CV error
    fit_k <- function(k) {
      suppressWarnings(tryCatch(
        cv_estimate_rt(
          dist_gamma      = c(si_shape, si_rate),
          observed_counts = incidence,
          korder          = k,
          nsol            = 200,
          maxiter         = 1e6
        ),
        error = function(e) NULL
      ))
    }
    
    fit2 <- fit_k(2)
    fit3 <- fit_k(3)
    
    # Pick korder with lower CV error; default to korder=2 if CV unavailable
    # (korder=2 is safer at series edges, especially for low R0)
    cvm_ok <- function(fit) {
      !is.null(fit) && length(fit$cvm) > 0 && any(is.finite(fit$cvm))
    }
    
    best_fit <- if (cvm_ok(fit2) && cvm_ok(fit3)) {
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
    
    if (is.null(best_fit)) return(NULL)
    
    # (f) Confidence band — try lambda.min first, fall back to lambda.1se
    #     (same options as advisor's code)
    get_cb <- function(lam) {
      cb <- suppressWarnings(tryCatch(
        confband(best_fit, lambda = lam),
        error = function(e) NULL
      ))
      if (is.null(cb) || any(!is.finite(cb$fit))) return(NULL)
      cb
    }
    
    cb <- get_cb("lambda.min")
    if (is.null(cb)) cb <- get_cb("lambda.1se")
    if (is.null(cb)) return(NULL)
    
    # (g) Apply burn-in of 2x mean GI, cap Rt at 20
    burn_in <- ceiling(2 * mean_gi)
    
    data.frame(
      sim_id    = sim_id,
      date      = days,
      median_rt = cb$fit,
      q025_rt   = cb$`2.5%`,
      q975_rt   = cb$`97.5%`
    ) %>%
      filter(date > burn_in) %>%
      mutate(across(c(median_rt, q025_rt, q975_rt), ~ pmin(.x, 20)))
    
  }, error = function(e) NULL)
}


# ── STEP 3: RUN IN PARALLEL ACROSS SIMULATIONS ───────────────────────────────

run_parallel <- function(saved_data, how_infected, ncores) {
  
  sim_nums <- unique(saved_data$transitions$id)
  cat("  Simulations:", length(sim_nums), "| Cores:", ncores, "\n")
  
  results <- parallel::mclapply(sim_nums, function(sid) {
    select <- dplyr::select   # re-declare in forked process
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


# ── STEP 4: VALIDATION vs ABM GROUND TRUTH ───────────────────────────────────

validate <- function(results, abm_data, date_min = 12) {
  
  abm_clean <- abm_data %>%
    rename(date = source_exposure_date) %>%
    filter(date >= date_min) %>%
    select(date, mean_rt)
  
  results %>%
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
}

print_validation <- function(val, label) {
  cat("\n  Validation:", label, "\n")
  cat("  N sims     :", nrow(val), "\n")
  cat("  Mean RMSE  :", round(mean(val$RMSE),   4), "\n")
  cat("  Median RMSE:", round(median(val$RMSE), 4), "\n")
  cat("  Mean MAE   :", round(mean(val$MAE),    4), "\n")
  cat("  Median MAE :", round(median(val$MAE),  4), "\n\n")
}


# ── MAIN LOOP ─────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n", sep = "")
cat("  RtEstim Estimation — All Scenarios\n")
cat(strrep("=", 60), "\n\n", sep = "")

for (model_type in model_types) {
  for (R0_val in R0_values) {
    for (transition in transitions) {
      
      label  <- paste(model_type, "| R0 =", R0_val, "|", transition)
      prefix <- paste0("saved_data/", model_type, "_R0_", R0_val, "_n_1e+05_nsim_100")
      
      cat(strrep("-", 60), "\n")
      cat("Scenario:", label, "\n")
      
      # Check data exists
      if (!file.exists(paste0(prefix, "_metadata.rds"))) {
        cat("  WARNING: Data not found, skipping.\n\n")
        next
      }
      
      # Load data
      saved_data <- load_seir_data(prefix)
      
      # Estimate Rt
      rt_results <- run_parallel(saved_data, transition, ncores)
      
      if (is.null(rt_results)) {
        cat("  WARNING: No results produced.\n\n")
        next
      }
      
      # Save
      out_stem <- file.path(out_dir, paste0(model_type, "_R0_", R0_val, "_", transition))
      saveRDS(rt_results,  paste0(out_stem, ".rds"))
      write.csv(rt_results, paste0(out_stem, ".csv"), row.names = FALSE)
      cat("  Saved:", basename(out_stem), "\n")
      
      # Validate vs ABM ground truth
      abm_file <- paste0(prefix, "_rt_ci.rds")
      if (file.exists(abm_file)) {
        abm_data <- readRDS(abm_file)
        val      <- validate(rt_results, abm_data)
        print_validation(val, label)
      } else {
        cat("  WARNING: ABM ground-truth file not found:", abm_file, "\n\n")
      }
    }
  }
}

cat(strrep("=", 60), "\n", sep = "")
cat("  Done. Now run 02_rtestim_plots.R\n")
cat(strrep("=", 60), "\n\n", sep = "")
