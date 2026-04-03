# data.R - Save epiworldR SEIR Simulation Results
# This script runs SEIR simulations and saves the key data components
# for later use with Cori and Wallinga-Teunis methods

# Load Required Libraries
suppressWarnings(suppressMessages(invisible({
  library(epiworldR)   # For epidemiological modeling and simulations
  library(data.table)  # Efficient data manipulation and aggregation
  library(dplyr)       # Data wrangling and manipulation
  library(tidyverse)   # Collection of data manipulation and visualization tools
  library(igraph)      # For Erdős–Rényi graph generation
})))

# ==== CONFIGURABLE R0 VALUES ====
R0_values <- c(1.5, 2.0, 3.0, 5.0)
#R0_values <- 2

# Function to run simulation and save data
save_seir_simulation_data <- function(
    model_type = "full",  # "full" for SEIRCONN, "partial" for SEIR with ER network
    name = "Covid",
    n = 1e5,
    prevalence = 100/1e5,
    contact_rate = 20.0,
    R0 = R0_values,
    recovery_rate = 1.0/7.0,
    incubation_days = 4,
    ndays = 150,
    nsim = 100,
    seed = 1234,
    save_dir = "saved_data/",
    deg = 20  # Average degree for ER network
) {
  
  # Create save directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Calculate transmission rate
  transmission_rate <- R0 * recovery_rate / contact_rate
  
  # Initialize model based on type
  if (model_type == "full") {
    # Full contact model (SEIRCONN)
    model <- ModelSEIRCONN(
      name = name,
      n = n,
      prevalence = prevalence,
      contact_rate = contact_rate,
      transmission_rate = transmission_rate,
      recovery_rate = recovery_rate,
      incubation_days = incubation_days
    )
  } else if (model_type == "partial") {
    # Partial contact model (SEIR with Erdős–Rényi network)
    set.seed(seed)
    
    # Generate Erdős–Rényi graph using igraph
    p_er <- deg / (n - 1)
    g <- erdos.renyi.game(n = n, p.or.m = p_er, type = "gnp", directed = FALSE)
    net2 <- as_edgelist(g, names = FALSE)
    net2 <- apply(net2, 2, as.integer) - 1L  # convert to integer and 0-indexed
    transmission_rate <- R0 * recovery_rate / deg
    # Create model from edgelist
    model <- ModelSEIR(
      name = paste(name, "ER"),
      prevalence = prevalence,
      transmission_rate = transmission_rate,
      incubation_days = incubation_days,
      recovery_rate = recovery_rate
    )
    
    # Create agents from edgelist
    agents_from_edgelist(
      model = model,
      source = net2[, 1],
      target = net2[, 2],
      size = as.integer(n),
      directed = FALSE
    )
    
    cat("Created Erdős–Rényi network with average degree:", deg, "\n")
    cat("Network has", nrow(net2), "edges\n")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Create saver object
  saver <- make_saver(
    "total_hist",
    "transmission",
    "transition",
    "reproductive",
    "generation"
  )
  
  # Run simulation
  cat("Running", nsim, "simulations for", ndays, "days with R0 =", R0, "...\n")
  run_multiple(model, ndays = ndays, nsim = nsim, seed = seed, saver = saver, nthreads = 18)
  
  # Get results
  results <- run_multiple_get_results(model, nthreads = 18)
  
  # Prepare transition data (the main dataframe for Cori/Wallinga-Teunis)
  df_transitions <- results$transition %>%
    rename(id = sim_num) %>%
    unite("transition", from, to, sep = " to ") %>%
    pivot_wider(names_from = transition, values_from = counts, values_fill = list(counts = 0))
  
  # Clean column names
  colnames(df_transitions) <- gsub(" ", "_", colnames(df_transitions))
  colnames(df_transitions) <- gsub("to", "to", colnames(df_transitions))
  colnames(df_transitions) <- tolower(colnames(df_transitions))
  
  # Floor the date and transition counts
  df_transitions$date <- floor(df_transitions$date)
  
  # Floor transition counts if they exist
  transition_cols <- grep("_to_", colnames(df_transitions), value = TRUE)
  for (col in transition_cols) {
    if (col %in% colnames(df_transitions)) {
      df_transitions[[col]] <- floor(df_transitions[[col]])
    }
  }
  
  # Get reproductive numbers
  reproductive_data <- results$reproductive
  reproductive_data <- subset(reproductive_data, source >= 0)
  
  # Calculate confidence intervals for Rt (ABM method)
  ci_results_rt <- reproductive_data %>%
    group_by(sim_num, source_exposure_date) %>%
    summarise(rt = mean(rt), .groups = "drop") %>%
    group_by(source_exposure_date) %>%
    summarise(
      mean_rt = mean(rt),
      ci_lower = quantile(rt, probs = 0.025),
      ci_upper = quantile(rt, probs = 0.975),
      n = n(),
      .groups = "drop"
    )
  
  # Get generation interval data
  generation_data <- results$generation %>%
    filter(gentime > 0)  # Remove zero generation times
  
  # Calculate generation interval statistics by simulation
  gi_stats <- generation_data %>%
    group_by(sim_num) %>%
    summarise(
      mean_gi = mean(gentime),
      sd_gi = sd(gentime),
      median_gi = median(gentime),
      n_gi = n(),
      .groups = "drop"
    )
  
  # Overall generation interval statistics
  overall_gi_stats <- list(
    mean_gi = mean(generation_data$gentime),
    sd_gi = sd(generation_data$gentime),
    median_gi = median(generation_data$gentime),
    n_total = nrow(generation_data)
  )
  
  # Create metadata
  metadata <- list(
    model_type = model_type,
    name = name,
    n = n,
    prevalence = prevalence,
    contact_rate = contact_rate,
    R0 = R0,
    transmission_rate = transmission_rate,
    recovery_rate = recovery_rate,
    incubation_days = incubation_days,
    ndays = ndays,
    nsim = nsim,
    seed = seed,
    run_date = Sys.time(),
    overall_gi_stats = overall_gi_stats,
    network_type = ifelse(model_type == "partial", "Erdos-Renyi", "Full Contact"),
    avg_degree = ifelse(model_type == "partial", deg, contact_rate)
  )
  
  # Create file prefix
  file_prefix <- paste0(save_dir, model_type, "_R0_", R0, "_n_", formatC(n, format = "e", digits = 0), "_nsim_", nsim)
  
  # Save all data components
  cat("Saving data to:", file_prefix, "...\n")
  
  # Save main transition dataframe (most important for Cori/Wallinga-Teunis)
  saveRDS(df_transitions, paste0(file_prefix, "_transitions.rds"))
  
  # Save reproductive number data
  saveRDS(reproductive_data, paste0(file_prefix, "_reproductive.rds"))
  saveRDS(ci_results_rt, paste0(file_prefix, "_rt_ci.rds"))
  
  # Save generation interval data
  saveRDS(generation_data, paste0(file_prefix, "_generation.rds"))
  saveRDS(gi_stats, paste0(file_prefix, "_gi_stats.rds"))
  
  # Save other results
  saveRDS(results$transmission, paste0(file_prefix, "_transmission.rds"))
  saveRDS(results$total_hist, paste0(file_prefix, "_total_hist.rds"))
  
  # Save metadata
  saveRDS(metadata, paste0(file_prefix, "_metadata.rds"))
  
  # Also save as CSV for easy inspection
  write.csv(df_transitions, paste0(file_prefix, "_transitions.csv"), row.names = FALSE)
  write.csv(ci_results_rt, paste0(file_prefix, "_rt_ci.csv"), row.names = FALSE)
  write.csv(gi_stats, paste0(file_prefix, "_gi_stats.csv"), row.names = FALSE)
  
  cat("Data saved successfully!\n")
  cat("Files saved with prefix:", file_prefix, "\n")
  
  # Return summary
  return(list(
    file_prefix = file_prefix,
    metadata = metadata,
    n_transitions = nrow(df_transitions),
    n_reproductive = nrow(reproductive_data),
    n_generation = nrow(generation_data),
    transition_columns = colnames(df_transitions)
  ))
}

# Function to load saved data
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

# Example usage: Save data for different R0 values and model types
if (FALSE) {  # Set to TRUE to run
  
  # Parameters
  n <- 1e5
  prevalence <- 100/n
  contact_rate <- 20.0
  recovery_rate <- 1.0/7.0
  incubation_days <- 4
  ndays <- 150
  nsim <- 100
  seed <- 1234
  
  # R0 values to test
  R0_values <- c(1.5, 2.0, 3.0, 5.0)
  
  # Save data for full contact model
  for (R0_val in R0_values) {
    save_seir_simulation_data(
      model_type = "full",
      R0 = R0_val,
      n = n,
      prevalence = prevalence,
      contact_rate = contact_rate,
      recovery_rate = recovery_rate,
      incubation_days = incubation_days,
      ndays = ndays,
      nsim = nsim,
      seed = seed
    )
  }
  
  # Save data for partial contact model (Erdős–Rényi)
  for (R0_val in R0_values) {
    save_seir_simulation_data(
      model_type = "partial",
      R0 = R0_val,
      n = n,
      prevalence = prevalence,
      contact_rate = 100,
      recovery_rate = recovery_rate,
      incubation_days = incubation_days,
      ndays = ndays,
      nsim = nsim,
      seed = seed,
      deg = 20  # Average degree for ER network
    )
  }
  
  # Example of loading data
  # data_full_R0_2 <- load_seir_data("saved_data/full_R0_2_n_1e+05_nsim_1000")
  # data_partial_R0_2 <- load_seir_data("saved_data/partial_R0_2_n_1e+05_nsim_1000")
}

# Helper function to inspect saved data
inspect_saved_data <- function(file_prefix) {
  data <- load_seir_data(file_prefix)
  
  cat("\n=== DATA INSPECTION ===\n")
  cat("Model type:", data$metadata$model_type, "\n")
  cat("Network type:", data$metadata$network_type, "\n")
  cat("R0:", data$metadata$R0, "\n")
  cat("Population size:", data$metadata$n, "\n")
  cat("Number of simulations:", data$metadata$nsim, "\n")
  cat("Simulation days:", data$metadata$ndays, "\n")
  
  if (data$metadata$model_type == "partial") {
    cat("Average degree:", data$metadata$avg_degree, "\n")
  }
  
  cat("\nTransition data dimensions:", dim(data$transitions), "\n")
  cat("Available transition columns:\n")
  transition_cols <- grep("_to_", colnames(data$transitions), value = TRUE)
  for (col in transition_cols) {
    cat(" -", col, "\n")
  }
  
  cat("\nGeneration interval stats:\n")
  cat(" - Overall mean GI:", round(data$metadata$overall_gi_stats$mean_gi, 2), "\n")
  cat(" - Overall SD GI:", round(data$metadata$overall_gi_stats$sd_gi, 2), "\n")
  
  cat("\nReproductive number data:", nrow(data$reproductive), "rows\n")
  cat("Date range:", min(data$transitions$date), "to", max(data$transitions$date), "\n")
  
  return(invisible(data))
}
