# =============================================================================
#  FILE 2 OF 2 — RtEstim vs ABM Ground Truth: Plots
#
#  Purpose : Read the Rt estimates produced by 01_rtestim_estimation.R
#            and generate interactive HTML plots comparing RtEstim to
#            the ABM ground truth.
#
#  Run AFTER 01_rtestim_estimation.R.
#
#  Output layout:
#    One HTML file per transition type, each containing a 2×4 grid:
#      rows    → model type (full, partial)
#      columns → R0 (1.5, 2.0, 3.0, 5.0)
#
#    rt_plots/
#      grid_susceptible_to_exposed.html
#      grid_exposed_to_infected.html
#      index.html                         ← open this in a browser
#
#  Each panel shows:
#    ● ABM ground truth   — dark line + shaded 95% CI band
#    ● RtEstim estimate   — green line + shaded 95% CI band
#    ── Rt = 1 reference dashed line
#
#  Dependencies:
#    install.packages(c("dplyr", "plotly", "htmlwidgets"))
#
#  Author : [Your name]
#  Date   : April 2026
# =============================================================================


# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

library(dplyr)
library(plotly)
library(htmlwidgets)

select <- dplyr::select


# ── 1. CONFIG ─────────────────────────────────────────────────────────────────

CFG <- list(
  R0_values   = c(1.5, 2.0, 3.0, 5.0),
  model_types = c("full", "partial"),
  transitions = c("susceptible_to_exposed", "exposed_to_infected"),
  date_min    = 12,     # first day shown in plots (post-burn-in)
  date_max    = 75,     # last day shown
  in_dir_rt   = "rtestim_results",
  in_dir_abm  = "saved_data",
  out_dir     = "rt_plots"
)

dir.create(CFG$out_dir, showWarnings = FALSE)

# Colour scheme — only two series needed
COLORS <- list(
  ABM     = list(line = "#2C3E50",  ribbon = "rgba(44,  62,  80,  0.13)"),
  RtEstim = list(line = "#27AE60",  ribbon = "rgba(39,  174, 96,  0.20)")
)


# ── 2. DATA LOADERS ───────────────────────────────────────────────────────────

load_rtestim <- function(model_type, R0_val, transition) {
  stem <- file.path(
    CFG$in_dir_rt,
    sprintf("%s_R0_%s_%s", model_type, R0_val, transition)
  )
  rds <- paste0(stem, ".rds")
  csv <- paste0(stem, ".csv")
  if (file.exists(rds)) return(readRDS(rds))
  if (file.exists(csv)) return(read.csv(csv))
  NULL
}

load_abm <- function(model_type, R0_val) {
  path <- file.path(
    CFG$in_dir_abm,
    sprintf("%s_R0_%s_n_1e+05_nsim_100_rt_ci.rds", model_type, R0_val)
  )
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}


# ── 3. SUMMARISE DATA FOR ONE PANEL ───────────────────────────────────────────
#
#  For RtEstim: average the point estimate across simulations;
#  use the 2.5/97.5 percentile of per-simulation medians as the
#  uncertainty band (simulation-to-simulation variability).
#
#  For ABM: the loaded object already contains mean_rt, ci_lower, ci_upper
#  computed across the ensemble of ABM runs.

prepare_panel_data <- function(model_type, R0_val, transition) {
  
  rt_raw <- load_rtestim(model_type, R0_val, transition)
  abm_raw <- load_abm(model_type, R0_val)
  
  # RtEstim summary
  rt_summary <- NULL
  if (!is.null(rt_raw) && nrow(rt_raw) > 0) {
    rt_summary <- rt_raw %>%
      filter(date >= CFG$date_min, date <= CFG$date_max) %>%
      group_by(date) %>%
      summarise(
        mean_rt  = mean(median_rt,            na.rm = TRUE),
        ci_lower = quantile(median_rt, 0.025, na.rm = TRUE),
        ci_upper = quantile(median_rt, 0.975, na.rm = TRUE),
        .groups  = "drop"
      )
  }
  
  # ABM summary
  abm_summary <- NULL
  if (!is.null(abm_raw) && nrow(abm_raw) > 0) {
    abm_summary <- abm_raw %>%
      rename(date = source_exposure_date) %>%
      filter(date >= CFG$date_min, date <= CFG$date_max) %>%
      select(date, mean_rt, ci_lower, ci_upper)
  }
  
  list(rt = rt_summary, abm = abm_summary)
}


# ── 4. BUILD ONE PANEL ────────────────────────────────────────────────────────
#
#  A single plotly trace set for one (model_type, R0) cell of the grid.
#  `show_legend` is TRUE only for the top-left panel so the legend
#  appears exactly once in the combined subplot.

build_panel <- function(data, R0_val, model_type, show_legend = FALSE) {
  
  p <- plot_ly()
  
  # ── ABM ground truth ──────────────────────────────────────────────────────
  if (!is.null(data$abm) && nrow(data$abm) > 0) {
    p <- p %>%
      add_ribbons(
        data      = data$abm,
        x         = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = COLORS$ABM$ribbon,
        line      = list(width = 0),
        legendgroup = "ABM", showlegend = FALSE, hoverinfo = "skip"
      ) %>%
      add_lines(
        data = data$abm,
        x = ~date, y = ~mean_rt,
        line        = list(color = COLORS$ABM$line, width = 3.5),
        name        = "ABM (True R\u209c)",
        legendgroup = "ABM",
        showlegend  = show_legend,
        hovertemplate = paste0(
          "<b>ABM (True R\u209c)</b><br>",
          "Day: %{x}<br>R\u209c: %{y:.3f}<extra></extra>"
        )
      )
  }
  
  # ── RtEstim ───────────────────────────────────────────────────────────────
  if (!is.null(data$rt) && nrow(data$rt) > 0) {
    p <- p %>%
      add_ribbons(
        data      = data$rt,
        x         = ~date, ymin = ~ci_lower, ymax = ~ci_upper,
        fillcolor = COLORS$RtEstim$ribbon,
        line      = list(width = 0),
        legendgroup = "RtEstim", showlegend = FALSE, hoverinfo = "skip"
      ) %>%
      add_lines(
        data = data$rt,
        x = ~date, y = ~mean_rt,
        line        = list(color = COLORS$RtEstim$line, width = 2.8),
        name        = "RtEstim",
        legendgroup = "RtEstim",
        showlegend  = show_legend,
        hovertemplate = paste0(
          "<b>RtEstim</b><br>",
          "Day: %{x}<br>R\u209c: %{y:.3f}<extra></extra>"
        )
      )
  }
  
  # ── Rt = 1 reference ──────────────────────────────────────────────────────
  p <- p %>%
    add_segments(
      x = CFG$date_min, xend = CFG$date_max, y = 1, yend = 1,
      line      = list(color = "#95A5A6", width = 1.4, dash = "dash"),
      showlegend = FALSE, hoverinfo = "skip"
    )
  
  # ── Axis layout ───────────────────────────────────────────────────────────
  axis_style <- list(
    gridcolor = "#E0E4EA", gridwidth = 1,
    showgrid  = TRUE,
    showline  = TRUE, linewidth = 1.5, linecolor = "#2C3E50", mirror = TRUE,
    tickfont  = list(size = 11)
  )
  
  p %>% layout(
    xaxis = c(
      axis_style,
      list(
        title = list(text = "<b>Day</b>", font = list(size = 12)),
        range = c(CFG$date_min, CFG$date_max)
      )
    ),
    yaxis = c(
      axis_style,
      list(
        title = list(text = "<b>R<sub>t</sub></b>", font = list(size = 12))
      )
    ),
    plot_bgcolor  = "#F8F9FA",
    paper_bgcolor = "white",
    hovermode     = "x unified",
    # Panel annotation: model type label (top-left) + R0 label (top-right)
    annotations = list(
      list(
        x = 0.02, y = 1.06, xref = "paper", yref = "paper",
        text      = sprintf("<b>%s</b>", toupper(model_type)),
        showarrow = FALSE, xanchor = "left",
        font      = list(size = 11, color = "#555")
      ),
      list(
        x = 0.98, y = 1.06, xref = "paper", yref = "paper",
        text      = sprintf("<b>R\u2080 = %s</b>", R0_val),
        showarrow = FALSE, xanchor = "right",
        font      = list(size = 13, color = "#2C3E50")
      )
    ),
    margin = list(t = 45, b = 40, l = 50, r = 10)
  )
}


# ── 5. BUILD 2×4 GRID FOR ONE TRANSITION ─────────────────────────────────────
#
#  Rows = model types (full, partial)
#  Columns = R0 values (1.5, 2.0, 3.0, 5.0)

build_grid <- function(transition) {
  
  trans_label <- if (transition == "susceptible_to_exposed") {
    "S \u2192 E (Susceptible to Exposed)"
  } else {
    "E \u2192 I (Exposed to Infected)"
  }
  
  cat("\n  Building grid for:", trans_label, "\n")
  
  # Collect all panels row-by-row (subplot fills left→right, top→bottom)
  panels     <- list()
  panel_idx  <- 1
  
  for (model_type in CFG$model_types) {
    for (R0_val in CFG$R0_values) {
      
      cat(sprintf("    Panel %d: %s | R0 = %s ... ", panel_idx, model_type, R0_val))
      
      data <- prepare_panel_data(model_type, R0_val, transition)
      
      # Show legend only on the very first (top-left) panel
      show_leg <- (panel_idx == 1)
      
      panels[[panel_idx]] <- build_panel(data, R0_val, model_type, show_leg)
      panel_idx <- panel_idx + 1
      cat("done\n")
    }
  }
  
  # Combine into subplot: 2 rows × 4 columns
  fig <- do.call(
    subplot,
    c(panels, list(nrows = 2, margin = 0.055, shareX = FALSE, shareY = FALSE))
  )
  
  # Column headers for R0 values
  col_positions <- c(0.115, 0.365, 0.615, 0.865)
  col_annotations <- lapply(seq_along(CFG$R0_values), function(i) {
    list(
      x         = col_positions[i],
      y         = 1.03,
      xref      = "paper", yref = "paper",
      text      = sprintf("<b>R\u2080 = %s</b>", CFG$R0_values[i]),
      showarrow = FALSE, xanchor = "center",
      font      = list(size = 14, color = "#2C3E50")
    )
  })
  
  # Row labels for model types
  row_positions   <- c(0.78, 0.22)
  row_annotations <- lapply(seq_along(CFG$model_types), function(i) {
    list(
      x         = -0.04,
      y         = row_positions[i],
      xref      = "paper", yref = "paper",
      text      = sprintf("<b>%s</b>", toupper(CFG$model_types[i])),
      showarrow = FALSE, xanchor = "right", textangle = -90,
      font      = list(size = 13, color = "#555")
    )
  })
  
  all_annotations <- c(col_annotations, row_annotations)
  
  fig <- fig %>%
    layout(
      title = list(
        text = paste0(
          "<b style='font-size:22px'>RtEstim vs ABM Ground Truth</b><br>",
          "<span style='font-size:15px; color:#666'>",
          trans_label,
          " | Shaded bands = 95% simulation interval",
          "</span>"
        ),
        x = 0.5, xanchor = "center", y = 0.99
      ),
      showlegend = TRUE,
      legend = list(
        title       = list(text = "<b>Method</b>", font = list(size = 12)),
        orientation = "v",
        x = 1.01, y = 0.5,
        bgcolor     = "rgba(255,255,255,0.97)",
        bordercolor = "#BDC3C7", borderwidth = 1.5,
        font        = list(size = 13),
        tracegroupgap = 12
      ),
      annotations   = all_annotations,
      paper_bgcolor = "white",
      margin        = list(t = 110, r = 160, b = 60, l = 80)
    ) %>%
    config(
      displayModeBar = TRUE,
      displaylogo    = FALSE,
      modeBarButtonsToRemove = c("select2d", "lasso2d"),
      toImageButtonOptions = list(
        format   = "png",
        filename = paste0("rtestim_vs_abm_", transition),
        width    = 1800, height = 900, scale = 2
      )
    )
  
  fig
}


# ── 6. INDEX PAGE ─────────────────────────────────────────────────────────────
#
#  A minimal landing page that links to the two HTML grids.

write_index <- function(files_written) {
  
  cards <- paste(sapply(files_written, function(f) {
    fname  <- basename(f)
    label  <- if (grepl("susceptible", fname)) "S → E Transition" else "E → I Transition"
    sub    <- if (grepl("susceptible", fname)) {
      "Susceptible to Exposed incidence"
    } else {
      "Exposed to Infected incidence"
    }
    icon <- if (grepl("susceptible", fname)) "🟠" else "🟢"
    sprintf('
    <a href="%s" class="card" target="_blank">
      <div class="icon">%s</div>
      <h3>%s</h3>
      <p>%s</p>
      <span class="tag">2 models × 4 R₀ values</span>
    </a>', fname, icon, label, sub)
  }), collapse = "\n")
  
  html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>RtEstim vs ABM — Results</title>
  <style>
    * { margin: 0; padding: 0; box-sizing: border-box; }
    body {
      font-family: "Georgia", serif;
      background: #F0F2F5;
      min-height: 100vh;
      display: flex; align-items: center; justify-content: center;
      padding: 40px 20px;
    }
    .container {
      max-width: 860px; width: 100%%;
      background: white;
      border-radius: 16px;
      box-shadow: 0 8px 40px rgba(0,0,0,0.10);
      overflow: hidden;
    }
    .header {
      background: #2C3E50;
      color: white;
      padding: 48px 48px 36px;
    }
    .header h1 { font-size: 28px; font-weight: normal; letter-spacing: -0.5px; }
    .header p  { margin-top: 10px; font-size: 15px; color: #A0AEC0; font-family: monospace; }
    .meta {
      padding: 24px 48px;
      background: #FAFBFC;
      border-bottom: 1px solid #E8ECF0;
      font-size: 14px; color: #555;
      display: flex; gap: 32px; flex-wrap: wrap;
    }
    .meta span b { color: #2C3E50; }
    .legend {
      padding: 20px 48px;
      background: #FAFBFC;
      border-bottom: 1px solid #E8ECF0;
      display: flex; gap: 24px; align-items: center; flex-wrap: wrap;
    }
    .dot { width: 14px; height: 14px; border-radius: 50%; display: inline-block; margin-right: 6px; }
    .legend-item { display: flex; align-items: center; font-size: 13px; color: #333; }
    .cards {
      padding: 40px 48px;
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 20px;
    }
    .card {
      display: block; text-decoration: none;
      border: 1.5px solid #E0E4EA;
      border-radius: 12px; padding: 28px 24px;
      transition: all 0.2s ease;
      color: #2C3E50;
    }
    .card:hover {
      border-color: #27AE60;
      box-shadow: 0 4px 20px rgba(39,174,96,0.15);
      transform: translateY(-2px);
    }
    .card .icon  { font-size: 32px; margin-bottom: 12px; }
    .card h3     { font-size: 17px; margin-bottom: 6px; }
    .card p      { font-size: 13px; color: #777; margin-bottom: 12px; }
    .card .tag   {
      display: inline-block; padding: 3px 10px;
      background: #EAF7F0; color: #27AE60;
      border-radius: 20px; font-size: 12px; font-weight: bold;
    }
    .footer {
      padding: 20px 48px;
      border-top: 1px solid #E8ECF0;
      font-size: 12px; color: #AAA; text-align: right;
    }
  </style>
</head>
<body>
<div class="container">
  <div class="header">
    <h1>RtEstim vs ABM Ground Truth</h1>
    <p>01_rtestim_estimation.R → 02_rtestim_plots.R</p>
  </div>
  <div class="meta">
    <span><b>Models:</b> Full, Partial</span>
    <span><b>R₀:</b> 1.5, 2.0, 3.0, 5.0</span>
    <span><b>Simulations per scenario:</b> 100</span>
    <span><b>Methods:</b> ABM ground truth, RtEstim</span>
  </div>
  <div class="legend">
    <div class="legend-item">
      <span class="dot" style="background:#2C3E50"></span> ABM (True Rₜ)
    </div>
    <div class="legend-item">
      <span class="dot" style="background:#27AE60"></span> RtEstim
    </div>
    <div class="legend-item" style="color:#999; font-size:12px;">
      Shaded bands = 95%% simulation interval across 100 runs
    </div>
  </div>
  <div class="cards">
    %s
  </div>
  <div class="footer">Generated April 2026</div>
</div>
</body>
</html>', cards)
  
  index_path <- file.path(CFG$out_dir, "index.html")
  writeLines(html, index_path)
  cat("  ✓ Index page →", index_path, "\n")
}


# ── 7. MAIN ───────────────────────────────────────────────────────────────────

cat("\n", strrep("═", 62), "\n", sep = "")
cat("  RtEstim Plots — Ground Truth vs RtEstim\n")
cat("  Reading from : ", CFG$in_dir_rt,  "/\n", sep = "")
cat("  Writing to   : ", CFG$out_dir, "/\n",    sep = "")
cat(strrep("═", 62), "\n", sep = "")

files_written <- c()

for (transition in CFG$transitions) {
  
  cat("\n", strrep("─", 62), "\n", sep = "")
  cat("Transition:", transition, "\n")
  
  fig <- build_grid(transition)
  
  out_file <- file.path(
    CFG$out_dir,
    paste0("grid_", transition, ".html")
  )
  
  saveWidget(fig, out_file, selfcontained = TRUE)
  files_written <- c(files_written, out_file)
  cat("  ✓ Saved →", out_file, "\n")
}

write_index(files_written)

cat("\n", strrep("═", 62), "\n", sep = "")
cat("  Done! Open rt_plots/index.html in your browser.\n")
cat(strrep("═", 62), "\n\n", sep = "")