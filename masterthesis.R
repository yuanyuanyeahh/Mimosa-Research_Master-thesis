# =========================
# Mimosa sjSDM — Full script 
# =========================

# 1) Packages
suppressPackageStartupMessages({
  library(sjSDM)
  library(tidyverse)
  library(forcats)
  library(ggplot2)
  library(stringr)
  library(patchwork)  
})

# 2) File paths
Y_FILE    <- "/Users/gabriel/Desktop/master/defense/processed data/y_matrix_now.csv"
X_FILE    <- "/Users/gabriel/Desktop/master/defense/processed data/X_matrix_now.csv"
NB_FILE   <- "/Users/gabriel/Desktop/master/defense/processed data/neighbor_matrix_now.csv"
SITE_FILE <- "/Users/gabriel/Desktop/master/defense/processed data/site_level.csv"

# 3) Helpers
impute_scale <- function(M){
  M <- as.matrix(M)
  keep <- colSums(!is.na(M)) > 0
  M <- M[, keep, drop = FALSE]
  for(j in seq_len(ncol(M))){
    v <- M[, j]
    med <- suppressWarnings(stats::median(v, na.rm = TRUE))
    if(!is.finite(med)) med <- 0
    v[!is.finite(v)] <- NA
    v[is.na(v)] <- med
    M[, j] <- v
  }
  vvar <- apply(M, 2, stats::var, na.rm = TRUE); vvar[is.na(vvar)] <- 0
  M <- M[, vvar > 0, drop = FALSE]
  scale(M)
}
norm_plotid <- function(x){
  id <- stringr::str_extract(x, "(?i)plot\\s*\\d+")
  tolower(gsub("\\s+","", id))
}

# 4) Read data
y_raw    <- readr::read_csv(Y_FILE,    show_col_types = FALSE)
x_raw    <- readr::read_csv(X_FILE,    show_col_types = FALSE)
nb_raw   <- readr::read_csv(NB_FILE,   show_col_types = FALSE)
site_raw <- readr::read_csv(SITE_FILE, show_col_types = FALSE)

# 5) Response & coords
Y_mat <- y_raw |>
  rename(Mimosa = pre) |>
  mutate(Mimosa = as.integer(Mimosa %in% c(1,"1","yes","present",TRUE))) |>
  filter(Mimosa_ID %in% x_raw$Mimosa_ID) |>
  column_to_rownames("Mimosa_ID") |>
  as.matrix()

coords <- x_raw |>
  select(Mimosa_ID, Longitude, Latitude) |>
  column_to_rownames("Mimosa_ID") |>
  as.matrix()

# 6) Circle-level factors (A/B circles): light / road / water — scaled
X_circle <- x_raw |>
  transmute(
    Mimosa_ID,
    Circle_Light             = suppressWarnings(readr::parse_number(as.character(Light))),
    Circle_Distance_to_Road  = suppressWarnings(readr::parse_number(as.character(Distance_to_Road))),
    Circle_Distance_to_Water = suppressWarnings(readr::parse_number(as.character(Distance_to_Water)))
  ) |>
  column_to_rownames("Mimosa_ID") |>
  impute_scale()
circle_names <- colnames(X_circle)

# 7) Site-level factors (for FULL)
x_aug  <- x_raw    |> mutate(PlotID = norm_plotid(Mimosa_ID))
site2  <- site_raw |> mutate(PlotID = norm_plotid(site_id))

X_site_full <- x_aug |>
  select(Mimosa_ID, PlotID) |>
  left_join(site2, by = "PlotID") |>
  transmute(
    Mimosa_ID,
    Site_Distance_to_Road   = suppressWarnings(readr::parse_number(as.character(Distance_to_Road))),
    Site_Distance_to_Water  = suppressWarnings(readr::parse_number(as.character(Distance_to_Water))),
    elevation_m             = suppressWarnings(readr::parse_number(as.character(elevation_m))),
    slope_deg               = suppressWarnings(readr::parse_number(as.character(slope_deg))),
    canopy_cover_percent    = suppressWarnings(readr::parse_number(as.character(canopy_cover_percent))),
    soil_moisture_percent   = suppressWarnings(readr::parse_number(as.character(soil_moisture_percent))),
    temperature_avg_c       = suppressWarnings(readr::parse_number(as.character(temperature_avg_c))),
    precipitation_mm        = suppressWarnings(readr::parse_number(as.character(precipitation_mm)))
  ) |>
  column_to_rownames("Mimosa_ID") |>
  impute_scale()
site_names <- colnames(X_site_full)

# 8) Neighbour factors (wide 0/1 -> scale)
nb_clean <- nb_raw |>
  mutate(Neighbor_Species = trimws(Neighbor_Species)) |>
  filter(!is.na(Neighbor_Species),
         Neighbor_Species != "",
         !Neighbor_Species %in% c("NA","N.A.","NA.","na","NaN",
                                  "Unknown","Unidentified","Unidentified_sp")) |>
  mutate(value = 1L) |>
  distinct(Mimosa_ID, Neighbor_Species, .keep_all = TRUE)

X_nb <- tidyr::pivot_wider(
  nb_clean,
  id_cols     = Mimosa_ID,
  names_from  = Neighbor_Species,
  values_from = value,
  values_fill = 0L,
  values_fn   = max,
  names_repair= "unique"
) |>
  column_to_rownames("Mimosa_ID") |>
  as.matrix()
X_nb <- X_nb[, colSums(X_nb != 0) > 0, drop = FALSE]
X_nb <- impute_scale(X_nb)
nb_names <- colnames(X_nb)

# 9) Align & FULL design
ids <- Reduce(intersect, list(rownames(Y_mat), rownames(coords),
                              rownames(X_circle), rownames(X_site_full), rownames(X_nb)))
Y_use   <- Y_mat    [ids,, drop = FALSE]
SP_use  <- coords   [ids,, drop = FALSE]
X_circle<- X_circle [ids,, drop = FALSE]
X_siteF <- X_site_full[ids,, drop = FALSE]
X_nb    <- X_nb     [ids,, drop = FALSE]
X_all   <- cbind(X_circle, X_siteF, X_nb)

# 10) Fit FULL sjSDM
set.seed(42)
mod_full <- sjSDM(
  Y       = Y_use,
  env     = linear(X_all,  ~ ., lambda = 1e-3),
  spatial = linear(SP_use, ~ Longitude + Latitude),
  iter    = 2000,
  family  = binomial("probit"),
  device  = "cpu"
)

# 11) Extract FULL betas (env)
get_beta_matrix <- function(fit, Xref){
  b <- try(coef(fit, "env")$env[[1]], silent = TRUE)
  if (inherits(b, "try-error") || is.null(b)) b <- try(coef(fit)$env[[1]], silent = TRUE)
  if (is.null(b)) stop("Cannot extract env coefficients.")
  if (nrow(b) == 1) b <- t(b)
  rownames(b) <- colnames(Xref)[seq_len(nrow(b))]
  b
}
B_full   <- get_beta_matrix(mod_full, X_all)
beta_all <- as.numeric(B_full[, 1]); names(beta_all) <- rownames(B_full)

# 12) Site-only sjSDM（site + spatial）
site_cols_needed <- c("Distance_to_Water","Distance_to_Road",
                      "elevation_m","slope_deg","canopy_cover_percent",
                      "soil_moisture_percent","temperature_avg_c","precipitation_mm")

site_exp_only <- x_aug |>
  select(Mimosa_ID, PlotID) |>
  left_join(
    site_raw |>
      mutate(PlotID = norm_plotid(site_id)) |>
      select(PlotID, any_of(site_cols_needed)),
    by = "PlotID"
  ) |>
  rename(Site_Distance_to_Water = Distance_to_Water,
         Site_Distance_to_Road  = Distance_to_Road) |>
  mutate(across(-Mimosa_ID, ~ suppressWarnings(readr::parse_number(as.character(.))))) |>
  distinct(Mimosa_ID, .keep_all = TRUE) |>
  column_to_rownames("Mimosa_ID") |>
  as.matrix()

ids_site <- intersect(rownames(site_exp_only), rownames(coords))

Y_site <- y_raw |>
  transmute(Mimosa_ID, Mimosa = as.integer(pre %in% c(1,"1","yes","present",TRUE))) |>
  tidyr::complete(Mimosa_ID = ids_site, fill = list(Mimosa = 0L)) |>
  distinct(Mimosa_ID, .keep_all = TRUE) |>
  column_to_rownames("Mimosa_ID") |>
  as.matrix()

X_site_std <- site_exp_only[ids_site, , drop = FALSE]
for(j in seq_len(ncol(X_site_std))){
  v <- X_site_std[, j]
  med <- suppressWarnings(stats::median(v, na.rm = TRUE))
  if(!is.finite(med)) med <- 0
  v[!is.finite(v)] <- NA
  v[is.na(v)] <- med
  X_site_std[, j] <- v
}
vv <- apply(X_site_std, 2, stats::var, na.rm = TRUE); vv[is.na(vv)] <- 0
X_site_std <- X_site_std[, vv > 0, drop = FALSE]
X_site_std <- impute_scale(X_site_std)

SP_site <- coords[ids_site, , drop = FALSE]

set.seed(42)
mod_site_only <- sjSDM(
  Y       = Y_site,
  env     = linear(X_site_std, ~ .),
  spatial = linear(SP_site, ~ Longitude + Latitude),
  iter    = 1000,
  family  = binomial("probit"),
  device  = "cpu"
)

# 12b) robust env-beta extractor for Site-only
get_env_betas_aligned <- function(fit, Xref){
  xnames <- colnames(Xref)
  co <- try(suppressWarnings(sjSDM::coefficients(fit, "env")), silent = TRUE)
  if (!inherits(co, "try-error") && !is.null(co) && !is.null(co$env)) {
    b <- co$env; if (is.list(b)) b <- b[[1]]
    b <- as.matrix(b); if (nrow(b) == 1) b <- t(b)
    rn <- rownames(b)
    if (!is.null(rn) && any(grepl("^(intercept|bias|\\(Intercept\\))$", rn, ignore.case = TRUE))) {
      b <- b[!grepl("^(intercept|bias|\\(Intercept\\))$", rn, ignore.case = TRUE), , drop = FALSE]
      rn <- rownames(b)
    }
    if (nrow(b) == length(xnames) + 1) { b <- b[-1, , drop = FALSE]; rn <- rownames(b) }
    est <- setNames(rep(NA_real_, length(xnames)), xnames)
    if (!is.null(rn)) {
      hit <- intersect(rn, xnames); est[hit] <- as.numeric(b[hit, 1])
    } else {
      k <- min(nrow(b), length(xnames)); est[seq_len(k)] <- as.numeric(b[seq_len(k), 1])
    }
    return(est)
  }
  # fallback: parse summary
  sm <- utils::capture.output(summary(fit))
  hdr <- grep("^\\s*Coefficients\\s*\\(beta\\)\\s*:\\s*$", sm)
  if (length(hdr) > 0) {
    i <- hdr[length(hdr)] + 1L
    while (i <= length(sm) && grepl("^\\s*$", sm[i])) i <- i + 1L
    rows <- character()
    while (i <= length(sm) && nzchar(trimws(sm[i]))) { rows <- c(rows, sm[i]); i <- i + 1L }
    if (length(rows) > 0) {
      parse_line <- function(ln){
        ln2 <- gsub("\\t", " ", gsub(" +", " ", trimws(ln)))
        toks <- strsplit(ln2, " +")[[1]]
        beta <- suppressWarnings(as.numeric(toks[length(toks)]))
        name <- paste(toks[-length(toks)], collapse = " "); tibble(name = name, beta = beta)
      }
      df <- dplyr::bind_rows(lapply(rows, parse_line))
      est <- setNames(rep(NA_real_, length(xnames)), xnames)
      nm  <- intersect(df$name, xnames)
      if (length(nm)) est[nm] <- df$beta[match(nm, df$name)]
      return(est)
    }
  }
  stop("Cannot extract env coefficients.")
}
b_site_vec <- get_env_betas_aligned(mod_site_only, X_site_std)

# ---- Plotting utilities (return ggplot + optional save) ----
plot_block_beta <- function(names_vec, beta_vec, title_txt, file_png = NULL,
                            top_n = 30, always_include = character(),
                            wrap_width = 28, base_height = 3.5, per_row = 0.32,
                            y_text_size = 10, width_in = 7.2, save_plot = TRUE){
  avail  <- intersect(names_vec, names(beta_vec))
  forced <- intersect(always_include, names(beta_vec))
  pick0  <- tibble(Factor = avail, Beta = beta_vec[avail]) |>
    filter(is.finite(Beta)) |>
    arrange(desc(abs(Beta)))
  head_keep <- head(pick0$Factor, top_n)
  keep <- unique(c(forced, head_keep))
  df <- pick0 |> filter(Factor %in% keep)
  if(nrow(df) == 0){
    warning("No factors to plot for: ", title_txt)
    return(invisible(NULL))
  }
  df <- df |>
    mutate(PlotLabel = stringr::str_wrap(Factor, width = wrap_width)) |>
    mutate(PlotLabel = forcats::fct_reorder(PlotLabel, abs(Beta)))
  lim <- max(abs(df$Beta)); if(!is.finite(lim) || lim == 0) lim <- 1
  p <- ggplot(df, aes(PlotLabel, Beta)) +
    geom_col(width = 0.6, fill = "grey30") +
    coord_flip() +
    scale_y_continuous(limits = c(-1.15*lim, 1.15*lim),
                       expand = expansion(mult = 0.02)) +
    labs(title = title_txt, x = NULL, y = "Standardized \u03B2 (probit)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = y_text_size))
  if (isTRUE(save_plot) && !is.null(file_png)) {
    h <- base_height + per_row * nrow(df)
    ggsave(file_png, p, width = width_in, height = h, dpi = 300)
  }
  return(p)
}

plot_neighbour_split <- function(beta_vec, file_png = NULL,
                                 n_each = 10, wrap_width = 18, y_text = 9,
                                 width_in = 10.5, save_plot = TRUE) {
  df <- tibble(Species = names(beta_vec), Beta = as.numeric(beta_vec)) |>
    filter(is.finite(Beta))
  top_pos <- df |> arrange(desc(Beta)) |> slice_head(n = n_each) |> mutate(Sign = "Positive")
  top_neg <- df |> arrange(Beta)       |> slice_head(n = n_each) |> mutate(Sign = "Negative")
  dd <- bind_rows(top_neg, top_pos) |>
    mutate(Label = stringr::str_wrap(Species, width = wrap_width)) |>
    group_by(Sign) |>
    mutate(Label = forcats::fct_reorder(Label, abs(Beta))) |>
    ungroup()
  n_max <- max(table(dd$Sign))
  h <- 4 + 0.34 * n_max
  p <- ggplot(dd, aes(Label, Beta)) +
    geom_col(width = 0.6, fill = "grey30") +
    coord_flip() +
    facet_grid(Sign ~ ., scales = "free_y", switch = "y") +
    labs(title = "(C) Effect of neighbour factors on Mimosa occurrence",
         x = NULL, y = "Standardized \u03B2 (probit)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = y_text),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold"))
  if (isTRUE(save_plot) && !is.null(file_png)) {
    ggsave(file_png, p, width = width_in, height = h, dpi = 300)
  }
  return(p)
}

# 14) Make individual figures
# Circle A/B factors — light/road/water
p_circle <- plot_block_beta(
  names_vec = circle_names, beta_vec = beta_all,
  title_txt = "(A) Effects of circle-level factors (circles A and B) on Mimosa occurrence",
  file_png  = "beta_circlesAB.png",
  top_n = 30,
  always_include = c("Circle_Distance_to_Road","Circle_Distance_to_Water","Circle_Light"),
  wrap_width = 28, y_text_size = 10, width_in = 7.2, save_plot = TRUE
)

# Site-level factors — include road & water etc.
site_keep <- c("Site_Distance_to_Road","Site_Distance_to_Water",
               "elevation_m","slope_deg","canopy_cover_percent",
               "soil_moisture_percent","temperature_avg_c","precipitation_mm")

p_site <- plot_block_beta(
  names_vec = site_keep, beta_vec = b_site_vec,
  title_txt = "(B) Effects of site-level factors on Mimosa occurrence",
  file_png  = "beta_site.png",
  top_n = length(site_keep),
  always_include = c("Site_Distance_to_Road"),
  wrap_width = 28, y_text_size = 10, width_in = 7.2, save_plot = TRUE
)

# Neighbour factors — split by sign
p_nb <- plot_neighbour_split(
  beta_vec  = beta_all[nb_names],
  file_png  = "beta_neighbour.png",
  n_each    = 5,
  wrap_width= 18,
  y_text    = 9,
  width_in  = 10.5,
  save_plot = TRUE
)

# 1) Circle-level factors (A/B circles)
p_circle <- plot_block_beta(
  names_vec = circle_names, beta_vec = beta_all,
  title_txt = "(A) Effects of circle-level factors (circles A and B) on Mimosa occurrence",
  file_png  = "beta_circlesAB.png",
  top_n = 30,
  always_include = intersect(c("Circle_Distance_to_Road","Circle_Distance_to_Water","Circle_Light"), circle_names),
  wrap_width = 28, y_text_size = 10, width_in = 7.2, save_plot = TRUE
)
print(p_circle)

# 2) Site-level factors（保持原样，不修改）
p_site <- plot_block_beta(
  names_vec = site_keep, beta_vec = b_site_vec,
  title_txt = "(B) Effects of site-level factors on Mimosa occurrence",
  file_png  = "beta_site.png",
  top_n = length(site_keep),
  always_include = c("Site_Distance_to_Road"),
  wrap_width = 28, y_text_size = 10, width_in = 7.2, save_plot = TRUE
)
print(p_site)

# 3) Neighbour factors
p_nb <- plot_neighbour_split(
  beta_vec  = beta_all[nb_names],
  file_png  = "beta_neighbour.png",
  n_each    = 5, wrap_width = 18, y_text = 9, width_in = 10.5, save_plot = TRUE
)
print(p_nb)

cat("Saved files in:", getwd(), "\n",
    "- beta_circlesAB.png\n",
    "- beta_site.png\n",
    "- beta_neighbour.png\n")

