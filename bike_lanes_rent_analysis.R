# =========================================================
# Bike Lanes + Housing Project
# Multi-year panel analysis — Cook County census tracts
# Years: 2012, 2016, 2022
#
# Strategy: tract-level pseudo-panel using multiple bike
# lane shapefiles. Compute long-difference (2012→2022)
# and run a delta regression: Δrent ~ Δbike_km + controls.
#
# Based on: Woolley (2018), "The Effect of Bike Lane
# Infrastructure on Urban Housing Markets"
# =========================================================


# ----------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------
# install.packages(c("tidycensus", "tidyverse", "janitor",
#                    "sf", "broom", "modelsummary", "fixest",
#                    "spdep"))

library(tidycensus)
library(tidyverse)
library(janitor)
library(sf)
library(broom)
library(modelsummary)
library(fixest)    # for panel models with fixed effects
library(spdep)     # for spatial autocorrelation diagnostics


# ----------------------------------------------------------
# 1. Config
# ----------------------------------------------------------

# Years
YEARS <- c(2012, 2016, 2022)

# Paths to bike lane shapefiles, one per year.
BIKE_SHAPEFILES <- list(
  "2016" = "~/Documents/DataSciencefiles/Data Viz Society Mentorship/Bike Routes (Deprecated February 2020)_20260329/geo_export_447a3ddf-10f3-4129-9497-31d73756732b.shp",
  "2020" = "~/Documents/DataSciencefiles/Data Viz Society Mentorship/CDOT_Bikeways_2016_0311_20260329/geo_export_290d4563-a7ab-4263-b749-7be9ba5b259f.shp",
  "2022" = "~/Documents/DataSciencefiles/Data Viz Society Mentorship/Chicago Bike Routes/geo_export_1573f37d-b117-47a5-88c9-4a41800fc421.shp"
)

# Projected CRS for Illinois — used for area/length calculations
# EPSG 26971 = NAD83 / Illinois East (meters)
IL_CRS <- 26971

# census_api_key("YOUR_KEY_HERE", install = TRUE)  # run once


# ----------------------------------------------------------
# 2. Pull ACS data for all years
# ----------------------------------------------------------
# NOTE: ACS 5-year estimates are labeled by their end year.
# "2012" = 2008–2012 pooled; "2022" = 2018–2022 pooled.
# This is expected and fine for a long-difference design.

acs_vars <- c(
  median_rent        = "B25064_001",   # Median gross rent
  median_income      = "B19013_001",   # Median household income
  total_population   = "B01003_001",   # Total population
  total_housing_units = "B25001_001",  # Total housing units
  white_nonhisp      = "B03002_003",   # White non-Hispanic
  poverty_denom      = "B17001_001",   # Poverty universe
  poverty_num        = "B17001_002"    # Below poverty line
)

get_cook_tract_acs <- function(year) {
  get_acs(
    geography = "tract",
    variables = acs_vars,
    state     = "IL",
    county    = "Cook",
    year      = year,
    survey    = "acs5",
    geometry  = FALSE
  ) %>%
    clean_names() %>%
    mutate(year = year)
}

acs_long <- map_dfr(YEARS, get_cook_tract_acs)

acs_wide <- acs_long %>%
  select(geoid, name, year, variable, estimate) %>%
  pivot_wider(names_from = variable, values_from = estimate)

acs_panel <- acs_wide %>%
  mutate(
    pct_white_nonhisp = if_else(total_population > 0,
                                white_nonhisp / total_population, NA_real_),
    pct_poverty       = if_else(poverty_denom > 0,
                                poverty_num / poverty_denom, NA_real_)
  ) %>%
  select(geoid, name, year, median_rent, median_income,
         total_population, total_housing_units,
         pct_white_nonhisp, pct_poverty) %>%
  arrange(geoid, year)

write_csv(acs_panel, "acs_panel_cook_tracts.csv")
cat("ACS panel rows:", nrow(acs_panel), "\n")


# ----------------------------------------------------------
# 3. Pull Cook County tract geometry (2022 boundaries)
# ----------------------------------------------------------
# We use a single geometry vintage throughout. Tract boundaries
# do shift between decennial censuses but are stable within a
# 5-year ACS window — this is an accepted trade-off for a
# tract-level pseudo-panel.

cook_tracts_sf <- get_acs(
  geography = "tract",
  variables = "B01003_001",
  state     = "IL",
  county    = "Cook",
  year      = 2022,
  survey    = "acs5",
  geometry  = TRUE
) %>%
  clean_names() %>%
  select(geoid, name, geometry) %>%
  st_transform(IL_CRS) %>%
  st_make_valid()

st_write(cook_tracts_sf, "cook_tracts_2022.gpkg", delete_dsn = TRUE)
cat("Tract count:", nrow(cook_tracts_sf), "\n")


# ----------------------------------------------------------
# 4. Compute tract-level bike lane length for each year
# ----------------------------------------------------------
# We aggregate total lane km per tract without trying to
# match individual segments across years (Woolley used
# installation-date data; we don't have that for Chicago).
# The tract-level aggregate is the unit of analysis throughout.

compute_bike_length <- function(shp_path, year_label, tracts_sf) {
  
  cat("Processing bike lanes for", year_label, "...\n")
  
  bike_sf <- st_read(shp_path, quiet = TRUE) %>%
    st_transform(IL_CRS) %>%
    st_make_valid()
  
  intersected <- st_intersection(bike_sf, tracts_sf)
  
  lengths_df <- intersected %>%
    mutate(seg_length_m = as.numeric(st_length(geometry))) %>%
    st_drop_geometry() %>%
    group_by(geoid) %>%
    summarise(total_bike_length_m = sum(seg_length_m, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(year = as.integer(year_label))
  
  lengths_df
}

bike_lengths_list <- imap(
  BIKE_SHAPEFILES,
  ~ compute_bike_length(.x, .y, cook_tracts_sf)
)

bike_lengths_raw <- bind_rows(bike_lengths_list)

# Fill zeros: tracts absent from intersection had no bike lanes
# NOTE: very small values (< 1m) are likely GIS edge artefacts
# from a segment clipping slightly into a neighbouring tract.
# We zero these out to avoid noise — consistent with Woolley (2018).
all_geoids  <- tibble(geoid = cook_tracts_sf$geoid)
all_years   <- tibble(year  = as.integer(YEARS))

bike_lengths_complete <- crossing(all_geoids, all_years) %>%
  left_join(bike_lengths_raw, by = c("geoid", "year")) %>%
  replace_na(list(total_bike_length_m = 0)) %>%
  mutate(
    total_bike_length_m = if_else(total_bike_length_m < 1, 0,
                                  total_bike_length_m),
    bike_length_km = total_bike_length_m / 1000
  )

write_csv(bike_lengths_complete, "bike_lengths_by_tract_year.csv")
cat("Bike length rows:", nrow(bike_lengths_complete), "\n")


# ----------------------------------------------------------
# 5. Build the full stacked panel
# ----------------------------------------------------------

tract_area <- cook_tracts_sf %>%
  mutate(tract_area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(geoid, tract_area_km2)

panel <- acs_panel %>%
  left_join(bike_lengths_complete, by = c("geoid", "year")) %>%
  left_join(tract_area,            by = "geoid") %>%
  mutate(pop_density = total_population / tract_area_km2) %>%
  filter(
    total_population > 100,   # drop near-empty tracts
    !is.na(median_rent),
    !is.na(median_income)
  )

write_csv(panel, "panel_cook_tracts.csv")
cat("Panel rows:", nrow(panel), "\n")
glimpse(panel)


# ----------------------------------------------------------
# 6. Cross-sectional model (2022 only) — replicates original
# ----------------------------------------------------------
# Kept here for direct comparison with the earlier analysis.
# R² ≈ 0.51 from the original run is the baseline to beat.

panel_2022 <- panel %>% filter(year == 2022)

model_cross <- lm(
  median_rent ~ bike_length_km + median_income + pct_white_nonhisp +
    pop_density + total_housing_units,
  data = panel_2022
)

summary(model_cross)
tidy(model_cross) %>% write_csv("model_cross_sectional_2022.csv")


# ----------------------------------------------------------
# 7. Long-difference dataset: 2012 → 2022
# ----------------------------------------------------------
# This is the primary identification strategy: within-tract
# changes absorb all time-invariant confounders (geography,
# transit access, stable neighborhood character).
#
#
# Also keep baseline 2012 levels as controls for
# mean-reversion: high-rent tracts in 2012 may have
# less room to grow (Woolley notes similar ceiling effects).

delta_df <- panel %>%
  filter(year %in% c(2012, 2022)) %>%
  arrange(geoid, year) %>%
  group_by(geoid) %>%
  filter(n() == 2) %>%   # keep only tracts observed in both years
  summarise(
    # Changes (outcome and key predictor)
    d_rent          = median_rent[year == 2022]      - median_rent[year == 2012],
    d_bike_km       = bike_length_km[year == 2022]   - bike_length_km[year == 2012],
    # Changes in controls
    d_income        = median_income[year == 2022]    - median_income[year == 2012],
    d_pct_white     = pct_white_nonhisp[year == 2022]- pct_white_nonhisp[year == 2012],
    d_pct_poverty   = pct_poverty[year == 2022]      - pct_poverty[year == 2012],
    d_pop_density   = pop_density[year == 2022]      - pop_density[year == 2012],
    # Baseline levels (control for mean-reversion)
    base_rent       = median_rent[year == 2012],
    base_bike_km    = bike_length_km[year == 2012],
    base_income     = median_income[year == 2012],
    .groups = "drop"
  )

write_csv(delta_df, "delta_dataset_2012_2022.csv")
cat("Long-difference rows:", nrow(delta_df), "\n")
glimpse(delta_df)

# Sanity checks on deltas
summary(delta_df$d_rent)
summary(delta_df$d_bike_km)
cat("Tracts that gained bike lanes:", sum(delta_df$d_bike_km > 0), "\n")
cat("Tracts that lost bike lanes:  ", sum(delta_df$d_bike_km < 0), "\n")
cat("Tracts with no change:        ", sum(delta_df$d_bike_km == 0), "\n")


# ----------------------------------------------------------
# 8. Delta regression (primary model)
# ----------------------------------------------------------
# Δrent ~ Δbike_km + Δcontrols + baseline_rent
#
# Interpretation of d_bike_km:
#   Positive + significant → bike lane additions associated
#     with rent growth (gentrification signal)
#   Negative + significant → additions associated with rent
#     suppression (possibly congestion/disruption, as Woolley
#     found in NYC)
#   Insignificant → consistent with cross-sectional null;
#     time-invariant neighbourhood traits drove the
#     cross-sectional association

model_delta <- lm(
  d_rent ~ d_bike_km + d_income + d_pct_white + d_pct_poverty +
    d_pop_density + base_rent,
  data = delta_df
)

summary(model_delta)
tidy(model_delta) %>% write_csv("model_delta_2012_2022.csv")


# ----------------------------------------------------------
# 9. Two-period pseudo-panel (2012→2016, 2016→2022)
# ----------------------------------------------------------
# Gives more observations but requires a period fixed effect.
# Use this as a robustness check, not the primary model.
# Clustered SEs by GEOID account for within-tract correlation.

delta_panel <- panel %>%
  arrange(geoid, year) %>%
  group_by(geoid) %>%
  mutate(
    d_rent        = median_rent       - lag(median_rent),
    d_bike_km     = bike_length_km    - lag(bike_length_km),
    d_income      = median_income     - lag(median_income),
    d_pct_white   = pct_white_nonhisp - lag(pct_white_nonhisp),
    d_pct_poverty = pct_poverty       - lag(pct_poverty),
    base_rent     = lag(median_rent),
    period        = paste0(lag(year), "_", year)
  ) %>%
  ungroup() %>%
  filter(!is.na(d_rent))

model_panel_fe <- feols(
  d_rent ~ d_bike_km + d_income + d_pct_white + d_pct_poverty + base_rent |
    period,                          # period fixed effect
  data    = delta_panel,
  cluster = ~geoid                   # cluster SEs within tract
)

summary(model_panel_fe)


# ----------------------------------------------------------
# 10. Spatial autocorrelation diagnostics
# ----------------------------------------------------------
# Adjacent tracts are not independent. Moran's I on residuals
# tells us whether spatial clustering is a problem.

# Build spatial weights from tract centroids
tract_centroids <- cook_tracts_sf %>%
  filter(geoid %in% delta_df$geoid) %>%
  arrange(geoid)

coords <- st_coordinates(st_centroid(tract_centroids))
nb     <- knn2nb(knearneigh(coords, k = 5))
w      <- nb2listw(nb, style = "W")

# Match residual order to weight matrix order
delta_for_moran <- delta_df %>%
  filter(geoid %in% tract_centroids$geoid) %>%
  arrange(geoid)

resids <- residuals(
  lm(d_rent ~ d_bike_km + d_income + d_pct_white + d_pct_poverty +
       d_pop_density + base_rent,
     data = delta_for_moran)
)

moran_result <- moran.test(resids, w)
print(moran_result)


# ----------------------------------------------------------
# 11. Model comparison table
# ----------------------------------------------------------

modelsummary(
  list(
    "Cross-sectional (2022)"  = model_cross,
    "Long-difference (Δ)"     = model_delta,
    "Panel FE (two periods)"  = model_panel_fe
  ),
  stars  = c("*" = 0.05, "**" = 0.01, "***" = 0.001),
  output = "model_comparison.html"
)

cat("\nDone. Output files written:\n")
cat("  acs_panel_cook_tracts.csv\n")
cat("  bike_lengths_by_tract_year.csv\n")
cat("  panel_cook_tracts.csv\n")
cat("  delta_dataset_2012_2022.csv\n")
cat("  model_cross_sectional_2022.csv\n")
cat("  model_delta_2012_2022.csv\n")
cat("  model_comparison.html\n")