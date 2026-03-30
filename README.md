# Bike Lanes & Housing Costs in Chicago

An analysis of the relationship between bike lane infrastructure
and median rent at the census tract level in Cook County, IL.

Built as part of the Data Visualization Society mentorship program.

---

## Research Question

Is bike lane infrastructure associated with housing costs
(specifically median rent) at the neighborhood level?

---

## Methodology

- Census tract-level panel across three years: 2012, 2016, 2022
- ACS 5-year estimates pulled via `tidycensus`
- Bike lane shapefiles from the Chicago Data Portal
- Primary model: long-difference regression (Δrent ~ Δbike_km + controls)
- Cross-sectional 2022 model included for comparison

Loosely modeled after Woolley (2018), "The Effect of Bike Lane
Infrastructure on Urban Housing Markets," adapted for Chicago
and a temporal identification strategy.

---

## Repository Structure
```
├── bike_lanes_housing_panel.R   # Main analysis script
├── .gitignore
└── README.md
```

---

## How to Reproduce

1. Clone this repo
2. Get a free Census API key at https://api.census.gov/data/key_signup.html
3. In R, run:
```r
   tidycensus::census_api_key("YOUR_KEY_HERE", install = TRUE)
```
4. Download bike lane shapefiles for 2012, 2016, and 2022
   from the Chicago Data Portal and update the `BIKE_SHAPEFILES`
   paths in the config section of the script
5. Run `bike_lanes_housing_panel.R` top to bottom

---

## Data Sources

- **ACS 5-year estimates**: U.S. Census Bureau via `tidycensus`
  - Years: 2012, 2016, 2022
  - Geography: Census tracts, Cook County IL
  - Variables: median rent, median income, population,
    housing units, race/ethnicity, poverty
- **Bike lane shapefiles**: Chicago Data Portal
  - https://data.cityofchicago.org

---

## Key Packages

- `tidycensus` — Census data
- `sf` — spatial operations
- `tidyverse` — data wrangling
- `fixest` — panel regression with fixed effects
- `spdep` — spatial autocorrelation diagnostics
- `modelsummary` — regression output tables

---

## Status

Work in progress. Visualization dashboard coming soon.
