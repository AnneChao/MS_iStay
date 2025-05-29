# ========================================================================================================== #

# Install and upload libraries ----

library(tidyverse)
library(lmerTest)

# ========================================================================================================== #
# 
# This code includes four parts:
# 
# (1) Figure 1. Gamma, alpha and synchrony profiles for q between 0 and 2 within each plot. 
# (2) Figure 2. Biodiversity–stability and biodiversity–synchrony relationships based on 76 plots.
# (3) Figure 3. Temporal effects of species richness on stability and synchrony based on 12 consecutive overlapping 10-year moving window.
# (4) Figure 4. Biodiversity–stability and biodiversity–synchrony relationships based on 20 sets
#
# See "Brief guide" for details. 
# 
# ========================================================================================================== #


library(devtools)
install_github("AnneChao/iStay")    # Press 'Enter' to skip update options
library(iStay)

data("Jena_plot_biomass_data")
data("Jena_species_biomass_data")

source("Source R code.txt")


# ========================================================================================================== #
# Figure 1. Gamma, alpha and synchrony profiles for q between 0 and 2 within each plot.
df1 <- list(B4A14_2 = Jena_species_biomass_data$B4A14_B4_2,
            B1A04_4 = Jena_species_biomass_data$B1A04_B1_4)
df1 <- Map(function(x, nm) {
  rownames(x) <- paste0(nm, rownames(x))
  x
}, df1, names(df1))

output_fig_1 <- Stay_Multiple(df1, order.q = seq(0.01, 2, 0.01))
output_fig_1b <- list(
  B1A04_4 = Stay_Single(df1$B1A04_4, order.q = seq(0.01, 2, 0.01)),
  B4A14_2 = Stay_Single(df1$B4A14_2, order.q = seq(0.01, 2, 0.01))
) |> 
  bind_rows(.id = "Site")

## Figure 1 (a)

fig_1a(output_fig_1)

## Figure 1 (b)

fig_1b(output_fig_1, output_fig_1b)

## Figure 1 (c)

fig_1c(output_fig_1)


# ========================================================================================================== #
# Figure 2. Biodiversity–stability and biodiversity–synchrony relationships based on 76 plots.

split_names2 <- str_split(names(Jena_species_biomass_data), "_", simplify = TRUE)

structure2 <- data.frame(
  block = split_names2[, 2],
  log2_sowndiv = log2(as.numeric(split_names2[, 3]))
)

structure2c <- structure2 |> filter(log2_sowndiv != 0)

output2 <- Stay_Multiple(Jena_species_biomass_data, order.q = c(0.5, 1, 2))

## Figure 2 (a)

output_fig_2a <- LMM_2_to_4(output2, structure = structure2, metric_name = "Gamma")

fig2_or_4(output_fig_2a, metric_name = "Gamma")

## Figure 2 (b)

output_fig_2b <- LMM_2_to_4(output2, structure = structure2, metric_name = "Alpha")

fig2_or_4(output_fig_2b, metric_name = "Alpha")

## Figure 2 (c)

output_fig_2c <- LMM_2_to_4(output2 |> filter(Synchrony != 1), structure = structure2c, metric_name = "Synchrony")

fig2_or_4(output_fig_2c, metric_name = "Synchrony")


# ========================================================================================================== #
# Figure 3. Temporal effects of species richness on stability and synchrony based on 12 consecutive overlapping 10-year moving window.

### Left

split_names3 <- str_split(names(Jena_species_biomass_data), "_", simplify = TRUE)

# Create 10-year moving windows, excluding those containing 2004
year_windows <- lapply(2003:2015, function(start) {
  yrs <- if (start == 2003) c(2003, 2005:2013) else start:(start + 9)
  if (2004 %in% yrs) return(NULL)
  as.character(yrs)
}) |> compact()

names(year_windows) <- as.character(c(2003, 2005:2015))

output_fig_3_left <- lapply(names(year_windows), function(start_year) {
  window_years <- year_windows[[start_year]]
  
  biomass_data <- lapply(Jena_species_biomass_data, \(df) df[, window_years, drop = FALSE])
  
  output3_left <- Stay_Multiple(biomass_data, order.q = c(0.5, 1, 2))
  
  output3_left |>
    mutate(
      Start_year = as.numeric(start_year),
      log2_sowndiv = rep(as.numeric(split_names3[, 3]), 3)
    )
})

Summary_fig_3_left <- bind_rows(output_fig_3_left) |>
  mutate(Order_q = paste0("q = ", Order_q)) |>
  group_by(Start_year, Order_q, log2_sowndiv) |>
  summarise(
    mean_gamma = mean(Gamma),
    mean_alpha = mean(Alpha),
    mean_syn = mean(Synchrony),
    .groups = "drop"
  )

## Figure 3 (a) left

fig_3_left(Summary_fig_3_left)$Gamma_plot

## Figure 3 (b) left

fig_3_left(Summary_fig_3_left)$Alpha_plot

## Figure 3 (c) left

fig_3_left(Summary_fig_3_left)$Synchrony_plot


### Right

structure3 <- data.frame(
  block = split_names3[, 2],
  log2_sowndiv = log2(as.numeric(split_names3[, 3]))
)

## Figure 3 (a) right

output_fig_3a_right <- slope_3(metric_name = "Gamma")

fig_3_right(output_fig_3a_right)

## Figure 3 (b) right

output_fig_3b_right <- slope_3(metric_name = "Alpha")

fig_3_right(output_fig_3b_right)

## Figure 3 (c) right

output_fig_3c_right <- slope_3(metric_name = "Synchrony")

fig_3_right(output_fig_3c_right)

# ========================================================================================================== #
# Figure 4. Biodiversity–stability and biodiversity–synchrony relationships based on 20 sets

split_names4 <- str_split(names(Jena_plot_biomass_data), "_", simplify = TRUE)

structure4 <- data.frame(
  block = split_names4[, 1],
  log2_sowndiv = log2(as.numeric(split_names4[, 2]))
)

output4 <- Stay_Multiple(Jena_plot_biomass_data, order.q = c(0.5, 1, 2))

## Figure 4 (a)

output_fig_4a <- LMM_2_to_4(output4, structure = structure4, metric_name = "Gamma")

fig2_or_4(output_fig_4a, metric_name = "Gamma")

## Figure 4 (b)

output_fig_4b <- LMM_2_to_4(output4, structure = structure4, metric_name = "Alpha")

fig2_or_4(output_fig_4b, metric_name = "Alpha")

## Figure 4 (c)

output_fig_4c <- LMM_2_to_4(output4, structure = structure4, metric_name = "Synchrony")

fig2_or_4(output_fig_4c, metric_name = "Synchrony")

