# ---------------------------
# Load Required Libraries
# ---------------------------
library(tidyverse)
library(readxl)
library(ggExtra)
library(patchwork)
library(lubridate)
library(moments)

# ---------------------------
# Load and Preprocess COVID-19 Data
# ---------------------------
covid <- read_excel("Data sets/covid-vaccinations-vs-covid-death-rate.xlsx") %>%
  rename(Deaths = 5, Vaccination = 6) %>%
  mutate(Day = as.Date(Day, format = "%m/%d/%Y"))

# Aggregate data by country
covid_agg <- covid %>%
  group_by(Entity) %>%
  summarize(
    total_deaths = sum(Deaths, na.rm = TRUE),
    total_vaccination = sum(Vaccination, na.rm = TRUE)
  ) %>%
  ungroup()

# Perform k-means clustering on the two aggregated variables
set.seed(123)
covid_clusters <- kmeans(covid_agg %>% dplyr::select(total_deaths, total_vaccination), centers = 3)
covid_agg$cluster <- as.factor(covid_clusters$cluster)

# Compute correlation for annotation within each cluster
covid_corr <- covid_agg %>% 
  group_by(cluster) %>% 
  summarize(corr = round(cor(total_deaths, total_vaccination, use = "complete.obs"),2))

# ---------------------------
# Create Individual Plots for Each Cluster with Marginal Densities
# ---------------------------
# Split data by cluster
covid_split <- split(covid_agg, covid_agg$cluster)

# Create a plot for each cluster
covid_plots <- map(names(covid_split), function(cl) {
  df <- covid_split[[cl]]
  # Compute correlation (if there are enough points, otherwise use NA)
  corr_val <- round(cor(df$total_deaths, df$total_vaccination, use = "complete.obs"),2)
  
  p <- ggplot(df, aes(x = total_vaccination, y = total_deaths)) +
    geom_point(color = "darkblue", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linetype = "dashed") +
    labs(title = paste("COVID-19 Cluster", cl),
         x = "Total Vaccinations",
         y = "Total Deaths",
         subtitle = paste("?? =", corr_val)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  # Add marginal density plots
  ggMarginal(p, type = "density", fill = "gray", alpha = 0.4)
})

# Combine the individual cluster plots into a grid (2 columns)
combined_covid <- wrap_plots(covid_plots, ncol = 2) +
  plot_annotation(
    title = "COVID-19: Vaccination vs Deaths Stratified by Cluster",
    subtitle = "Each panel shows a different stratum with marginal density plots and correlation annotation"
  )

# Display the combined plot
print(combined_covid)
#####################################Data 2 (Ebola Data: Generation vs. Days Since Start)
# Load necessary libraries
library(tidyverse)
library(ggExtra)
library(patchwork)
library(lubridate)
library(outbreaks)

# Load and preprocess Ebola data
data("ebola_sim")
ebola_df <- if ("linelist" %in% names(ebola_sim)) ebola_sim$linelist else ebola_sim

ebola_df <- ebola_df %>%
  mutate(
    date_of_onset    = as.Date(date_of_onset),
    days_since_start = as.numeric(date_of_onset - min(date_of_onset, na.rm = TRUE))
  ) %>%
  filter(!is.na(days_since_start))

# Cluster on the auxiliary variable: days_since_start
set.seed(123)
ebola_clusters <- kmeans(ebola_df %>% dplyr::select(days_since_start), centers = 4)
ebola_df$cluster <- as.factor(ebola_clusters$cluster)

# Split the data by cluster and generate individual plots
ebola_split <- split(ebola_df, ebola_df$cluster)

ebola_plots <- map(names(ebola_split), function(cl) {
  df <- ebola_split[[cl]]
  # Compute correlation within this cluster
  corr_val <- round(cor(df$days_since_start, df$generation, use = "complete.obs"), 2)
  
  p <- ggplot(df, aes(x = days_since_start, y = generation)) +
    geom_point(color = "forestgreen", size = 2.5, alpha = 0.8) +
    geom_smooth(method = "loess", color = "darkorange", se = TRUE) +
    labs(title = paste("Ebola Cluster", cl),
         x = "Days Since Start",
         y = "Generation",
         subtitle = paste("?? =", corr_val)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  # Add marginal density plots
  ggMarginal(p, type = "density", fill = "gray", alpha = 0.4)
})

# Combine the cluster plots into a 2x2 grid
ebola_combined <- wrap_plots(ebola_plots, ncol = 2) +
  plot_annotation(
    title = "Ebola: Generation vs Days Since Start by Cluster",
    subtitle = "Each panel shows a different stratum with marginal density plots and correlation annotations"
  )

print(ebola_combined)
##############################PBC Data: Bili vs. Age
# Load necessary libraries
library(tidyverse)
library(ggExtra)
library(patchwork)
library(survival)

# Load and preprocess PBC data
data("pbc", package = "survival")
pbc_df <- pbc %>%
  mutate(
    age  = as.numeric(age),
    bili = as.numeric(bili)
  ) %>%
  filter(!is.na(age))

# Cluster on the auxiliary variable: age
set.seed(123)
pbc_clusters <- kmeans(pbc_df %>% dplyr::select(age), centers = 4)
pbc_df$cluster <- as.factor(pbc_clusters$cluster)

# Split the data by cluster and generate individual plots
pbc_split <- split(pbc_df, pbc_df$cluster)

pbc_plots <- map(names(pbc_split), function(cl) {
  df <- pbc_split[[cl]]
  # Compute correlation within this cluster
  corr_val <- round(cor(df$age, df$bili, use = "complete.obs"), 2)
  
  p <- ggplot(df, aes(x = age, y = bili)) +
    geom_point(color = "purple", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
    labs(title = paste("PBC Cluster", cl),
         x = "Age",
         y = "Bili (Serum Bilirubin)",
         subtitle = paste("?? =", corr_val)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggMarginal(p, type = "density", fill = "gray", alpha = 0.4)
})

# Combine the cluster plots into a grid
pbc_combined <- wrap_plots(pbc_plots, ncol = 2) +
  plot_annotation(
    title = "PBC: Bili vs Age by Cluster",
    subtitle = "Each panel shows a different stratum with marginal density plots and correlation annotations"
  )

print(pbc_combined)
###############################Birth Weight Data: Birth Weight vs. Mother's Age
# Load necessary libraries
library(tidyverse)
library(ggExtra)
library(patchwork)
library(MASS)

# Load and preprocess Birth Weight data
data("birthwt", package = "MASS")
birthwt_df <- birthwt %>%
  mutate(
    age = as.numeric(age),
    bwt = as.numeric(bwt)
  ) %>%
  filter(!is.na(age))

# Cluster on the auxiliary variable: age (mother's age)
set.seed(123)
birthwt_clusters <- kmeans(birthwt_df %>% dplyr::select(age), centers = 4)
birthwt_df$cluster <- as.factor(birthwt_clusters$cluster)

# Split the data by cluster and generate individual plots
birthwt_split <- split(birthwt_df, birthwt_df$cluster)

birthwt_plots <- map(names(birthwt_split), function(cl) {
  df <- birthwt_split[[cl]]
  # Compute correlation within this cluster
  corr_val <- round(cor(df$age, df$bwt, use = "complete.obs"), 2)
  
  p <- ggplot(df, aes(x = age, y = bwt)) +
    geom_point(color = "darkred", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, linetype = "dashed") +
    labs(title = paste("Birthwt Cluster", cl),
         x = "Mother's Age",
         y = "Birth Weight (grams)",
         subtitle = paste("?? =", corr_val)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggMarginal(p, type = "density", fill = "gray", alpha = 0.4)
})

# Combine the cluster plots into a grid
birthwt_combined <- wrap_plots(birthwt_plots, ncol = 2) +
  plot_annotation(
    title = "Birth Weight: BWT vs Mother's Age by Cluster",
    subtitle = "Each panel shows a different stratum with marginal density plots and correlation annotations"
  )

print(birthwt_combined)
