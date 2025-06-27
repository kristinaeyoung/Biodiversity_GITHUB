# 01_Data_formatting

library(dplyr)
library(tidyr)
library(vegan)

# Combine species identity into a single column
df_prepped <- biodiversity_cover_compiled %>%
  mutate(species_combo = paste(sp1, sp2, sp3, sp4, sep = "_"))

# Filter for specific target dates
target_dates <- c("2017_07_04")
df_prepped <- df_prepped %>%
  filter(date %in% target_dates)

# Pivot the raw cover data to long format
df_long <- df_prepped %>%
  pivot_longer(
    cols = c(B1, B2, L1, L2),
    names_to = "species",
    values_to = "cover"
  ) %>%
  filter(!is.na(cover))

# Scale cover by the number of frames per unit (assumed 5)
df_scaled <- df_long %>%
  mutate(scaled_cover = cover * (1 / 5))

# Summarize total cover per species per unit, per date
unit_species_cover <- df_scaled %>%
  group_by(date, bench, water, origin, unit, species, species_combo, spp_num) %>%
  summarise(
    unit_cover = sum(scaled_cover, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot wider so species become columns
unit_cover_wide <- unit_species_cover %>%
  pivot_wider(
    names_from = species,
    values_from = unit_cover,
    values_fill = 0
  )

# Identify species columns (everything else is metadata)
species_cols <- setdiff(
  names(unit_cover_wide),
  c("date", "bench", "water", "origin", "B1", "B2", "L1", "L2", "fxn", "species_combo", "spp_num", "unit")
)

# Compute biodiversity metrics using only B1, B2, L1, and L2 for total_cover
biodiversity_metrics <- unit_cover_wide %>%
  rowwise() %>%
  mutate(
    richness = sum(c_across(all_of(species_cols)) > 0),
    total_cover = sum(c_across(c("B1", "B2", "L1", "L2")), na.rm = TRUE)
    # Optional:
    # shannon = diversity(c_across(all_of(species_cols)), index = "shannon"),
    # evenness = ifelse(richness > 1, shannon / log(richness), NA)
  ) %>%
  ungroup()


# Count number of replicates (rows) per combination *per date*
biodiversity_metrics_with_reps <- biodiversity_metrics %>%
  group_by(water, origin, species_combo) %>%
  mutate(reps = row_number()) %>%
  ungroup()

# Summary table — how many reps per combo per date (long view)
biodiversity_metrics_with_reps %>%
  select(date, water, origin, species_combo, reps) %>%
  arrange(date, water, origin, species_combo, reps) %>%
  print(n = 300)

# Clean summary: count number of replicates per treatment combo per date
replicate_summary <- biodiversity_metrics %>%
  group_by(date, water, origin, species_combo) %>%
  summarise(num_reps = n(), .groups = "drop") %>%
  arrange(date, water, origin, species_combo)

# View the table
print(replicate_summary, n = 300)

# Adding functional data
# Select only the desired columns from Functional_Data_compiled
functional_selected <- Functional_Data_compiled %>%
  select(unit, slakevalue, torvane, ugPO4, TotN, TotC, ugNH4, ugNO3)

# Join with biodiversity_metrics_with_reps by "unit"
combined_data <- biodiversity_metrics_with_reps %>%
  left_join(functional_selected, by = "unit")

# View the result to check the merge
print(head(combined_data))

# Remove rows with any NA values
combined_data_clean <- combined_data %>%
  drop_na()

# View cleaned data

new_replicate_summary <- combined_data_clean %>%
  group_by(date, water, origin, species_combo) %>%
  summarise(num_reps = n(), .groups = "drop") %>%
  arrange(date, water, origin, species_combo)))

# Calculate 1.5*IQR cutoff
q1 <- quantile(combined_data_clean$ugNH4, 0.25, na.rm = TRUE)
q3 <- quantile(combined_data_clean$ugNH4, 0.75, na.rm = TRUE)
iqr <- q3 - q1
upper_bound <- q3 + 1.5 * iqr

# Filter and show rows with ugNH4 above the upper bound
outliers <- combined_data_clean %>%
  filter(ugNH4 > upper_bound)

combined_data_clean <- combined_data_clean %>%
  filter(unit != 124)

# Export the replicate summary table to CSV
write.csv(combined_data_clean, "biodiversity_clean_final.csv", row.names = FALSE)