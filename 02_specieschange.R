library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(purrr)
library(mgcv)
library(writexl)

# ----- CONSTANTS -----
target_dates_all <- c("2017_02_20", "2017_03_27", "2017_04_24", "2017_05_30", "2017_07_04")
target_dates_first_last <- c("2017_02_20", "2017_07_04")
frames_per_unit <- 5

cb_palette <- c(
  "B1" = "#E69F00",   # orange
  "B2" = "#56B4E9",   # sky blue
  "L1" = "#009E73",   # bluish green
  "L2" = "#A6761D"    # red-brown / ochre
)

# ----- DATA PREPARATION -----

# Filter early for target dates, then combine species identity
df_prepped <- biodiversity_cover_compiled %>%
  filter(date %in% target_dates_all) %>%
  mutate(species_combo = paste(sp1, sp2, sp3, sp4, sep = "_"))

# Pivot longer and scale cover by frames per unit
df_long <- df_prepped %>%
  pivot_longer(cols = c(B1, B2, L1, L2),
               names_to = "species",
               values_to = "cover") %>%
  filter(!is.na(cover)) %>%
  mutate(scaled_cover = cover / frames_per_unit)

# Summarize total scaled cover per unit/species/date
unit_species_cover <- df_long %>%
  group_by(date, bench, water, origin, fxn, unit, species, species_combo, spp_num) %>%
  summarise(unit_cover = sum(scaled_cover, na.rm = TRUE), .groups = "drop")

# Pivot wider for biodiversity calculations
unit_cover_wide <- unit_species_cover %>%
  pivot_wider(names_from = species, values_from = unit_cover, values_fill = 0)

# Identify species columns for biodiversity metrics
species_cols <- setdiff(names(unit_cover_wide), c("date", "bench", "water", "origin", "fxn", "species_combo", "spp_num", "unit"))

# Compute biodiversity metrics per unit/date
biodiversity_metrics <- unit_cover_wide %>%
  rowwise() %>%
  mutate(
    richness = sum(c_across(all_of(species_cols)) > 0),
    total_cover = sum(c_across(all_of(species_cols)))
    # shannon = diversity(c_across(all_of(species_cols)), index = "shannon"),
    # evenness = ifelse(richness > 1, shannon / log(richness), NA)
  ) %>%
  ungroup()

# Add replicate counts per grouping for plotting
biodiversity_metrics_with_reps <- biodiversity_metrics %>%
  group_by(water, origin, fxn, species_combo, date) %>%
  mutate(reps = row_number()) %>%
  ungroup()

# ----- HELPER FUNCTION FOR COVER SUMMARY -----
cover_summary_calc <- function(df, group_vars) {
  df %>%
    pivot_longer(cols = c(B1, B2, L1, L2), names_to = "species", values_to = "percent_cover") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      mean_cover = mean(percent_cover, na.rm = TRUE),
      se_cover = sd(percent_cover, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

# ----- PLOTS -----
cover_summary_full <- cover_summary_calc(biodiversity_metrics_with_reps, c("date", "species_combo", "origin", "species"))

plot_origin <- function(data, origin_code, title_suffix) {
  ggplot(filter(data, origin == origin_code), aes(x = species_combo, y = mean_cover, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = mean_cover - se_cover, ymax = mean_cover + se_cover),
                  position = position_dodge(width = 0.9), width = 0.2) +
    facet_grid(. ~ date, scales = "free") +
    scale_fill_manual(values = cb_palette) +
    labs(
      title = paste(title_suffix, "(Origin =", origin_code, ")"),
      x = "Species Combo",
      y = "Mean Percent Cover ± SE",
      fill = "Species"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8),
      panel.spacing = unit(0.8, "lines")
    )
}

plot_origin(cover_summary_full, "O", "Organ Pipe")
plot_origin(cover_summary_full, "C", "Colorado Plateau")

# ----- FIRST & LAST DATE ANALYSIS -----

df_first_last <- biodiversity_metrics_with_reps %>%
  filter(date %in% target_dates_first_last)

cover_summary_first_last <- cover_summary_calc(df_first_last, c("date", "water", "species_combo", "origin", "species"))

plot_origin_first_last <- function(data, origin_code, title_suffix) {
  ggplot(filter(data, origin == origin_code), aes(x = species_combo, y = mean_cover, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = mean_cover - se_cover, ymax = mean_cover + se_cover),
                  position = position_dodge(width = 0.9), width = 0.2) +
    facet_grid(. ~ date, scales = "free_x") +
    scale_fill_manual(values = cb_palette) +
    labs(
      title = paste(title_suffix, "(Origin =", origin_code, ")"),
      x = "Species Combo",
      y = "Mean Percent Cover ± SE",
      fill = "Species"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.8, "lines")
    )
}

plot_origin_first_last(cover_summary_first_last, "O", "Organ Pipe (First and Last Dates)")
plot_origin_first_last(cover_summary_first_last, "C", "Colorado Plateau (First and Last Dates)")

# ----- RELATIVE COVER PLOT -----

cover_long <- biodiversity_metrics_with_reps %>%
  filter(date %in% target_dates_all) %>%
  select(unit, origin, water, date, B1, B2, L1, L2, reps) %>%
  pivot_longer(cols = c(B1, B2, L1, L2),
               names_to = "species",
               values_to = "cover") %>%
  group_by(unit, date) %>%
  mutate(total_cover = sum(cover, na.rm = TRUE),
         prop_cover = ifelse(total_cover > 0, cover / total_cover, 0)) %>%
  ungroup()

ggplot(cover_long, aes(x = date, y = prop_cover, fill = species)) +
  geom_bar(stat = "identity") +
  facet_grid(origin ~ date, scales = "free") +
  scale_fill_manual(values = cb_palette) +
  labs(
    x = "Date",
    y = "Relative Cover",
    fill = "Species",
    title = "Species Composition and Relative Cover per Sample by Origin and Date"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "lines")
  )

## PERMANOVA ON THE IMPACTS OF COVER

cover_wide_prop <- cover_long %>%
  select(unit, origin, water, date, species, prop_cover) %>%
  pivot_wider(
    names_from = species,
    values_from = prop_cover,
    values_fill = 0
  )

# Get species columns
species_cols <- c("B1", "B2", "L1", "L2")

# Remove rows where all species values are NA or 0
species_data_clean <- cover_wide_prop %>%
  filter(rowSums(across(all_of(species_cols), ~ is.na(.) | . == 0)) < length(species_cols))

# Separate species matrix and grouping info
species_mat <- species_data_clean %>%
  select(all_of(species_cols)) %>%
  replace_na(list(B1 = 0, B2 = 0, L1 = 0, L2 = 0))  # Optional: Replace NAs with 0

group_info <- species_data_clean %>%
  select(origin, date)

# Run PERMANOVA
adonis_result <- adonis2(species_mat ~ origin * date, data = group_info, method = "bray", permutations = 999)

print(adonis_result)

# Convert adonis output to data frame
adonis_df <- as.data.frame(adonis_result)

# Add row names as a new column
adonis_df <- tibble::rownames_to_column(adonis_df, var = "Source")

# Clean and format the table
permanova_table <- adonis_df %>%
  dplyr::select(Source, Df, SumOfSqs, R2, F = F, `Pr(>F)`) %>%
  dplyr::mutate(
    across(where(is.numeric), ~ signif(., 4)),
    `Pr(>F)` = ifelse(Source == "Residual" | Source == "Total", NA, 
                      ifelse(`Pr(>F)` < 0.001, "<0.001", signif(`Pr(>F)`, 3)))
  )

dispersions <- betadisper(vegdist(species_mat), group_info$origin)
anova(dispersions)

plot(dispersions)
# Export to CSV (which Excel can open)
write.csv(permanova_table, file = "PERMANOVA_results.csv", row.names = FALSE)

#-------------GENERALIZED ADDITIVE MODELS-------------------

biodiversity_metrics_with_reps <- biodiversity_metrics_with_reps %>%
  mutate(
    species_combo = as.character(species_combo),
    bench = as.character(bench),
    origin = as.character(origin),
    reps = as.character(reps),  # converting from integer to character
    plot_id = interaction(species_combo, bench, origin, reps, drop = TRUE)
  )

# Filter biodiversity_metrics by target dates
biodiv_filtered <- biodiversity_metrics_with_reps %>%
  filter(date %in% target_dates_all)

# Calculate total cover per group and relative cover per cover type
biodiv_relative <- biodiv_filtered %>%
  group_by(species_combo, date, origin) %>%
  mutate(total_cover = sum(B1, B2, L1, L2, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_longer(cols = c(B1, B2, L1, L2), names_to = "cover_type", values_to = "cover_value") %>%
  mutate(relative_cover = cover_value / total_cover)

# Convert date strings to Date and numeric for GAMM
biodiv_relative <- biodiv_relative %>%
  mutate(
    date = lubridate::ymd(gsub("_", "-", date)),
    date_num = as.numeric(date)
  )

# Summarize relative cover by cover_type, origin, date (average over species_combo)
summary_by_group <- biodiv_relative %>%
  group_by(cover_type, origin, date, date_num) %>%
  summarise(
    mean_rel_cover = mean(relative_cover, na.rm = TRUE),
    se_rel_cover = sd(relative_cover, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Fit GAMMs per cover_type and origin (without random effects)
gamm_results <- summary_by_group %>%
  group_by(origin, cover_type) %>%
  group_map(~ {
    dat <- .
    
    gamm_mod <- tryCatch(
      mgcv::gamm(mean_rel_cover ~ s(date_num, k = 3), data = dat),
      error = function(e) NULL
    )
    
    if (is.null(gamm_mod)) {
      dat %>% mutate(pred = NA_real_)
    } else {
      dat %>% mutate(pred = predict(gamm_mod$gam, newdata = dat))
    }
  }, .keep = TRUE) %>%
  bind_rows()

# Plot results separated by origin ------------------------------------

plot_gam <- function(origin_code) {
  ggplot() +
    geom_point(data = summary_by_group %>% filter(origin == origin_code), 
               aes(x = date, y = mean_rel_cover, color = cover_type), 
               size = 2, alpha = 0.6) +
    geom_errorbar(data = summary_by_group %>% filter(origin == origin_code), 
                  aes(x = date, ymin = mean_rel_cover - se_rel_cover, ymax = mean_rel_cover + se_rel_cover, color = cover_type),
                  width = 3, alpha = 0.6) +
    geom_line(data = gamm_results %>% filter(origin == origin_code), 
              aes(x = date, y = pred, color = cover_type), linewidth = 1) +
    facet_wrap(~ cover_type, scales = "free_y") +
    scale_color_manual(values = cb_palette) +  # Apply your color palette here
    scale_x_date(date_labels = "%b %d, %Y", date_breaks = "1 month") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste0("Origin ", origin_code, ": Observed and Predicted Relative Cover Over Time"),
      x = "Date",
      y = "Relative Cover",
      color = "Cover Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_O <- plot_gam("O")
plot_C <- plot_gam("C")

print(plot_O)
print(plot_C)

# Extract GAM summaries including p-value, edf, and trend ---------------------

# Get group keys and data split by origin and cover_type
groups <- summary_by_group %>%
  group_by(origin, cover_type) %>%
  group_keys()

data_splits <- summary_by_group %>%
  group_split(origin, cover_type)

# Extract model summaries per group
gam_summaries <- map2_dfr(data_splits, seq_along(data_splits), function(dat, i) {
  current_origin <- groups$origin[i]
  current_cover_type <- groups$cover_type[i]
  
  mod <- tryCatch(
    mgcv::gam(mean_rel_cover ~ s(date_num, k = 3), data = dat),
    error = function(e) NULL
  )
  
  if (is.null(mod)) {
    tibble(
      origin = current_origin,
      cover_type = current_cover_type,
      p_value = NA_real_,
      edf = NA_real_,
      r_squared = NA_real_,
      dev_explained = NA_real_,
      trend = NA_character_
    )
  } else {
    # Extract p-value and edf for smooth term
    p_val <- summary(mod)$s.table[1, "p-value"]
    edf_val <- summary(mod)$s.table[1, "edf"]
    
    # Extract R-squared and deviance explained
    r_sq <- summary(mod)$r.sq
    dev_exp <- summary(mod)$dev.expl
    
    # Predict fitted values for trend
    date_seq <- seq(min(dat$date_num), max(dat$date_num), length.out = 100)
    preds <- predict(mod, newdata = data.frame(date_num = date_seq))
    
    trend_dir <- if ((preds[length(preds)] - preds[1]) > 0) {
      "increasing"
    } else if ((preds[length(preds)] - preds[1]) < 0) {
      "decreasing"
    } else {
      "stable"
    }
    
    tibble(
      origin = current_origin,
      cover_type = current_cover_type,
      p_value = p_val,
      edf = edf_val,
      r_squared = r_sq,
      dev_explained = dev_exp,
      trend = trend_dir
    )
  }
})

# Export to Excel (adjust path as needed)
write_xlsx(gam_summaries, "gam_summaries_GlobalChangeBiology.xlsx")

# Export summarized predictions if needed
output_table <- gamm_results %>%
  arrange(origin, cover_type, date)
output_table
# write.csv(output_table, "relative_cover_predictions.csv", row.names = FALSE)