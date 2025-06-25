# Load required libraries
library(dplyr)
library(vegan)
library(tidyr)
library(ggplot2)
library(tibble) 

# Functional response variables
func_vars <- c("ugPO4", "TotN", "TotC", "ugNH4", "ugNO3")

# Replace negative values with zero
biodiversity_clean[func_vars] <- lapply(biodiversity_clean[func_vars], function(x) {
  x[x < 0] <- 0
  return(x)
})

# Convert to factors
biodiversity_clean$origin <- factor(biodiversity_clean$origin)
biodiversity_clean$water <- factor(biodiversity_clean$water)

# Create relative cover variables
biodiversity_clean <- biodiversity_clean %>%
  mutate(across(c(B1, B2, L1, L2), ~ .x / sum(c_across(B1:L2)), .names = "rel_{col}")) %>%
  mutate(across(starts_with("rel_"), ~ replace_na(.x, 0)))

# Define environmental variables
env_vars <- c("rel_B1", "rel_B2", "rel_L1", "rel_L2")

# Prepare functional matrix and transform
func_matrix <- biodiversity_clean %>% select(all_of(func_vars))
func_hel <- decostand(func_matrix, method = "hellinger")

# RDA model
rda_model <- rda(func_hel ~ (rel_B1 + rel_B2 + rel_L1 + rel_L2) * water * origin, data = biodiversity_clean)

# Run permutation tests
anova_full <- anova(rda_model, permutations = 999)
anova_terms <- anova(rda_model, by = "terms", permutations = 999)
anova_axes <- anova(rda_model, by = "axis", permutations = 999)

# Axis table
axis_table <- data.frame(
  Axis = rownames(anova_axes),
  Df = anova_axes[, "Df"],
  Variance = anova_axes[, "Variance"],
  F_value = anova_axes[, "F"],
  P_value = anova_axes[, "Pr(>F)"]
)
axis_table$Axis <- gsub("RDA", "Axis ", axis_table$Axis)
write.csv(axis_table, "RDA_axis_permutation_results.csv", row.names = FALSE)

# Prepare term-level table
anova_terms_df <- as.data.frame(anova_terms)
anova_terms_df$Significance <- cut(
  anova_terms_df$`Pr(>F)`, 
  breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", ".", "")
)
anova_terms_df$Term <- rownames(anova_terms_df)
anova_terms_df <- anova_terms_df[, c("Term", "Df", "Variance", "F", "Pr(>F)", "Significance")]
write.csv(anova_terms_df, "RDA_terms_results.csv", row.names = FALSE)

# Extract significant environmental variables (p < 0.05)
sig_env_vars <- anova_terms_df %>%
  filter(Term %in% env_vars & `Pr(>F)` < 0.05) %>%
  pull(Term)

# Extract scores
site_scores <- scores(rda_model, display = "sites") %>% as.data.frame()
species_scores <- scores(rda_model, display = "species") %>% as.data.frame()
env_scores <- scores(rda_model, display = "bp") %>% as.data.frame()

# Label rows
site_scores$Site <- rownames(site_scores)
species_scores$Variable <- rownames(species_scores)
env_scores$EnvVar <- rownames(env_scores)

# Filter only significant env arrows
env_scores_sig <- env_scores %>% filter(EnvVar %in% sig_env_vars)

# Add origin for plotting
site_scores$Origin <- biodiversity_clean$origin

# Arrow scaling factor
arrow_multiplier <- 2

site_scores$Water <- biodiversity_clean$water  # Add water factor
scores(rda_model, display = "cn")

# Get centroids of factor levels
cn_scores <- scores(rda_model, display = "cn") %>% as.data.frame()
cn_scores$Factor <- rownames(cn_scores)

# RDA biplot
p <- ggplot() +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, fill = Origin), shape = 21, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("C" = "steelblue", "O" = "forestgreen")) +
  labs(
    x = "RDA1", y = "RDA2",
    fill = "Origin"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  # Species arrows
  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = RDA1 * arrow_multiplier, yend = RDA2 * arrow_multiplier),
               arrow = arrow(length = unit(0.25, "cm")), color = "red", size = 1) +
  geom_text(data = species_scores,
            aes(x = RDA1 * arrow_multiplier * 1.1, y = RDA2 * arrow_multiplier * 1.1, label = Variable),
            color = "red", size = 4, fontface = "bold", hjust = 0.5) +
  # Significant environmental arrows
  geom_segment(data = env_scores_sig,
               aes(x = 0, y = 0, xend = RDA1 * arrow_multiplier, yend = RDA2 * arrow_multiplier),
               arrow = arrow(length = unit(0.25, "cm")), color = "darkgreen", size = 1, linetype = "dashed") +
  geom_text(data = env_scores_sig,
            aes(x = RDA1 * arrow_multiplier * 1.1, y = RDA2 * arrow_multiplier * 1.1, label = EnvVar),
            color = "darkgreen", size = 4, fontface = "italic", hjust = 0.5)

# Display plot
print(p)

