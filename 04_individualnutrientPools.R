library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)

cb_palette <- c(
  "B1" = "#E69F00",   # orange
  "B2" = "#56B4E9",   # sky blue
  "L1" = "#009E73",   # bluish green
  "L2" = "#A6761D"    # red-brown / ochre
)

biodiversity_clean <- biodiversity_clean_final

biodiversity_clean <- biodiversity_clean %>%
  filter(!(unit %in% c(124, 63)))

# Split data by origin
biodiversity_C <- biodiversity_clean %>% filter(origin == "C")
biodiversity_O <- biodiversity_clean %>% filter(origin == "O")

# Create filtered datasets for NO3 models only
biodiversity_C_no3 <- biodiversity_C %>% filter(ugNO3 > 0)
biodiversity_O_no3 <- biodiversity_O %>% filter(ugNO3 > 0)


# Fit models (reusing C and O from before for NH4)
models <- list(
  NH4_C = lm(log(ugNH4) ~ L1, data = biodiversity_C),
  NH4_O = lm(log(ugNH4) ~ L1, data = biodiversity_O)
)

model_results <- lapply(names(models), function(name) {
  model <- models[[name]]
  shapiro_p <- shapiro.test(residuals(model))$p.value
  summary_model <- summary(model)
  
  data.frame(
    Model = name,
    Estimate_L1 = summary_model$coefficients["L1", "Estimate"],
    P_value_L1 = summary_model$coefficients["L1", "Pr(>|t|)"],
    Adj_R2 = summary_model$adj.r.squared,
    Shapiro_p = shapiro_p
  )
}) %>%
  bind_rows()

print(model_results)

## CREATING GRAPHS
# Model for origin C
model_nh4_C <- lm(log(ugNH4) ~ L1, data = biodiversity_C)

# Create new data for prediction
pred_C <- data.frame(
  L1 = seq(min(biodiversity_C$L1, na.rm = TRUE), max(biodiversity_C$L1, na.rm = TRUE), length.out = 100)
)

# Predict with confidence intervals on log scale
pred_fit_C <- predict(model_nh4_C, newdata = pred_C, interval = "confidence")
pred_C$fit <- exp(pred_fit_C[, "fit"])
pred_C$lwr <- exp(pred_fit_C[, "lwr"])
pred_C$upr <- exp(pred_fit_C[, "upr"])

# Plot
p_nh4_C <- ggplot(biodiversity_C, aes(x = L1, y = ugNH4)) +
  geom_point(color = cb_palette["L1"], size = 2.5, alpha = 0.8) +
  geom_line(data = pred_C, aes(x = L1, y = fit), color = cb_palette["L1"], size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_C, aes(x = L1, ymin = lwr, ymax = upr), alpha = 0.2, fill = cb_palette["L1"], inherit.aes = FALSE) +
  labs(title = "Origin C", x = expression(L[1]~"cover (%)"), y = expression(NH[4]^+~(µg~g^{-1}))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_nh4_C)


# ORGAN PIPE
# Prediction data frame
pred_O <- data.frame(
  L1 = seq(min(biodiversity_O$L1, na.rm = TRUE), max(biodiversity_O$L1, na.rm = TRUE), length.out = 100)
)

# Predict on log scale with confidence intervals (returns a matrix)
pred_mat <- predict(model_nh4_O, newdata = pred_O, interval = "confidence")

# Convert to data frame and rename columns to avoid duplication
pred_df <- as.data.frame(pred_mat)
colnames(pred_df) <- c("fit_log", "lwr_log", "upr_log")

# Bind with prediction data
pred_O <- cbind(pred_O, pred_df)

# Back-transform predictions and confidence intervals
pred_O <- pred_O %>%
  mutate(
    fit = exp(fit_log),
    lwr = exp(lwr_log),
    upr = exp(upr_log)
  )

p_nh4_O <- ggplot(biodiversity_O, aes(x = L1, y = ugNH4)) +
  geom_point(color = cb_palette["L1"], size = 2.5, alpha = 0.8) +
  geom_line(data = pred_O, aes(x = L1, y = fit), color = cb_palette["L1"], size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_O, aes(x = L1, ymin = lwr, ymax = upr), alpha = 0.2, fill = cb_palette["L1"], inherit.aes = FALSE) +
  labs(title = "Origin O", x = expression(L[1]~"cover (%)"), y = expression(NH[4]^+~(µg~g^{-1}))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_nh4_O)


