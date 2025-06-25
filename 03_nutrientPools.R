library(dplyr)
library(tidyr)
library(lme4)

biodiversity_clean <- biodiversity_clean_final

biodiversity_clean <- biodiversity_clean %>%
  filter(unit != 124)

# Assuming biodiversity_long is already in long format
biodiversity_C <- biodiversity_clean %>% filter(origin == "C")
biodiversity_O <- biodiversity_clean %>% filter(origin == "O")

biodiversity_C <- biodiversity_C %>%
  mutate(across(c(B1, B2, L1, L2), ~ . + 0.01))

model_glm <- lm(log(ugPO4) ~ B2,
                 data = biodiversity_O)

summary(model_glm)


ggplot(biodiversity_clean, aes(y = ugPO4, x = B2)) +
  geom_point(alpha = 0.6) +
  facet_grid(origin ~ ., scale = "free")

ggplot(biodiversity_clean, aes(y = (1 + log(TotC)), x = B2)) +
  geom_point(alpha = 0.6) +
  facet_grid(origin ~ ., scale = "free")

+
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")),
              se = FALSE, color = "blue") +
  facet_grid(water ~ cover_species, scales = "free") +
  labs(
    x = "Cover Value (+0.01)",
    y = "ugNH4",
    title = "GLM Fits of Cover vs ugNH4 by Water and Species"
  ) +
  theme_minimal()




# Add fitted values and residuals to the data
biodiversity_C$fitted <- fitted(model_long)
biodiversity_C$residuals <- resid(model_long)

# Pivot B1, B2, L1, L2 into long format
biodiversity_long <- biodiversity_C %>%
  pivot_longer(
    cols = c(B1, B2, L1, L2),
    names_to = "cover_species",
    values_to = "cover_value"
  )

ggplot(biodiversity_long, aes(x = cover_value, y = ugNH4, color = cover_species)) +
  facet_grid(water ~ cover_species, scale = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "log(Cover Value)", y = "log(ugNH4)", title = "Cover Value vs log(ugNH4) by Species")

model_interact <- lm((1 + log(ugNH4)) ~ B1 + L1, data = biodiversity_C)

summary(model_interact)

# Residuals vs Fitted
plot(model_interact, which = 1)

# Normal Q-Q
plot(model_interact, which = 2)

# Scale-Location (Homoscedasticity)
plot(model_interact, which = 3)

# Cook's distance
plot(model_interact, which = 4)

