# ============================================================
# DETERMINANTS OF NUS SPECIES ADOPTION
# Jalq'a and Mojocoya, Bolivia
# Reproducible analysis: data cleaning -> descriptive statistics -> EFA ->
# GLMM -> region × species interaction -> adoption breadth
# -> figures and tables
# ============================================================

# --------------------------------
# 0. PACKAGES
# -------------------------------
packages <- c(
  "tidyverse", "janitor", "readr", "lme4", "broom.mixed",
  "emmeans", "performance", "DHARMa", "psych", "GPArotation",
  "ordinal", "glmmTMB", "patchwork", "scales", "viridis",
  "openxlsx", "parameters"
)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
    }
  }
}

install_if_missing(packages)

library(tidyverse)
library(janitor)
library(readr)
library(lme4)
library(broom.mixed)
library(emmeans)
library(performance)
library(DHARMa)
library(psych)
library(GPArotation)
library(ordinal)
library(glmmTMB)
library(patchwork)
library(scales)
library(viridis)
library(openxlsx)
library(parameters)

options(stringsAsFactors = FALSE)
set.seed(20260808)

# ---------------------------
# 1. DATA PATHS
# ---------------------------
# Change this path if NUS.csv is located in another folder.
data_path <- "NUS.csv"

dir.create("NUS_results", showWarnings = FALSE)
dir.create("NUS_results/tables", showWarnings = FALSE)
dir.create("NUS_results/figures", showWarnings = FALSE)

# ---------------------------
# 2. DATA IMPORT
# ---------------------------
# The provided dataset is semicolon-separated and encoded in Latin-1.
raw_data <- read.csv(
  data_path,
  sep = ";",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

data <- raw_data %>%
  janitor::clean_names()

cat("Initial rows:", nrow(data), "\n")
cat("Columns:", ncol(data), "\n")

# ---------------------------
# 3. TEXT NORMALIZATION
# ---------------------------
normalize_yes_no <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- stringr::str_to_lower(x)

  case_when(
    x %in% c("si", "sí", "yes", "y", "1") ~ "Yes",
    x %in% c("no", "n", "0") ~ "No",
    TRUE ~ NA_character_
  )
}

yes_no_vars <- c(
  "has_garden",
  "knows_plant",
  "willing_to_cultivate"
)

for (v in yes_no_vars) {
  if (v %in% names(data)) {
    data[[v]] <- normalize_yes_no(data[[v]])
  }
}

data <- data %>%
  mutate(
    name = as.factor(name),
    bioculture = factor(bioculture),
    community = factor(community),
    education = factor(education),
    sex = factor(sex),
    plant_name = factor(plant_name),
    species = factor(species),
    status = factor(status),

    willing_bin = if_else(
      willing_to_cultivate == "Yes", 1L,
      if_else(willing_to_cultivate == "No", 0L, NA_integer_)
    ),

    knows_bin = if_else(
      knows_plant == "Yes", 1L,
      if_else(knows_plant == "No", 0L, NA_integer_)
    ),

    has_garden_bin = if_else(
      has_garden == "Yes", 1L,
      if_else(has_garden == "No", 0L, NA_integer_)
    )
  )

# If has_garden contains NA values because of different labels:
data <- data %>%
  mutate(
    has_garden_bin = case_when(
      str_to_lower(as.character(has_garden)) %in%
        c("yes", "si", "sí") ~ 1L,

      str_to_lower(as.character(has_garden)) == "no" ~ 0L,

      TRUE ~ has_garden_bin
    )
  )

# ---------------------------
# 4. DATA QUALITY AUDIT
# ---------------------------

# 4.1 Exact duplicates
exact_duplicates <- data %>%
  filter(duplicated(.))

write.csv(
  exact_duplicates,
  "NUS_results/tables/exact_duplicates.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "Exact duplicates:",
  nrow(exact_duplicates),
  "\n"
)

# Remove ONLY exact duplicates.
data <- data %>%
  distinct()

# 4.2 Consistency of farmer-level variables
farmer_vars <- c(
  "bioculture",
  "community",
  "education",
  "age",
  "sex",
  "family_size",
  "has_garden",
  "garden_caregiver"
)

consistency_audit <- map_dfr(
  farmer_vars,
  function(v) {

    data %>%
      group_by(name) %>%
      summarise(
        n_values = n_distinct(
          .data[[v]],
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      summarise(
        variable = v,
        farmers_with_inconsistency =
          sum(n_values > 1, na.rm = TRUE)
      )
  }
)

write.csv(
  consistency_audit,
  "NUS_results/tables/farmer_consistency_audit.csv",
  row.names = FALSE
)

# Inconsistent age values: not automatically corrected.
inconsistent_age <- data %>%
  group_by(name) %>%
  summarise(
    n_ages = n_distinct(age, na.rm = TRUE),
    ages = paste(
      sort(unique(age[!is.na(age)])),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  filter(n_ages > 1)

write.csv(
  inconsistent_age,
  "NUS_results/tables/inconsistent_age.csv",
  row.names = FALSE
)

cat(
  "Farmers with inconsistent ages:",
  nrow(inconsistent_age),
  "\n"
)

# ---------------------------
# 5. FARMER × SPECIES REPETITIONS
# ---------------------------
farmer_species_repetitions <- data %>%
  count(
    name,
    plant_name,
    sort = TRUE
  ) %>%
  filter(n > 1)

write.csv(
  farmer_species_repetitions,
  "NUS_results/tables/farmer_species_repetitions.csv",
  row.names = FALSE
)

cat(
  "Repeated farmer × species combinations:",
  nrow(farmer_species_repetitions),
  "\n"
)

# IMPORTANT:
# Non-identical repetitions are NOT automatically removed.
# In particular, Isaño should be reviewed before the final analysis.

# ---------------------------
# 6. SAMPLE SUMMARY
# ---------------------------
sample_summary <- tibble(
  n_records = nrow(data),
  n_farmers = n_distinct(data$name),
  n_communities = n_distinct(data$community),
  n_regions = n_distinct(data$bioculture),
  n_species = n_distinct(data$plant_name),
  mean_age = mean(data$age, na.rm = TRUE),
  sd_age = sd(data$age, na.rm = TRUE),
  mean_family_size = mean(
    data$family_size,
    na.rm = TRUE
  ),
  sd_family_size = sd(
    data$family_size,
    na.rm = TRUE
  )
)

write.csv(
  sample_summary,
  "NUS_results/tables/sample_summary.csv",
  row.names = FALSE
)

print(sample_summary)

# ---------------------------
# 7. DESCRIPTIVE TABLE BY REGION
# ---------------------------
region_table <- data %>%
  distinct(name, .keep_all = TRUE) %>%
  group_by(bioculture) %>%
  summarise(
    n = n(),
    mean_age = mean(age, na.rm = TRUE),
    sd_age = sd(age, na.rm = TRUE),
    mean_family_size = mean(
      family_size,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  region_table,
  "NUS_results/tables/regional_descriptive_table.csv",
  row.names = FALSE
)

# ---------------------------
# 8. WILLINGNESS BY SPECIES × REGION
# ---------------------------
species_region_table <- data %>%
  group_by(
    bioculture,
    plant_name
  ) %>%
  summarise(
    n = n(),
    willing = sum(
      willing_bin == 1,
      na.rm = TRUE
    ),
    proportion = mean(
      willing_bin,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = 100 * proportion
  )

write.csv(
  species_region_table,
  "NUS_results/tables/willingness_species_region.csv",
  row.names = FALSE
)

# ---------------------------
# 9. FIGURE 1: HEATMAP
# ---------------------------
p_heat <- ggplot(
  species_region_table,
  aes(
    x = bioculture,
    y = fct_reorder(
      plant_name,
      percentage
    ),
    fill = percentage
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        percentage
      )
    ),
    size = 3.3,
    fontface = "plain"
  ) +
  scale_fill_viridis(
    option = "D",
    limits = c(0, 100),
    name = "Willingness to cultivate (%)"
  ) +
  labs(
    x = "Biocultural region",
    y = "NUS species",
    title = ""
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      face = "italic"
    ),
    plot.title = element_text(
      face = "plain"
    )
  )

ggsave(
  "NUS_results/figures/Figure_1_adoption_heatmap.png",
  p_heat,
  width = 8,
  height = 8,
  dpi = 400
)
# ---------------------------
# 10. BARRIERS AND FACILITATING FACTORS
# ---------------------------
barriers <- c(
  "perceived_as_weed",
  "lack_of_knowledge_benefits",
  "seed_scarcity",
  "loss_of_knowledge",
  "loss_of_culinary_tradition",
  "sensory_rejection",
  "climate_limitation",
  "water_scarcity",
  "space_limitation",
  "limited_availability"
)

facilitating_factors <- c(
  "nutritional_value",
  "medicinal_value",
  "climate_advantage",
  "economic_savings",
  "commercial_opportunity",
  "biological_pest_control",
  "soil_improvement",
  "ease_of_cultivation"
)

barriers <- barriers[barriers %in% names(data)]
facilitating_factors <- facilitating_factors[
  facilitating_factors %in% names(data)
]

# ---------------------------
# 11. ITEM DISTRIBUTION
# ---------------------------
item_summary <- bind_rows(
  map_dfr(barriers, ~tibble(
    variable = .x,
    group = "Barrier",
    n_unique = n_distinct(data[[.x]], na.rm = TRUE),
    mean = mean(data[[.x]], na.rm = TRUE),
    sd = sd(data[[.x]], na.rm = TRUE)
  )),

  map_dfr(facilitating_factors, ~tibble(
    variable = .x,
    group = "Facilitating factor",
    n_unique = n_distinct(data[[.x]], na.rm = TRUE),
    mean = mean(data[[.x]], na.rm = TRUE),
    sd = sd(data[[.x]], na.rm = TRUE)
  ))
)

write.csv(
  item_summary,
  "NUS_results/tables/summary_barriers_facilitating_factors.csv",
  row.names = FALSE
)

# ---------------------------
# 12. ORDINAL EFA:
# EXPLORATORY FACTOR ANALYSIS
# ---------------------------
# A polychoric correlation matrix is used because
# the items are ordinal.
# The three-factor solution is an initial solution
# that should be evaluated using parallel analysis.

run_efa <- function(data, vars, name, nf = 3) {

  if (length(vars) < 3) return(NULL)

  X <- data %>%
    select(all_of(vars)) %>%
    mutate(
      across(
        everything(),
        ~as.numeric(.x)
      )
    ) %>%
    select(
      where(
        ~n_distinct(.x, na.rm = TRUE) >= 2
      )
    )

  vars_ok <- names(X)

  if (length(vars_ok) < 3) return(NULL)

  poly_cor <- psych::polychoric(X)$rho

  png(
    paste0(
      "NUS_results/figures/Parallel_",
      name,
      ".png"
    ),
    width = 1800,
    height = 1400,
    res = 220
  )

  parallel_analysis <- psych::fa.parallel(
    poly_cor,
    n.obs = nrow(X),
    fa = "fa",
    fm = "minres",
    main = paste(
      "Parallel analysis:",
      name
    )
  )

  dev.off()

  # Initial factor solution
  efa <- psych::fa(
    poly_cor,
    nfactors = nf,
    n.obs = nrow(X),
    fm = "minres",
    rotate = "oblimin"
  )

  write.csv(
    as.data.frame(
      unclass(efa$loadings)
    ),
    paste0(
      "NUS_results/tables/loadings_",
      name,
      ".csv"
    )
  )

  return(
    list(
      data = X,
      rho = poly_cor,
      parallel = parallel_analysis,
      model = efa
    )
  )
}

efa_barriers <- run_efa(
  data,
  barriers,
  "barriers",
  nf = 3
)

efa_facilitating_factors <- run_efa(
  data,
  facilitating_factors,
  "facilitating_factors",
  nf = 3
)

# ---------------------------
# 13. ROBUST BARRIER AND
# FACILITATING FACTOR INDICES
# ---------------------------
# While the factor structure is being evaluated,
# standardized indices are created as a sensitivity analysis.
#
# For barriers:
# higher values = greater perceived barriers.
#
# For facilitating factors:
# higher values = stronger perceived facilitating factors.

data <- data %>%
  mutate(
    barrier_index = if (
      length(barriers) > 0
    )
      rowMeans(
        across(all_of(barriers)),
        na.rm = TRUE
      )
    else NA_real_,

    facilitating_factor_index = if (
      length(facilitating_factors) > 0
    )
      rowMeans(
        across(all_of(facilitating_factors)),
        na.rm = TRUE
      )
    else NA_real_
  )

# ---------------------------
# 14. STANDARDIZATION OF
# CONTINUOUS PREDICTORS
# ---------------------------

# Check required variables
variables_to_standardize <- c(
  "age",
  "family_size",
  "barrier_index",
  "facilitating_factor_index",
  "cultural_revaluation",
  "institutional_promotion",
  "has_interest_nus"
)

missing_variables <- setdiff(
  variables_to_standardize,
  names(data)
)

if (length(missing_variables) > 0) {
  stop(
    paste(
      "The following variables are missing from the dataset:",
      paste(
        missing_variables,
        collapse = ", "
      )
    )
  )
}

# Z-score standardization
data <- data %>%
  mutate(
    age_z = as.numeric(scale(age)),
    family_size_z = as.numeric(
      scale(family_size)
    ),
    barrier_index_z = as.numeric(
      scale(barrier_index)
    ),
    facilitating_factor_index_z = as.numeric(
      scale(facilitating_factor_index)
    ),
    cultural_revaluation_z = as.numeric(
      scale(cultural_revaluation)
    ),
    institutional_promotion_z = as.numeric(
      scale(institutional_promotion)
    ),
    has_interest_nus_z = as.numeric(
      scale(has_interest_nus)
    )
  )

# ---------------------------
# 15. DATASET FOR GLMM
# ---------------------------
model_data <- data %>%
  filter(
    !is.na(willing_bin),
    !is.na(name),
    !is.na(plant_name),
    !is.na(bioculture)
  ) %>%
  droplevels()

cat(
  "GLMM observations:",
  nrow(model_data),
  "\n"
)

cat(
  "Farmers:",
  n_distinct(model_data$name),
  "\n"
)

cat(
  "Species:",
  n_distinct(model_data$plant_name),
  "\n"
)

# ---------------------------
# 16. MAIN GLMM
# ---------------------------
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(broom.mixed)

# ---------------------------
# GLMM WITHOUT 'COMMUNITY'
# ---------------------------
glmm_model <- glmmTMB(
  willing_bin ~
    bioculture +
    knows_bin +
    age_z +
    sex +
    education +
    has_garden_bin +
    has_interest_nus_z +
    cultural_revaluation_z +
    institutional_promotion_z +
    barrier_index_z +
    facilitating_factor_index_z +
    (1 | plant_name),

  data = model_data,

  family = binomial(
    link = "logit"
  )
)

summary(glmm_model)

# Nakagawa's R2
performance::r2_nakagawa(
  glmm_model
)

# ---------------------------
# 17. ODDS RATIOS AND 95% CI
# ---------------------------
or_glmm <- broom.mixed::tidy(
  glmm_model,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

write.csv(
  or_glmm,
  "NUS_results/tables/GLMM_odds_ratios.csv",
  row.names = FALSE
)

print(or_glmm)

# ---------------------------
# 18. GLMM DIAGNOSTICS
# ---------------------------
print(
  performance::check_model(
    glmm_model
  )
)

print(
  performance::check_singularity(
    glmm_model
  )
)

# DHARMa residual diagnostics
sim_glmm <- DHARMa::simulateResiduals(
  glmm_model,
  n = 1000
)

png(
  "NUS_results/figures/DHARMa_diagnostics.png",
  width = 1800,
  height = 1400,
  res = 220
)

plot(sim_glmm)

while (!is.null(dev.list())) {
  dev.off()
}

# ---------------------------
# 19. FIGURE 2: FOREST PLOT
# ---------------------------
or_plot <- or_glmm %>%
  filter(
    term != "(Intercept)"
  ) %>%
  mutate(
    term = dplyr::recode(
      term,

      "biocultureMojocoya" =
        "Biocultural region: Mojocoya vs. Jalq'a",

      "knows_bin" =
        "Species knowledge",

      "age_z" =
        "Age",

      "sexMale" =
        "Sex: male",

      "educationPrimary" =
        "Education: primary",

      "educationProfessional" =
        "Education: professional",

      "educationSecondary" =
        "Education: secondary",

      "has_garden_bin" =
        "Home garden ownership",

      "has_interest_nus_z" =
        "Interest in NUS",

      "cultural_revaluation_z" =
        "Cultural revaluation",

      "institutional_promotion_z" =
        "Institutional promotion",

      "barrier_index_z" =
        "Perceived barriers",

      "facilitating_factor_index_z" =
        "Perceived facilitating factors"
    )
  )

p_forest <- ggplot(
  or_plot,
  aes(
    x = estimate,
    y = reorder(
      term,
      estimate
    )
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    color = "grey40"
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3,
    color = "#006D77"
  ) +
  scale_x_log10() +
  labs(
    x = "Odds Ratio (log scale)",
    y = NULL,
    title =
      "Determinants of willingness to cultivate NUS species"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "NUS_results/figures/Figure_2_GLMM_forest_plot.png",
  p_forest,
  width = 9,
  height = 6,
  dpi = 400
)
# ---------------------------
# ---------------------------
# 20. REGION × SPECIES MODEL
# ---------------------------
# This model assesses whether the potential adoption of each species
# depends on the biocultural context.

region_species_model <- glmer(
  willing_bin ~
    bioculture * plant_name +
    knows_bin +
    age_z +
    sex +
    has_garden_bin +
    has_interest_nus_z +
    barrier_index_z +
    facilitating_factor_index_z +
    (1 | name),

  data = model_data,

  family = binomial(
    link = "logit"
  ),

  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 2e5
    )
  )
)

summary(region_species_model)


# ---------------------------
# 21. MODEL COMPARISON
# ---------------------------
library(lmtest)

class(glmm_model)
class(region_species_model)

formula(glmm_model)
formula(region_species_model)

lrtest(
  glmm_model,
  region_species_model
)


# ---------------------------
# 22. MARGINAL PREDICTIONS
# ---------------------------
emm <- emmeans(
  region_species_model,
  ~ bioculture | plant_name,
  type = "response"
)

emm_df <- as.data.frame(emm)

write.csv(
  emm_df,
  "NUS_results/tables/region_species_predictions.csv",
  row.names = FALSE
)


# ---------------------------
# 23. FIGURE 3:
# REGION × SPECIES PREDICTIONS
# ---------------------------
p_interaction <- ggplot(
  emm_df,
  aes(
    x = prob,
    y = reorder(
      plant_name,
      prob
    ),
    color = bioculture
  )
) +
  geom_point(
    position = position_dodge(
      width = 0.5
    ),
    size = 3
  ) +
  geom_errorbar(
    aes(
      xmin = asymp.LCL,
      xmax = asymp.UCL
    ),
    position = position_dodge(
      width = 0.5
    ),
    width = 0.15
  ) +
  scale_x_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    limits = c(0, 1)
  ) +
  labs(
    x = "Adjusted probability of willingness to adopt",
    y = "NUS species",
    color = "Biocultural context",
    title = "Adjusted probability of adoption by biocultural context and species"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    axis.text.y = element_text(
      face = "italic"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "NUS_results/figures/Figure_3_biocultural_context_species_interaction.png",
  p_interaction,
  width = 9,
  height = 8,
  dpi = 400
)


# ---------------------------
# 24. NUMBER OF ACCEPTED SPECIES
# PER FARMER
# ---------------------------
adoption_breadth <- data %>%
  group_by(name) %>%
  summarise(
    bioculture = first(bioculture),
    community = first(community),
    age = first(age),
    sex = first(sex),
    education = first(education),
    family_size = first(family_size),
    has_garden_bin = first(has_garden_bin),
    has_interest_nus = first(has_interest_nus),
    cultural_revaluation = first(
      cultural_revaluation
    ),
    institutional_promotion = first(
      institutional_promotion
    ),

    n_species_evaluated = n_distinct(
      plant_name
    ),

    n_species_accepted = sum(
      willing_bin == 1,
      na.rm = TRUE
    ),

    acceptance_proportion = mean(
      willing_bin,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

write.csv(
  adoption_breadth,
  "NUS_results/tables/farmer_adoption_breadth.csv",
  row.names = FALSE
)


# ---------------------------
# 25. BETA-BINOMIAL MODEL
# FOR ADOPTION BREADTH
# ---------------------------
# Use only if the number of evaluated species is reasonably
# comparable across farmers. The model directly uses
# the number of successes and failures.

adoption_breadth_data <- adoption_breadth %>%
  filter(
    n_species_evaluated > 0
  ) %>%
  mutate(
    age_z = as.numeric(
      scale(age)
    ),

    family_size_z = as.numeric(
      scale(family_size)
    ),

    cultural_revaluation_z = as.numeric(
      scale(cultural_revaluation)
    ),

    institutional_promotion_z = as.numeric(
      scale(institutional_promotion)
    ),

    has_interest_nus_z = as.numeric(
      scale(has_interest_nus
    ))
  )

head(
  adoption_breadth_data
)

write.csv(
  adoption_breadth_data,
  "NUS_results/tables/adoption_breadth_data.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


beta_binomial_model <- glmmTMB(
  cbind(
    n_species_accepted,
    n_species_evaluated -
      n_species_accepted
  ) ~
    bioculture +
    age_z +
    sex +
    education +
    has_garden_bin +
    has_interest_nus_z +
    cultural_revaluation_z +
    institutional_promotion_z,

  family = betabinomial(
    link = "logit"
  ),

  data = adoption_breadth_data
)

summary(
  beta_binomial_model
)

write.csv(
  broom.mixed::tidy(
    beta_binomial_model,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ),
  "NUS_results/tables/beta_binomial_adoption_breadth_OR.csv",
  row.names = FALSE
)


# ---------------------------
# 26. FIGURE 4:
# ADOPTION BREADTH
# ---------------------------
p_adoption_breadth <- ggplot(
  adoption_breadth,
  aes(
    x = bioculture,
    y = n_species_accepted,
    fill = bioculture
  )
) +
  geom_boxplot(
    alpha = 0.75,
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    alpha = 0.45,
    size = 2
  ) +
  scale_y_continuous(
    breaks = 0:max(
      adoption_breadth$n_species_evaluated,
      na.rm = TRUE
    )
  ) +
  labs(
    x = "Biocultural region",
    y = "Number of accepted species",
    title = "Adoption breadth of NUS species by biocultural region"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "NUS_results/figures/Figure_4_adoption_breadth.png",
  p_adoption_breadth,
  width = 7,
  height = 6,
  dpi = 400
)


# ---------------------------
# 27. ORDINAL MODEL:
# INTEREST IN NUS
# ---------------------------
# The Has_interest_NUS variable is treated as ordinal
# if the original -3 to +3 scale is retained.

ordinal_data <- data %>%
  filter(
    !is.na(has_interest_nus)
  ) %>%
  mutate(
    interest_ordinal = ordered(
      has_interest_nus,
      levels = sort(
        unique(has_interest_nus)
      )
    ),

    age_z = as.numeric(
      scale(age)
    ),

    family_size_z = as.numeric(
      scale(family_size)
    )
  )


# Fixed-effects ordinal model.
# A random effect is not included because the effective
# number of farmers is small. A CLMM can be considered
# if the ordinal scale is properly validated.

ordinal_model <- ordinal::clm(
  interest_ordinal ~
    bioculture +
    age_z +
    sex +
    education +
    has_garden_bin +
    cultural_revaluation +
    institutional_promotion,

  data = ordinal_data,

  link = "logit",

  Hess = TRUE
)

summary(
  ordinal_model
)


# ---------------------------
# 28. EXPORT RESULTS TO EXCEL
# ---------------------------
wb <- createWorkbook()

addWorksheet(
  wb,
  "Sample"
)

writeData(
  wb,
  "Sample",
  sample_summary
)


addWorksheet(
  wb,
  "Species_region"
)

writeData(
  wb,
  "Species_region",
  species_region_table
)


addWorksheet(
  wb,
  "GLMM_OR"
)

writeData(
  wb,
  "GLMM_OR",
  or_glmm
)


addWorksheet(
  wb,
  "Interaction"
)

writeData(
  wb,
  "Interaction",
  emm_df
)


addWorksheet(
  wb,
  "Adoption_breadth"
)

writeData(
  wb,
  "Adoption_breadth",
  adoption_breadth
)


addWorksheet(
  wb,
  "Items"
)

writeData(
  wb,
  "Items",
  item_summary
)


saveWorkbook(
  wb,
  "NUS_results/NUS_results.xlsx",
  overwrite = TRUE
)


# ---------------------------
# 29. SESSION INFORMATION
# AND REPRODUCIBILITY
# ---------------------------
writeLines(
  capture.output(
    sessionInfo()
  ),
  "NUS_results/sessionInfo.txt"
)


# ---------------------------
# 30. FINAL MESSAGE
# ---------------------------
cat(
  "\n====================================================\n"
)

cat(
  "ANALYSIS COMPLETED\n"
)

cat(
  "Results saved in: NUS_results/\n"
)

cat(
  "Main figures:\n"
)

cat(
  "  Figure_1_adoption_heatmap.png\n"
)

cat(
  "  Figure_2_GLMM_forest_plot.png\n"
)

cat(
  "  Figure_3_biocultural_context_species_interaction.png\n"
)

cat(
  "  Figure_4_adoption_breadth.png\n"
)

cat(
  "IMPORTANT: Review Isaño and inconsistent age values\n"
)

cat(
  "before considering the results as final.\n"
)

cat(
  "====================================================\n"
)