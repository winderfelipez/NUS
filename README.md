# NUS
# ============================================================
# DETERMINANTES DE LA ADOPCIÓN DE ESPECIES NUS
# Jalq'a y Mojocoya, Bolivia
# Análisis reproducible: limpieza -> descriptivos -> EFA ->
# GLMM -> interacción región × especie -> amplitud de adopción
# -> figuras y tablas
# ============================================================

# --------------------------------
# 0. PAQUETES
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
# 1. RUTAS
# ---------------------------
# Cambie esta ruta si NUS.csv está en otra carpeta.
ruta_datos <- "NUS.csv"

dir.create("NUS_resultados", showWarnings = FALSE)
dir.create("NUS_resultados/tablas", showWarnings = FALSE)
dir.create("NUS_resultados/figuras", showWarnings = FALSE)

# ---------------------------
# 2. IMPORTACIÓN
# ---------------------------
# La base proporcionada está separada por ';' y codificada en Latin-1.
datos_raw <- read.csv(
  ruta_datos,
  sep = ";",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

datos <- datos_raw %>%
  janitor::clean_names()

cat("Filas iniciales:", nrow(datos), "\n")
cat("Columnas:", ncol(datos), "\n")

# ---------------------------
# 3. NORMALIZACIÓN DE TEXTO
# ---------------------------
normalizar_si_no <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- stringr::str_to_lower(x)
  case_when(
    x %in% c("si", "sí", "yes", "y", "1") ~ "Si",
    x %in% c("no", "n", "0") ~ "No",
    TRUE ~ NA_character_
  )
}

vars_si_no <- c(
  "has_garden",
  "knows_plant",
  "willing_to_cultivate"
)

for (v in vars_si_no) {
  if (v %in% names(datos)) datos[[v]] <- normalizar_si_no(datos[[v]])
}

datos <- datos %>%
  mutate(
    name = as.factor(name),
    bioculture = factor(bioculture),
    community = factor(community),
    education = factor(education),
    sex = factor(sex),
    plant_name = factor(plant_name),
    species = factor(species),
    status = factor(status),
    willing_bin = if_else(willing_to_cultivate == "Si", 1L,
                          if_else(willing_to_cultivate == "No", 0L, NA_integer_)),
    knows_bin = if_else(knows_plant == "Si", 1L,
                        if_else(knows_plant == "No", 0L, NA_integer_)),
    has_garden_bin = if_else(has_garden == "Yes", 1L,
                             if_else(has_garden == "No", 0L, NA_integer_))
  )

# Si Has_garden quedó como NA por etiquetas distintas:
datos <- datos %>%
  mutate(
    has_garden_bin = case_when(
      str_to_lower(as.character(has_garden)) %in% c("yes", "si", "sí") ~ 1L,
      str_to_lower(as.character(has_garden)) == "no" ~ 0L,
      TRUE ~ has_garden_bin
    )
  )

# ---------------------------
# 4. AUDITORÍA DE CALIDAD
# ---------------------------

# 4.1 Duplicados exactos
duplicados_exactos <- datos %>%
  filter(duplicated(.))

write.csv(
  duplicados_exactos,
  "NUS_resultados/tablas/duplicados_exactos.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Duplicados exactos:", nrow(duplicados_exactos), "\n")

# Eliminar SOLO duplicados exactos.
datos <- datos %>% distinct()

# 4.2 Consistencia de variables del agricultor
vars_agricultor <- c(
  "bioculture", "community", "education", "age", "sex",
  "family_size", "has_garden", "garden_caregiver"
)

consistencia <- map_dfr(vars_agricultor, function(v) {
  datos %>%
    group_by(name) %>%
    summarise(n_valores = n_distinct(.data[[v]], na.rm = TRUE),
              .groups = "drop") %>%
    summarise(
      variable = v,
      agricultores_con_inconsistencia = sum(n_valores > 1, na.rm = TRUE)
    )
})

write.csv(
  consistencia,
  "NUS_resultados/tablas/auditoria_consistencia_agricultor.csv",
  row.names = FALSE
)

# Edad inconsistente: no se corrige automáticamente.
edad_inconsistente <- datos %>%
  group_by(name) %>%
  summarise(
    n_edades = n_distinct(age, na.rm = TRUE),
    edades = paste(sort(unique(age[!is.na(age)])), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_edades > 1)

write.csv(
  edad_inconsistente,
  "NUS_resultados/tablas/edad_inconsistente.csv",
  row.names = FALSE
)

cat("Agricultores con edades inconsistentes:",
    nrow(edad_inconsistente), "\n")

# ---------------------------
# 5. REVISIÓN DE REPETICIONES AGRICULTOR × ESPECIE
# ---------------------------
rep_agri_especie <- datos %>%
  count(name, plant_name, sort = TRUE) %>%
  filter(n > 1)

write.csv(
  rep_agri_especie,
  "NUS_resultados/tablas/repeticiones_agricultor_especie.csv",
  row.names = FALSE
)

cat("Combinaciones agricultor × especie repetidas:", nrow(rep_agri_especie), "\n")

# IMPORTANTE:
# No se eliminan automáticamente las repeticiones no idénticas.
# En particular, revisar Isaño antes del análisis final.

# ---------------------------
# 6. RESUMEN DE LA MUESTRA
# ---------------------------
resumen_muestra <- tibble(
  n_registros = nrow(datos),
  n_agricultores = n_distinct(datos$name),
  n_comunidades = n_distinct(datos$community),
  n_regiones = n_distinct(datos$bioculture),
  n_especies = n_distinct(datos$plant_name),
  edad_media = mean(datos$age, na.rm = TRUE),
  edad_sd = sd(datos$age, na.rm = TRUE),
  tam_familiar_media = mean(datos$family_size, na.rm = TRUE),
  tam_familiar_sd = sd(datos$family_size, na.rm = TRUE)
)

write.csv(
  resumen_muestra,
  "NUS_resultados/tablas/resumen_muestra.csv",
  row.names = FALSE
)

print(resumen_muestra)

# ---------------------------
# 7. TABLA DESCRIPTIVA POR REGIÓN
# ---------------------------
tabla_region <- datos %>%
  distinct(name, .keep_all = TRUE) %>%
  group_by(bioculture) %>%
  summarise(
    n = n(),
    edad_media = mean(age, na.rm = TRUE),
    edad_sd = sd(age, na.rm = TRUE),
    familia_media = mean(family_size, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  tabla_region,
  "NUS_resultados/tablas/tabla_descriptiva_region.csv",
  row.names = FALSE
)

# ---------------------------
# 8. DISPOSICIÓN POR ESPECIE × REGIÓN
# ---------------------------
tabla_especie_region <- datos %>%
  group_by(bioculture, plant_name) %>%
  summarise(
    n = n(),
    dispuestos = sum(willing_bin == 1, na.rm = TRUE),
    proporcion = mean(willing_bin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    porcentaje = 100 * proporcion
  )

write.csv(
  tabla_especie_region,
  "NUS_resultados/tablas/disposicion_especie_region.csv",
  row.names = FALSE
)

# ---------------------------
# 9. FIGURA 1: HEATMAP
# ---------------------------
p_heat <- ggplot(
  tabla_especie_region,
  aes(x = bioculture, y = fct_reorder(plant_name, porcentaje),
      fill = porcentaje)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", porcentaje)),
            size = 3.3, fontface = "plain") +
  scale_fill_viridis(
    option = "D",
    limits = c(0, 100),
    name = "Willingness to cultivate (%)"
  ) +
  labs(
    x = "Biocultural region",
    y = "Species NUS",
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "plain")
  )

ggsave(
  "NUS_resultados/figuras/Figura_1_heatmap_adopcion.png",
  p_heat, width = 8, height = 8, dpi = 400
)

# ---------------------------
# 10. VARIABLES DE BARRERAS Y FACILITADORES
# ---------------------------
barreras <- c(
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

facilitadores <- c(
  "nutritional_value",
  "medicinal_value",
  "climate_advantage",
  "economic_savings",
  "commercial_opportunity",
  "biological_pest_control",
  "soil_improvement",
  "ease_of_cultivation"
)

barreras <- barreras[barreras %in% names(datos)]
facilitadores <- facilitadores[facilitadores %in% names(datos)]

# ---------------------------
# 11. DISTRIBUCIÓN DE LOS ÍTEMS
# ---------------------------
resumen_items <- bind_rows(
  map_dfr(barreras, ~tibble(
    variable = .x,
    grupo = "Barrera",
    n_unicos = n_distinct(datos[[.x]], na.rm = TRUE),
    media = mean(datos[[.x]], na.rm = TRUE),
    sd = sd(datos[[.x]], na.rm = TRUE)
  )),
  map_dfr(facilitadores, ~tibble(
    variable = .x,
    grupo = "Facilitador",
    n_unicos = n_distinct(datos[[.x]], na.rm = TRUE),
    media = mean(datos[[.x]], na.rm = TRUE),
    sd = sd(datos[[.x]], na.rm = TRUE)
  ))
)

write.csv(
  resumen_items,
  "NUS_resultados/tablas/resumen_barreras_facilitadores.csv",
  row.names = FALSE
)

# ---------------------------
# 12. EFA ORDINAL: ANÁLISIS FACTORIAL EXPLORATORIO
# ---------------------------
# Se usa matriz policórica porque los ítems son ordinales.
# La solución de 3 factores es una solución inicial que debe
# contrastarse con análisis paralelo.

hacer_efa <- function(data, vars, nombre, nf = 3) {

  if (length(vars) < 3) return(NULL)

  X <- data %>%
    select(all_of(vars)) %>%
    mutate(across(everything(), ~as.numeric(.x))) %>%
    select(where(~n_distinct(.x, na.rm = TRUE) >= 2))

  vars_ok <- names(X)

  if (length(vars_ok) < 3) return(NULL)

  cor_poly <- psych::polychoric(X)$rho

  png(
    paste0("NUS_resultados/figuras/Parallel_", nombre, ".png"),
    width = 1800, height = 1400, res = 220
  )
  pa <- psych::fa.parallel(
    cor_poly,
    n.obs = nrow(X),
    fa = "fa",
    fm = "minres",
    main = paste("Análisis paralelo:", nombre)
  )
  dev.off()

  # Solución inicial
  efa <- psych::fa(
    cor_poly,
    nfactors = nf,
    n.obs = nrow(X),
    fm = "minres",
    rotate = "oblimin"
  )

  write.csv(
    as.data.frame(unclass(efa$loadings)),
    paste0("NUS_resultados/tablas/cargas_", nombre, ".csv")
  )

  return(list(
    datos = X,
    rho = cor_poly,
    paralelo = pa,
    modelo = efa
  ))
}

efa_barreras <- hacer_efa(
  datos, barreras, "barreras", nf = 3
)

efa_facilitadores <- hacer_efa(
  datos, facilitadores, "facilitadores", nf = 3
)

# ---------------------------
# 13. ÍNDICES ROBUSTOS DE BARRERAS/FACILITADORES
# ---------------------------
# Mientras se revisa la estructura factorial, se crean índices
# estandarizados como análisis de sensibilidad.
# Para barreras: valores mayores = mayor barrera.
# Para facilitadores: valores mayores = mayor facilitador.

datos <- datos %>%
  mutate(
    indice_barreras = if (length(barreras) > 0)
      rowMeans(across(all_of(barreras)), na.rm = TRUE) else NA_real_,
    indice_facilitadores = if (length(facilitadores) > 0)
      rowMeans(across(all_of(facilitadores)), na.rm = TRUE) else NA_real_
  )
# ---------------------------
# 14. ESTANDARIZACIÓN DE PREDICTORES CONTINUOS
# ---------------------------
# Verificar variables requeridas
variables_estandarizar <- c(
  "age",
  "family_size",
  "indice_barreras",
  "indice_facilitadores",
  "cultural_revaluation",
  "institutional_promotion",
  "has_interest_nus"
)

faltantes <- setdiff(variables_estandarizar, names(datos))

if (length(faltantes) > 0) {
  stop(
    paste(
      "Las siguientes variables no existen en la base:",
      paste(faltantes, collapse = ", ")
    )
  )
}
# Estandarización Z
datos <- datos %>%
  mutate(
    age_z = as.numeric(scale(age)),
    family_size_z = as.numeric(scale(family_size)),
    indice_barreras_z = as.numeric(scale(indice_barreras)),
    indice_facilitadores_z = as.numeric(scale(indice_facilitadores)),
    cultural_revaluation_z = as.numeric(scale(cultural_revaluation)),
    institutional_promotion_z = as.numeric(scale(institutional_promotion)),
    has_interest_nus_z = as.numeric(scale(has_interest_nus))
  )
# ---------------------------
# 15. DATASET PARA EL GLMM
# ---------------------------
datos_modelo <- datos %>%
  filter(
    !is.na(willing_bin),
    !is.na(name),
    !is.na(plant_name),
    !is.na(bioculture)
  ) %>%
  droplevels()

cat("N observaciones GLMM:", nrow(datos_modelo), "\n")
cat("Agricultores:", n_distinct(datos_modelo$name), "\n")
cat("Especies:", n_distinct(datos_modelo$plant_name), "\n")

# ---------------------------
# 16. GLMM PRINCIPAL
# ---------------------------
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(broom.mixed)
# ---------------------------
# 1. GLMM MODEL WITHOUT 'community'
# ---------------------------
modelo_glmm_tmb <- glmmTMB(
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
    indice_barreras_z +
    indice_facilitadores_z +
    (1 | plant_name),
  data = datos_modelo,
  family = binomial(link = "logit")
)

summary(modelo_glmm_tmb)

# Check Nakagawa's R2
performance::r2_nakagawa(modelo_glmm_tmb)

# ---------------------------
# 2. ODDS RATIOS AND 95% CI
# ---------------------------
or_glmm <- broom.mixed::tidy(
  modelo_glmm_tmb,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

write.csv(
  or_glmm,
  "NUS_resultados/tablas/GLMM_odds_ratios.csv",
  row.names = FALSE
)
print(or_glmm)

# ---------------------------
# 3. GLMM DIAGNOSTICS
# ---------------------------
print(performance::check_model(modelo_glmm_tmb))
print(performance::check_singularity(modelo_glmm_tmb))

# DHARMa residuals
sim_glmm <- DHARMa::simulateResiduals(modelo_glmm_tmb, n = 1000)
png(
  "NUS_resultados/figuras/diagnostico_DHARMa_.png",
  width = 1800, height = 1400, res = 220
)
plot(sim_glmm)

while (!is.null(dev.list())) {
  dev.off()
}

# ---------------------------
# 4. FIGURE 2: FOREST PLOT
# ---------------------------
or_plot <- or_glmm %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = dplyr::recode(
      term,
      "biocultureMojocoya" = "Biocultural region: Mojocoya vs. Jalq'a",
      "knows_bin" = "Species knowledge",
      "age_z" = "Age",
      "sexMale" = "Sex: male",
      "educationPrimary" = "Education: Primary",
      "educationProfessional" = "Education: Professional",
      "educationSecondary" = "Education: Secondary",
      "has_garden_bin" = "Home garden ownership",
      "has_interest_nus_z" = "Interest in NUS",
      "cultural_revaluation_z" = "Cultural revaluation",
      "institutional_promotion_z" = "Institutional promotion",
      "indice_barreras_z" = "Perceived barriers",
      "indice_facilitadores_z" = "Perceived facilitators"
    )
  )

p_forest <- ggplot(
  or_plot,
  aes(
    x = estimate,
    y = reorder(term, estimate)
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
    title = "Determinants of willingness to cultivate NUS species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  "NUS_resultados/figuras/Figura_2_forestplot_GLMM.png",
  p_forest,
  width = 9,
  height = 6,
  dpi = 400
)
# ---------------------------
# 20. MODELO REGIÓN × ESPECIE
# ---------------------------
# Este modelo responde si la adopción potencial de cada especie
# depende del contexto biocultural.

modelo_interaccion <- glmer(
  willing_bin ~
    bioculture * plant_name +
    knows_bin +
    age_z +
    sex +
    has_garden_bin +
    has_interest_nus_z +
    indice_barreras_z +
    indice_facilitadores_z +
    (1 | name),
  data = datos_modelo,
  family = binomial(link = "logit"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
summary(modelo_interaccion)
# ---------------------------
# 21. COMPARACIÓN DE MODELOS
# ---------------------------
library(lmtest)
class(modelo_glmm_tmb)
class(modelo_interaccion)

formula(modelo_glmm_tmb)

formula(modelo_interaccion)

lrtest(modelo_glmm_tmb, modelo_interaccion)
# ---------------------------
# 22. PREDICCIONES MARGINALES
# ---------------------------
emm <- emmeans(
  modelo_interaccion,
  ~ bioculture | plant_name,
  type = "response"
)

emm_df <- as.data.frame(emm)

write.csv(
  emm_df,
  "NUS_resultados/tablas/predicciones_region_especie.csv",
  row.names = FALSE
)

# ---------------------------
# 23. FIGURA 3: PREDICCIONES REGIÓN × ESPECIE
# ---------------------------
p_interaccion <- ggplot(
  emm_df,
  aes(
    x = prob,
    y = reorder(plant_name, prob),
    color = bioculture
  )
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_errorbar(
    aes(xmin = asymp.LCL, xmax = asymp.UCL),
    position = position_dodge(width = 0.5),
    width = 0.15
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Adjusted probability of willingness to adopt",
    y = "NUS species",
    color = "Biocultural context",
    title = "Adjusted probability of adoption by biocultural context and species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  "NUS_resultados/figuras/Figure_3_biocultural_context_species_interaction.png",
  p_interaccion,
  width = 9,
  height = 8,
  dpi = 400
)

# ---------------------------
# 24. NÚMERO DE ESPECIES ACEPTADAS POR AGRICULTOR
# ---------------------------
amplitud <- datos %>%
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
    cultural_revaluation = first(cultural_revaluation),
    institutional_promotion = first(institutional_promotion),
    n_especies_evaluadas = n_distinct(plant_name),
    n_especies_aceptadas = sum(willing_bin == 1, na.rm = TRUE),
    proporcion_aceptacion = mean(willing_bin, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(
  amplitud,
  "NUS_resultados/tablas/amplitud_adopcion_agricultor.csv",
  row.names = FALSE
)

# ---------------------------
# 25. MODELO BETA-BINOMIAL PARA AMPLITUD
# ---------------------------
# Solo usar si el número de especies evaluadas es comparable
# entre agricultores. El modelo utiliza éxitos/fracasos directamente.

datos_amplitud <- amplitud %>%
  filter(n_especies_evaluadas > 0) %>%
  mutate(
    age_z = as.numeric(scale(age)),
    family_size_z = as.numeric(scale(family_size)),
    cultural_revaluation_z = as.numeric(scale(cultural_revaluation)),
    institutional_promotion_z = as.numeric(scale(institutional_promotion)),
    has_interest_nus_z = as.numeric(scale(has_interest_nus))
  )

head(datos_amplitud)
write.csv(
  datos_amplitud,
  "datos_amplitud.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

modelo_betabin <- glmmTMB(
  cbind(
    n_especies_aceptadas,
    n_especies_evaluadas - n_especies_aceptadas
  ) ~
    bioculture +
    age_z +
    sex +
    education +
    has_garden_bin +
    has_interest_nus_z +
    cultural_revaluation_z +
    institutional_promotion_z,
  family = betabinomial(link = "logit"),
  data = datos_amplitud
)

summary(modelo_betabin)

write.csv(
  broom.mixed::tidy(
    modelo_betabin,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ),
  "NUS_resultados/tablas/betabinomial_amplitud_OR.csv",
  row.names = FALSE
)

# ---------------------------
# 26. FIGURA 4: AMPLITUD DE ADOPCIÓN
# ---------------------------
p_amplitud <- ggplot(
  amplitud,
  aes(
    x = bioculture,
    y = n_especies_aceptadas,
    fill = bioculture
  )
) +
  geom_boxplot(alpha = 0.75, width = 0.55,
               outlier.shape = NA) +
  geom_jitter(
    width = 0.08,
    alpha = 0.45,
    size = 2
  ) +
  scale_y_continuous(
    breaks = 0:max(amplitud$n_especies_evaluadas, na.rm = TRUE)
  ) +
  labs(
    x = "Biocultural region", # Traducido
    y = "Number of accepted species", # Traducido
    title = "" # Traducido
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(
  "NUS_resultados/figuras/Figure_4_adoption_breadth.png", # Nombre de archivo traducido
  p_amplitud, width = 7, height = 6, dpi = 400
)

# ---------------------------
# 27. MODELO ORDINAL: INTERÉS EN NUS
# ---------------------------
# Variable Has_interest_NUS se trata como ordinal si conserva
# la escala original de -3 a +3.

datos_ordinal <- datos %>%
  filter(!is.na(has_interest_nus)) %>%
  mutate(
    interes_ord = ordered(
      has_interest_nus,
      levels = sort(unique(has_interest_nus))
    ),
    age_z = as.numeric(scale(age)),
    family_size_z = as.numeric(scale(family_size))
  )

# Modelo ordinal fijo. Se evita aleatorio porque el tamaño
# efectivo de agricultores es pequeño; puede ampliarse con clmm()
# si la codificación de la escala está validada.
modelo_ordinal <- ordinal::clm(
  interes_ord ~
    bioculture +
    age_z +
    sex +
    education +
    has_garden_bin +
    cultural_revaluation +
    institutional_promotion,
  data = datos_ordinal,
  link = "logit",
  Hess = TRUE
)

summary(modelo_ordinal)

# ---------------------------
# 28. EXPORTACIÓN A EXCEL
# ---------------------------
wb <- createWorkbook()

addWorksheet(wb, "Muestra")
writeData(wb, "Muestra", resumen_muestra)

addWorksheet(wb, "Especie_region")
writeData(wb, "Especie_region", tabla_especie_region)

addWorksheet(wb, "GLMM_OR")
writeData(wb, "GLMM_OR", or_glmm)

addWorksheet(wb, "Interaccion")
writeData(wb, "Interaccion", emm_df)

addWorksheet(wb, "Amplitud")
writeData(wb, "Amplitud", amplitud)

addWorksheet(wb, "Items")
writeData(wb, "Items", resumen_items)

saveWorkbook(
  wb,
  "NUS_resultados/Resultados_NUS.xlsx",
  overwrite = TRUE
)

# ---------------------------
# 29. SESIÓN Y REPRODUCIBILIDAD
# ---------------------------
writeLines(
  capture.output(sessionInfo()),
  "NUS_resultados/sessionInfo.txt"
)

# ---------------------------
# 30. MENSAJE FINAL
# ---------------------------
cat("\n====================================================\n")
cat("ANÁLISIS FINALIZADO\n")
cat("Resultados guardados en: NUS_resultados/\n")
cat("Figuras principales:\n")
cat("  Figura_1_heatmap_adopcion.png\n")
cat("  Figura_2_forestplot_GLMM.png\n")
cat("  Figura_3_interaccion_region_especie.png\n")
cat("  Figura_4_amplitud_adopcion.png\n")
cat("IMPORTANTE: revisar Isaño y edades inconsistentes antes\n")
cat("de considerar los resultados como definitivos.\n")
cat("====================================================\n")

