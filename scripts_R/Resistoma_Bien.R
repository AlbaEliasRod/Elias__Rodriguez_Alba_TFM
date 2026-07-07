library(tidyverse)
library(pheatmap)

# 1. Cargar la definición de genes y resistencias
# Leer el archivo sin encabezados para no perder la primera fila
genes_def_raw <- read_delim("genes_y_familias.txt", delim = "\t", col_names = FALSE)
# Limpiar: eliminar la fila que dice "GENE/RESISTANCE" si aparece, 
# y asignar nombres reales.
genes_def <- genes_def_raw %>%
  rename(GENE = X1, RESISTANCE = X2) %>%
  filter(GENE != "GENE") # Elimina la fila de encabezado que está en medio de los datos

# 2. Cargar los metadatos y la matriz de abundancia (usando delimitador de tabulación)
metadatos_raw <- read_csv("metadatos_analisis.csv")
matriz <- read_delim("matriz_resistoma.csv", delim = "\t")

# 3. Transformar la matriz a formato largo (Tidy data)
# Esto pasa los nombres de los genes de columnas a una sola columna llamada 'GENE'
df_largo_familias <- matriz %>%
  pivot_longer(cols = -c(`#FILE`, NUM_FOUND), names_to = "GENE", values_to = "Valor_Raw") %>%
  filter(Valor_Raw != ".") %>% 
  mutate(Abundancia = map_dbl(str_split(Valor_Raw, ";"), ~mean(as.numeric(.x), na.rm = TRUE))) %>%
  inner_join(genes_def, by = "GENE") %>%
  separate_rows(RESISTANCE, sep = ";") %>%
  mutate(RESISTANCE = trimws(RESISTANCE)) %>%
  # Clasificación en macrocategorías
  mutate(Familia = case_when(
    str_detect(RESISTANCE, "Amoxicillin|Ampicillin|Cef|bla|Ticarcillin|Piperacillin") ~ "Beta-lactams",
    str_detect(RESISTANCE, "Tetracycline|Doxycycline|Minocycline|tet|Tigecycline") ~ "Tetracyclines",
    str_detect(RESISTANCE, "Gentamicin|Kanamycin|Streptomycin|Amikacin|Tobramycin|Spectinomycin|ant\\(|aph\\(|aac\\(") ~ "Aminoglycosides",
    str_detect(RESISTANCE, "Chloramphenicol|Florfenicol|cat") ~ "Phenicols",
    str_detect(RESISTANCE, "Erythromycin|Clindamycin|Azithromycin|Lincomycin|lnu|erm|mef|mph|msr") ~ "MLSB",
    str_detect(RESISTANCE, "Trimethoprim|Sulfamethoxazole|sul|dfr") ~ "Sulfonamides/Trimethoprim",
    TRUE ~ "Others"
  )) %>%
  # Calculamos el promedio del gen/familia por muestra
  group_by(`#FILE`, Familia) %>%
  summarize(Abundancia_Media = mean(Abundancia, na.rm = TRUE), .groups = 'drop') %>%
  # Eliminamos el ".tab" y creamos la columna limpia de cruce
  mutate(run_accession = str_remove(`#FILE`, "\\.tab"))

# 5. Generar la matriz numérica pura (Pacientes en filas, Familias en columnas)
matriz_familias <- df_largo_familias %>%
  select(run_accession, Familia, Abundancia_Media) %>%
  distinct(run_accession, Familia, .keep_all = TRUE) %>% 
  pivot_wider(names_from = Familia, values_from = Abundancia_Media, values_fill = 0) %>% 
  column_to_rownames("run_accession") %>% 
  as.matrix()

# 6. Generar la tabla de metadatos para la barra lateral (Filas)
# Nos aseguramos de que tenga exactamente las mismas filas que la matriz
metadatos_rows <- tibble(run_accession = rownames(matriz_familias)) %>% 
  left_join(metadatos_raw, by = "run_accession") %>% 
  select(run_accession, diagnosis) %>% 
  rename(Grupo = diagnosis) %>% 
  column_to_rownames("run_accession")

# 7. Graficar el Heatmap definitivo
png("heatmap_familias_log_perfecto.png", width = 1800, height = 1500, res = 150)

pheatmap(log10(matriz_familias + 1),
         annotation_row = metadatos_rows,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Resistoma por Familias de Antibióticos")

dev.off()

# ==============================================================================
# 2. GRÁFICO TODOS LOS GENES
# ==============================================================================
# 1. Procesar, promediar y dar formato LARGO a los genes individuales
df_largo_genes <- matriz %>%
  pivot_longer(cols = -c(`#FILE`, NUM_FOUND), names_to = "GENE", values_to = "Valor_Raw") %>%
  filter(Valor_Raw != ".") %>%
  # Limpiamos los valores compuestos extrayendo la media
  mutate(Abundancia = map_dbl(str_split(Valor_Raw, ";"), ~mean(as.numeric(.x), na.rm = TRUE))) %>%
  select(`#FILE`, GENE, Abundancia) %>%
  # Eliminamos el ".tab" y creamos la columna limpia de cruce
  mutate(run_accession = str_remove(`#FILE`, "\\.tab"))

# 2. Generar la matriz numérica pura (Pacientes en filas, Genes individuales en columnas)
matriz_final_genes <- df_largo_genes %>%
  select(run_accession, GENE, Abundancia) %>%
  distinct(run_accession, GENE, .keep_all = TRUE) %>% 
  pivot_wider(names_from = GENE, values_from = Abundancia, values_fill = 0) %>% 
  column_to_rownames("run_accession") %>% 
  as.matrix()

# 3. Generar la tabla de metadatos para la barra lateral (Filas)
# Vincula de forma exacta los pacientes presentes en la matriz matemática
metadatos_genes <- tibble(run_accession = rownames(matriz_final_genes)) %>% 
  left_join(metadatos_raw, by = "run_accession") %>% 
  select(run_accession, diagnosis) %>% 
  rename(Grupo = diagnosis) %>% 
  column_to_rownames("run_accession")

# 4. Graficar el Heatmap de alta resolución
png("heatmap_todos_los_genes_log_perfecto.png", width = 4000, height = 1800, res = 150)

pheatmap(log10(matriz_final_genes + 1),
         annotation_row = metadatos_genes,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_colnames = TRUE,   
         fontsize_col = 4,       
         fontsize_row = 8,
         main = "Resistoma de Alta Resolución: Genes Individuales")

dev.off()

# ==============================================================================
# 3. ESTADISTICAS FAMILIAS
# ==============================================================================
# 1. Calculamos el test estadístico familia por familia
tabla_estadistica_familias <- df_largo_familias %>%
  left_join(metadatos_raw, by = "run_accession") %>%
  filter(!is.na(diagnosis)) %>%
  group_by(Familia) %>%
  summarise(
    Media_Abundancia_CD = mean(Abundancia_Media[diagnosis == "CD"], na.rm = TRUE),
    Media_Abundancia_nonIBD = mean(Abundancia_Media[diagnosis == "nonIBD"], na.rm = TRUE),
    # Aplicamos el test de Wilcoxon a las macrocategorías
    P_Valor = tryCatch(
      wilcox.test(Abundancia_Media ~ diagnosis, data = pick(everything()))$p.value,
      error = function(e) NA
    ),
    .groups = 'drop'
  ) %>%
  # Clasificamos según dónde es mayor la abundancia
  mutate(Asociado_A = if_else(Media_Abundancia_CD > Media_Abundancia_nonIBD, "Crohn (CD)", "Control (nonIBD)")) %>%
  arrange(P_Valor)

# 2. Guardamos el resultado en un archivo CSV en tu carpeta
write_csv(tabla_estadistica_familias, "TFM_familias_estadistica.csv")

# 3. Mostramos el resultado en la consola para verlo ya mismo
print(tabla_estadistica_familias)

# ==============================================================================
# 4. ESTADISTICAS GENES
# ==============================================================================
# 1. Aseguramos que tenemos los metadatos cruzados con el formato largo
df_estadistica_prep <- df_largo_genes %>%
  left_join(metadatos_raw, by = "run_accession") %>%
  filter(!is.na(diagnosis))

# 2. Calculamos las métricas y el test estadístico gen por gen
tabla_genes_asociados <- df_estadistica_prep %>%
  group_by(GENE) %>%
  summarise(
    Muestras_Con_El_Gen = sum(Abundancia > 0),
    Media_Abundancia_CD = mean(Abundancia[diagnosis == "CD"], na.rm = TRUE),
    Media_Abundancia_nonIBD = mean(Abundancia[diagnosis == "nonIBD"], na.rm = TRUE),
    # Ejecutamos el test no paramétrico de Wilcoxon
    P_Valor = tryCatch(
      wilcox.test(Abundancia ~ diagnosis, data = pick(everything()))$p.value,
      error = function(e) NA
    ),
    .groups = 'drop'
  ) %>%
  # Filtramos los genes que no pudieron calcularse (por estar ausentes)
  filter(!is.na(P_Valor)) %>%
  # Asociamos el gen al grupo donde su presencia/media es mayor
  mutate(Asociado_A = if_else(Media_Abundancia_CD > Media_Abundancia_nonIBD, "Crohn (CD)", "Control (nonIBD)")) %>%
  # Ordenamos de los más significativos (p-valor más bajo) a los menos
  arrange(P_Valor)

# 3. Guardamos el resultado en un archivo CSV en tu carpeta de trabajo
write_csv(tabla_genes_asociados, "TFM_genes_asociados_grupos.csv")