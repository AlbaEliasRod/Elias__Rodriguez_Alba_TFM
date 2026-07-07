library(dplyr)
library(readr)

# 1. Cargamos el reporte de ENA (el que tiene los SRR)
# Busca la columna 'run_accession' (SRR...) y 'sample_alias' (HSM...)
ena_data <- read_tsv("read_run_PRJNA389280.tsv")

# 2. Cargamos los metadatos de IBDMDB (el que tiene diagnóstico, edad, etc.)
ibdmdb_meta <- read.csv("hmp2_metadata.csv")

# Unimos por los códigos que empiezan por 'G'
meta_unida <- ena_data %>%
  inner_join(ibdmdb_meta, by = c("sample_alias" = "Project"))

# Filtramos por diagnóstico y tipo de muestra
df_estudio <- meta_unida %>%
  filter(diagnosis %in% c("CD", "nonIBD"),
         data_type == "metagenomics") %>%
  filter(grepl("Stool", IntervalName))

# Nos quedamos con la mejor muestra por paciente (usando reads_filtered)
df_pacientes_unicos <- df_estudio %>%
  group_by(Participant.ID) %>%
  arrange(desc(reads_filtered)) %>%
  slice(1) %>%
  ungroup()

# Comprobamos cuántos pacientes tenemos de cada grupo
table(df_pacientes_unicos$diagnosis)

library(MatchIt)
# 1. Limpieza: Aseguramos que usamos age_at_diagnosis y no hay NAs
df_limpio_total <- df_pacientes_unicos %>%
  # Nos aseguramos de incluir ambos grupos
  filter(diagnosis %in% c("CD", "nonIBD")) %>%
  # Usamos consent_age que sí está disponible para los sanos
  filter(!is.na(consent_age), !is.na(sex), !is.na(Antibiotics)) %>%
  # Limpiamos niveles y creamos el grupo
  mutate(diagnosis = droplevels(as.factor(diagnosis)),
         group = ifelse(diagnosis == "CD", 1, 0))

print("Conteo por grupo:")
print(table(df_limpio_total$group))

# 2. Creamos la variable 'group' como un FACTOR de dos niveles (0 y 1)
# Esto es lo que MatchIt prefiere para evitar errores de tipo de dato
df_limpio_total$group <- factor(ifelse(df_limpio_total$diagnosis == "CD", 1, 0))

# 3. Verificación de seguridad antes del MatchIt
print(table(df_limpio_total$group)) # Debería salir: 0: 22, 1: 57 (o similar)
print(levels(df_limpio_total$group)) # Debería salir: "0" "1"

# 2. Ejecución del Matching mejorado
set.seed(1234)
library(optmatch)
mod_match_optimo <- matchit(group ~ consent_age + sex + Antibiotics, 
                            data = df_limpio_total, 
                            method = "optimal") 

# 3. Extraer el dataset final
datos_tfm_balanceados <- match.data(mod_match_optimo)
print("Tabla final del TFM:")
table(datos_tfm_balanceados$diagnosis)

# 1. Mira si las medias de edad ahora son casi iguales
summary(mod_match_optimo)


library(cobalt)

# El Love Plot mostrará que los puntos ahora están casi en la línea del cero
love.plot(mod_match_optimo, 
          thresholds = c(m = .1), 
          var.order = "unadjusted",
          abs = TRUE,
          title = "Balance de Covariables (Método Óptimo)")


# Extraer los datos finales
datos_tfm_balanceados <- match.data(mod_match_optimo)

# Verificar el diagnóstico
table(datos_tfm_balanceados$diagnosis)

# Guardar en CSV para no perder el trabajo
write.csv(datos_tfm_balanceados, "metadatos_finales_TFM_balanceados.csv", row.names = FALSE)


# Esto te dará los nombres de las 44 muestras que "ganaron" el matching
muestras_a_procesar <- datos_tfm_balanceados$run_accession 

# Muestra los primeros para confirmar
head(muestras_a_procesar)

# Guarda esta lista en un archivo de texto para usarlo fuera de R
write.table(muestras_a_procesar, "muestras_fastp.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Guardamos un csv sólo con los nombres, diagnostico, sexo, edad y antibióticos
metadatos <- datos_tfm_balanceados[, c("run_accession", "diagnosis", "consent_age", "sex", "Antibiotics")]
write.table(metadatos, "metadatos.csv", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
