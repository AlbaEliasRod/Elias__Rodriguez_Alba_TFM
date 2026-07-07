# ==============================================================================
# 1. CARGA DE LIBRERÍAS (Instala las que falten con install.packages() o BiocManager)
# ==============================================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("phyloseq", "microbiome", "DESeq2"))
ainstall.packages("tidyverse")

library(phyloseq)
library(tidyverse)
library(microbiome)
library(vegan)
library(scales)

# ==============================================================================
# 2. CARGA DE METADATOS (CSV del MatchIt)
# ==============================================================================
# 1. Establece tu carpeta de trabajo (donde tienes los .report y el CSV)
setwd("~/TFM/TFM_R_Analysis") 

# 2. Cargar tu CSV de metadatos
metadata <- read.csv("metadatos_analisis.csv", header = TRUE)
# Aseguramos que la columna 'run_accession' sea el ID de las filas
rownames(metadata) <- metadata$run_accession

# ==============================================================================
# 3. IMPORTACIÓN DE REPORTES DE KRAKEN2
# ==============================================================================
# Buscamos todos los archivos .report
archivos <- list.files(pattern = "\\.report$")

# Función para extraer solo los Géneros (rank == 'G')
leer_kraken <- function(archivo) {
  df <- read.delim(archivo, header = FALSE)
  colnames(df) <- c("pct", "clade_counts", "tax_counts", "rank", "taxid", "name")
  
  df %>% 
    filter(rank == "G") %>% 
    mutate(name = str_trim(name)) %>%
    select(name, clade_counts)
}

# Leemos y unimos todos los archivos en una sola tabla
lista_datos <- map(archivos, ~{
  muestra_id <- str_remove(.x, "\\.report")
  d <- leer_kraken(.x)
  colnames(d)[2] <- muestra_id
  return(d)
})

tabla_abundancia <- lista_datos %>% 
  reduce(full_join, by = "name") %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("name")

# ==============================================================================
# 4. CREACIÓN DEL OBJETO PHYLOSEQ
# ==============================================================================
OTU <- otu_table(as.matrix(tabla_abundancia), taxa_are_rows = TRUE)
SAMP <- sample_data(metadata)
ps <- phyloseq(OTU, SAMP)

# Limpieza: Eliminar géneros con abundancia 0 en todas las muestras
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# ==============================================================================
# ABUNDANCIA RELATIVA A NIVEL DE FILO
# ==============================================================================
# 1. Obtenemos todos los nombres de taxones que tienes
# 1. Creamos una tabla con 2 columnas para asegurar que tenga dimensiones positivas
taxones_nombres <- taxa_names(ps)
tax_mat <- matrix(NA, nrow = length(taxones_nombres), ncol = 2)
rownames(tax_mat) <- taxones_nombres
colnames(tax_mat) <- c("Kingdom", "Phylum")

# 2. Asignamos los reinos y filos de forma masiva
tax_mat[, "Kingdom"] <- "Bacteria"

# Asignación de Filos basada en tus géneros detectados
tax_mat[rownames(tax_mat) %in% c("Bacteroides", "Phocaeicola", "Parabacteroides", "Tannerella", "Alistipes", "Paraprevotella"), "Phylum"] <- "Bacteroidota"
tax_mat[rownames(tax_mat) %in% c("Veillonella", "Pediococcus", "Enterocloster", "Faecalibacterium"), "Phylum"] <- "Bacillota"
tax_mat[rownames(tax_mat) %in% c("Escherichia", "Shigella", "Leclercia", "Alkalicella"), "Phylum"] <- "Pseudomonadota"
tax_mat[rownames(tax_mat) %in% c("Akkermansia"), "Phylum"] <- "Verrucomicrobiota"
tax_mat[rownames(tax_mat) %in% c("Cohcovirus"), "Phylum"] <- "Vira"

# Los que falten los agrupamos en "Otros"
tax_mat[is.na(tax_mat[, "Phylum"]), "Phylum"] <- "Otros"

# 3. Reasignamos al objeto phyloseq
tax_table(ps) <- tax_table(tax_mat)

# 4. Intentamos el gráfico de nuevo (esta vez funcionará porque hay >1 columna)
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

# 1. Preparar los datos: calcular la media por grupo
data_plot <- psmelt(ps_phylum_rel) %>%
  group_by(diagnosis, Phylum) %>%
  summarize(Abundancia = mean(Abundance), .groups = 'drop') %>%
  mutate(Abundancia = Abundancia * 100) # Convertir a porcentaje

# 2. Reordenar los Filos para que los más abundantes queden abajo (como en tu imagen)
# Esto hace que Bacteroidetes y Firmicutes sean la base del gráfico
data_plot$Phylum <- reorder(data_plot$Phylum, data_plot$Abundancia)

# 3. Crear el gráfico
p_final_filos <- ggplot(data_plot, aes(x = diagnosis, y = Abundancia, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white") +
  scale_fill_brewer(palette = "Paired") + # Paleta de colores similar a la imagen
  theme_minimal() +
  labs(title = "Composición taxonómica a nivel de filo",
       subtitle = "Comparación de abundancia relativa media (%)",
       x = "", 
       y = "Abundancia relativa (%)",
       fill = "Taxa") +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12)
  )

# 4. Mostrar y guardar
print(p_final_filos)
ggsave("Composicion_Filos_Grupos.png", p_final_filos, width = 7, height = 8, dpi = 300)


# 1. Calcular las medias exactas (en porcentaje)
tabla_datos_filos <- psmelt(ps_phylum_rel) %>%
  group_by(diagnosis, Phylum) %>%
  summarize(Media_Porcentaje = mean(Abundance) * 100, .groups = 'drop') %>%
  pivot_wider(names_from = diagnosis, values_from = Media_Porcentaje) %>%
  # Ordenar por abundancia en el grupo CD
  arrange(desc(CD)) %>%
  # Redondear a dos decimales para que sea legible
  mutate(across(where(is.numeric), ~ round(.x, 2)))
print(tabla_datos_filos)
write.csv(tabla_datos_filos, "Tabla_Filos_TFM_Final.csv", row.names = FALSE)

# ==============================================================================
# NUEVO APARTADO: COMPOSICIÓN TAXONÓMICA A NIVEL DE GÉNERO
# ==============================================================================
# 1. Extraemos la matriz de cuentas directas del objeto phyloseq
counts_mat <- as.data.frame(as(otu_table(ps), "matrix"))
if (taxa_are_rows(ps) == FALSE) counts_mat <- as.data.frame(t(counts_mat))

# 2. Porcentajes por muestra
prop_mat <- as.data.frame(apply(counts_mat, 2, function(x) (x / sum(x)) * 100))
prop_mat$OTU <- rownames(prop_mat)

# 3. Datos clínicos
metadata_df <- data.frame(sample_data(ps))
metadata_df$SampleID <- rownames(metadata_df)

# 4. PASO CLAVE: Forzamos el filtrado estricto antes de calcular medias
data_genero <- prop_mat %>%
  pivot_longer(cols = -OTU, names_to = "SampleID", values_to = "Abundancia") %>%
  left_join(metadata_df, by = "SampleID") %>%
  filter(OTU != "Homo") %>%            # <-- Borrado físico de la fila
  filter(!grepl("Homo", OTU))          # <-- Filtro de seguridad extra por si acaso

# 5. Calculamos las medias microbiológicas reales
data_genero_medias <- data_genero %>%
  group_by(diagnosis, OTU) %>%
  summarize(Abundancia = mean(Abundancia), .groups = 'drop')

# 6. Seleccionamos el Top 14 microbiano + Akkermansia
top_generos <- data_genero_medias %>%
  group_by(OTU) %>%
  summarize(Total = sum(Abundancia)) %>%
  arrange(desc(Total)) %>%
  slice_head(n = 14) %>%
  pull(OTU)

top_generos <- unique(c(top_generos, "Akkermansia"))

# 7. Agrupamos y recolocamos "Otros"
data_plot_genero <- data_genero_medias %>%
  mutate(Genero = ifelse(OTU %in% top_generos, OTU, "Otros")) %>%
  group_by(diagnosis, Genero) %>%
  summarize(Abundancia = sum(Abundancia), .groups = 'drop')

data_plot_genero$Genero <- factor(data_plot_genero$Genero, 
                                  levels = c(setdiff(top_generos, "Otros"), "Otros"))

# 8. GUARDADO FORZADO SOBRE EL ARCHIVO ORIGINAL
tabla_generos_valores <- data_plot_genero %>%
  pivot_wider(names_from = diagnosis, values_from = Abundancia) %>%
  arrange(desc(CD)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write.csv(tabla_generos_valores, "Tabla_Generos_TFM_Final.csv", row.names = FALSE)

# 6. Forzamos el orden de la leyenda
data_plot_genero$Genero <- factor(data_plot_genero$Genero, 
                                  levels = c(setdiff(top_generos, "Otros"), "Otros"))

# Esto estira la barra de CD y la de nonIBD en los datos reales
data_plot_genero_100 <- data_plot_genero %>%
  group_by(diagnosis, Genero) %>%
  summarise(Abundancia = sum(Abundancia, na.rm = TRUE), .groups = 'drop') %>%
  group_by(diagnosis) %>%
  mutate(Abundancia = (Abundancia / sum(Abundancia)) * 100) %>%
  ungroup()

# 7. Dibujamos el gráfico (Sintaxis moderna sin avisos)
p_generos_grupos <- ggplot(data_plot_genero_100, aes(x = diagnosis, y = Abundancia, fill = Genero)) +
  # Usamos tu geom_col normal que no rompe los factores ni los colores
  geom_col(width = 0.55, color = "white", linewidth = 0.2) + 
  
  # Tu paleta de colores original (aseguramos los niveles con top_generos)
  scale_fill_manual(values = c(hue_pal()(length(unique(data_plot_genero_100$Genero)) - 1), "grey85")) + 
  
  # Eje Y limpio: va de 0 a 100 reales porque los datos ya suman 100
  scale_y_continuous(
    limits = c(0, 100), 
    breaks = seq(0, 100, by = 25),
    expand = c(0, 0) # Pega las barras perfectamente al techo del 100% y al suelo del 0%
  ) +
  
  theme_minimal() +
  labs(title = "Composición Taxonómica General a nivel de Género",
       subtitle = "Comparación de abundancias relativas medias por grupo",
       x = "Grupo de Diagnóstico", 
       y = "Abundancia Relativa Media (%)",
       fill = "Género") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12),
    legend.text = element_text(face = "italic")
  )

print(p_generos_grupos)

# 8. Guardamos el gráfico final limpio en tu carpeta
ggsave("Composicion_Generos_Por_Grupo_100.png", p_generos_grupos, width = 8, height = 8, dpi = 300)

# 9. Exportar la tabla numérica exacta en CSV para tus anexos o tablas de soporte
tabla_datos_generos <- data_plot_genus %>%
  pivot_wider(names_from = diagnosis, values_from = Abundancia) %>%
  arrange(desc(CD)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write.csv(tabla_datos_generos, "Tabla_Generos_TFM_Final.csv", row.names = FALSE)

# ==============================================================================
# 5. ANÁLISIS DE DIVERSIDAD ALFA (Objetivo 3)
# ==============================================================================
# Calculamos las métricas con la librería microbiome
tab_alfa <- microbiome::alpha(ps, index = c("shannon", "observed", "evenness_pielou"))

# Unimos con los metadatos (diagnóstico)
tab_alfa$diagnosis <- sample_data(ps)$diagnosis
tab_alfa$run_accession <- rownames(tab_alfa)

# Guardamos la tabla para el TFM
write.csv(tab_alfa, "resultados_diversidad_alfa.csv")

library(ggpubr)
# --- GRÁFICO 1: RIQUEZA (OBSERVED) ---
# Mide cuántos géneros diferentes hay.
p_observed_pro <- ggplot(tab_alfa, aes(x = diagnosis, y = observed, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8, size = 0.7, width = 0.6) +
  geom_jitter(color = "black", width = 0.2, alpha = 0.4, size = 2) +
  scale_fill_manual(values = c("CD" = "#F8766D", "nonIBD" = "#00BFC4")) +
  theme_light() +
  labs(x = "Grupo", y = "Riqueza (Número de Géneros)") +
  # Ajusta el límite superior según tus datos (ej. si el máximo es 120, pon 150)
  expand_limits(y = 0) + 
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.4) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

print(p_observed_pro)

# --- GRÁFICO 2: SHANNON ---
# Mide abundancia y uniformidad. Es el más importante.
p_shannon_pro <- ggplot(tab_alfa, aes(x = diagnosis, y = diversity_shannon, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8, size = 0.7, width = 0.6) +
  geom_jitter(color = "black", width = 0.2, alpha = 0.4, size = 2) +
  scale_fill_manual(values = c("CD" = "#F8766D", "nonIBD" = "#00BFC4")) +
  theme_light() +
  labs(x = "Grupo", y = "Índice de Shannon") +
  scale_y_continuous(breaks = seq(0, 3.5, 0.5), limits = c(0, 3.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.4, label.y = 3.3) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

print(p_shannon_pro)

# --- GRÁFICO 3: EVENNESS (PIELOU) ---
# Mide cómo de equilibradas están las poblaciones.
p_evenness_pro <- ggplot(tab_alfa, aes(x = diagnosis, y = evenness_pielou, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8, size = 0.7, width = 0.6) +
  geom_jitter(color = "black", width = 0.2, alpha = 0.4, size = 2) +
  scale_fill_manual(values = c("CD" = "#F8766D", "nonIBD" = "#00BFC4")) +
  theme_light() +
  labs(x = "Grupo", y = "Índice de Equidad (Evenness)") +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), limits = c(0, 0.6)) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.4, label.y = 0.55) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

print(p_evenness_pro)

# Guardar
ggsave("Grafico_Observed_Final.png", p_observed_pro, width = 6, height = 5)
ggsave("Grafico_Shannon_Final.png", p_shannon_pro, width = 6, height = 5)
ggsave("Grafico_Evenness_Final.png", p_evenness_pro, width = 6, height = 5)


# Estadística descriptiva
library(dplyr)
library(tidyr)

# Calculamos estadística descriptiva completa (Mediana, IQR, Mínimo y Máximo)
tabla_estadisticos <- tab_alfa %>%
  group_by(diagnosis) %>%
  summarise(
    # Riqueza (Observed)
    Obs_Mediana = median(observed),
    Obs_IQR = IQR(observed),
    # Shannon
    Sha_Mediana = median(diversity_shannon),
    Sha_IQR = IQR(diversity_shannon),
    # Evenness
    Eve_Mediana = median(evenness_pielou),
    Eve_IQR = IQR(evenness_pielou)
  )

# Transponemos la tabla para que sea más fácil de leer en el TFM
tabla_tfm <- pivot_longer(tabla_estadisticos, cols = -diagnosis, 
                          names_to = "Metrica", values_to = "Valor")

print(tabla_tfm)
# Guardar en CSV para no perder los datos
write.csv(tabla_tfm, "tablas_medianas_diversidad.csv")

# ==============================================================================
# 6. ANÁLISIS DE DIVERSIDAD BETA (PCoA) (Objetivo 3)
# ==============================================================================
# Normalizamos (transformamos a abundancia relativa) para comparar
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Función optimizada para el estilo del panel
graficar_pcoa_final <- function(ps_obj, dist_method, titulo) {
  dist_mat <- phyloseq::distance(ps_obj, method = dist_method)
  ord <- ordinate(ps_obj, method = "PCoA", distance = dist_mat)
  
  # Cálculo de varianza para los ejes
  evals <- ord$values$Eigenvalues
  pc1 <- round(evals[1] / sum(evals) * 100, 1)
  pc2 <- round(evals[2] / sum(evals) * 100, 1)
  
  plot_ordination(ps_obj, ord, color = "diagnosis") +
    geom_point(size = 2.5, alpha = 0.7) +
    # Añadimos elipses de confianza (95%)
    stat_ellipse(aes(group = diagnosis), linetype = 2, size = 0.5, alpha = 0.5) +
    scale_color_manual(values = c("CD" = "#F8766D", "nonIBD" = "#00BFC4")) +
    theme_minimal() +
    labs(title = titulo, 
         subtitle = paste0("PCoA1: ", pc1, "% | PCoA2: ", pc2, "%"),
         x = "PCoA 1", y = "PCoA 2") +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey80", fill = NA, size = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
}

# Generar los dos gráficos principales
p_bray <- graficar_pcoa_final(ps_rel, "bray", "Bray-Curtis")
p_jaccard <- graficar_pcoa_final(ps_rel, "jaccard", "Jaccard")

# Unir con espacio y leyenda
panel_pro <- (p_bray + plot_spacer() + p_jaccard) + 
  plot_layout(widths = c(1, 0.1, 1), guides = 'collect') & 
  theme(legend.position = "right") &
  labs(color = "Diagnóstico")

# --- GUARDADO MAESTRO ---
# Esta proporción (2:1) evita que los gráficos se compriman
ggsave("Panel_Beta_TFM_Final.png", 
       panel_pro, 
       width = 12, 
       height = 5.5, 
       dpi = 300)


# 1. Aseguramos que las distancias existan (por si se borraron)
dist_bray <- phyloseq::distance(ps_rel, method = "bray")
dist_jaccard <- phyloseq::distance(ps_rel, method = "jaccard", binary = TRUE)

# 2. Extraemos la metadata
metadata <- data.frame(sample_data(ps_rel))

# 3. Lanzamos PERMANOVA para Bray-Curtis
set.seed(123)
permanova_bray <- adonis2(dist_bray ~ diagnosis, data = metadata)
print("--- RESULTADOS BRAY-CURTIS ---")
print(permanova_bray)

# 4. Lanzamos PERMANOVA para Jaccard
set.seed(123)
permanova_jaccard <- adonis2(dist_jaccard ~ diagnosis, data = metadata)
print("--- RESULTADOS JACCARD ---")
print(permanova_jaccard)

# ==============================================================================
# 7. IDENTIDFICACIÓN DE BIOMARCADORES
# ==============================================================================
# Filtramos el objeto 'ps' para quedarnos con todos los taxones EXCEPTO "Homo"
todos_los_taxones <- taxa_names(ps)
taxones_filtrados <- todos_los_taxones[todos_los_taxones != "Homo"]
ps <- prune_taxa(taxones_filtrados, ps)

# 1. Convertir objeto phyloseq a DESeq2
diag_ds2 <- phyloseq_to_deseq2(ps, ~ diagnosis)

# 2. Ejecutar el test de abundancia diferencial
# Usamos el test de Wald y un ajuste para datos de microbioma
ds2_run <- DESeq(diag_ds2, test = "Wald", fitType = "local")

# 3. Especificamos que queremos comparar CD (Grupo 1) frente a nonIBD (Referencia)
res <- results(ds2_run, contrast = c("diagnosis", "CD", "nonIBD"), cooksCutoff = FALSE)
res_df <- as.data.frame(res)

# 4. Filtrar por significación (p-adj < 0.05)
# El 'padj' es el p-valor corregido para evitar falsos positivos
biomarcadores <- res_df[which(res_df$padj < 0.05), ]
biomarcadores$taxón <- rownames(biomarcadores)

# Ordenar por importancia (log2FoldChange)
biomarcadores <- biomarcadores[order(biomarcadores$log2FoldChange), ]
print(paste("Se han identificado", nrow(biomarcadores), "géneros biomarcadores."))

# 5. Creamos una tabla limpia con los 36 biomarcadores
tabla_anexo <- biomarcadores %>%
  select(taxón, log2FoldChange, pvalue, padj) %>%
  # Añadimos una columna descriptiva para que sea más legible
  mutate(Estado = ifelse(log2FoldChange > 0, "Aumentado en CD", "Disminuido en CD")) %>%
  # Redondeamos para que quede bonito
  mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  arrange(desc(log2FoldChange))

# 6. Guardamos en un archivo que puedas abrir en Excel y copiar a Word
write.csv(tabla_anexo, "Anexo_Biomarcadores_DESeq2.csv", row.names = FALSE)

# 7. Seleccionamos los 20 con mayor cambio (ahora el hueco de Homo lo tendrá un microbio)
top_biomarcadores <- biomarcadores %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)

# 8. Creamos el gráfico de barras con la leyenda corregida
p_biomarcadores_tfm <- ggplot(top_biomarcadores, 
                              aes(x = reorder(taxón, log2FoldChange), 
                                  y = log2FoldChange, 
                                  fill = log2FoldChange > 0)) +
  geom_col(color = "black", alpha = 0.8, width = 0.7, linewidth = 0.2) + 
  coord_flip() + 
  scale_fill_manual(values = c("TRUE" = "#F8766D", "FALSE" = "#00BFC4"), 
                    breaks = c("FALSE", "TRUE"),
                    labels = c("Enriquecido en nonIBD", "Enriquecido en CD"),
                    name = "Tendencia de Abundancia") +
  theme_minimal() +
  labs(title = "Identificación de Biomarcadores Diferenciales",
       subtitle = "Análisis de Abundancia Diferencial (DESeq2, p-adj < 0.05)",
       x = "Género Bacteriano / Viral",
       y = "Log2 Fold Change (CD / nonIBD)") +
  theme(axis.text.y = element_text(face = "italic", size = 10),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5),
        panel.grid.minor = element_blank())

# 9. Mostrar y guardar el archivo final impecable
print(p_biomarcadores_tfm)
ggsave("Biomarcadores_Final_TFM_bien.png", p_biomarcadores_tfm, width = 9, height = 8, dpi = 300)

# ==============================================================================
# 7.1. IDENTIDFICACIÓN DE BIOMARCADORES
# ==============================================================================
# Función robusta para crear boxplots
crear_boxplot_biomarcador <- function(ps_obj, taxon_name) {
  
# 1. Convertimos a dataframe y nos aseguramos de que los nombres de taxones sean accesibles
df_plot <- psmelt(ps_obj)
  
# 2. Buscamos en qué columna está nuestro taxón (por si no se llama 'Genus')
# Este paso filtra el dataframe buscando el nombre del taxón en cualquier columna taxonómica
df_especifico <- df_plot %>%
  filter(if_any(where(is.character), ~ .x == taxon_name))
  
# Si el dataframe está vacío, avisamos
if(nrow(df_especifico) == 0) {
  stop(paste("No se encontró el taxón", taxon_name, "en el objeto phyloseq."))
}
  
# 3. Creamos el gráfico
ggplot(df_especifico, aes(x = diagnosis, y = Abundance, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = diagnosis)) +
  scale_fill_manual(values = c("CD" = "#F8766D", "nonIBD" = "#00BFC4")) +
  scale_color_manual(values = c("CD" = "#F8766D", "nonIBD" = "#00BFC4")) +
  theme_minimal() +
  labs(title = paste("Abundancia de", taxon_name),
        x = "", y = "Abundancia Relativa") +
  scale_y_log10() +
  theme(legend.position = "none",
        plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
        axis.title.y = element_text(size = 10),
        panel.grid.major.x = element_blank())
}

# Intentamos generar los gráficos de nuevo
# (Asegúrate de que los nombres coincidan con los de tu gráfico de barras)
try({
  p_akk <- crear_boxplot_biomarcador(ps_rel, "Akkermansia")
  p_vei <- crear_boxplot_biomarcador(ps_rel, "Veillonella")
  p_coh <- crear_boxplot_biomarcador(ps_rel, "Cohcovirus")
  
  panel_boxplots <- p_akk | p_vei | p_coh
  print(panel_boxplots)
  ggsave("Boxplots_Biomarcadores_Corregido.png", panel_boxplots, width = 12, height = 4)
})
