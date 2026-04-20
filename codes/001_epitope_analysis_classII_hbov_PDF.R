install.packages(c("tidyverse", "pheatmap", "ggseqlogo", "factoextra"))

library(tidyverse)
library(readr)

# Import
df_wide   <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class2/consolidated_peptide_table_de2f9162.csv",na.string=c("","NA"),header=T, sep =",")
df_2allele <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class2/MHCII_consolidated_peptide_table_5182db8f_sort.csv",na.string=c("","NA"),header=T, sep =",")
df_long   <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class2/peptide_table_de2f9162.csv",na.string=c("","NA"),header=T, sep =",")


## 1 Data Import & Tidying
# Tidy wide → long for heatmap analyses
df_el <- df_wide %>%
  select(peptide, starts_with("binding.netmhciipan_el.percentile")) %>%
  pivot_longer(-peptide,
               names_to  = "allele",
               names_prefix = "binding.netmhciipan_el.percentile.",
               values_to = "el_percentile")

# STEP 0: Always check exact column names after import
colnames(df_wide)
# Look for the column containing "overall" or "median"
grep("overall|median", colnames(df_wide), value = TRUE)

# Safe reference — assign the exact name to a variable
median_col <- grep("overall", colnames(df_wide), value = TRUE)[1]
cat("Using column:", median_col, "\n")


## 2 Top Binder Identification

# Overall top binders by global median
top_binders <- df_wide %>%
  arrange(.data[[median_col]]) %>%
  slice_head(n = 20)

print(top_binders[, c("peptide", median_col)])

# Per-allele top binders
top_per_allele <- df_el %>%
  filter(el_percentile <= 2) %>%
  count(allele, sort = TRUE)


## 3. Heatmap — All Alleles vs. All Peptides


library(pheatmap)

# Build matrix
heat_mat <- df_wide %>%
  select(peptide, starts_with("binding.netmhciipan_el.percentile")) %>%
  column_to_rownames("peptide") %>%
  as.matrix()

colnames(heat_mat) <- gsub("binding.netmhciipan_el.percentile.", "", colnames(heat_mat))
heat_mat

hm <- pheatmap(
  heat_mat,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  color           = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
  breaks          = seq(0, 100, length.out = 101),
  fontsize_row    = 6,
  fontsize_col    = 7,
  main            = "MHC-II EL Binding Percentile Heatmap",
  filename        = "heatmap_mhcii.png",
  width           = 14, height = 12
)



pheatmap(
  heat_mat,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  color           = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
  breaks          = seq(0, 100, length.out = 101),
  fontsize_row    = 6,
  fontsize_col    = 7,
  main            = "MHC-II EL Binding Percentile Heatmap",
  filename        = "heatmap_mhcii.png",
  width           = 14, height = 12
)
hm

## 4. Broad-Coverage Peptide Selection

# Count how many alleles each peptide binds at ≤ 10th percentile
coverage <- df_el %>%
  filter(el_percentile <= 10) %>%
  group_by(peptide) %>%
  summarise(n_alleles_covered = n_distinct(allele)) %>%
  arrange(desc(n_alleles_covered))

print(head(coverage, 20))

# Visualise
coverage_fig <- ggplot(coverage, aes(x = reorder(peptide, n_alleles_covered), y = n_alleles_covered)) +
  geom_col(fill = "#2166ac") +
  coord_flip() +
  labs(title  = "Promiscuous MHC-II Binders (EL ≤ 10%)",
       x = "Peptide", y = "No. of Alleles Covered") +
  theme_minimal(base_size = 10)
coverage_fig

## 5. Allele-Specific Analysis (DRB1*12:02 vs DRB1*15:02)
# From the 2-allele sorted table
# Check exact column names in df_2allele
grep("global|median", colnames(df_2allele), value = TRUE)

# Safe references
global_col  <- grep("global", colnames(df_2allele), value = TRUE)[1]
drb12_col   <- grep("DRB1.12", colnames(df_2allele), value = TRUE)[1]
drb15_col   <- grep("DRB1.15", colnames(df_2allele), value = TRUE)[1]

cat("Global median col:", global_col, "\n")
cat("DRB1*12:02 col:   ", drb12_col,  "\n")
cat("DRB1*15:02 col:   ", drb15_col,  "\n")

# Classify binders
df_2allele <- df_2allele %>%
  mutate(binder_class = case_when(
    .data[[global_col]] <= 10 ~ "Strong",
    .data[[global_col]] <= 30 ~ "Moderate",
    TRUE                      ~ "Weak"
  ))

# Scatter plot
df2_allelle <- ggplot(df_2allele,
       aes(x = .data[[drb12_col]],
           y = .data[[drb15_col]],
           label = peptide,
           color = binder_class)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_color_manual(values = c(Strong = "#d73027", Moderate = "#fc8d59", Weak = "grey60")) +
  labs(title = "DRB1*12:02 vs DRB1*15:02 Binding Comparison",
       x     = paste("EL Percentile —", drb12_col),
       y     = paste("EL Percentile —", drb15_col),
       color = "Binder Class") +
  theme_minimal()
df2_allelle

library(factoextra)

# Scale and cluster
mat_scaled <- scale(heat_mat)
set.seed(42)
km <- kmeans(mat_scaled, centers = 4, nstart = 25)

df_wide$cluster <- km$cluster

# Visualise
fviz_cluster(km, data = mat_scaled,
             geom = "point", ellipse.type = "convex",
             ggtheme = theme_minimal(),
             main = "Peptide Clusters by MHC-II Binding Profile")


## 6. Position & Core Motif Analysis
library(ggseqlogo)  # install.packages("ggseqlogo")

# ── Check exact column names in df_long ──────────────────────────────────────
grep("core|score|percentile", colnames(df_long), value = TRUE)

# Assign safe references (adjust if your column names differ)
core_col       <- grep("core",       colnames(df_long), value = TRUE)[1]
score_col      <- grep("score",      colnames(df_long), value = TRUE)[1]
percentile_col <- grep("percentile", colnames(df_long), value = TRUE)[1]

# ── 6a. Sequence logo with anchor-position annotations ──────────────────────
cores_top <- df_long %>%
  filter(.data[[percentile_col]] <= 5) %>%
  pull(all_of(core_col))

n_cores     <- length(cores_top)
n_peptides  <- df_long %>% filter(.data[[percentile_col]] <= 5) %>% pull(peptide) %>% n_distinct()
n_alleles_r <- df_long %>% filter(.data[[percentile_col]] <= 5) %>% pull(allele)  %>% n_distinct()

# MHC-II anchor positions: P1, P4, P6, P9 (primary anchors highlighted)
anchor_df <- data.frame(
  x     = c(1, 4, 6, 9),
  label = c("P1\n(anchor)", "P4\n(anchor)", "P6\n(anchor)", "P9\n(anchor)")
)


ggseqlogo(cores_top, method = "bits") +
  # Shade anchor columns
  annotate("rect", xmin = c(0.5,3.5,5.5,8.5), xmax = c(1.5,4.5,6.5,9.5),
           ymin = -Inf, ymax = Inf, fill = "#fdae61", alpha = 0.2) +
  annotate("text", x = c(1,4,6,9), y = -0.1,
           label = "▲", color = "#d73027", size = 3) +
  labs(
    title    = "MHC-II Core Binding Motif (EL ≤ 5th percentile)",
    subtitle = sprintf("%d cores | %d unique peptides | %d alleles  |  ▲ = MHC-II anchor positions (P1/P4/P6/P9)",
                       n_cores, n_peptides, n_alleles_r),
    x        = "Core Position",
    y        = "Information Content (bits)"
  ) +
  theme(
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

# ── 6b. Per-allele logos (faceted) ───────────────────────────────────────────
cores_by_allele <- df_long %>%
  filter(.data[[percentile_col]] <= 5) %>%
  select(allele, all_of(core_col)) %>%
  split(.[["allele"]]) %>%
  lapply(function(x) x[[core_col]])

# Only keep alleles with ≥ 3 cores (enough for a meaningful logo)
cores_by_allele <- Filter(function(x) length(x) >= 3, cores_by_allele)

ggseqlogo(cores_by_allele, method = "bits", ncol = 4) +
  labs(
    title    = "Per-Allele Core Motifs (EL ≤ 5%)",
    subtitle = sprintf("%d alleles shown (≥ 3 binders each)", length(cores_by_allele)),
    x        = "Core Position",
    y        = "Bits"
  )

# ── 6c. Position distribution with peptide labels ────────────────────────────
pos_df <- df_long %>%
  filter(.data[[percentile_col]] <= 10) %>%
  distinct(peptide, start, end) %>%
  mutate(mid = (start + end) / 2)

ggplot(pos_df, aes(x = start)) +
  geom_histogram(binwidth = 20, fill = "#4393c3", color = "white", alpha = 0.85) +
  # Rug + labels for the strongest individual binders
  geom_rug(aes(x = start), color = "#2166ac", alpha = 0.5) +
  geom_text(
    data  = pos_df %>% slice_min(start, n = 5, with_ties = FALSE),
    aes(x = start, y = 0, label = peptide),
    angle = 90, vjust = -0.3, hjust = 0, size = 2.5, color = "#d73027"
  ) +
  labs(
    title    = "Positional Distribution of Top MHC-II Binding Peptides (EL ≤ 10%)",
    subtitle = sprintf("%d unique peptides across protein sequence", nrow(pos_df)),
    x        = "Start Position in Antigen Sequence (aa)",
    y        = "Number of Peptides",
    caption  = "Red labels: 5 N-terminal-most peptides. Rug = individual peptide positions."
  ) +
  theme_minimal(base_size = 11)



library(factoextra)
library(ggrepel)   # install.packages("ggrepel") — non-overlapping labels

# ── 8a. Find optimal k with elbow + silhouette ─────────────────────────────
mat_scaled <- scale(heat_mat)

# Elbow plot (within-cluster sum of squares)
fviz_nbclust(mat_scaled, kmeans, method = "wss", k.max = 10) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "#d73027") +
  labs(
    title    = "Elbow Method — Optimal Number of Clusters",
    subtitle = "Red dashed line marks selected k",
    x        = "Number of Clusters (k)",
    y        = "Total Within-Cluster SS"
  )

# Silhouette plot (confirms cluster cohesion)
fviz_nbclust(mat_scaled, kmeans, method = "silhouette", k.max = 10) +
  labs(
    title    = "Silhouette Method — Cluster Quality",
    subtitle = "Higher average silhouette = better-defined clusters",
    x        = "Number of Clusters (k)",
    y        = "Average Silhouette Width"
  )

# ── 8b. Fit k-means with chosen k ────────────────────────────────────────────────
set.seed(42)
k_chosen <- 4   # adjust based on elbow/silhouette above
km <- kmeans(mat_scaled, centers = k_chosen, nstart = 25)

df_wide$cluster <- factor(km$cluster)  # factor for cleaner ggplot legends

# ── 8c. PCA scatter — label top binders per cluster ─────────────────────────
pca_res  <- prcomp(mat_scaled)
pca_df   <- as.data.frame(pca_res$x[, 1:2])
pca_df$peptide <- rownames(mat_scaled)
pca_df$cluster <- df_wide$cluster
pca_df$median_el <- df_wide[[median_col]]

# Variance explained for axis labels
var_exp  <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

# Flag top binder per cluster for labelling
label_df <- pca_df %>%
  group_by(cluster) %>%
  slice_min(median_el, n = 3, with_ties = FALSE)  # 3 best per cluster

cluster_epitope <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(aes(size = -median_el), alpha = 0.8) +     # bigger = stronger binder
  stat_ellipse(level = 0.85, linetype = "dashed") +
  geom_label_repel(
    data        = label_df,
    aes(label   = peptide),
    size        = 2.5,
    box.padding = 0.4,
    show.legend = FALSE
  ) +
  scale_size_continuous(
    name   = "Binding strength",
    labels = function(x) paste0(-round(x, 1), "th pctile")
  ) +
  labs(
    title    = "Peptide Clusters by MHC-II Binding Profile (PCA)",
    subtitle = sprintf("k = %d  |  Labels = top 3 binders per cluster  |  Point size ∝ binding strength",
                       k_chosen),
    x        = sprintf("PC1 (%s%% variance)", var_exp[1]),
    y        = sprintf("PC2 (%s%% variance)", var_exp[2]),
    color    = "Cluster",
    caption  = "Dashed ellipses = 85% confidence regions"
  ) +
  theme_minimal(base_size = 11)
cluster_epitope 

# ── 8d. Cluster binding-profile boxplot ─────────────────────────────────────────
# Shows the spread of global median EL percentile within each cluster
df_wide %>%
  ggplot(aes(x = cluster, y = .data[[median_col]], fill = cluster)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  geom_hline(yintercept = c(2, 10), linetype = "dashed",
             color = c("#d73027", "#fc8d59"), linewidth = 0.7) +
  annotate("text", x = 0.55, y = 2,  label = "Strong (≤2%)",    hjust = 0, size = 3, color = "#d73027") +
  annotate("text", x = 0.55, y = 10, label = "Moderate (≤10%)", hjust = 0, size = 3, color = "#fc8d59") +
  labs(
    title    = "EL Percentile Distribution by Cluster",
    subtitle = "Lower = stronger MHC-II binder",
    x        = "Cluster",
    y        = "Global Median EL Percentile",
    caption  = "Dashed lines: strong (≤2%) and moderate (≤10%) binder thresholds"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# ── 8e. Cluster summary table ────────────────────────────────────────────────────
df_wide %>%
  group_by(cluster) %>%
  summarise(
    n_peptides      = n(),
    median_el       = round(median(.data[[median_col]]), 2),
    min_el          = round(min(.data[[median_col]]),    2),
    strong_binders  = sum(.data[[median_col]] <= 2),
    moderate_binders = sum(.data[[median_col]] <= 10),
    example_peptides = paste(head(peptide, 3), collapse = ", ")
  ) %>%
  arrange(median_el) %>%
  print(width = Inf)








# Use median_col variable (defined in Section 2) for safe column reference
final_candidates <- df_wide %>%
  filter(
    .data[[median_col]] <= 15,
    rowSums(select(., starts_with("binding.netmhciipan_el.percentile")) <= 10) >= 10
  ) %>%
  arrange(.data[[median_col]]) %>%
  select(peptide, all_of(median_col))

print(final_candidates)
write_csv(final_candidates, "final_vaccine_candidates.csv")
