# ```{r packages}
library(tidyverse)
library(readr)
library(pheatmap)
library(factoextra)
library(ggrepel)
library(ggseqlogo)
library(knitr)
library(kableExtra)
library(cowplot)
library(corrplot)   # method correlation matrix
library(patchwork)  # multi-panel composition
# ```

# ```{r import}
df_wide    <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class2/consolidated_peptide_table_696d662e.csv",na.string=c("","NA"),header=T, sep =",")
df_2allele <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class2/MHCII_consolidated_peptide_table_5182db8f_sort.csv",na.string=c("","NA"),header=T, sep =",")
df_long    <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class2/peptide_table_de2f9162.csv",na.string=c("","NA"),header=T, sep =",")

# ---- Method registry: name, regex pattern, plot color ---

method_registry <- tibble::tribble(
  ~method,            ~pattern,                                              ~color,
  "NetMHCIIpan_EL",   "binding\\.netmhciipan_el\\.percentile\\.",            "#d73027",
  "NetMHCIIpan_BA",   "binding\\.netmhciipan_ba\\.percentile\\.",            "#fc8d59",
  "SMM_align",        "binding\\.smm_align\\.(adjusted_)?percentile\\.",     "#4393c3",
  "Sturniolo",        "binding\\.tepitope\\.(adjusted_)?percentile\\.",      "#762a83"
)

method_registry <- method_registry %>%
  rowwise() %>%
  mutate(n_cols = sum(grepl(pattern, colnames(df_wide)))) %>%
  ungroup() %>%
  filter(n_cols > 0)

print(method_registry %>% select(method, n_cols))
# Expected:
#   NetMHCIIpan_EL   29
#   SMM_align        24
#   Sturniolo        11
# ```

# ```{r tidy}
# Build a tidy long table with one row per (peptide, method, allele)
df_methods <- map_dfr(seq_len(nrow(method_registry)), function(i) {
  m   <- method_registry$method[i]
  pat <- method_registry$pattern[i]
  cols <- grep(pat, colnames(df_wide), value = TRUE)
  if (length(cols) == 0) return(NULL)
  df_wide %>%
    select(peptide, all_of(cols)) %>%
    pivot_longer(
      -peptide,
      names_to     = "allele",
      names_prefix = pat,
      values_to    = "percentile"
    ) %>%
    mutate(method = m, .before = allele)
}) %>%
  filter(!is.na(percentile))

# Per-peptide summary across alleles, by method
per_method_summary <- df_methods %>%
  group_by(peptide, method) %>%
  summarise(
    median_pct  = median(percentile, na.rm = TRUE),
    min_pct     = min(percentile,    na.rm = TRUE),
    n_strong    = sum(percentile <= 2,  na.rm = TRUE),
    n_moderate  = sum(percentile <= 10, na.rm = TRUE),
    n_alleles   = n(),
    .groups = "drop"
  )
# ```
# Top Binder Identification (per method)

# ```{r top-binders-per-method}
top_per_method <- per_method_summary %>%
  group_by(method) %>%
  arrange(median_pct, .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup()

top_per_method %>%
  select(method, peptide, median_pct, min_pct, n_strong, n_moderate, n_alleles) %>%
  kable(
    caption   = "Top 15 peptides by median percentile, per method (lower = stronger)",
    col.names = c("Method", "Peptide", "Median %", "Min %", "# Strong (<=2%)", "# Moderate (<=10%)", "# Alleles"),
    digits    = 2
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
  collapse_rows(columns = 1, valign = "top")
# ```


# ```{r consensus-binders}
# Consensus binder = a peptide that is a *promiscuous* binder in >= 2 methods.
# Definition (literature standard for MHC-II vaccine selection):
#   method-level pass = peptide binds >= n_alleles_min alleles at <= 10% in that method
#   consensus         = method-level pass holds in >= 2 methods
# This is far more sensible than "median percentile across all alleles" because
# most peptides are non-binders for most alleles in any panel - taking a median
# across the full panel washes out true promiscuous binders.

n_alleles_min <- 3   # tunable: try 3 (exploratory), 5 (tight), 8 (paper-grade)

build_consensus <- function(summary_df, n_min) {
  summary_df %>%
    mutate(passes = n_moderate >= n_min) %>%
    group_by(peptide) %>%
    summarise(
      methods_pass        = sum(passes, na.rm = TRUE),
      methods_tested      = sum(!is.na(n_moderate)),
      methods_passed_list = paste(method[which(passes)], collapse = ", "),
      total_allele_hits   = sum(n_moderate, na.rm = TRUE),  # method x allele moderate-binder count
      best_min_pct        = min(min_pct,    na.rm = TRUE),  # strongest single allele binding
      .groups = "drop"
    ) %>%
    filter(methods_tested > 0) %>%
    mutate(consensus_class = case_when(
      methods_pass == methods_tested & methods_tested >= 2 ~ "Full consensus",
      methods_pass >= 2                                    ~ "Majority consensus",
      methods_pass == 1                                    ~ "Single-method only",
      TRUE                                                 ~ "Non-binder"
    )) %>%
    arrange(desc(methods_pass), desc(total_allele_hits), best_min_pct)
}

# Try n_alleles_min first; auto-relax if empty
threshold_used <- n_alleles_min
consensus <- build_consensus(per_method_summary, threshold_used)

if (sum(consensus$methods_pass >= 2) == 0 && threshold_used > 1) {
  threshold_used <- max(1, threshold_used - 2)
  consensus <- build_consensus(per_method_summary, threshold_used)
  message(sprintf("No peptides cleared >=%d alleles in >=2 methods. Relaxed to >=%d.",
                  n_alleles_min, threshold_used))
}

top_consensus <- consensus %>% filter(methods_pass >= 2) %>% slice_head(n = 25)

if (nrow(top_consensus) == 0) {
  cat("**No promiscuous consensus binders even at the relaxed threshold.** ",
      "Try `consensus %>% filter(methods_pass >= 1) %>% slice_head(n = 25)` ",
      "or inspect `per_method_summary %>% arrange(min_pct)` for the best ",
      "single-allele binders.\n")
} else {
  top_consensus %>%
    kable(
      caption   = sprintf("Consensus binders - bind >=%d alleles at <=10%% in >=2 methods (%d shown)",
                          threshold_used, nrow(top_consensus)),
      col.names = c("Peptide", "# Methods Passed", "# Methods Tested", "Methods",
                    "Total Allele Hits", "Best Min %", "Class"),
      digits    = 2
    ) %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
    row_spec(which(top_consensus$consensus_class == "Full consensus"),
             background = "#fde0dc", bold = TRUE)
}
# ```




# ```{r heatmap-multi, fig.height=14, fig.width=16}
# Build one heatmap per detected method, then arrange side-by-side
build_heat <- function(method_name, pattern) {
  cols <- grep(pattern, colnames(df_wide), value = TRUE)
  if (length(cols) == 0) return(NULL)
  mat <- df_wide %>%
    select(peptide, all_of(cols)) %>%
    column_to_rownames("peptide") %>%
    as.matrix()
  colnames(mat) <- gsub(pattern, "", colnames(mat))
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color        = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
    breaks       = seq(0, 100, length.out = 101),
    fontsize_row = 5,
    fontsize_col = 7,
    main         = paste0(method_name, " - EL/Score Percentile"),
    silent       = TRUE
  )
}

heat_list <- pmap(
  list(method_registry$method, method_registry$pattern),
  build_heat
)

# Display each heatmap (one per page block)
for (i in seq_along(heat_list)) {
  if (!is.null(heat_list[[i]])) grid::grid.newpage(); grid::grid.draw(heat_list[[i]]$gtable)
}
# ```

# ```{r method-correlation, fig.width=6, fig.height=6}
# Wide matrix: rows = (peptide, allele), cols = methods, values = percentile
method_wide <- df_methods %>%
  pivot_wider(names_from = method, values_from = percentile) %>%
  select(-peptide, -allele) %>%
  drop_na()

cor_mat <- cor(method_wide, method = "spearman", use = "pairwise.complete.obs")

corrplot(
  cor_mat,
  method     = "color",
  type       = "upper",
  addCoef.col = "black",
  number.cex  = 0.9,
  tl.col      = "black",
  tl.srt      = 45,
  col         = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
  title       = "Spearman correlation between MHC-II methods",
  mar         = c(0, 0, 2, 0)
)
# ```


# ```{r method-correlation, fig.width=6, fig.height=6}
# Wide matrix: rows = (peptide, allele), cols = methods, values = percentile
method_wide <- df_methods %>%
  pivot_wider(names_from = method, values_from = percentile) %>%
  select(-peptide, -allele) %>%
  drop_na()

cor_mat <- cor(method_wide, method = "spearman", use = "pairwise.complete.obs")

corrplot(
  cor_mat,
  method     = "color",
  type       = "upper",
  addCoef.col = "black",
  number.cex  = 0.9,
  tl.col      = "black",
  tl.srt      = 45,
  col         = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
  title       = "Spearman correlation between MHC-II methods",
  mar         = c(0, 0, 2, 0)
)
# ```


# ```{r promiscuous-table}
coverage_per_method <- df_methods %>%
  filter(percentile <= 10) %>%
  group_by(peptide, method) %>%
  summarise(n_alleles = n_distinct(allele), .groups = "drop")

coverage_combined <- coverage_per_method %>%
  pivot_wider(names_from = method, values_from = n_alleles, values_fill = 0) %>%
  mutate(total_method_hits = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total_method_hits))

coverage_combined %>%
  slice_head(n = 15) %>%
  kable(
    caption = "Top 15 promiscuous binders - alleles covered at <=10% per method"
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
# ```

# ```{r promiscuous-plot, fig.height=8}
top_pep_ids <- coverage_combined %>% slice_head(n = 20) %>% pull(peptide)

coverage_per_method %>%
  filter(peptide %in% top_pep_ids) %>%
  ggplot(aes(x = reorder(peptide, n_alleles, sum), y = n_alleles, fill = method)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = setNames(method_registry$color, method_registry$method)) +
  coord_flip() +
  labs(
    title    = "Promiscuity by method (top 20 combined peptides)",
    subtitle = "Bars = # alleles bound at percentile <=10% per method",
    x        = "Peptide",
    y        = "# Alleles Covered",
    fill     = "Method"
  ) +
  theme_minimal(base_size = 10)
# ```

# ```{r seqlogo-per-method, fig.height=3.5, fig.width=12}
# Requires df_long to have a 'method' column. If not, falls back to single panel.
if (!is.na(method_col_long)) {
  cores_by_method <- df_long %>%
    filter(.data[[percentile_col]] <= 5) %>%
    select(method = all_of(method_col_long), core = all_of(core_col)) %>%
    drop_na() %>%
    split(.$method) %>%
    map(~ pull(.x, core))
  ggseqlogo(cores_by_method, method = "bits", ncol = length(cores_by_method)) +
    labs(title = "Core binding motif (percentile <=5%) - per method")
} else {
  cores_top <- df_long %>% filter(.data[[percentile_col]] <= 5) %>% pull(all_of(core_col))
  ggseqlogo(cores_top, method = "bits") +
    labs(title = "Core binding motif (percentile <=5%) - all methods combined")
}
# ```

# ```{r pos-dist}
if (!is.na(method_col_long)) {
  pos_df <- df_long %>%
    filter(.data[[percentile_col]] <= 10) %>%
    distinct(peptide, start, end, method = .data[[method_col_long]])
} else {
  pos_df <- df_long %>%
    filter(.data[[percentile_col]] <= 10) %>%
    distinct(peptide, start, end) %>%
    mutate(method = "NetMHCIIpan_EL")
}

ggplot(pos_df, aes(x = start, fill = method)) +
  geom_histogram(binwidth = 20, color = "white", alpha = 0.85, position = "stack") +
  scale_fill_manual(values = setNames(method_registry$color, method_registry$method)) +
  labs(
    title    = "Positional distribution of top binders (<=10%)",
    subtitle = sprintf("%d peptide-method records", nrow(pos_df)),
    x        = "Start position in antigen (aa)",
    y        = "Count",
    fill     = "Method"
  ) +
  theme_minimal(base_size = 11)
# ```

# ```{r optimal-k, fig.width=10, fig.height=4}
# Build a consensus matrix: rows = peptides, cols = (method × allele)
consensus_mat <- df_methods %>%
  unite("feature", method, allele, sep = "__") %>%
  pivot_wider(names_from = feature, values_from = percentile) %>%
  column_to_rownames("peptide") %>%
  as.matrix()

# Mean-impute remaining NAs (e.g., Sturniolo missing for DQ/DP) before scaling
consensus_mat[is.na(consensus_mat)] <- mean(consensus_mat, na.rm = TRUE)
mat_scaled <- scale(consensus_mat)

p_elbow <- fviz_nbclust(mat_scaled, kmeans, method = "wss", k.max = 10) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "#d73027") +
  labs(title = "Elbow Method (consensus features)")
p_sil   <- fviz_nbclust(mat_scaled, kmeans, method = "silhouette", k.max = 10) +
  labs(title = "Silhouette Method (consensus features)")
plot_grid(p_elbow, p_sil, ncol = 2)
# ```

# ```{r kmeans-consensus}
set.seed(42)
k_chosen <- 4
km       <- kmeans(mat_scaled, centers = k_chosen, nstart = 25)

pca_res        <- prcomp(mat_scaled)
pca_df         <- as.data.frame(pca_res$x[, 1:2])
pca_df$peptide <- rownames(mat_scaled)
pca_df$cluster <- factor(km$cluster)

# Use median of per-method medians as size aesthetic
median_lookup <- per_method_summary %>%
  group_by(peptide) %>%
  summarise(median_of_medians = median(median_pct, na.rm = TRUE), .groups = "drop")
pca_df <- left_join(pca_df, median_lookup, by = "peptide")
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

label_df <- pca_df %>% group_by(cluster) %>% slice_min(median_of_medians, n = 3, with_ties = FALSE)

ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(aes(size = -median_of_medians), alpha = 0.8) +
  stat_ellipse(level = 0.85, linetype = "dashed") +
  geom_label_repel(data = label_df, aes(label = peptide), size = 2.5, box.padding = 0.4, show.legend = FALSE) +
  labs(
    title    = "Peptide clusters - consensus MHC-II profile (PCA)",
    subtitle = sprintf("k = %d, features = method x allele percentiles", k_chosen),
    x = sprintf("PC1 (%s%% var)", var_exp[1]),
    y = sprintf("PC2 (%s%% var)", var_exp[2])
  ) +
  theme_minimal(base_size = 11)
# ```

# ```{r cluster-summary}
# Defensive: if `km` is stale (length mismatch with current mat_scaled),
# re-run kmeans so cluster_assign is always consistent. This protects against
# the common case where Chunk 11 was re-sourced after Chunk 11b.
if (!exists("km") || length(km$cluster) != nrow(mat_scaled)) {
  message("km is stale or missing (",
          if (exists("km")) length(km$cluster) else "NA",
          " clusters vs ", nrow(mat_scaled),
          " peptides). Re-running kmeans.")
  set.seed(42)
  k_chosen <- if (exists("k_chosen")) k_chosen else 4
  km <- kmeans(mat_scaled, centers = k_chosen, nstart = 25)
}

cluster_assign <- tibble(peptide = rownames(mat_scaled), cluster = factor(km$cluster))

per_method_summary %>%
  left_join(cluster_assign, by = "peptide") %>%
  group_by(cluster, method) %>%
  summarise(
    n_peptides       = n_distinct(peptide),
    median_of_median = round(median(median_pct), 2),
    n_strong_total   = sum(n_strong),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = method,
              values_from = c(median_of_median, n_strong_total),
              names_glue  = "{method}__{.value}") %>%
  kable(caption = "Cluster summary by method - median percentile and total strong-binder hits") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
# ```


# ```{r shortlist}
# Selection rules (promiscuity-based, anchored on NetMHCIIpan_EL):
#   R1. consensus pass: methods_pass >= 2 (uses `consensus` from Chunk 5b)
#   R2. EL evidence:    NetMHCIIpan_EL n_moderate >= 1 (>=1 allele at <=10%)
#   R3. coverage:       total_allele_hits >= 5 across methods
# If empty, R2 and R3 are progressively relaxed and a message is printed.

# Defensive: rebuild `consensus` if Chunk 5b wasn't sourced this session
if (!exists("consensus")) {
  message("`consensus` not found - rebuilding inline (re-source Chunk 5b for tunable params).")
  consensus <- per_method_summary %>%
    mutate(passes = n_moderate >= 3) %>%
    group_by(peptide) %>%
    summarise(
      methods_pass        = sum(passes, na.rm = TRUE),
      methods_tested      = sum(!is.na(n_moderate)),
      methods_passed_list = paste(method[which(passes)], collapse = ", "),
      total_allele_hits   = sum(n_moderate, na.rm = TRUE),
      best_min_pct        = min(min_pct,    na.rm = TRUE),
      .groups = "drop"
    )
}

# Pull EL-specific evidence to anchor on the IEDB-recommended primary method
el_anchor <- per_method_summary %>%
  filter(method == "NetMHCIIpan_EL") %>%
  transmute(peptide,
            el_n_moderate = n_moderate,
            el_min_pct    = min_pct,
            el_median     = median_pct)

build_shortlist <- function(min_methods_pass, min_total_hits, require_el) {
  out <- consensus %>%
    filter(methods_pass >= min_methods_pass,
           total_allele_hits >= min_total_hits) %>%
    left_join(el_anchor,      by = "peptide") %>%
    left_join(cluster_assign, by = "peptide")
  if (require_el) out <- out %>% filter(!is.na(el_n_moderate), el_n_moderate >= 1)
  out %>% arrange(desc(methods_pass), desc(total_allele_hits), best_min_pct)
}

# Try strict, then progressively relax
shortlist <- build_shortlist(min_methods_pass = 2, min_total_hits = 5, require_el = TRUE)
rule_used <- "strict (>=2 methods, >=5 total hits, EL evidence required)"

if (nrow(shortlist) == 0) {
  shortlist <- build_shortlist(2, 3, TRUE)
  rule_used <- "relaxed (>=2 methods, >=3 total hits, EL evidence required)"
}
if (nrow(shortlist) == 0) {
  shortlist <- build_shortlist(2, 1, FALSE)
  rule_used <- "exploratory (>=2 methods, no EL/coverage anchor)"
}

if (nrow(shortlist) == 0) {
  cat("**No final candidates at any threshold.** ",
      "Inspect `consensus %>% filter(methods_pass >= 1) %>% arrange(desc(total_allele_hits))` ",
      "for the most promiscuous single-method binders.\n")
} else {
  message("Shortlist rule used: ", rule_used)
  write_csv(shortlist, "final_vaccine_candidates_multimethod.csv")
  shortlist %>%
    slice_head(n = 30) %>%
    select(peptide, methods_pass, methods_tested, methods_passed_list,
           total_allele_hits, best_min_pct, el_n_moderate, el_min_pct, cluster) %>%
    kable(
      caption   = sprintf("%d final candidates - %s", nrow(shortlist), rule_used),
      col.names = c("Peptide", "# Methods Pass", "# Methods Tested", "Methods",
                    "Total Allele Hits", "Best Min %", "EL # Moderate", "EL Min %", "Cluster"),
      digits    = 2
    ) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "bordered"), full_width = FALSE) %>%
    row_spec(0, bold = TRUE, background = "#2166ac", color = "white")
}
# ```
