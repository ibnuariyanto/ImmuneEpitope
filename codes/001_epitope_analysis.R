# ============================================================
# MHC Class I Epitope Analysis Pipeline
# Publication-grade | NetMHCpan Output
# Author: [Your Name] | Date: 2026
# ============================================================

# Install required packages (run once)
pkgs <- c("tidyverse", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
          "viridis", "gridExtra", "knitr", "kableExtra", "writexl",
          "VennDiagram", "UpSetR", "corrplot", "scales", "patchwork")
install.packages(setdiff(pkgs, rownames(installed.packages())))

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(knitr)
library(kableExtra)
library(writexl)
library(VennDiagram)
library(UpSetR)
library(corrplot)
library(scales)
library(patchwork)

# ── Load Data ─────────────────────────────────────────────────
df_all    <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class1/peptide_table_40dba94e.csv",na.string=c("","NA"),header=T, sep =",")
df_pat    <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class1/MHCI_peptide_table_e04f7760_sort.csv",na.string=c("","NA"),header=T, sep =",")
dist_all  <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class1/allele_distances_table_40dba94e.csv",na.string=c("","NA"),header=T, sep =",")
dist_pat  <- read.csv("D:/06_RCoding/ImmuneEpitope/data/hbov_class1/MHCI_allele_distances_table_e04f7760_sort.csv",na.string=c("","NA"),header=T, sep =",")

# Standardize column names
df_all <- df_all %>% rename_with(tolower) %>%
  rename(peptide = peptide,
         allele  = allele,
         el_pct  = `netmhcpan_el.percentile`,
         el_score= `netmhcpan_el.score`,
         median_rank = `median.binding.percentile`)

df_pat <- df_pat %>% rename_with(tolower) %>%
  rename(peptide = peptide,
         allele  = allele,
         el_pct  = `netmhcpan_el.percentile`,
         el_score= `netmhcpan_el.score`,
         median_rank = `median.binding.percentile`)

# ============================================================
# 2. BINDING CLASSIFICATION (NetMHCpan Standard)
# ============================================================

classify_binding <- function(df) {
  df %>%
    mutate(binding_class = case_when(
      el_pct <= 0.5  ~ "Strong Binder (SB)",
      el_pct <= 2.0  ~ "Weak Binder (WB)",
      TRUE           ~ "Non-Binder (NB)"
    ),
    binding_class = factor(binding_class,
                           levels = c("Strong Binder (SB)", "Weak Binder (WB)", "Non-Binder (NB)")))
}

df_all <- classify_binding(df_all)
df_pat <- classify_binding(df_pat)

# Summary table
bind_summary <- df_all %>%
  count(binding_class) %>%
  mutate(pct = round(n / sum(n) * 100, 2))
print(bind_summary)

# ============================================================
# 3. PROMISCUITY SCORING
# SB = 2pts, WB = 1pt | Promiscuous if >= 3 alleles (SB+WB)
# ============================================================

promiscuity <- df_all %>%
  filter(binding_class != "Non-Binder (NB)") %>%
  group_by(peptide) %>%
  summarise(
    n_alleles_total   = n_distinct(allele),
    n_SB              = sum(binding_class == "Strong Binder (SB)"),
    n_WB              = sum(binding_class == "Weak Binder (WB)"),
    promiscuity_score = n_SB * 2 + n_WB,
    best_el_pct       = min(el_pct),
    median_el_pct     = median(el_pct),
    alleles_bound     = paste(sort(unique(allele)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    is_promiscuous = n_alleles_total >= 3,
    tier = case_when(
      n_SB >= 3                         ~ "Tier 1 - Highly Promiscuous SB",
      n_SB >= 1 & n_alleles_total >= 3  ~ "Tier 1 - Promiscuous",
      n_alleles_total >= 2              ~ "Tier 2 - Moderately Promiscuous",
      TRUE                              ~ "Tier 3 - Single Allele"
    )
  ) %>%
  arrange(desc(promiscuity_score), best_el_pct)

# Top candidates
top_candidates <- promiscuity %>% filter(tier %in% c("Tier 1 - Highly Promiscuous SB",
                                                     "Tier 1 - Promiscuous"))
cat("Top epitope candidates:", nrow(top_candidates), "\n")


p1 <- df_all %>%
  count(allele, binding_class) %>%
  ggplot(aes(x = reorder(allele, -n), y = n, fill = binding_class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Strong Binder (SB)" = "#d73027",
                               "Weak Binder (WB)"   = "#fc8d59",
                               "Non-Binder (NB)"    = "#e0e0e0")) +
  labs(title = "MHC Class I Binding Distribution per Allele",
       x = "HLA Allele", y = "Number of Peptides",
       fill = "Binding Class") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold"))
p1
ggsave("Fig1_binding_distribution.pdf", p1, width = 12, height = 6, dpi = 300)
ggsave("Fig1_binding_distribution.png", p1, width = 12, height = 6, dpi = 300)

p2 <- promiscuity %>%
  ggplot(aes(x = promiscuity_score, fill = tier)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.85) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Peptide Promiscuity Score Distribution",
       x = "Promiscuity Score (SB=2, WB=1)",
       y = "Number of Peptides", fill = "Tier") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
p2
ggsave("Fig2_promiscuity_distribution.pdf", p2, width = 8, height = 5, dpi = 300)

top30 <- top_candidates$peptide[1:min(30, nrow(top_candidates))]

heat_mat <- df_all %>%
  filter(peptide %in% top30) %>%
  select(peptide, allele, el_pct) %>%
  pivot_wider(names_from = allele, values_from = el_pct,
              values_fn = min, values_fill = 100) %>%
  column_to_rownames("peptide") %>%
  as.matrix()

# Cap at 5% for visualization clarity
heat_mat_capped <- pmin(heat_mat, 5)

pheatmap(heat_mat_capped,
         color = colorRampPalette(c("#d73027","#fc8d59","#fee090","#e0f3f8","#f0f0f0"))(100),
         breaks = seq(0, 5, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize_row = 7, fontsize_col = 8,
         main = "EL%ile Heatmap: Top Epitope Candidates x HLA Alleles (capped at 5%)",
         filename = "Fig3_binding_heatmap.pdf",
         width = 14, height = 10)

p4 <- promiscuity %>%
  filter(n_alleles_total >= 2) %>%
  ggplot(aes(x = best_el_pct, y = n_alleles_total,
             size = promiscuity_score, color = tier)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = . %>% filter(tier %in% c("Tier 1 - Highly Promiscuous SB",
                                    "Tier 1 - Promiscuous")) %>%
      slice_min(best_el_pct, n = 15),
    aes(label = peptide), size = 3, max.overlaps = 20
  ) +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_color_brewer(palette = "Set1") +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Epitope Candidates: Binding Affinity vs. Promiscuity",
       x = "Best EL Percentile (lower = stronger)",
       y = "Number of Alleles Bound",
       size = "Promiscuity Score", color = "Tier") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
p4

ggsave("Fig4_bubble_promiscuity.pdf", p4, width = 10, height = 7, dpi = 300)
