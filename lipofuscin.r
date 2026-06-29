library(ggplot2)
library(patchwork)
library(irr)

raw <- read.delim("Dog Lipofuscin Quantifications - Dog.tsv",
                  sep = "\t", header = TRUE, check.names = FALSE,
                  na.strings = c("", "NA"), strip.white = TRUE)

# Clean column names: strip whitespace and control characters
colnames(raw) <- trimws(gsub("[[:cntrl:]]", "", colnames(raw)))

# Helper: find column by partial match
col_match <- function(pattern) {
  idx <- grep(pattern, colnames(raw), ignore.case = TRUE)
  if (length(idx) == 0) stop(paste("No column matching:", pattern))
  colnames(raw)[idx[1]]
}

# Keep rows with a Donor ID
id_col <- col_match("^Donor ID$")
raw <- raw[!is.na(raw[[id_col]]) & raw[[id_col]] != "", ]

# Parse age to numeric years
parse_age <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\+", "", x)
  x <- gsub("\\s*(years?|yr)\\s*$", "", x, ignore.case = TRUE)
  suppressWarnings(as.numeric(x))
}

# Parse lipofuscin % (strip "%" sign)
parse_pct <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("%", "", x)
  suppressWarnings(as.numeric(x))
}

df <- data.frame(
  donor_id     = raw[[col_match("^Donor ID$")]],
  age          = parse_age(raw[[col_match("^Age$")]]),
  cds          = trimws(as.character(raw[[col_match("Signs of confusion")]])),
  lipo_pct_T   = parse_pct(raw[[col_match("Lipofuscin %.*T")]]),
  lipo_pct_H   = parse_pct(raw[[col_match("Lipofuscin %.*H")]]),
  lipo_count_T = suppressWarnings(as.numeric(raw[[col_match("annotations.*T")]])),
  lipo_count_H = suppressWarnings(as.numeric(raw[[col_match("annotations.*H")]]))
)

# Binary CDS variable (exclude Unknown and NA)
df$cds_binary <- ifelse(grepl("Yes", df$cds, ignore.case = TRUE), "Yes",
                        ifelse(df$cds == "No", "No", NA))
df$cds_binary <- factor(df$cds_binary, levels = c("No", "Yes"))

# Colors
col_no  <- "#4A90D9"  # blue for No CDS
col_yes <- "#D94A4A"  # red for Yes CDS
fill_no  <- "#A8CBF0"
fill_yes <- "#F0A8A8"


# --- Helpers ---

# Format p-value: scientific notation if < 0.001, otherwise 3 decimals
fmt_p <- function(p) {
  if (p < 0.001) {
    formatC(p, format = "e", digits = 1)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

# Spearman annotation string
spearman_label <- function(x, y) {
  ok <- complete.cases(x, y)
  test <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  paste0("rho == ", formatC(test$estimate, format = "f", digits = 2),
         "~~p == ", fmt_p(test$p.value),
         "~~n == ", sum(ok))
}

# Mann-Whitney annotation string
mw_label <- function(data, col) {
  no  <- data[[col]][data$cds_binary == "No"]
  yes <- data[[col]][data$cds_binary == "Yes"]
  test <- wilcox.test(no, yes, exact = FALSE)
  paste0("U == ", formatC(test$statistic, format = "f", digits = 0),
         "~~p == ", fmt_p(test$p.value))
}


# ============================================================
# FIGURE 1: Age vs Lipofuscin (4 panels, 2x2)
#            Red dashed regression line, blue CI band
# ============================================================

p_age_pctT <- ggplot(df, aes(x = age, y = lipo_pct_T)) +
  geom_smooth(method = "lm", se = TRUE, fill = fill_no, color = col_yes,
              linetype = "dashed", linewidth = 0.8, alpha = 0.3) +
  geom_point(size = 2.5, alpha = 0.8, color = col_no) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
           label = spearman_label(df$age, df$lipo_pct_T), parse = TRUE) +
  labs(x = "Age (years)", y = "Lipofuscin density (thalamus)") +
  theme_bw(base_size = 11)

p_age_pctH <- ggplot(df, aes(x = age, y = lipo_pct_H)) +
  geom_smooth(method = "lm", se = TRUE, fill = fill_no, color = col_yes,
              linetype = "dashed", linewidth = 0.8, alpha = 0.3) +
  geom_point(size = 2.5, alpha = 0.8, color = col_no) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
           label = spearman_label(df$age, df$lipo_pct_H), parse = TRUE) +
  labs(x = "Age (years)", y = "Lipofuscin density (hippocampal)") +
  theme_bw(base_size = 11)

p_age_cntT <- ggplot(df, aes(x = age, y = lipo_count_T)) +
  geom_smooth(method = "lm", se = TRUE, fill = fill_no, color = col_yes,
              linetype = "dashed", linewidth = 0.8, alpha = 0.3) +
  geom_point(size = 2.5, alpha = 0.8, color = col_no) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
           label = spearman_label(df$age, df$lipo_count_T), parse = TRUE) +
  labs(x = "Age (years)", y = "Lipofuscin objects (thalamus)") +
  theme_bw(base_size = 11)

p_age_cntH <- ggplot(df, aes(x = age, y = lipo_count_H)) +
  geom_smooth(method = "lm", se = TRUE, fill = fill_no, color = col_yes,
              linetype = "dashed", linewidth = 0.8, alpha = 0.3) +
  geom_point(size = 2.5, alpha = 0.8, color = col_no) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
           label = spearman_label(df$age, df$lipo_count_H), parse = TRUE) +
  labs(x = "Age (years)", y = "Lipofuscin objects (hippocampal)") +
  theme_bw(base_size = 11)

fig1 <- (p_age_pctT | p_age_pctH) / (p_age_cntT | p_age_cntH) +
  plot_annotation(tag_levels = "a")

print(fig1)

ggsave("age_lipofuscin.png", plot = fig1, width = 8, height = 6, dpi = 600)

# ============================================================
# FIGURE 2: CDS (binary) vs Lipofuscin (4 panels, 2x2)
#            Blue = No CDS, Red = Yes CDS
#            Mean crossbar + jitter, colored by group
# ============================================================

df_cds <- df[!is.na(df$cds_binary), ]

make_cds_plot <- function(data, col, ylab) {
  ggplot(data, aes(x = cds_binary, y = .data[[col]])) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.5,
                 aes(color = cds_binary)) +
    geom_jitter(width = 0.08, size = 2.5, alpha = 0.8,
                aes(color = cds_binary)) +
    scale_color_manual(values = c("No" = col_no, "Yes" = col_yes)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
             label = mw_label(data, col), parse = TRUE) +
    labs(x = "Cognitive dysfunction", y = ylab) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
}

p_cds_pctT <- make_cds_plot(df_cds, "lipo_pct_T", "Lipofuscin density (thalamus)")
p_cds_pctH <- make_cds_plot(df_cds, "lipo_pct_H", "Lipofuscin density (hippocampal)")
p_cds_cntT <- make_cds_plot(df_cds, "lipo_count_T", "Lipofuscin objects (thalamus)")
p_cds_cntH <- make_cds_plot(df_cds, "lipo_count_H", "Lipofuscin objects (hippocampal)")

fig2 <- (p_cds_pctT | p_cds_pctH) / (p_cds_cntT | p_cds_cntH) +
  plot_annotation(tag_levels = "a")

print(fig2)

ggsave("cds_lipofuscin.png", plot = fig2, width = 8, height = 6, dpi = 600)

# ============================================================
# INTER-RATER RELIABILITY: Blood clearance grades (Region H, T)
#            Quadratic-weighted Cohen's kappa (ordinal, 2 raters)
# ============================================================

r1_H <- as.numeric(raw[[col_match("Rater 1 Blood Clearance - Region H")]])
r2_H <- as.numeric(raw[[col_match("Rater 2 Blood Clearance - Region H")]])
r1_T <- as.numeric(raw[[col_match("Rater 1 Blood Clearance - Region T")]])
r2_T <- as.numeric(raw[[col_match("Rater 2 Blood Clearance - Region T")]])

# Pooled across both regions
pooled <- data.frame(rater1 = c(r1_H, r1_T), rater2 = c(r2_H, r2_T))
pooled <- pooled[complete.cases(pooled), ]
k_pooled <- kappa2(pooled, weight = "squared")

cat("Pooled:\n"); print(k_pooled)
