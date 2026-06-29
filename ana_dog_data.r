library(readr)
library(dplyr)
library(stringr)

raw <- read_tsv("Dog ratings - Sheet1.tsv", col_types = cols(.default = col_character()))

# the first gross column header carries the scoring legend; trim it to its prefix
names(raw) <- sub("^(ACA R Gross).*", "\\1", names(raw))
names(raw) <- trimws(names(raw))   # strip stray trailing CR/space (e.g. Region T)

mins <- function(x) {
  x <- str_to_lower(str_trim(x))
  if (is.na(x) || str_detect(x, "tod|record")) return(NA_real_)
  m <- as.numeric(str_match(x, "(\\d+)\\s*m")[, 2])
  if (is.na(m) || m == 0) NA_real_ else m
}
press <- function(x) {
  x <- str_to_lower(str_trim(x))
  if (is.na(x) || x == "nr" || x == "") return(NA_real_)
  suppressWarnings(as.numeric(x))
}
litres <- function(x) {
  x <- str_to_lower(str_trim(x))
  if (is.na(x) || x == "nr" || x == "") return(NA_real_)
  suppressWarnings(as.numeric(str_replace(x, "\\s*l", "")))
}

df <- raw %>%
  rename(
    age = `Age (years)`, weight = `Weight (lbs)`,
    p_low_raw = `Lowest pressure recorded`, p_high_raw = `Highest pressure recorded`,
    dur_raw = `Duration of perfusion`, vol_raw = `Volume used (Liters)`,
    pmi_raw = `PMI (minutes)`
  ) %>%
  mutate(
    age       = suppressWarnings(as.numeric(age)),
    weight    = suppressWarnings(as.numeric(weight)),
    pmi_min   = vapply(pmi_raw,  mins,   numeric(1)),
    dur_min   = vapply(dur_raw,  mins,   numeric(1)),
    p_low     = vapply(p_low_raw,  press, numeric(1)),
    p_high    = vapply(p_high_raw, press, numeric(1)),
    volume_L  = vapply(vol_raw,  litres, numeric(1)),
    flow_per_lb = (volume_L * 1000) / dur_min / weight   # mL/min/lb
  )

# non-perfused donors: exclude from all perfusion parameters and quality scores,
# but keep their age, weight, and PMI
non_perfused <- str_detect(str_to_lower(coalesce(df$p_high_raw, "")), "not perfused|no perfusion") |
  str_detect(str_to_lower(coalesce(df$p_low_raw,  "")), "not perfused|no perfusion")

df <- df %>%
  mutate(
    p_low       = ifelse(non_perfused, NA_real_, p_low),
    p_high      = ifelse(non_perfused, NA_real_, p_high),
    dur_min     = ifelse(non_perfused, NA_real_, dur_min),
    volume_L    = ifelse(non_perfused, NA_real_, volume_L),
    flow_per_lb = ifelse(non_perfused, NA_real_, flow_per_lb)
  )

gross_cols <- c("ACA R Gross","MCA R Gross","PCA R Gross","CB R Gross",
                "ACA L Gross","MCA L Gross","PCA L Gross","CB L Gross")
ct_cols    <- c("ACA R CT","MCA R CT","PCA R CT","CB R CT",
                "ACA L CT","MCA L CT","PCA L CT","CB L CT")
hist_cols  <- c("Consensus Blood Clearance (0-3 Scale) - Region H",
                "Consensus Blood Clearance (0-3 Scale) - Region T")
df <- df %>% mutate(across(all_of(c(gross_cols, ct_cols, hist_cols)),
                           ~ suppressWarnings(as.numeric(.))))

# null out quality scores for non-perfused donors before computing composites
df <- df %>% mutate(across(all_of(c(gross_cols, ct_cols, hist_cols)),
                           ~ ifelse(non_perfused, NA_real_, .)))

row_mean_min <- function(m, min_n = 4)
  apply(m, 1, function(r) if (sum(!is.na(r)) >= min_n) mean(r, na.rm = TRUE) else NA_real_)

df <- df %>%
  mutate(
    gross_composite = row_mean_min(as.matrix(across(all_of(gross_cols)))),
    ct_composite    = row_mean_min(as.matrix(across(all_of(ct_cols)))),
    # histology score = mean of the two consensus region scores (0-3 scale,
    # like gross/CT); both regions are scored together so this needs both
    hist_composite  = row_mean_min(as.matrix(across(all_of(hist_cols))), min_n = 2)
  )

num_vars <- df %>%
  select(age, weight, pmi_min, dur_min, p_low, p_high, volume_L, flow_per_lb,
         gross_composite, ct_composite, hist_composite)

vlab <- c("Age","Weight","PMI (min)","Duration (min)","Lowest Pressure (PSI)",
          "Highest Pressure (PSI)","Volume (L)","Flow (mL/min/lb)",
          "Gross Score","CT Score","Histology Score")
spm_vars <- num_vars
colnames(spm_vars) <- vlab

stars <- function(p) if (is.na(p)) "" else if (p<.001) "***" else if (p<.01) "**" else if (p<.05) "*" else ""

panel_lower <- function(x, y, ...) {
  points(x, y, pch = 19, cex = 0.6, col = adjustcolor("#444444", alpha.f = 0.5))
  ok <- is.finite(x) & is.finite(y)
  p <- NA
  if (sum(ok) >= 4) { ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman")); p <- ct$p.value }
  if (sum(ok) >= 3) {
    f <- lm(y[ok] ~ x[ok]); xs <- range(x[ok]); ys <- coef(f)[1] + coef(f)[2]*xs
    s <- !is.na(p) && p < .05
    lines(xs, ys, col = if (s) "red" else "grey60", lwd = if (s) 1.6 else 1.1, lty = if (s) 1 else 2)
  }
  txt <- if (sum(ok) >= 4)
    paste0("rho=", format(ct$estimate, digits=2), stars(p),
           "\np=", format.pval(p, digits=2, eps=1e-4), "  n=", sum(ok))
  else paste0("n=", sum(ok))
  u <- par("usr"); tx <- u[1]+0.05*(u[2]-u[1]); ty <- u[4]-0.06*(u[4]-u[3])
  rect(tx, ty - strheight(txt, cex=0.85),
       tx + max(strwidth(strsplit(txt,"\n")[[1]], cex=0.85)), ty,
       col = adjustcolor("white", alpha.f=0.75), border = NA)
  text(tx, ty, txt, cex = 0.85, adj = c(0,1))
}

panel_upper <- function(x, y, ...) {
  ok <- is.finite(x) & is.finite(y); u <- par("usr")
  if (sum(ok) >= 4) {
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman")); r <- as.numeric(ct$estimate)
    it <- min(abs(r), 1)
    fill <- if (r >= 0) rgb(1-it, 1-it, 1) else rgb(1, 1-it, 1-it)
    rect(u[1], u[3], u[2], u[4], col = fill, border = NA)
    col <- if (it > 0.5) "white" else "black"
    text(mean(u[1:2]), mean(u[3:4]), paste0(format(r, digits=2), stars(ct$p.value)),
         cex = 1.8, font = 2, col = col)
    text(mean(u[1:2]), u[3]+0.12*(u[4]-u[3]), paste0("n=", sum(ok)), cex = 1.3, col = col)
  } else text(mean(u[1:2]), mean(u[3:4]), paste0("n=", sum(ok)), cex = 0.8)
}

panel_diag <- function(x, ...) {
  u <- par("usr"); par(usr = c(u[1:2], 0, 1.5)); on.exit(par(usr = u))
  xx <- x[is.finite(x)]; h <- hist(xx, plot = FALSE, breaks = 10)
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts/max(h$counts),
       col = "grey85", border = "white")
  d <- density(xx); lines(d$x, d$y/max(d$y), col = "grey30")
}

png("perfusion_scatterplot_matrix.png", width = 4600, height = 4600, res = 300)
pairs(spm_vars, labels = vlab,
      lower.panel = panel_lower, upper.panel = panel_upper, diag.panel = panel_diag,
      gap = 0.3, cex.labels = 0.9, cex.axis = 0.95, oma = c(4,4,4,4))
dev.off()

get_cor <- function(a, b) {
  x <- spm_vars[[a]]; y <- spm_vars[[b]]
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  cat(sprintf("%-22s vs %-22s rho=%.2f  p=%.4g  n=%d\n",
              a, b, ct$estimate, ct$p.value, sum(ok)))
}

get_cor("Gross Score", "CT Score")
get_cor("Histology Score", "Gross Score")
get_cor("Histology Score", "CT Score")
get_cor("Weight", "Gross Score")
get_cor("Flow (mL/min/lb)", "Gross Score")
get_cor("Age", "Weight")

