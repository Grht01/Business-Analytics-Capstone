###############################################################################
#  Stance Analysis of Public Climate and Sustainability Discussion
#  Complete estimation script for Chapter IV
#
#  Author : Golam Rahat (ID 1112310093)
#  Course : BSA 4401 Business Analytics Capstone
#  R      : tested against R 4.5.2, marginaleffects 0.32.0, lme4 1.1-38
#
#  INPUT   Stance_sample_data.csv
#          This file IS the 250,000-tweet stratified draw described in
#          Section 3.3. The seed-42 draw from the full CCTD happened when the
#          file was created and is NOT repeated here. This script starts
#          downstream, applying only the three complete-case restrictions that
#          take the draw to the 196,406-tweet analytic sample.
#
#  OUTPUT  results/  fourteen CSV files plus a session log
#
#  HOW TO RUN
#          Put this file and Stance_sample_data.csv in the same folder, set
#          that folder as the working directory, then:  source("analysis.R")
#          Runtime is roughly 15 to 45 minutes; Equation (4) is the slow step.
#
#  FIXES INCORPORATED (each of these broke an earlier run)
#    1. broom and broom.mixed removed. They pull a purrr version dependency
#       that fails with "namespace purrr 1.2.0 is already loaded, but
#       >= 1.2.1 is required". All model output is extracted with base R.
#    2. McFadden computed against an explicitly fitted null model. Using
#       update(m, . ~ 1) inside a helper function raises a scoping error.
#    3. Equation (4) uses nAGQ = 0 and is wrapped in try(), so a failure or a
#       long fit cannot kill the rest of the run.
#    4. Average marginal effect contrasts use hypothesis = "b2 - b1 = 0".
#       The string "pairwise" is rejected by marginaleffects 0.32.0, and the
#       ~pairwise shorthand returns differently named columns across versions
#       and can silently produce zero rows.
#    5. Every write is preceded by a row-count check, so the script stops with
#       a clear message rather than writing an empty file.
#    6. Table 4.4 is written to disk. An earlier version only printed it.
#    7. Debug traps are cleared at the top, so a stray breakpoint or
#       options(error = recover) cannot trap the run in the browser.
###############################################################################

## ==========================================================================
## 0. Environment hygiene
## ==========================================================================
options(error = NULL)            # no recover()/browser() error handler
suppressWarnings(try(undebug(data.frame), silent = TRUE))
suppressWarnings(try(undebug(lapply),     silent = TRUE))

pkgs <- c("tidyverse", "lme4", "sandwich", "lmtest", "car", "pROC",
          "nnet", "marginaleffects")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# Seed governs only the permutation placebo in Section 9. It does not
# reproduce the original draw, which is already fixed in the input file.
set.seed(42)
dir.create("results", showWarnings = FALSE)

ok   <- function(msg) cat("[ok] ", msg, "\n", sep = "")
step <- function(msg) cat("\n[..] ", msg, "\n", sep = "")

## Base-R replacement for broom::tidy(). Handles glm, glmerMod and coeftest.
tidy_base <- function(m) {
  cf <- if (inherits(m, "glmerMod")) summary(m)$coefficients
        else if (inherits(m, "coeftest")) unclass(m)
        else summary(m)$coefficients
  tibble(term = rownames(cf), estimate = cf[, 1], std_error = cf[, 2],
         statistic = cf[, 3], p_value = cf[, 4], odds_ratio = exp(cf[, 1]))
}

## Guarded write: refuses to create an empty file.
safe_write <- function(x, path) {
  if (is.null(x) || nrow(x) == 0)
    stop("Refusing to write an empty file: ", path,
         ". Inspect the object above.", call. = FALSE)
  write_csv(x, path); ok(paste0(path, "  (", nrow(x), " rows)"))
}

## ==========================================================================
## 1. Load, and construct the analytic sample (Section 3.3)
## ==========================================================================
step("Section 1: loading data and constructing the analytic sample")

raw <- read_csv("Stance_sample_data.csv", show_col_types = FALSE)
stopifnot(nrow(raw) == 250000)          # guards against loading a different file

cat("Provided stratified draw N =", nrow(raw), "\n")
cat("  neutral stance       :", sum(raw$stance  == "neutral"),     "\n")
cat("  undefined gender     :", sum(raw$gender  == "undefined"),   "\n")
cat("  unclassified economy :", sum(raw$economy == "Unclassified"), "\n")

df <- raw %>%
  filter(stance  %in% c("believer", "denier"),
         gender  %in% c("male", "female"),
         economy %in% c("Developed (High income)", "Developing (LMIC)")) %>%
  mutate(
    Y  = as.integer(stance == "believer"),                  # believer = 1
    # Centering is on the ANALYTIC sample because filter() precedes mutate().
    # Centering shifts the intercept only; slopes and interactions are
    # unaffected by the choice of centering base.
    X1 = sentiment - mean(sentiment, na.rm = TRUE),          # mean-centered
    X2 = as.integer(aggressiveness == "aggressive"),         # aggressive = 1
    X3 = temperature_avg - mean(temperature_avg, na.rm = TRUE),
    Z1 = as.integer(gender  == "male"),                      # male = 1
    M  = as.integer(economy == "Developed (High income)"),   # developed = 1
    country = factor(country),
    year    = as.integer(substr(as.character(created_at), 1, 4))
  )

stopifnot(nrow(df) == 196406)           # must match Section 3.3
sd_X1 <- sd(df$X1); sd_X3 <- sd(df$X3)

cat("\nAnalytic N =", nrow(df),
    "| countries =", nlevels(droplevels(df$country)),
    "| developed =", sum(df$M == 1), "| developing =", sum(df$M == 0), "\n")
cat("SD(X1) =", round(sd_X1, 5), " SD(X3) =", round(sd_X3, 5), "\n")

safe_write(
  tibble(step = c("Stratified draw provided", "Less neutral stance",
                  "Less undefined gender", "Less unclassified economy",
                  "Analytic sample"),
         n = c(250000L,
               250000L - sum(raw$stance == "neutral"),
               NA_integer_, NA_integer_, nrow(df)),
         note = c("input file", "believer or denier only",
                  "male or female only", "valid World Bank income group",
                  "complete cases on all modelled variables")),
  "results/table_4_0_sample_construction.csv")
ok("Section 1 complete")

## ==========================================================================
## 2. Table 4.1 descriptive statistics
## ==========================================================================
step("Section 2: descriptive statistics")

desc <- df %>%
  select(Y, X1, X2, X3, Z1, M, sentiment, temperature_avg) %>%
  pivot_longer(everything(), names_to = "var", values_to = "v") %>%
  group_by(var) %>%
  summarise(N = sum(!is.na(v)), mean = mean(v, na.rm = TRUE),
            sd = sd(v, na.rm = TRUE), min = min(v, na.rm = TRUE),
            p25 = quantile(v, .25, na.rm = TRUE),
            median = median(v, na.rm = TRUE),
            p75 = quantile(v, .75, na.rm = TRUE),
            max = max(v, na.rm = TRUE), .groups = "drop")
print(as.data.frame(desc), digits = 4)
safe_write(desc, "results/table_4_1_descriptives.csv")

## ==========================================================================
## 3. Table 4.2 correlations and variance inflation factors
## ==========================================================================
step("Section 3: correlations and VIF")

cormat <- cor(df %>% select(Y, X1, X2, X3, Z1, M), use = "complete.obs")
print(round(cormat, 4))
write.csv(round(cormat, 4), "results/table_4_2_correlations.csv")
ok("results/table_4_2_correlations.csv")

vif_mod <- glm(Y ~ X1 + X2 + X3 + Z1 + M, data = df, family = binomial())
vif_vals <- car::vif(vif_mod)
print(round(vif_vals, 4))
safe_write(tibble(term = names(vif_vals), VIF = as.numeric(vif_vals)),
           "results/table_4_2_vif.csv")

## ==========================================================================
## 4. Table 4.4 nested estimation, Equations (1) to (3)
## ==========================================================================
step("Section 4: nested logistic estimation")

f1 <- Y ~ X1 + X2 + X3
f2 <- Y ~ X1 + X2 + X3 + Z1 + M
f3 <- Y ~ X1 + X2 + X3 + Z1 + M + X1:M + X2:M + X3:M

eq1 <- glm(f1, data = df, family = binomial())
eq2 <- glm(f2, data = df, family = binomial())
eq3 <- glm(f3, data = df, family = binomial())

# Country-clustered standard errors (Cameron & Miller, 2015; Abadie et al., 2023)
clust  <- function(m, d = df) coeftest(m, vcov = vcovCL(m, cluster = d$country))
eq3_cl <- clust(eq3)
print(eq3_cl)

# McFadden against an explicitly fitted null model. Do not use update() here:
# calling update(m, . ~ 1) inside a helper raises an environment-scoping error.
null_mod <- glm(Y ~ 1, data = df, family = binomial())
ll0      <- as.numeric(logLik(null_mod))
mcfadden <- function(m) 1 - as.numeric(logLik(m)) / ll0

cat("\nMcFadden pseudo R2:  Eq1 =", round(mcfadden(eq1), 4),
    " Eq2 =", round(mcfadden(eq2), 4),
    " Eq3 =", round(mcfadden(eq3), 4), "\n")

tab44 <- bind_rows(
  tidy_base(eq1)    %>% mutate(column = "(1) Eq. (1)"),
  tidy_base(eq2)    %>% mutate(column = "(2) Eq. (2)"),
  tidy_base(eq3)    %>% mutate(column = "(3) Eq. (3) naive"),
  tidy_base(eq3_cl) %>% mutate(column = "(4) Eq. (3) clustered")
) %>% select(column, term, estimate, std_error, odds_ratio, p_value)
safe_write(tab44, "results/table_4_4_nested.csv")

safe_write(tibble(spec = c("Eq1", "Eq2", "Eq3"),
                  mcfadden = c(mcfadden(eq1), mcfadden(eq2), mcfadden(eq3)),
                  n = c(nobs(eq1), nobs(eq2), nobs(eq3))),
           "results/table_4_4_fit.csv")

## ==========================================================================
## 5. Equation (4) country random intercept, the primary specification
## ==========================================================================
step("Section 5: Equation (4) random intercept. This is the slow step.")

# Set to FALSE to skip. The clustered and fixed effects specifications already
# provide inference that respects the country structure.
RUN_GLMER <- TRUE

if (RUN_GLMER) {
  eq4 <- try(glmer(Y ~ X1 + X2 + X3 + Z1 + M + X1:M + X2:M + X3:M + (1 | country),
                   data = df, family = binomial(), nAGQ = 0,
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 2e5))),
             silent = TRUE)

  if (inherits(eq4, "try-error")) {
    cat("[!!] glmer failed. Continuing with clustered and FE specifications.\n")
    cat(as.character(eq4), "\n")
  } else {
    print(summary(eq4))
    eq4_tab <- tidy_base(eq4) %>%
      mutate(ci_low  = exp(estimate - 1.96 * std_error),
             ci_high = exp(estimate + 1.96 * std_error))
    safe_write(eq4_tab, "results/table_4_4_eq4_random_intercept.csv")

    vc  <- as.data.frame(lme4::VarCorr(eq4))
    s2u <- vc$vcov[1]
    safe_write(tibble(sigma2_u = s2u,
                      icc = s2u / (s2u + pi^2 / 3),   # latent-scale ICC
                      n_countries = as.integer(lme4::ngrps(eq4)),
                      n_obs = nobs(eq4)),
               "results/table_4_4_eq4_variance.csv")
    cat("Country intercept variance =", round(s2u, 4),
        " latent-scale ICC =", round(s2u / (s2u + pi^2 / 3), 4), "\n")
  }
}

## ==========================================================================
## 6. Fit and discrimination (Section 3.7)
## ==========================================================================
step("Section 6: fit and discrimination")

p_hat   <- predict(eq3, type = "response")
roc_obj <- pROC::roc(df$Y, p_hat, quiet = TRUE)
yhat    <- as.integer(p_hat >= 0.5)

# Denier is the minority class, so precision and recall are reported for it.
tp <- sum(yhat == 0 & df$Y == 0)
fp <- sum(yhat == 0 & df$Y == 1)
fn <- sum(yhat == 1 & df$Y == 0)

fit_tab <- tibble(
  auc                = as.numeric(pROC::auc(roc_obj)),
  accuracy           = mean(yhat == df$Y),
  denier_precision   = ifelse(tp + fp > 0, tp / (tp + fp), 0),
  denier_recall      = ifelse(tp + fn > 0, tp / (tp + fn), 0),
  believer_base_rate = mean(df$Y)
)
print(as.data.frame(fit_tab), digits = 4)
safe_write(fit_tab, "results/table_4_3_fit_discrimination.csv")

## ==========================================================================
## 7. Average marginal effects by economy, and the moderation contrasts
## ==========================================================================
step("Section 7: average marginal effects (H4 to H6)")

ame <- avg_slopes(eq3, variables = c("X1", "X2", "X3"), by = "M",
                  vcov = ~country)
print(ame)
safe_write(as_tibble(ame), "results/table_4_ame_by_economy.csv")

# Formal moderation test: developed minus developing, one cue at a time.
#
# hypothesis = "b2 - b1 = 0" is the explicit string form. With by = "M" there
# are exactly two estimates, b1 (M = 0, developing) and b2 (M = 1, developed).
# Do NOT use "pairwise" (rejected outright) or ~pairwise (returns columns whose
# names vary by version and can silently yield zero rows). Requesting all three
# variables at once would also return 15 cross-comparisons, most meaningless.
labs     <- c(X1 = "Sentiment", X2 = "Aggressiveness", X3 = "Temperature anomaly")
hyp_id   <- c(X1 = "H4",        X2 = "H5",             X3 = "H6")
sd_scale <- c(X1 = sd_X1,       X2 = 1,                X3 = sd_X3)

ame_diff <- bind_rows(lapply(c("X1", "X2", "X3"), function(v) {
  d <- avg_slopes(eq3, variables = v, by = "M",
                  hypothesis = "b2 - b1 = 0", vcov = ~country)
  if (nrow(d) == 0)
    stop("avg_slopes returned no rows for ", v, ". Inspect the object.",
         call. = FALSE)
  s <- unname(sd_scale[v])
  tibble(hypothesis = unname(hyp_id[v]),
         cue        = unname(labs[v]),
         scaling    = if (v == "X2") "discrete 0 to 1" else "per 1 SD",
         diff_pp    = d$estimate[1]  * s * 100,   # developed minus developing
         se_pp      = d$std.error[1] * s * 100,
         statistic  = d$statistic[1],
         p_value    = d$p.value[1],
         ci_low_pp  = d$conf.low[1]  * s * 100,
         ci_high_pp = d$conf.high[1] * s * 100)
}))
print(as.data.frame(ame_diff), digits = 4)
safe_write(ame_diff, "results/table_4_ame_difference.csv")

## ==========================================================================
## 8. Table 4.5 subgroup estimation (Section 4.5)
## ==========================================================================
step("Section 8: subgroup estimation by economy")

sub_fit <- function(g) {
  d  <- df %>% filter(M == g)
  m  <- glm(Y ~ X1 + X2 + X3 + Z1, data = d, family = binomial())
  ct <- coeftest(m, vcov = vcovCL(m, cluster = d$country))
  tibble(group = if (g == 1) "Developed" else "Developing",
         term = rownames(ct), estimate = ct[, 1], se = ct[, 2],
         p = ct[, 4], odds_ratio = exp(ct[, 1]),
         n = nrow(d), n_countries = nlevels(droplevels(d$country)),
         mcfadden = 1 - as.numeric(logLik(m)) /
                        as.numeric(logLik(glm(Y ~ 1, data = d, family = binomial()))))
}
subgroups <- bind_rows(sub_fit(1), sub_fit(0))
print(as.data.frame(subgroups), digits = 4)
safe_write(subgroups, "results/table_4_5_subgroups.csv")

## ==========================================================================
## 9. Section 4.6 robustness, in the order declared in Section 3.7
## ==========================================================================
step("Section 9: robustness checks R1 to R5")

## R1 alternative outcome coding: multinomial retaining the neutral category
mn <- raw %>%
  filter(gender  %in% c("male", "female"),
         economy %in% c("Developed (High income)", "Developing (LMIC)")) %>%
  mutate(X1 = sentiment - mean(sentiment),
         X2 = as.integer(aggressiveness == "aggressive"),
         X3 = temperature_avg - mean(temperature_avg),
         Z1 = as.integer(gender == "male"),
         M  = as.integer(economy == "Developed (High income)"),
         S  = relevel(factor(stance), ref = "denier"))
r1 <- nnet::multinom(S ~ X1 + X2 + X3 + Z1 + M + X1:M + X2:M + X3:M,
                     data = mn, trace = FALSE)
r1_cf <- summary(r1)$coefficients; r1_se <- summary(r1)$standard.errors
r1_tab <- bind_rows(lapply(rownames(r1_cf), function(lv)
  tibble(outcome = lv, term = colnames(r1_cf),
         estimate = r1_cf[lv, ], std_error = r1_se[lv, ],
         odds_ratio = exp(r1_cf[lv, ]),
         p_value = 2 * pnorm(-abs(r1_cf[lv, ] / r1_se[lv, ])))))
safe_write(r1_tab, "results/table_4_6_r1_multinomial.csv")
cat("R1 multinomial N =", nrow(mn), "\n")

## R2 alternative estimator: country fixed effects, full sample
eq_fe <- glm(Y ~ X1 + X2 + X3 + Z1 + X1:M + X2:M + X3:M + country,
             data = df, family = binomial())
safe_write(bind_rows(
  tidy_base(eq3)   %>% mutate(spec = "Eq3 pooled"),
  tidy_base(eq_fe) %>% filter(!grepl("^country", term)) %>%
                       mutate(spec = "Country FE")
), "results/table_4_6_r2_estimators.csv")

## R3 alternative sample window: temporal subsamples
r3 <- bind_rows(lapply(list(c(2006, 2012), c(2013, 2016), c(2017, 2019)),
  function(w) {
    d  <- df %>% filter(year >= w[1], year <= w[2])
    m  <- glm(f3, data = d, family = binomial())
    ct <- coeftest(m, vcov = vcovCL(m, cluster = d$country))
    tibble(window = paste0(w[1], "-", w[2]), term = rownames(ct),
           odds_ratio = exp(ct[, 1]), p = ct[, 4], n = nrow(d))
  }))
safe_write(r3, "results/table_4_6_r3_temporal.csv")

## R4 alternative moderator coding: high plus upper-middle against the rest
df <- df %>% mutate(M_alt = as.integer(grepl("High|Upper", income_group,
                                             ignore.case = TRUE)))
r4 <- glm(Y ~ X1 + X2 + X3 + Z1 + M_alt + X1:M_alt + X2:M_alt + X3:M_alt,
          data = df, family = binomial())
safe_write(tidy_base(clust(r4)), "results/table_4_6_r4_alt_moderator.csv")

## R5 falsification: randomly permuted moderator, interactions should be null
df_p <- df %>% mutate(M = sample(M))
r5   <- glm(f3, data = df_p, family = binomial())
r5_ct <- coeftest(r5, vcov = vcovCL(r5, cluster = df_p$country))
keep  <- intersect(c("X1:M", "X2:M", "X3:M"), rownames(r5_ct))
print(r5_ct[keep, , drop = FALSE])
safe_write(tibble(term = keep, estimate = r5_ct[keep, 1],
                  se = r5_ct[keep, 2], p = r5_ct[keep, 4]),
           "results/table_4_6_r5_placebo.csv")

## ==========================================================================
## 10. Verification against the values reported in Chapter IV
## ==========================================================================
step("Section 10: verification against Chapter IV")

get_or <- function(tab, col, tm) {
  v <- tab$odds_ratio[tab$column == col & tab$term == tm]
  if (length(v)) round(v, 3) else NA_real_
}
check <- tibble(
  quantity = c("Analytic N", "Countries", "Eq3 X1 OR", "Eq3 X2 OR",
               "Eq3 X1:M OR", "Eq3 X2:M OR", "Eq3 X3:M OR", "McFadden Eq3"),
  obtained = c(nrow(df), nlevels(droplevels(df$country)),
               get_or(tab44, "(3) Eq. (3) naive", "X1"),
               get_or(tab44, "(3) Eq. (3) naive", "X2"),
               get_or(tab44, "(3) Eq. (3) naive", "X1:M"),
               get_or(tab44, "(3) Eq. (3) naive", "X2:M"),
               get_or(tab44, "(3) Eq. (3) naive", "X3:M"),
               round(mcfadden(eq3), 4)),
  chapter  = c(196406, 142, 6.632, 0.414, 0.534, 1.517, 1.023, 0.0675)
) %>% mutate(match = ifelse(abs(obtained - chapter) < 0.005, "OK", "CHECK"))
print(as.data.frame(check))
safe_write(check, "results/table_4_verification.csv")

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

cat("\n==========================================================\n")
cat("All sections complete. Files written to results/\n")
print(list.files("results"))
cat("==========================================================\n")
