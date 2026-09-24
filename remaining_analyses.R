###############################################################################
#  REMAINING ANALYSES
#  Stance Analysis of Public Climate and Sustainability Discussion
#  Golam Rahat, ID 1112310093, BSA 4401
#
#  This script finishes the three things still outstanding after the main run.
#  It is STANDALONE: it rebuilds everything from the CSV, so it does not matter
#  what is or is not left in your R session. Start a fresh session and run it.
#
#     PART 1  Sample construction with real sequential counts
#             (replaces the two NA rows in table_4_0_sample_construction.csv)
#
#     PART 2  Placebo test, permuted at COUNTRY level
#             (replaces the broken observation-level output, which reported
#              X3:M at p = 1.7e-15, a spurious result)
#
#     PART 3  Sentiment-artifact sensitivity check
#             (replaces the "high confidence subsample" check promised in
#              Chapter III but impossible here: the dataset carries no
#              classifier confidence scores to threshold on)
#
#  HOW TO RUN
#     Put this file and Stance_sample_data.csv in the same folder, set that
#     folder as the working directory, then:   source("remaining_analyses.R")
#     Runtime is about two to four minutes. No mixed model is refitted.
###############################################################################

options(error = NULL)                      # no debugger traps
suppressWarnings(try(undebug(data.frame), silent = TRUE))

pkgs <- c("tidyverse", "sandwich", "lmtest")
miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(42)
dir.create("results", showWarnings = FALSE)
ok <- function(m) cat("[ok] ", m, "\n", sep = "")

safe_write <- function(x, path) {
  if (is.null(x) || nrow(x) == 0)
    stop("Refusing to write an empty file: ", path, call. = FALSE)
  write_csv(x, path); ok(paste0(path, "  (", nrow(x), " rows)"))
}

## Base-R tidier. broom is avoided deliberately: it pulls a purrr version
## dependency that has already broken one run in this project.
tidy_base <- function(m) {
  cf <- if (inherits(m, "coeftest")) unclass(m) else summary(m)$coefficients
  tibble(term = rownames(cf), estimate = cf[, 1], std_error = cf[, 2],
         statistic = cf[, 3], p_value = cf[, 4], odds_ratio = exp(cf[, 1]))
}

## ==========================================================================
## Load and rebuild the analytic sample
## ==========================================================================
cat("\n[..] Loading data\n")
raw <- read_csv("Stance_sample_data.csv", show_col_types = FALSE)
stopifnot(nrow(raw) == 250000)

df <- raw %>%
  filter(stance  %in% c("believer", "denier"),
         gender  %in% c("male", "female"),
         economy %in% c("Developed (High income)", "Developing (LMIC)")) %>%
  mutate(
    Y  = as.integer(stance == "believer"),
    X1 = sentiment - mean(sentiment, na.rm = TRUE),
    X2 = as.integer(aggressiveness == "aggressive"),
    X3 = temperature_avg - mean(temperature_avg, na.rm = TRUE),
    Z1 = as.integer(gender  == "male"),
    M  = as.integer(economy == "Developed (High income)"),
    country = factor(country)
  )
stopifnot(nrow(df) == 196406)
ok(paste0("Analytic sample rebuilt: ", nrow(df), " tweets, ",
          nlevels(droplevels(df$country)), " countries"))

FORM  <- Y ~ X1 + X2 + X3 + Z1 + M + X1:M + X2:M + X3:M
clust <- function(m, d) coeftest(m, vcov = vcovCL(m, cluster = d$country))

## ==========================================================================
## PART 1  Sample construction, real sequential counts
## ==========================================================================
cat("\n[..] PART 1: sample construction\n")

s1 <- raw %>% filter(stance %in% c("believer", "denier"))
s2 <- s1  %>% filter(gender %in% c("male", "female"))
s3 <- s2  %>% filter(economy %in% c("Developed (High income)", "Developing (LMIC)"))

sample_tab <- tibble(
  step = c("Stratified draw provided",
           "Less tweets labeled neutral on stance",
           "Less tweets with no confident gender label",
           "Less tweets with no World Bank income group",
           "Analytic sample"),
  n       = c(nrow(raw), nrow(s1), nrow(s2), nrow(s3), nrow(df)),
  removed = c(NA_integer_, nrow(raw) - nrow(s1), nrow(s1) - nrow(s2),
              nrow(s2) - nrow(s3), NA_integer_),
  note = c("input file", "believer or denier only", "male or female only",
           "valid World Bank income group",
           "complete cases on all modelled variables")
)
print(as.data.frame(sample_tab))
safe_write(sample_tab, "results/table_4_0_sample_construction.csv")

## country composition, for Chapter III Table 3.3
comp <- df %>%
  count(economy, country, name = "tweets") %>%
  group_by(economy) %>% arrange(desc(tweets), .by_group = TRUE) %>%
  mutate(rank = row_number(), share_of_group = tweets / sum(tweets)) %>%
  ungroup()
safe_write(comp, "results/table_3_3_country_composition.csv")

## ==========================================================================
## PART 2  Placebo, permuted at COUNTRY level
## ==========================================================================
cat("\n[..] PART 2: placebo (country-level permutation)\n")

# WHY COUNTRY LEVEL. M is a country-level variable: every tweet from one
# country shares its value, which is exactly why the clustered standard errors
# on the interactions are large. Permuting M across individual tweets destroys
# that structure, collapses the clustered SE and manufactures significance
# under the null. That is what produced X3:M at p = 1.7e-15 in the earlier run.
cmap <- df %>% distinct(country, M)
cmap$M_perm <- sample(cmap$M)
cat("    countries relabelled: ", sum(cmap$M != cmap$M_perm), " of ",
    nrow(cmap), "\n", sep = "")

df_p <- df %>% select(-M) %>%
  left_join(cmap %>% select(country, M = M_perm), by = "country")

stopifnot(nrow(df_p) == nrow(df))
# M must still be constant within country, otherwise the permutation is wrong
stopifnot(all((df_p %>% group_by(country) %>%
                 summarise(k = n_distinct(M), .groups = "drop"))$k == 1))

r5    <- glm(FORM, data = df_p, family = binomial())
r5_ct <- clust(r5, df_p)
keep  <- intersect(c("X1:M", "X2:M", "X3:M"), rownames(r5_ct))
print(r5_ct[keep, , drop = FALSE])

placebo <- tidy_base(r5_ct) %>%
  filter(term %in% keep) %>%
  mutate(permutation = "country level", n = nrow(df_p))
safe_write(placebo, "results/table_4_6_r5_placebo.csv")

if (all(placebo$p_value > .05)) {
  ok("All three placebo interactions are null, as the falsification test requires.")
} else {
  cat("[!!] At least one placebo interaction is significant. Do NOT report the\n",
      "     moderation results as surviving falsification until this is understood.\n", sep = "")
}

## ==========================================================================
## PART 3  Sentiment-artifact sensitivity check
## ==========================================================================
cat("\n[..] PART 3: sentiment-artifact sensitivity\n")

# WHY THIS REPLACES THE HIGH-CONFIDENCE SUBSAMPLE CHECK.
# Chapter III promised a subsample of tweets on which the classifiers report
# high confidence. The dataset has no confidence, probability or margin field,
# so that check cannot be run. The concern it was meant to address is that
# sentiment and stance are produced by classifiers reading the same text, so
# their association may be partly mechanical.
#
# The logic here targets the same concern with data that does exist. If the
# association were mostly classifier agreement, it should be concentrated in
# tweets whose tone is obvious, and should weaken or vanish once those are
# removed. Two tests:
#   3A  re-estimate on progressively narrower middle bands of sentiment,
#       discarding the tweets where tone is least ambiguous. Note that the
#       narrowest band is a severe test by construction: removing 70 percent
#       of the spread in sentiment leaves the slope weakly identified, so
#       attenuation there is expected and should be reported as such rather
#       than read as evidence of artifact.
#   3B  believer rate across sentiment deciles; a smooth monotone gradient is
#       harder to produce by a mechanical labelling artifact than a jump

## ---- 3A  middle-band restriction ----
bands <- list(
  list(lab = "Full sample",                 lo = 0.00, hi = 1.00),
  list(lab = "Middle 80 percent",           lo = 0.10, hi = 0.90),
  list(lab = "Middle 50 percent",           lo = 0.25, hi = 0.75),
  list(lab = "Middle 30 percent",           lo = 0.35, hi = 0.65)
)

band_rows <- bind_rows(lapply(bands, function(b) {
  d <- if (b$lo == 0) df else {
    q <- quantile(df$sentiment, c(b$lo, b$hi), na.rm = TRUE)
    df %>% filter(sentiment > q[1], sentiment < q[2])
  }
  m  <- glm(FORM, data = d, family = binomial())
  ct <- clust(m, d)
  tibble(
    band          = b$lab,
    n             = nrow(d),
    X1_odds_ratio = exp(ct["X1", 1]),
    X1_p          = ct["X1", 4],
    X1M_odds_ratio = exp(ct["X1:M", 1]),
    X1M_p          = ct["X1:M", 4]
  )
}))
print(as.data.frame(band_rows), digits = 4)
safe_write(band_rows, "results/table_4_6_r6_sentiment_bands.csv")

## ---- 3B  believer rate by sentiment decile ----
dec <- df %>%
  mutate(decile = ntile(sentiment, 10)) %>%
  group_by(decile) %>%
  summarise(n = n(),
            mean_sentiment = mean(sentiment),
            believer_rate  = mean(Y), .groups = "drop") %>%
  mutate(believer_pct = 100 * believer_rate)
print(as.data.frame(dec), digits = 4)
safe_write(dec, "results/table_4_6_r6_sentiment_deciles.csv")

mono <- all(diff(dec$believer_rate) > 0)
ct_test <- cor.test(dec$decile, dec$believer_rate, method = "spearman")
cat("    strictly increasing across all deciles: ", mono, "\n",
    "    Spearman rho = ", round(unname(ct_test$estimate), 4), "\n", sep = "")

## ---- 3C  specification without the weakest proxy ----
# Chapter III, section 3.4 flags aggressiveness as the least strong proxy.
noagg <- glm(Y ~ X1 + X3 + Z1 + M + X1:M + X3:M, data = df, family = binomial())
noagg_ct <- clust(noagg, df)
safe_write(tidy_base(noagg_ct) %>% mutate(spec = "aggressiveness excluded"),
           "results/table_4_6_r6_no_aggressiveness.csv")
cat("    X1 odds ratio without aggressiveness in the model: ",
    round(exp(noagg_ct["X1", 1]), 3), "\n", sep = "")

## ==========================================================================
## Session log and summary
## ==========================================================================
writeLines(capture.output(sessionInfo()), "results/sessionInfo_remaining.txt")

cat("\n==========================================================\n")
cat("DONE. Files written to results/\n")
print(list.files("results", pattern = "table_4_0|table_3_3|r5_placebo|r6_"))
cat("\nHow to read these results:\n")
cat("  1. PLACEBO. All three interaction p-values should be above .05. If any\n")
cat("     is significant, do not claim the moderation survives falsification.\n")
cat("  2. BANDS. Expect the sentiment odds ratio to hold up through the 80 and\n")
cat("     50 percent bands. The 30 percent band is a deliberately severe test:\n")
cat("     it leaves so little spread in sentiment that the slope is weakly\n")
cat("     identified, so attenuation there is expected and is NOT by itself\n")
cat("     evidence of artifact. Report all four bands, not a favourable subset.\n")
cat("  3. DECILES. A strictly increasing believer rate with Spearman rho near 1\n")
cat("     is the cleanest evidence against a mechanical labelling artifact.\n")
cat("==========================================================\n")
