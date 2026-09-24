###############################################################################
#  FALSIFICATION TEST, CORRECTED
#  Stance Analysis of Public Climate and Sustainability Discussion
#  Golam Rahat, ID 1112310093, BSA 4401
#
#  WHY THIS REPLACES THE SINGLE-DRAW PLACEBO
#
#  The earlier check relabelled the 142 countries once at random and read the
#  clustered p-values off that single draw. That design is uninformative here,
#  for a reason specific to this sample: country sizes are extremely unequal.
#  The United States alone is 56.9 percent of the analytic sample, the four
#  largest countries are 81.2 percent, and the median country has 35 tweets.
#
#  Under a single random relabelling, the share of tweets landing in the
#  permuted "developed" group ranges from about 2 percent to about 91 percent
#  across draws, with a median near 20 percent, against a true share of 91.6
#  percent. So one draw produces a contrast that need not resemble the real
#  one at all. Worse, because countries genuinely differ in believer rate and
#  in slope, almost ANY partition of countries generates a sizeable
#  interaction coefficient. A single draw will therefore often look
#  "significant" even when the labels carry no information, which is exactly
#  what happened: the one-draw run returned p = .002 and p = .018 on the two
#  linguistic interactions.
#
#  The correct procedure is a permutation TEST, not a single permutation. It
#  relabels the countries many times, builds the empirical null distribution
#  of the interaction coefficient, and asks where the observed coefficient
#  falls within it. That is calibrated for exactly the instability above.
#
#  HOW TO RUN
#     Put this file beside Stance_sample_data.csv, set the working directory,
#     then:  source("placebo_permutation_test.R")
#     Runtime is roughly four to eight minutes at B = 500.
###############################################################################

options(error = NULL)
pkgs <- c("tidyverse")
miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(42)
dir.create("results", showWarnings = FALSE)

B <- 500                      # number of permutations; 500 gives p to about .002

## ---- rebuild the analytic sample -----------------------------------------
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
    M  = as.integer(economy == "Developed (High income)")
  )
stopifnot(nrow(df) == 196406)

## ---- document the size imbalance that motivates this design ---------------
cs <- sort(table(df$country), decreasing = TRUE)
cat("\nCountry size imbalance\n")
cat("  countries              :", length(cs), "\n")
cat("  largest country share  :", sprintf("%.1f%% (%s)", 100*cs[1]/nrow(df), names(cs)[1]), "\n")
cat("  four largest share     :", sprintf("%.1f%%", 100*sum(cs[1:4])/nrow(df)), "\n")
cat("  median country size    :", median(as.numeric(cs)), "\n")

## ---- estimator: coefficients only, no clustering --------------------------
# The same statistic is used for the observed fit and for every permutation,
# so the comparison is internally consistent. Clustered standard errors are
# not needed here because the permutation distribution itself supplies the
# reference against which the observed coefficient is judged.
coefs <- function(Mv) {
  fit <- glm(df$Y ~ df$X1 + df$X2 + df$X3 + df$Z1 + Mv +
               df$X1:Mv + df$X2:Mv + df$X3:Mv, family = binomial())
  cf <- coef(fit)
  c(X1M = unname(cf["df$X1:Mv"]),
    X2M = unname(cf["df$X2:Mv"]),
    X3M = unname(cf["df$X3:Mv"]))
}

obs <- coefs(df$M)
cat("\nObserved interaction coefficients\n")
print(round(obs, 4))

## ---- permutation null: relabel COUNTRIES, not observations ---------------
countries <- sort(unique(df$country))
base_lab  <- df %>% group_by(country) %>% summarise(M = first(M), .groups = "drop") %>%
  arrange(country) %>% pull(M)
cidx <- match(df$country, countries)

cat("\nRunning", B, "permutations. This is the slow part.\n")
null <- matrix(NA_real_, nrow = B, ncol = 3,
               dimnames = list(NULL, c("X1M", "X2M", "X3M")))
t0 <- Sys.time()
for (b in seq_len(B)) {
  perm <- sample(base_lab)          # shuffle labels across countries
  null[b, ] <- coefs(perm[cidx])    # M stays constant within each country
  if (b %% 50 == 0)
    cat("   ", b, "/", B, "  (",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s)\n", sep = "")
}

## ---- two-sided permutation p-values --------------------------------------
res <- tibble(
  term            = c("X1:M", "X2:M", "X3:M"),
  observed_coef   = as.numeric(obs),
  observed_abs    = abs(as.numeric(obs)),
  null_median_abs = apply(abs(null), 2, median),
  null_p95_abs    = apply(abs(null), 2, quantile, 0.95),
  perm_p          = sapply(1:3, function(i)
                      (sum(abs(null[, i]) >= abs(obs[i])) + 1) / (B + 1)),
  B               = B,
  permutation     = "country level, B draws"
)
print(as.data.frame(res), digits = 4)
write_csv(res, "results/table_4_6_r5_placebo_permutation.csv")
write_csv(as_tibble(null), "results/table_4_6_r5_null_distribution.csv")

cat("\nInterpretation\n")
for (i in 1:3) {
  verdict <- if (res$perm_p[i] < .05)
    "observed effect exceeds the permutation null: survives falsification" else
    "observed effect lies inside the permutation null: does not survive"
  cat("  ", res$term[i], ": permutation p = ", sprintf("%.4f", res$perm_p[i]),
      "  ", verdict, "\n", sep = "")
}
cat("\nNote the null medians. Even under randomly assigned labels the\n")
cat("interaction coefficients are far from zero, which is precisely why a\n")
cat("single permutation draw could not have settled this question.\n")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_permutation.txt")
cat("\n[ok] Files written to results/\n")
