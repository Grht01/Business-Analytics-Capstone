# Stance Analysis of Public Climate and Sustainability Discussion

**Affective, Rhetorical, and Experiential Drivers of Climate Stance, Moderated by Country Economy**

Replication archive for a BBA Business Analytics capstone.

Golam Rahat (ID 1112310093) · BSA 4401, Summer 2026
Supervisor: Dr. Md. Qamruzzaman
School of Business and Economics, United International University

Project page: https://grht01.github.io/Business-Analytics-Capstone/

---

## What this is

This repository holds the data extract, estimation code, and output tables behind the capstone report. Running the two scripts on the supplied data reproduces every table and figure in Chapter IV and every appendix table.

The study asks whether three properties of a tweet predict its position on climate change: how positive it sounds, how aggressive its language is, and how unusual the local temperature was. It then asks whether the economic development level of the country the tweet comes from changes how strongly each of those works.

## Contents

| File | Description |
|---|---|
| `Stance_sample_data.zip` | Analytic data extract. A stratified random draw of 250,000 tweets from the Climate Change Twitter Dataset, seed 42. |
| `analysis.R` | Main estimation. Equations 1 to 4, diagnostics, average marginal effects, robustness checks R1 to R4, and a self-verification block. |
| `placebo_permutation_test.R` | Falsification test. 500 relabellings of the moderator across the 142 countries, building an empirical null distribution. |
| `results.zip` | Output tables as CSV. |

## Reproducing the analysis

1. Unzip `Stance_sample_data.zip` so `Stance_sample_data.csv` sits beside the scripts.
2. Set that folder as the R working directory.
3. Run the main estimation. Missing packages install automatically.

```r
source("analysis.R")
```

Expect 15 to 45 minutes. Equation 4, the country random-intercept model, is the slow step and can be skipped with `RUN_GLMER <- FALSE`.

4. Run the falsification test.

```r
source("placebo_permutation_test.R")
```

Expect 4 to 8 minutes at 500 permutations.

Output is written to `results/`. The last section of `analysis.R` prints a verification table comparing its own output against the values reported in the report. Every row should read OK.

Both scripts halt immediately if the input is not 250,000 rows or the analytic sample is not 196,406 rows, so loading the wrong file fails loudly instead of producing plausible but wrong numbers.

## Method

The outcome is binary, believer against denier, so estimation is by binary logistic regression. Country economy and the interactions built from it are country-level quantities applied to tweet-level observations, so tweets from one country are not independent with respect to the moderator. Inference uses country-clustered standard errors alongside a country random intercept, and moderation is interpreted through average marginal effects rather than raw interaction coefficients.

Analytic sample: 196,406 tweets, 142 countries, 2006 to 2019. 179,872 from developed economies, 16,534 from developing.

## Results

| Hypothesis | Effect | Outcome |
|---|---|---|
| H1 Sentiment | +4.38 pp per 1 SD | Supported |
| H2 Aggressiveness | -3.83 pp | Supported |
| H3 Temperature anomaly | +1.20 pp per 1 SD | Supported |
| H4 Economy moderates sentiment | OR 0.534 | Supported |
| H5 Economy moderates aggressiveness | OR 1.517 | Supported |
| H6 Economy moderates temperature | OR 1.023 | Not supported |

Marginal effects are for developed economies, in percentage points on the probability of a believer stance. All six decisions are identical under the pooled, clustered, fixed effects and random intercept specifications.

One result worth stating plainly: the model carries real information about the ordering of cases, with an area under the curve of 0.700, yet at the conventional 0.5 threshold it identifies no deniers at all and its accuracy of 0.9152 is exactly the believer base rate. Association at population scale and usable detection are different achievements.

## Software

R 4.5.2 on Windows 11 x64, with lme4 1.1-38, marginaleffects 0.32.0, sandwich 3.1-3, lmtest 0.9-40, car 3.1-3, pROC 1.19.0.1, nnet 7.3-20 and tidyverse 2.0.0.

The scripts avoid `broom` and `broom.mixed` deliberately, since they carry a `purrr` version dependency that can fail on an otherwise working installation. Model output is extracted with base R instead.

## Data provenance

The data here is a derived sample, not the original corpus. It comes from the Climate Change Twitter Dataset compiled and released by Effrosynidis, Karasakalidis, Sylaios and Arampatzis, documented in *Expert Systems with Applications* (2022) with the descriptive analysis in *PLOS ONE* (2022). All credit for the underlying collection and annotation belongs to them.

Anyone reusing this extract should cite the original dataset, check its licence terms, and treat every label as a machine classification carrying error rather than as verified ground truth.

## Citation

```
Rahat, G. (2026). Stance analysis of public climate and sustainability
discussion: Affective, rhetorical, and experiential drivers of climate
stance, moderated by country economy [Capstone report]. School of
Business and Economics, United International University.
```
