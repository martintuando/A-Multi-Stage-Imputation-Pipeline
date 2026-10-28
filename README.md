# MICPR: A Multi-Stage Imputation Pipeline for Suppressed Health Data

## The Problem: Data Suppression in Public Health Research

CDC WONDER suppresses death counts between 1 and 9 to protect individual privacy. Although this policy is necessary, it creates a challenge for statistical analysis: suppressed counts appear as missing values, even though the underlying number of deaths is greater than zero.

This results in a Missing Not at Random (MNAR) problem, since the probability of an observation being missing depends on the value itself. Excluding these observations can distort the analytical sample, while replacing them with a constant or zero can misrepresent the underlying mortality data.

The goal of MICPR was to develop a structured approach to handling these gaps. Rather than relying on a single imputation method, I combined temporal information, aggregated data and statistical imputation to make better use of the information available in the CDC WONDER records.

## The Solution: A Multi-Stage Imputation Pipeline

The pipeline combines three complementary components. Each addresses a different source of information: neighbouring years help estimate selected missing values, aggregated data provide a broader basis for estimating suppressed counts, and MICE handles missing values in the aggregate dataset.

### Stage 1: Rule-Based Temporal Substitution

The first stage uses observations from the preceding and following years within the same county. Its rules are inspired by the methodological approach discussed by Erdman et al. (2021), who examined methods for estimating suppressed values in protected public-health surveillance data.

The implementation applies three rules:

- **Upper-bound substitution:** If both adjacent years report at least 10 deaths, the missing value is set to 9, the highest count within the suppression range.
- **Neighbour-based estimation:** If one adjacent year reports at least 10 deaths and the other is either zero or missing, half of the available count is used as an estimate, rounded to an integer.
- **Unresolved gaps:** If neither condition is met, the value remains missing and is handled by subsequent stages.

These rules provide a simple way to use temporal information without assigning estimates to every missing observation. The implementation adapts the general methodological ideas discussed by Erdman et al.; it is not a direct replication of their procedure.

### Stage 2: Estimation from Aggregated Data

Temporal substitution alone does not resolve all suppressed observations. To address the remaining gaps, the pipeline uses a larger aggregated dataset.

The motivation is straightforward: a count suppressed in the annual data may be represented more clearly when observations are aggregated over a broader period. This provides additional information that can be used to estimate missing annual values, rather than relying solely on neighbouring years.

The aggregate estimates are converted into annual approximations and used to fill remaining missing observations. In the implementation, the estimates are divided by two and rounded to integers, with a minimum value of one.

This conversion is an approximation based on the project's aggregation scheme. It does not recover the original suppressed counts, and its suitability depends on the relationship between the aggregate data and the annual observations.

The approach is informed by the broader problem addressed by Erdman et al. (2021): finding practical ways to estimate suppressed public-health data when the original values are unavailable. Their work provides methodological context for this project, but does not validate the complete MICPR pipeline.

### Stage 3: Statistical Imputation with MICE

Missing death counts in the aggregate dataset are imputed using **Multiple Imputation by Chained Equations (MICE)** with **Predictive Mean Matching (PMM)**.

**Log transformation**

Death counts are transformed using `log(Deaths + 1)`. Adding one accommodates zero values, while the logarithm reduces the influence of large counts and can make the distribution more suitable for modelling.

**Population as a predictor**

Population is used as a predictor because county size provides relevant context for interpreting mortality counts. County codes are retained as identifiers rather than used as explanatory variables.

**Constrained post-processing**

After imputation, the estimates are transformed back to the original count scale. The pipeline then applies min-max rescaling to the interval from 1 to 9, rounds the results to integers and clips negative values to zero.

This step addresses a practical issue with the imputed estimates: their back-transformed values are not necessarily within the range of suppressed counts. Rescaling places the estimates within the intended interval, but it also changes their original scale and distribution. It should therefore be understood as a processing decision rather than a guarantee of statistical accuracy.

The MICE procedure generates five imputations over five iterations. The subsequent pipeline uses the first completed dataset rather than combining all five. Consequently, variation across imputations is not carried through to the final dataset.

## Data Processing and Output

The pipeline processes the two mortality outcomes separately before merging them at county-year level.

1. Load the CDC WONDER data for 2014–2019 and convert suppressed entries to missing values.
2. Apply the temporal substitution rules to eligible observations.
3. Impute missing aggregate counts using MICE and PMM.
4. Use the processed aggregate estimates to fill remaining gaps in the annual data.
5. Merge the outcomes and reshape the data into wide format.

The final output is `Master_Imputed_Outcomes_Wide.csv`. Each row represents a county, with separate variables for mental and behavioural disorder deaths, suicide deaths and population organised by year.

## Limitations

The true suppressed counts are unavailable, so the imputed values cannot be directly validated against the original observations. The aggregation-based estimates rely on assumptions about how the broader data relate to missing annual counts, while min-max rescaling changes the distribution of imputed values.

In addition, the final dataset uses only the first completed MICE dataset. The uncertainty across the five imputations is therefore not reflected in subsequent analyses. This should be considered when using the data for regression or other statistical modelling.

MICPR provides a structured way to handle suppressed observations using the information available in the data. It does not establish that every imputed count matches the true value.

## Implementation

The project is implemented in R using:

- `mice` for multiple imputation and Predictive Mean Matching
- `dplyr` for data manipulation and joins
- `tidyr` for reshaping the final dataset

A random seed of `42` is set for reproducibility. The file paths refer to the original development environment and must be updated before running the pipeline elsewhere.

## Reference

Erdman et al. (2021). *A Novel Imputation Approach for Sharing Protected Public Health Data.*

The temporal substitution rules and the broader motivation for estimating suppressed values are informed by this work. MICPR adapts these ideas to the available CDC WONDER data and combines them with aggregate-based estimation and MICE.
