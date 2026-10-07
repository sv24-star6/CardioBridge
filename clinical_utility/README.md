# Clinical-utility evaluation

This module is deliberately separated from the existing CardioBridge analyses because the current NHANES phenotype sample is **not linked** to the genetic/PGS data.

## Question

Does adding a CAD polygenic score provide incremental predictive and decision value beyond a clinical model in an independent linked cohort?

The framework compares:

1. **Clinical model** — age, sex and a supplied clinical risk score.
2. **PGS alone**.
3. **Clinical + PGS**.

It reports five-fold cross-validated AUC, Brier score, approximate calibration measures, change in AUC, categorical net reclassification improvement (5%, 10% and 20% risk boundaries), and decision-curve net benefit from 1% to 30% risk thresholds.

## Required validation data

A tab-delimited individual-level dataset containing:

`CAD, age, sex, clinical_score, PGS`

The outcome must be observed CAD status (0/1). The PGS should be calculated from the CardioBridge PGS module or an equivalently documented implementation. A real analysis should additionally account for ancestry/population structure, study design and appropriate clinical covariates.

Run:

```bash
python clinical_utility/evaluate_clinical_utility.py --data linked_validation_cohort.tsv.gz
```

## Interpretation

A higher AUC or positive NRI does **not** establish clinical utility. Decision-curve analysis asks whether a model has greater net benefit at clinically meaningful decision thresholds, but the thresholds and consequences must be justified for the intended use.

Evidence of clinical utility ultimately requires a defined clinical use case and evidence that implementing the score changes management in a beneficial way, ideally supported by prospective or interventional evaluation.

## Current CardioBridge status

**Framework implemented; independent clinical-utility analysis pending an appropriate linked genotype + phenotype validation cohort.**

No synthetic demonstration should be reported as patient evidence. No result from this module should be described as clinical utility unless it was generated from an appropriate external cohort and interpreted in the intended clinical context.
