# CardioBridge — NHANES biochemical prediction benchmark

Executed 7 October 2026. Version 0.1: phenotypic prediction module only.

## Question
Can six routine demographic and physical measurements predict measured fasting triglyceride concentration in a later NHANES survey cycle?

## Actual results
Training: NHANES 2013–2014 (2,553 adults) and 2015–2016 (2,254 adults), total 4,807. Held-out evaluation: 2017–2018, 2,393 adults. Total analysed: 7,200.

| Model | Weighted MAE, mg/dL | Weighted RMSE, mg/dL | Weighted R² |
|---|---:|---:|---:|
| Training-mean baseline | 55.49 | 96.23 | -0.0021 |
| Elastic net | 53.51 | 93.85 | 0.0469 |
| Random forest | 53.05 | 93.16 | 0.0608 |
| Gradient boosting | 52.49 | 93.03 | 0.0634 |

Gradient boosting had the lowest observed test errors, reducing RMSE by approximately 3.3% compared with the baseline. R² of 0.0634 represents approximately 6.3% explained variation in the weighted test outcome. These modest results do not establish clinical utility, and differences between models have not been tested for statistical significance. Model ranking here is descriptive; selecting a deployed model from these test results would require a new evaluation dataset.

## Reproduce
Python 3 with the packages listed in requirements.txt:

```bash
pip install -r requirements.txt
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python run_analysis.py
```

The script downloads 12 public CDC XPT files and caches them under data/. Raw files are excluded from the deliverable ZIP and should not be committed to GitHub. Internet access to wwwn.cdc.gov is required for the first run.

## Design and safeguards
- Adults aged at least 20 with a positive measured LBXTR and positive fasting-subsample weight WTSAF2YR.
- Predictors: age, recorded sex, BMI, waist circumference, mean systolic and mean diastolic blood pressure from the first three readings.
- Predicts measured concentration in mg/dL directly; no outcome clipping or removal of high triglyceride observations.
- Calendar-cycle holdout: no test participants in training. All preprocessing is fitted on training data only.
- Median imputation with missingness indicators and standard scaling precede each model.
- Fixed model configurations and seed 42. No tuning on the test set and no claim of an optimised model.
- Weighted test metrics use 2017–2018 fasting-subsample weights.
- Calculated LDL, triglyceride-derived ratios, and alternate-unit triglyceride variables are excluded to prevent target leakage.
- Source URLs, SHA-256 hashes and runtime versions are recorded in outputs/provenance.json.

## Outputs
- model_comparison.csv: actual test metrics and model training times.
- cohort_flow.csv: total, adult and eligible counts per cycle.
- missingness.csv: predictor missingness proportions by cycle.
- test_predictions.csv: predictions, measured outcomes, weights and public NHANES participant identifiers.
- subgroup_metrics.csv: descriptive metrics by recorded sex and age band.
- benchmark.png: held-out RMSE and observed-versus-predicted plot.
- provenance.json and run.log: reproducibility and execution records.

## Limits and next steps
Routine physical measurements alone provide limited predictive information. Medication use, dietary measures, smoking, other biomarkers, and genetic information were not included. Extreme triglyceride values remain in the analysis and can strongly influence RMSE. The later survey cycle is a temporal holdout within the same survey programme, not validation in a separate clinical cohort.

For the later genetic analyses, see CardioBridge_analysis_report.md.

## Sources
CDC NHANES component files: DEMO, BMX, BPX and TRIGLY for 2013–2014, 2015–2016 and 2017–2018.

These are secondary analyses of public de-identified data. Cite NHANES and relevant software in any manuscript. This repository is an independent portfolio prototype.
