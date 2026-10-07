# CardioBridge

Reproducible public-data cardiometabolic research prototype. Version 0.2, executed 7 October 2026.

Start with **CardioBridge_analysis_report.md** for methods, actual results and limitations. **NHANES_METHODS.md** documents the original prediction module.

## Completed analyses

- Triglyceride model benchmark in 7,200 NHANES adults, including 2,393 held-out participants.
- Cis-pQTL single-variant MR for six proteins, using INTERVAL and CAD GWAS summary statistics.
- Regional ABF colocalisation and prior sensitivity for all six proteins.
- Reactome annotation through UniProt cross-references.
- Batched inference on one million resampled rows (computational test, not new participants).

The strongest shared regional signal is APOE. ANGPTL1 is prior-sensitive. IL6R's MR association is not supported by strong colocalisation under the single-signal model used here. These findings are exploratory, not validated causal mechanisms or treatment recommendations.

## Run

Install requirements.txt and follow the commands in the report. Python scripts fetch public source files. Run run_analysis.py before run_scalability.py. Run genetics/run_genetics.py before the other genetics scripts.

## Repository contents

- run_analysis.py: original prediction benchmark.
- run_scalability.py: throughput benchmark.
- genetics/run_genetics.py: data retrieval, QC, allele matching and MR.
- genetics/remote_tabix.py: indexed public GWAS region retrieval.
- genetics/run_colocalisation.py: ABF colocalisation with prior sensitivity.
- genetics/run_pathways.py: Reactome cross-reference annotation.
- genetics/verify_results.py: source and numerical consistency checks.
- genetics/make_report.py: report and figure generation.
- outputs/ and genetics/outputs/: actual results, figures and provenance.

Downloaded raw data and caches are excluded. The analysis does not create a new GWAS, validate polygenic scores, link NHANES to genetic participants, or provide matched individual-level multi-omics. It does not establish clinical utility.