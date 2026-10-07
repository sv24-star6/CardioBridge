# CardioBridge polygenic-score module

This extension demonstrates a reproducible **polygenic-score calculation workflow** using a published coronary artery disease score from the PGS Catalog.

## Selected score

**PGS012581 / PRS169_CAD** predicts coronary artery disease. The PGS Catalog reports 169 variants, beta effect weights, and construction from genome-wide-significant SNPs with LD pruning at R² < 0.01. The Catalog released the score on 17 June 2026.

The source publication is Zheng J et al., *Journal of Internal Medicine* (2024). The Catalog reports external UK Biobank evaluation across multiple ancestry groups. 

## What this module does

`run_pgs.py` downloads the official PGS Catalog scoring file, records its SHA-256 hash and metadata, checks the scoring weights, harmonises a target dosage table by rsID and effect allele, calculates the weighted score as:

`PGS = sum(dosage_i × effect_weight_i)`

and reports variant coverage for each participant.

## Target genotype format

A tab-delimited file with:

`sample_id, rsID, effect_allele, other_allele, dosage`

where dosage is 0–2 copies of the score effect allele.

Run metadata/QC without individual-level genotypes:

```bash
python pgs/run_pgs.py
```

Calculate scores when an appropriate genotype dataset is available:

```bash
python pgs/run_pgs.py --genotypes target_genotypes.tsv.gz
```

## Important limitation

CardioBridge does **not** currently contain an independent individual-level genotype cohort with adjudicated CAD outcomes. Therefore this module demonstrates score acquisition, harmonisation and computation, but it does not claim independent discrimination, calibration, hazard ratios, odds ratios, or clinical validation.

NHANES participants used in the phenotype benchmark are not linked to the PGS data. The two modules must not be represented as matched multi-omic data.

## Recommended next validation

Apply the score to an accessible cohort containing individual-level genotype data, ancestry covariates and CAD phenotype. Report variant coverage, ancestry-specific score distribution, logistic/Cox association, discrimination, calibration and incremental performance over age/sex/clinical covariates.
