# Functional FTIR Analysis of Sheep-Milk Production Systems

## Reproducibility package accompanying the Journal of Dairy Science manuscript

This repository provides analytical R code, fully synthetic FTIR datasets, a data-generation script, and mathematical documentation for the study *Functional Analysis of FTIR Spectra for Classifying Sheep Milk from Altitude-Associated Production Systems*.

The original measured spectra and farm-level records are not publicly available due to data-sharing restrictions and are **not included**. The synthetic spectra were generated independently of the measured absorbances. They provide example inputs for examining the computational workflow, not an independent replication of the empirical results reported in the manuscript.

## Repository contents

- `R/FTIR_FUNCTIONS.R` — shared functions for spectral preprocessing, feature extraction, spectral-domain selection, and classification.
- `scripts/RUN_ALTITUDE.R` — primary sample-level nested cross-validation using SVM and multinomial elastic net (seed 69; 10 outer and 5 inner folds).
- `scripts/RUN_BREED.R` — exploratory breed-classification analysis.
- `scripts/RUN_OUTLIER.R` — functional outlier sensitivity analysis.
- `data_synthetic/FTIR_simulated_spectra_dataset_1.xlsx` — 460 synthetic spectra across 10 simulated farms.
- `data_synthetic/FTIR_simulated_spectra_dataset_2.xlsx` — 450 synthetic spectra across 9 simulated farms.
- `simulation/generate_synthetic_data.py` — synthetic-data generation script.
- `docs/SYNTHETIC_DATA_MATHEMATICS.md` — mathematical model, parameters, and limitations of the simulation.
- `tests/validate_package.py` — structural and generator reproducibility checks.
- `MANIFEST_SHA256.txt` — SHA-256 checksums of repository files.

The generic dataset filenames intentionally omit breed names. Dataset 1 follows the sampling design of the Sarda subset; Dataset 2 follows that of the Valle del Belice subset. Neither dataset contains measured spectra or original individual records.

## Software requirements

The analyses were developed in the R 4.5.1 software context. Required R packages include `readxl`, `dplyr`, `stringr`, `signal`, `pracma`, `MASS`, `caret`, `e1071`, `glmnet`, `future`, `future.apply`, and `roahd` (for the outlier sensitivity analysis). The simulation generator requires Python and NumPy.

Package installation and full R execution have **not yet been independently verified** for this release. Check dependencies and run the scripts in RStudio before archival deposition; computational requirements may be substantial.

## Running the R analyses

Start R or RStudio with the repository root as your working directory. Run each analysis in a fresh R session:

```r
source("scripts/RUN_ALTITUDE.R")
```

For the supplementary workflows, use a separate fresh R session for each:

```r
source("scripts/RUN_BREED.R")
```

```r
source("scripts/RUN_OUTLIER.R")
```

The scripts read the synthetic `.xlsx` inputs from `data_synthetic/` and compute their transformations directly. Any `.rds` files written during execution are outputs rather than required data inputs.

## Synthetic data generation

From the repository root, run:

```bash
python simulation/generate_synthetic_data.py
python tests/validate_package.py
```

The generation seed is **20261008**, distinct from the analytical cross-validation seed **69**. The two Excel inputs contain the columns `Azienda`, `Zona`, `Matricola`, and 1,060 simulated spectral columns (`VAR_240` through `VAR_1299`). The generation model combines smooth spectral components, Gaussian absorption features, simulated farm-level and individual effects, altitude-category effects, and correlated noise.

For the precise equations, fixed coefficients, and assumptions, see [`docs/SYNTHETIC_DATA_MATHEMATICS.md`](docs/SYNTHETIC_DATA_MATHEMATICS.md).

## Interpretation and data availability

Synthetic farm offsets, class effects, and spectral noise are artificial, **not estimates from the study data**. Performance measured on these synthetic files is illustrative and should not be compared numerically with the manuscript's reported results. The example files support examination of the computational procedures while the original experimental data remain restricted.

## Citation and archiving

The repository can be archived as a tagged release on Zenodo. Add the permanent DOI only after the corresponding release has been deposited.
