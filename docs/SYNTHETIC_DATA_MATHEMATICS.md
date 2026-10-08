# Mathematical construction of the synthetic FTIR datasets

## Purpose and provenance

These datasets are **algorithmically simulated**, not anonymised, perturbed, resampled, reconstructed, fitted, or inferred from the authors' measured FTIR absorbances. Fixed numerical parameters were chosen for demonstration, not estimated from milk. Their role is to provide realistically smooth *FTIR-like* inputs to the R workflow, not physically calibrated milk spectra and not reproducible estimates of reported prediction accuracy.

The **source of truth** is `simulation/generate_synthetic_data.py`. The equations below document the exact implemented generator, including the dependence caused by the use of `numpy.roll` at the spectral grid boundaries.

## Discretized spectral domain and fixed basis

Let j = 1,...,J with J = 1060 and define the wavenumber grid

\[\nu_j=925.92 + (j-1)\frac{5011.54-925.92}{1059}\quad (\mathrm{cm}^{-1}).\]

For k = 1,...,8 define Gaussian radial basis peaks

\[\phi_k(\nu)=\exp\left[-\frac12\left(\frac{\nu-\mu_k}{\sigma_k}\right)^2\right].\]

The fixed peak centers (cm^-1) are

\[(\mu_1,...,\mu_8)=(1050,1160,1545,1650,1745,2850,2920,3300),\]

and the fixed width parameters (cm^-1) are

\[(\sigma_1,...,\sigma_8)=(70,65,90,100,55,80,80,170).\]

**These widths are Gaussian standard-deviation parameters, not full widths at half maximum.**

A common synthetic reference curve is

\[m(\nu)=0.015+0.025\frac{\nu-925.92}{5011.54-925.92}+\sum_{k=1}^8 a_k\phi_k(\nu),\]

where

\[(a_1,...,a_8)=(0.09,0.12,0.14,0.22,0.26,0.10,0.12,0.06).\]

## Farm, class, individual and dataset effects

Let b in {1,2} index the dataset, f the simulated farm and i an observation. The observed synthetic intensity is

\[X_{bfi}(\nu_j)=m(\nu_j)+s_b+\sum_{k=1}^8 [F_{bfk}+C_{bfi,k}+h_{c(bf)}v_k]\phi_k(\nu_j)+N_{bfi,j}+U_{bfi}.\]

Here

- `s_1 = 0` and `s_2 = 0.012` are constant artificial shifts between synthetic datasets;
- `F_{bfk} ~ N(0, 0.015^2)` are random coefficients drawn **once per farm** for each peak k;
- `C_{bfi,k} ~ N(0, 0.009^2)` are independent individual coefficients for each peak;
- `h_plain = -0.010`, `h_hill = 0`, and `h_mountain = +0.010` specify an arbitrary simulated class signal;
- `v=(0.4,0.2,0.9,1,1,0.4,0.8,0.1)` specifies its fixed peak-specific weights;
- `U_{bfi} ~ N(0,0.004^2)` is a random **individual constant baseline offset**;
- all random values are generated using NumPy's `default_rng(20261008)` and in the order specified by the generator.

The altitude class effect is an **illustrative signal deliberately introduced** to exercise classification. It is **not** an estimate of a biological altitude effect. Farm and altitude labels are nested by construction.

## Correlated noise

For each individual, draw iid

\[E_{bfi,j}\sim N(0,0.003^2),\qquad j=1,\dots,1060.\]

The implemented noise is

\[N_{bfi,j}=\frac{E_{bfi,j}+E_{bfi,j-1}+E_{bfi,j+1}}{3},\]

with **circular indexing**, i.e., `E_0 = E_1060` and `E_1061 = E_1`, because the Python implementation uses `numpy.roll`. Neighboring intensities therefore have short-range noise correlation. The random offset U adds long-range within-spectrum correlation. The model does **not** use an empirically estimated noise covariance.

## Sample design, simulated IDs and file format

Dataset 1 contains 10 fictional farms and 460 rows, divided into plain / hill / mountain with farm counts 4 / 3 / 3. Individual sample counts by farm are:

`45, 44, 44, 44 | 47, 47, 46 | 48, 48, 47`.

Dataset 2 contains 9 fictional farms and 450 rows, divided into 3 / 3 / 3 farms, with 50 spectra per farm.

Rows are grouped by farm; `Azienda` and `Matricola` are generated identifiers `SYN_*`; `Zona` takes the Italian categories `Pianura`, `Collina`, and `Montagna`. Each Excel includes exactly three metadata columns followed by `VAR_240` through `VAR_1299`. Numeric intensities are saved to **seven digits after the decimal point** (`{v:.7f}`), so the packaged Excel values are rounded compared with the internal floating-point draws.

Datasets 1 and 2 are generated sequentially **using one continuous RNG stream**; regenerating only dataset 2 after resetting the seed will not yield the published dataset 2. The dataset tags `SW` and `VDB` in synthetic identifiers are purely programmatic labels, not released farm or animal identifiers. All generator outputs are artificial.

## Reproducibility and limitations

The generator is deterministic for the chosen seed and the same NumPy random-number implementation. Recorded R nested-CV seed is **69**, separate from the simulation seed **20261008**. Simulated results need not match real-data classification performance, feature selection, spectral chemistry or effect magnitudes. The supplied Excel workbooks support software demonstration and do **not** enable independent replication of the original numerical findings without restricted original data.
