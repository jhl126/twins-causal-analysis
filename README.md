# Causal Analysis of the Twin Birth Weight Effect

This project estimates the causal effect of being born the heavier twin on first-year infant mortality, and tests whether that effect differs for infants whose mothers used tobacco or alcohol during pregnancy. It uses the Twins dataset (Louizos et al.) and applies a Doubly Robust estimator for the average treatment effect and a T-Learner for the conditional (subgroup) treatment effect.

## Overview

- **Research question:** Does being the heavier twin reduce first-year mortality, and is this effect modified by maternal prenatal substance exposure?
- **Headline findings:** Being the heavier twin reduces first-year mortality by approximately 2.6 percentage points (95% CI [-0.0329, -0.0189]), robust across four estimation methods and three refutation tests. No statistically detectable evidence was found that prenatal substance exposure modifies this effect.

## Repository Structure

```
twins-causal-analysis/
├── README.md
├── pyproject.toml          # uv-managed dependencies
├── data/
│   ├── raw/                # downloaded source CSVs (gitignored)
│   └── processed/          # cleaned data (gitignored)
├── notebooks/
│   ├── 01_eda.ipynb                  # EDA, variable table, figures, causal DAG
│   ├── 02_ate_doubly_robust.ipynb    # ATE estimation (DR, regression, naive)
│   ├── 03_cate_t_learner.ipynb       # CATE / subgroup estimation (T-Learner)
│   └── 04_evaluation.ipynb           # DoWhy refutation tests and sensitivity analysis
├── src/
│   ├── __init__.py
│   ├── data_loader.py      # downloads and reshapes the Twins data
│   ├── preprocessing.py    # imputation, encoding, covariate selection
│   └── estimators.py       # DR, T-Learner, and cluster bootstrap functions
├── figures/                # generated figures used in the report
└── report/
    └── final_report.pdf
```

## Setup

This project uses `uv` for environment management.

```bash
# Clone the repository
git clone <YOUR_REPO_URL>
cd twins-causal-analysis

# Create the environment and install dependencies
uv sync

# Register the environment as a Jupyter kernel
uv run python -m ipykernel install --user --name=twins-causal --display-name="Python (twins-causal)"
```

When opening the notebooks, select the "Python (twins-causal)" kernel.

## Reproducing the Results

Each notebook is self-contained: it loads the data, applies preprocessing, and runs its analysis independently. The notebooks can be run in any order, though the numbering reflects the logical flow of the report.

The data is downloaded automatically from the source repository each time `load_twins_data()` is called, so no manual data download is required. An internet connection is needed.

| Notebook | Produces |
|:---|:---|
| `01_eda.ipynb` | Variable description table, distribution figures (Figures 1-3), and the causal DAG (Figure 4) |
| `02_ate_doubly_robust.ipynb` | The four ATE estimates, the doubly robust cluster-bootstrap confidence interval, and the propensity score diagnostic (Figure 5) |
| `03_cate_t_learner.ipynb` | The subgroup CATEs, the difference and its confidence interval, and the subgroup effect figure (Figure 6) |
| `04_evaluation.ipynb` | The DoWhy cross-check estimate, the three refutation test results, and the sensitivity analysis |

To regenerate all results, run each notebook top to bottom.

## Methods

- **Average Treatment Effect:** Doubly Robust estimator (Augmented Inverse Propensity Weighting), with logistic regression for both the propensity and outcome models. Confidence intervals are computed via a cluster bootstrap that resamples whole twin pairs to respect the within-pair correlation.
- **Conditional Average Treatment Effect:** T-Learner with logistic regression base learners, fit separately on the treated and control arms, with subgroup CATEs and their difference compared across the exposure subgroups.
- **Evaluation:** DoWhy refutation tests (placebo treatment, random common cause, data subset) and a sensitivity analysis for unobserved confounding.

## Data Source

The Twins dataset is sourced from the Counterfactual VAE repository associated with Louizos et al., and is filtered to same-sex twin pairs where both twins weigh under 2kg at birth. The data is downloaded directly from:
`https://github.com/AMLab-Amsterdam/CEVAE/tree/master/datasets/TWINS`

## Author

Josh Lim
ADSP 32029, Causal Models in Data Science
