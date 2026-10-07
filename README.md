# Active Credit-Portfolio Optimization

**Integrating Machine-Learning Underwriting, Endogenous Stiglitz–Weiss Risk Pricing, and Basel III/IV Economic Capital**

M.A. thesis in Financial Economics, Reichman University. Author: Niv Vaisman Zancho. Supervisor: Prof. Doron Avramov. Draft of September 2026; final submission expected late 2026.

📄 **[Read the thesis (PDF)](Vaisman_Zancho_2026_MA_Thesis_Active_Credit_Portfolio_Optimization.pdf)** · 💻 **[R implementation](active_credit_engine.R)** · 📊 **[Data](lendingclub_sample_10000_loans.csv)**

---

## In one paragraph

Most lenders use a machine-learning score only to accept or reject a loan whose terms were set elsewhere. This thesis builds a **model-agnostic, multi-stage underwriting engine** that instead co-optimizes the loan amount, maturity and interest rate. It also accounts for customer acquisition cost and the borrower's likelihood of accepting the offer, and it maximizes risk-adjusted return on capital. The engine is evaluated strictly out of sample against nine passive underwriting regimes on 10,000 LendingClub loans.

## Key results

| | Engine | Best passive regime |
|---|---|---|
| Expected net profit (base case) | **$12.55M** | $8.62M |
| Ratio | **1.456×** (paired-bootstrap 95% CI 1.40–1.51) | — |
| Lifetime RAROC | **99.9%** | 79.6% |

- **Robust.** The engine beats the best passive regime in every check. Those checks are:
  - six loss profiles (1.35–1.51×) and fourteen ±25% parameter shocks;
  - a ten-fold increase in demand elasticity;
  - zero acquisition cost and an alternative PD model;
  - twenty re-estimated random train/test splits.
- **Machine-Learning Underwriting Paradox.** Discrimination (AUC) does not order value. A model with a statistically superior AUC earns less. A calibrated and an uncalibrated Random Forest have identical AUC (0.650) but very different profit, because profit comes from calibration and pricing, not from ranking.
- **Mechanism.** About 94% of the gain comes from the rate axis. It concentrates in the middle of the credit distribution, not on the riskiest or lowest-income borrowers.

## Methods

The engine combines these components:
- **Default modelling:**
  - Random Forest with Platt calibration (5-fold out-of-fold) and logistic regression;
  - ROC/AUC with DeLong tests;
  - out-of-sample log-likelihood, AIC/BIC and Vuong tests.
- **Endogenous default (Stiglitz–Weiss).** Rate and payment-to-income elasticities of default are estimated from loan-level data, with a sub-grade fixed-effects identification check.
- **Contract optimizer:**
  - 1,178 candidate contracts per borrower, valued with discounted amortization and an explicit time-to-default and LGD;
  - a Basel ASRF capital charge inside the objective;
  - deployability constraints: maximum rate hike, rate multiple and PTI cap.
- **Screening mechanism.** A three-tier bargaining tree with heterogeneous, partially observed demand and Bayesian updating after a rejection.
- **Evaluation:**
  - expected, pseudo-realized and realized profit;
  - Vasicek economic capital by Monte Carlo;
  - paired bootstrap and a tornado sensitivity analysis.

## Reproducing

Requires R ≥ 4.5 with `dplyr`, `pROC` and `randomForest`. From the repository root, run:

```
Rscript active_credit_engine.R
```

The script and the data file are named `sidur_13.R` and `data_for_final.csv` in the thesis (Appendix E). Seeds are fixed, so a re-run reproduces every table in the thesis. The full run takes about 19 minutes. Outputs are `*_v13.csv` tables, `active_credit_engine_output.txt` and `active_credit_engine_plots.pdf`.

## Data

`lendingclub_sample_10000_loans.csv` is a sample of 10,000 loans from the public **LendingClub** loan dataset (source: LendingClub; also available on [Kaggle](https://www.kaggle.com/datasets/wordsforthewise/lending-club)), included so the results can be reproduced. It holds loan amount, rate, grade/sub-grade, employment length, revolving balance and utilization, FICO, home ownership, annual income, verification status and loan outcome. LendingClub serves only as a public, well-studied testbed for the framework. No proprietary or employer data is used anywhere in this work.

## Copyright

© 2026 Niv Vaisman Zancho. All rights reserved (the LendingClub data remain the property of their owner). The thesis and code are shared for review purposes only; no licence to copy, modify or redistribute them is granted. If you refer to this work, please cite it as:

> Vaisman Zancho, N. (2026). *Active Credit-Portfolio Optimization: Integrating Machine-Learning Underwriting, Endogenous Stiglitz–Weiss Risk Pricing, and Basel III/IV Economic Capital.* M.A. thesis, Reichman University.
