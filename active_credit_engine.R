# =========================================================================
# active_credit_engine.R  —  Active credit-underwriting engine, version 13 (2026-09-14)
# (named sidur_13.R in the thesis, Appendix E)
# 1:1 port of the Python reference implementation `engine13/` (validated).
# Changes vs sidur_12.R (see ROADMAP / ENGINE_IMPROVEMENTS_PROPOSAL):
#   1 loss profile inside the optimizer (tau_bar, LGD on balance-at-default)
#   2 Platt-calibrated RF probabilities (5-fold out-of-fold on train)
#   3 acceptance with explicit, heterogeneous APR elasticity (eta0, xi)
#   4 NAM gate on the full mechanism value (tree value) for Group 1
#   5 Basel ASRF capital charge inside the objective (hurdle h, rho(PD))
#   6 two-type borrowers with Bayesian updating after a Tier-1 rejection
#   7 deployability constraints (max hike, max rate multiple, PTI cap)
#   + pseudo-realized column, passive +NAM / +NAM+K rows, Basel other-retail
#     correlation in EC, extended tornado, breakeven elasticity, scenarios.
#   + ADDENDUM (2026-09-25, append-only): zero-CAC, paired bootstrap, axis decomposition, PD-source
#     robustness, sub-grade fixed-effects identification, heterogeneity, DeLong AUC tests.
#   + A8 (2026-09-26, append-only): likelihood comparison of PD models (OOS log-lik, AIC/BIC, Vuong), LR tests of the
#     Stiglitz-Weiss channel, optimality gap of the contract grid (2x / 5x finer grids).
# Run:  Rscript active_credit_engine.R   (working dir = repository root; ~19 min incl. the 20-split appendix and addendum)
# Outputs: *_v13.csv tables, active_credit_engine_output.txt (via sink), active_credit_engine_plots.pdf
# =========================================================================
# run from the repository root, where lendingclub_sample_10000_loans.csv lives (data_for_final.csv in the thesis)
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(randomForest) })
set.seed(42)
sink("active_credit_engine_output.txt", split = TRUE)
t_start <- Sys.time()

# ---- configuration (base case) -------------------------------------------
cfg <- list(
  seed_split = 42, seed_rf = 123, seed_oof = 124, seed_z = 101, seed_mc = 42, scale = 5,
  r_disc = 0.05, term_orig = 60,
  tau_bar = 12, lgd = 0.8,                          # fix 1
  pd_source = "rf_cal",                             # fix 2: 'rf_raw' | 'rf_cal' | 'logit'
  kappa = 0.02, phi = 0.20, gamma_cac = 0.15,       # CAC
  alpha = 1.0, lambda = 0.003, eta0 = 0.05, xi = 0.5, eta_cut = 0.02, eta_cap = 4.0,   # fix 3
  capital_in_objective = TRUE, hurdle = 0.10, rho_mode = "basel_retail", rho_fixed = 0.15,  # fix 5
  constrain = TRUE, max_hike = 5.0, max_rate_mult = 1.5, max_pti = 0.15,               # fix 7
  two_type = TRUE, type_mult_low = 0.5, type_mult_high = 1.5, type_prior_high = 0.5, attrition = 0.0,  # fix 6
  gate_mode = "tree", exclude_rejected_from_tier2 = TRUE,                              # fix 4
  n_mc = 10000, rho_ec_mode = "basel_retail", rho_ec_fixed = 0.15, var_level = 0.999
)
frozen_overrides <- list(tau_bar = 1, lgd = 1.0, pd_source = "rf_raw", xi = 0, capital_in_objective = FALSE, rho_mode = "fixed",
                         constrain = FALSE, two_type = FALSE, gate_mode = "best", exclude_rejected_from_tier2 = FALSE,
                         rho_ec_mode = "fixed")
with_cfg <- function(base, ...) { o <- list(...); for (n in names(o)) base[[n]] <- o[[n]]; base }

# ---- data -----------------------------------------------------------------
data_path <- "lendingclub_sample_10000_loans.csv"       # run from the repository root ...
if (!file.exists(data_path) && requireNamespace("here", quietly = TRUE)) data_path <- here::here("lendingclub_sample_10000_loans.csv")
stopifnot(file.exists(data_path))
raw <- read.csv(data_path, stringsAsFactors = FALSE)
raw$Default_numeric <- ifelse(raw$default == "Defaulted", 1, 0)
raw$emp_length[is.na(raw$emp_length) | raw$emp_length == ""] <- "Unknown"
raw$revol_util[is.na(raw$revol_util)] <- median(raw$revol_util, na.rm = TRUE)
# pool rare factor levels (<10 loans) so every random split sees every level (avoids "new levels" in predict)
sg <- as.character(raw$sub_grade); rare_sg <- names(which(table(sg) < 10)); sg[sg %in% rare_sg] <- "G_rare"; raw$sub_grade <- sg
ho <- as.character(raw$home_ownership); ho[ho %in% c("NONE", "OTHER")] <- "OTHER"; raw$home_ownership <- ho
cat("Pooled rare sub-grades:", paste(rare_sg, collapse = ", "), "-> G_rare (", sum(sg == "G_rare"), "loans )\n")
for (v in c("grade", "sub_grade", "home_ownership", "verification_status", "emp_length")) raw[[v]] <- as.factor(raw[[v]])
L_bar <- mean(raw$loan_amnt); R_bar <- mean(raw$int_rate)

# ---- cash-flow engine (Spitzer, discounted; loss profile; Basel) ------------
pmt_of <- function(L, R, n) { rc <- R / 100 / 12; L * rc / (1 - (1 + rc)^(-n)) }
annuity <- function(n, r) ifelse(n > 0, (1 - (1 + r)^(-pmax(n, 1e-12))) / r, 0)
balance_after <- function(L, R, n, t) { rc <- R / 100 / 12; p <- pmt_of(L, R, n); L * (1 + rc)^t - p * ((1 + rc)^t - 1) / rc }
# NPV_cond = -L + (1-PD)*pv_nd + PD*pv_d ; ead = LGD * balance before the missed payment
pv_components <- function(L, R, n, tau_bar, lgd, r_m) {
  p <- pmt_of(L, R, n); tau <- pmin(tau_bar, n); k <- tau - 1; B <- balance_after(L, R, n, k)
  list(pv_nd = p * annuity(n, r_m), pv_d = p * annuity(k, r_m) + (1 - lgd) * B * (1 + r_m)^(-tau), ead = lgd * B)
}
balance_pv_years <- function(L, R, n, r_m) {  # (1/12) * sum_{t<n} (B_t/L) v^t, closed form
  rc <- R / 100 / 12; p <- pmt_of(L, R, n); v <- 1 / (1 + r_m); q <- (1 + rc) * v
  q1 <- abs(q - 1) < 1e-12
  sum_q <- ifelse(q1, n, (1 - q^n) / ifelse(q1, 1, 1 - q)); sum_v <- (1 - v^n) / (1 - v)
  ((L - p / rc) * sum_q + (p / rc) * sum_v) / L / 12
}
basel_rho_retail <- function(pd) { w <- (1 - exp(-35 * pd)) / (1 - exp(-35)); 0.03 * w + 0.16 * (1 - w) }
rho_of <- function(pd, mode, fixed) if (mode == "basel_retail") basel_rho_retail(pd) else rep(fixed, length(pd))
basel_k <- function(pd, lgd, rho, level = 0.999) { pd <- pmin(pmax(pd, 1e-6), 1 - 1e-6)
  lgd * (pnorm((qnorm(pd) + sqrt(rho) * qnorm(level)) / sqrt(1 - rho)) - pd) }
clip01 <- function(p) pmin(1 - 1e-4, pmax(1e-4, p))

# ---- canonical split, RF, OOF Platt calibration, scorecards, fit_emp --------
rf_formula <- as.factor(Default_numeric) ~ loan_amnt + int_rate + grade + sub_grade + emp_length + revol_bal + revol_util +
  fico + home_ownership + annual_inc + verification_status
pt <- function(m, d) predict(m, newdata = d, type = "response")
platt_apply <- function(p) plogis(coef(platt)[1] + coef(platt)[2] * qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)))
# fit_split(): canonical 80/20 split + all train-only estimation; assigns the globals used downstream
# (train, test, y_tr, y_te, rf_model, platt, m_*, fit_emp, emp). Re-called in the repeated-splits appendix.
fit_split <- function(seed_split) {
  set.seed(seed_split); idx <- sample(1:nrow(raw), 0.8 * nrow(raw))
  train <<- raw[idx, ]; test <<- raw[-idx, ]; y_tr <<- train$Default_numeric; y_te <<- test$Default_numeric
  set.seed(cfg$seed_rf); rf_model <<- randomForest(rf_formula, data = train, ntree = 500, importance = TRUE)
  test$pd_rf_raw <<- predict(rf_model, newdata = test, type = "prob")[, 2]
  set.seed(cfg$seed_oof); folds <- sample(rep(1:5, length.out = nrow(train))); oof <- numeric(nrow(train))
  for (k in 1:5) { set.seed(cfg$seed_rf + k)
    m <- randomForest(rf_formula, data = train[folds != k, ], ntree = 500)
    oof[folds == k] <- predict(m, newdata = train[folds == k, ], type = "prob")[, 2] }
  train$pd_rf_oof <<- oof
  platt <<- glm(y_tr ~ qlogis(pmin(pmax(oof, 1e-6), 1 - 1e-6)), family = binomial)
  test$pd_rf_platt <<- platt_apply(test$pd_rf_raw)
  m_fico <<- glm(Default_numeric ~ fico, data = train, family = binomial)
  m_grade <<- glm(Default_numeric ~ grade, data = train, family = binomial)
  m_sub <<- glm(Default_numeric ~ sub_grade, data = train, family = binomial)
  m_logit <<- glm(Default_numeric ~ loan_amnt + int_rate + grade + sub_grade + emp_length + revol_bal + revol_util + fico +
                    home_ownership + annual_inc + verification_status, data = train, family = binomial)
  test$pd_fico <<- pt(m_fico, test); test$pd_grade <<- pt(m_grade, test); test$pd_sub <<- pt(m_sub, test); test$pd_logit <<- pt(m_logit, test)
  train$pmt_to_inc <<- pmt_of(train$loan_amnt, train$int_rate, 60) / (train$annual_inc / 12)
  train$dti <<- train$revol_bal / train$annual_inc * 100
  fit_emp <<- glm(Default_numeric ~ int_rate + pmt_to_inc + fico + dti + revol_util, data = train, family = binomial)
  emp <<- list(gamma_emp = unname(coef(fit_emp)["int_rate"]), delta_emp = unname(coef(fit_emp)["pmt_to_inc"]))
  invisible(NULL)
}
fit_split(cfg$seed_split)
cat("\n--- STIGLITZ-WEISS fit_emp (train) ---\n"); print(summary(fit_emp)$coefficients)
write.csv(summary(fit_emp)$coefficients, "Empirical_StiglitzWeiss_Regression_v13.csv")
cat(sprintf("Platt calibration: intercept %.4f slope %.4f\n", coef(platt)[1], coef(platt)[2]))
# calibration diagnostics on test
calib_diag <- function(p, y) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); z <- qlogis(p); m <- glm(y ~ z, family = binomial)
  c(brier = mean((p - y)^2), auc = as.numeric(auc(roc(y, p, quiet = TRUE))), cal_intercept = unname(coef(m)[1]),
    cal_slope = unname(coef(m)[2]), mean_pred = mean(p), obs_rate = mean(y)) }
cal_tab <- t(sapply(c("pd_rf_raw", "pd_rf_platt", "pd_logit", "pd_fico", "pd_grade", "pd_sub"), function(v) calib_diag(test[[v]], y_te)))
cat("\n--- PD CALIBRATION DIAGNOSTICS (test) ---\n"); print(round(cal_tab, 4)); write.csv(cal_tab, "Calibration_Diagnostics_v13.csv")
select_pd <- function(cf) switch(cf$pd_source, rf_raw = test$pd_rf_raw, rf_cal = test$pd_rf_platt, logit = test$pd_logit)

# ---- the engine ------------------------------------------------------------
GRID <- expand.grid(pct = seq(1.00, 0.10, by = -0.05), term = c(36, 60), shift = seq(-5, 10, by = 0.5))
cac_of <- function(L, R, cf) cf$kappa * L * (L / L_bar)^cf$phi * (R / R_bar)^cf$gamma_cac

optimize_contract <- function(l0, r0, pd0, D, inc, cf, em, pd_ref) {
  r_m <- cf$r_disc / 12; T0 <- cf$term_orig; G <- nrow(GRID)
  pmt0 <- pmt_of(l0, r0, T0); c0 <- pv_components(l0, r0, T0, cf$tau_bar, cf$lgd, r_m); inc_m <- inc / 12
  pd0 <- clip01(pd0); cond0 <- -l0 + (1 - pd0) * c0$pv_nd + pd0 * c0$pv_d; real0 <- -l0 + (1 - D) * c0$pv_nd + D * c0$pv_d
  pti0 <- pmt0 / inc_m
  Rn <- r0 + GRID$shift; Ln <- l0 * GRID$pct; Tn <- GRID$term; pmtn <- pmt_of(Ln, Rn, Tn)
  cc <- pv_components(Ln, Rn, Tn, cf$tau_bar, cf$lgd, r_m); ptin <- pmtn / inc_m
  pd_dyn <- clip01(plogis(qlogis(pd0) + em$gamma_emp * (Rn - r0) + em$delta_emp * (pmtn - pmt0) / inc_m))
  cond <- -Ln + (1 - pd_dyn) * cc$pv_nd + pd_dyn * cc$pv_d; real <- -Ln + (1 - D) * cc$pv_nd + D * cc$pv_d
  up <- pmax(0, (pd_dyn - pd0) / (1 - pd0)); down <- pmax(0, (pd0 - pd_dyn) / pd0); D_ps <- D * (1 - down) + (1 - D) * up
  real_ps <- -Ln + (1 - D_ps) * cc$pv_nd + D_ps * cc$pv_d
  if (cf$capital_in_objective) {
    capcost <- cf$hurdle * basel_k(pd_dyn, cf$lgd, rho_of(pd_dyn, cf$rho_mode, cf$rho_fixed), cf$var_level) * Ln * balance_pv_years(Ln, Rn, Tn, r_m)
    capcost0 <- cf$hurdle * basel_k(pd0, cf$lgd, rho_of(pd0, cf$rho_mode, cf$rho_fixed), cf$var_level) * l0 * balance_pv_years(l0, r0, T0, r_m)
  } else { capcost <- rep(0, G); capcost0 <- 0 }
  value <- cond - capcost; value0 <- cond0 - capcost0
  # acceptance: heterogeneous elasticity, K borrower types
  eta_i <- cf$eta0 * min(cf$eta_cap, (pd_ref / pd0)^cf$xi)
  mults <- if (cf$two_type) c(cf$type_mult_low, cf$type_mult_high) else 1
  priors <- if (cf$two_type) c(1 - cf$type_prior_high, cf$type_prior_high) else 1
  size_pmt <- (Ln / l0)^cf$alpha * exp(-cf$lambda * pmax(0, pmtn - pmt0)); dR <- Rn - r0
  p_k <- sapply(mults, function(m) pmin(1, pmax(0, size_pmt * ifelse(dR > 0, exp(-eta_i * m * dR), 1 + cf$eta_cut * (-dR)))))  # G x K
  p_k <- matrix(p_k, nrow = G); pbar <- as.vector(p_k %*% priors)
  valid <- Rn > 1
  if (cf$constrain) valid <- valid & (dR <= cf$max_hike + 1e-9) & (Rn <= cf$max_rate_mult * r0 + 1e-9) & (ptin <= max(cf$max_pti, pti0) + 1e-12)
  if (!is.null(cf$grid_mask) && cf$grid_mask != "all") {   # addendum A3: restrict the axes the engine may move
    on_p <- abs(GRID$pct - 1) < 1e-9; on_t <- GRID$term == T0; on_s <- abs(GRID$shift) < 1e-9
    valid <- valid & switch(cf$grid_mask, rate = on_p & on_t, rate_term = on_p, rate_amount = on_t, amount_term = on_s, none = on_p & on_t & on_s) }
  w <- ifelse(valid, pbar * value, -Inf); ib <- which.max(w); w_best <- w[ib]; pbar_best <- pbar[ib]; pk_best <- p_k[ib, ]
  a <- w_best > value0
  rej <- 1 - pbar_best
  post <- if (rej > 1e-12) priors * (1 - pk_best) / rej else priors
  p2_post <- as.vector(p_k %*% post); p2_prior <- pbar; surplus <- value - value0
  valid2 <- valid; if (cf$exclude_rejected_from_tier2) valid2[ib] <- FALSE
  S_post <- ifelse(valid2, p2_post * surplus, -Inf); S_prior <- ifelse(valid2, p2_prior * surplus, -Inf)
  if (a) { ii <- which.max(S_post); S_ii <- S_post[ii]; p2_used <- p2_post[ii] } else { ii <- which.max(S_prior); S_ii <- S_prior[ii]; p2_used <- p2_prior[ii] }
  found <- is.finite(S_ii) && S_ii > 0
  P1 <- if (a) pbar_best else 0; Pw <- if (a) (1 - pbar_best) * cf$attrition else 0
  reach2 <- if (a) (1 - pbar_best) * (1 - cf$attrition) else 1
  P2 <- reach2 * (if (found) p2_used else 0); P3 <- 1 - P1 - P2 - Pw
  cac_i <- cac_of(l0, r0, cf); g1 <- cond0 > 0
  tree_value <- P1 * value[ib] + P2 * value[ii] + P3 * value0
  gate <- if (cf$gate_mode == "tree" && g1) tree_value >= cac_i else w_best >= cac_i
  list(best_npv_w = w_best, best_value = value[ib], best_cond = cond[ib], best_real = real[ib], best_real_ps = real_ps[ib],
       best_l = Ln[ib], best_r = Rn[ib], best_term = Tn[ib], best_pd = pd_dyn[ib], best_ead = cc$ead[ib], best_pti = ptin[ib], best_p = pbar_best,
       best_pk1 = pk_best[1], best_pk2 = pk_best[length(pk_best)],
       found = found, ind_value = value[ii], ind_cond = cond[ii], ind_real = real[ii], ind_real_ps = real_ps[ii], ind_l = Ln[ii], ind_r = Rn[ii],
       ind_term = Tn[ii], ind_pd = pd_dyn[ii], ind_ead = cc$ead[ii], ind_p = p2_used, ind_pk1 = p_k[ii, 1], ind_pk2 = p_k[ii, ncol(p_k)],
       orig_value = value0, orig_cond = cond0, orig_real = real0, orig_pd = pd0, orig_ead = c0$ead, l_orig = l0, r_orig = r0, pti_orig = pti0,
       a = a, P1 = P1, P2 = P2, P3 = P3, Pw = Pw, tree_value = tree_value, cac = cac_i, gate = gate, g1 = g1, eta_i = eta_i)
}
build_contracts <- function(cf, em = emp, pdv = select_pd(cf)) {
  pd_ref <- median(clip01(pdv))   # rows -> numeric matrix (logicals become 0/1) -> data.frame (fast rbind)
  as.data.frame(do.call(rbind, lapply(seq_len(nrow(test)), function(i) unlist(optimize_contract(
    test$loan_amnt[i], test$int_rate[i], pdv[i], test$Default_numeric[i], test$annual_inc[i], cf, em, pd_ref)))))
}
expectation <- function(c) {
  g1 <- c$g1 == 1; gt <- c$gate == 1; s1 <- g1 & gt; s2 <- (!g1) & gt; n <- nrow(c); z <- numeric(n)
  ge <- z; gr <- z; gps <- z; gv <- z; cap <- z; pact <- z
  ge[s1] <- (c$P1 * c$best_cond + c$P2 * c$ind_cond + c$P3 * c$orig_cond)[s1]
  gr[s1] <- (c$P1 * c$best_real + c$P2 * c$ind_real + c$P3 * c$orig_real)[s1]
  gps[s1] <- (c$P1 * c$best_real_ps + c$P2 * c$ind_real_ps + c$P3 * c$orig_real)[s1]
  gv[s1] <- c$tree_value[s1]; cap[s1] <- (c$P1 * c$best_l + c$P2 * c$ind_l + c$P3 * c$l_orig)[s1]; pact[s1] <- (1 - c$Pw)[s1]
  ge[s2] <- (c$best_p * c$best_cond)[s2]; gr[s2] <- (c$best_p * c$best_real)[s2]; gps[s2] <- (c$best_p * c$best_real_ps)[s2]
  gv[s2] <- (c$best_p * c$best_value)[s2]; cap[s2] <- (c$best_p * c$best_l)[s2]; pact[s2] <- c$best_p[s2]
  list(gross_exp = sum(ge), gross_real = sum(gr), gross_ps = sum(gps), gross_value = sum(gv), cac = sum(ifelse(gt, c$cac, 0)),
       capital = sum(cap), approved = sum(pact), s1 = s1, s2 = s2, n_gate = sum(gt), n_g2_sourced = sum(s2), n_g2_rescued = sum(c$best_p[s2]),
       P1 = mean(c$P1[s1]), P2 = mean(c$P2[s1]), P3 = mean(c$P3[s1]))
}

# ---- risk: Vasicek on a directly instantiated 10k book ------------------------
risk_metrics <- function(losses, level = 0.999) { EL <- mean(losses); VaR <- quantile(losses, level, names = FALSE)
  tail <- losses[losses >= VaR]; list(EL = EL, VaR = VaR, CVaR = if (length(tail)) mean(tail) else VaR, EC = VaR - EL) }
vasicek_losses <- function(pd, ead, Z, cf, seed) { if (!length(pd)) return(rep(0, length(Z)))
  pd <- clip01(pd); rho <- rho_of(pd, cf$rho_ec_mode, cf$rho_ec_fixed); qn <- qnorm(pd); set.seed(seed)
  vapply(Z, function(z) sum(rbinom(length(pd), 1, pnorm((qn - sqrt(rho) * z) / sqrt(1 - rho))) * ead), numeric(1)) }
engine_mc <- function(c, ex, Z, cf) {
  idx <- rep(seq_len(nrow(c)), cf$scale); N <- length(idx); s1 <- ex$s1[idx]; s2 <- ex$s2[idx]; a <- c$a[idx] == 1; found <- c$found[idx] == 1
  pk1b <- c$best_pk1[idx]; pk2b <- c$best_pk2[idx]; pk1i <- c$ind_pk1[idx]; pk2i <- c$ind_pk2[idx]; ph <- if (cf$two_type) cf$type_prior_high else 0
  bc <- c$best_cond[idx]; ic <- c$ind_cond[idx]; oc <- c$orig_cond[idx]; br <- c$best_real[idx]; ir <- c$ind_real[idx]; orr <- c$orig_real[idx]
  bpd <- c$best_pd[idx]; ipd <- c$ind_pd[idx]; opd <- c$orig_pd[idx]; be <- c$best_ead[idx]; ie <- c$ind_ead[idx]; oe <- c$orig_ead[idx]
  cac10 <- ex$cac * cf$scale; n_mc <- length(Z); out <- matrix(NA_real_, n_mc, 4, dimnames = list(NULL, c("loss", "net_exp", "net_real", "active")))
  set.seed(cf$seed_mc)
  for (r in 1:n_mc) {
    high <- runif(N) < ph; p1 <- ifelse(high, pk2b, pk1b); p2 <- ifelse(high, pk2i, pk1i)
    u1 <- runif(N); u2 <- runif(N); u3 <- runif(N)
    t1 <- s1 & a & (u1 <= p1); walk <- s1 & a & !t1 & (u3 < cf$attrition); t2 <- s1 & !t1 & !walk & found & (u2 <= p2); t3 <- s1 & !t1 & !walk & !t2
    g2s <- s2 & (u1 <= p1); active <- t1 | t2 | t3 | g2s
    pdv <- ifelse(t1, bpd, ifelse(t2, ipd, ifelse(t3, opd, ifelse(g2s, bpd, 1e-4)))); eadv <- ifelse(t1, be, ifelse(t2, ie, ifelse(t3, oe, ifelse(g2s, be, 0))))
    rho <- rho_of(pdv, cf$rho_ec_mode, cf$rho_ec_fixed); pc <- pnorm((qnorm(clip01(pdv)) - sqrt(rho) * Z[r]) / sqrt(1 - rho))
    out[r, ] <- c(sum((rbinom(N, 1, pc) == 1 & active) * eadv), sum(t1 * bc + t2 * ic + t3 * oc + g2s * bc) - cac10,
                  sum(t1 * br + t2 * ir + t3 * orr + g2s * br) - cac10, sum(active))
  }
  out
}

# ---- evaluation of all regimes ---------------------------------------------
passive_row <- function(pd, cf, rule = "npv") {
  L <- test$loan_amnt; R <- test$int_rate; y <- y_te; r_m <- cf$r_disc / 12; cac <- cac_of(L, R, cf)
  cc <- pv_components(L, R, cf$term_orig, cf$tau_bar, cf$lgd, r_m); npv <- -L + (1 - pd) * cc$pv_nd + pd * cc$pv_d
  ap <- if (rule == "nam_k") (npv - cf$hurdle * basel_k(pd, cf$lgd, rho_of(pd, cf$rho_mode, cf$rho_fixed), cf$var_level) * L * balance_pv_years(L, R, cf$term_orig, r_m)) > cac
        else if (rule == "nam") npv > cac else npv > 0
  real <- -L + (1 - y) * cc$pv_nd + y * cc$pv_d; s <- cf$scale
  list(approved = sum(ap) * s, rate = mean(ap) * 100, capital = sum(L[ap]) * s, net_exp = (sum(npv[ap]) - sum(cac[ap])) * s,
       net_real = (sum(real[ap]) - sum(cac[ap])) * s, net_ps = NA, ap = ap, pd = pd, ead = cc$ead)
}
evaluate_all <- function(cf, do_mc = TRUE, n_mc = cf$n_mc, extra = list(), label = "") {
  set.seed(cf$seed_z); Z <- rnorm(n_mc); rows <- list(); dists <- list()
  pdsets <- list("FICO-Only Baseline" = test$pd_fico, "Bank Grade Model" = test$pd_grade, "Bank Sub-Grade Model" = test$pd_sub, "Random Forest ML" = select_pd(cf))
  sd <- 200
  for (nm in names(pdsets)) for (rule in c("npv", "nam", if (nm == "Random Forest ML") "nam_k")) {
    b <- passive_row(pdsets[[nm]], cf, rule); rm <- list(EL = NA, VaR = NA, CVaR = NA, EC = NA)
    if (do_mc) { sd <- sd + 1; rm <- risk_metrics(vasicek_losses(rep(b$pd[b$ap], cf$scale), rep(b$ead[b$ap], cf$scale), Z, cf, sd), cf$var_level) }
    rows[[length(rows) + 1]] <- data.frame(regime = paste0(nm, c(npv = "", nam = " +NAM", nam_k = " +NAM+K")[rule]), approved = b$approved, rate = b$rate,
      capital = b$capital, net_exp = b$net_exp, net_real = b$net_real, net_ps = NA, EL = rm$EL, VaR = rm$VaR, CVaR = rm$CVaR, EC = rm$EC,
      RAROC = if (!is.na(rm$EC) && rm$EC > 0) b$net_exp / rm$EC else NA, profit_per_dollar = b$net_exp / b$capital)
  }
  engines <- c(list("Multi-Stage Engine" = cf), extra)
  for (en in names(engines)) { ecf <- engines[[en]]; c <- build_contracts(ecf); ex <- expectation(c); s <- ecf$scale
    rm <- list(EL = NA, VaR = NA, CVaR = NA, EC = NA); ci <- c(NA, NA)
    if (do_mc) { mc <- engine_mc(c, ex, Z, ecf); rm <- risk_metrics(mc[, "loss"], ecf$var_level); ci <- quantile(mc[, "net_exp"], c(.025, .975)) }
    ne <- (ex$gross_exp - ex$cac) * s
    rows[[length(rows) + 1]] <- data.frame(regime = en, approved = ex$approved * s, rate = ex$approved / nrow(test) * 100, capital = ex$capital * s,
      net_exp = ne, net_real = (ex$gross_real - ex$cac) * s, net_ps = (ex$gross_ps - ex$cac) * s, EL = rm$EL, VaR = rm$VaR, CVaR = rm$CVaR, EC = rm$EC,
      RAROC = if (!is.na(rm$EC) && rm$EC > 0) ne / rm$EC else NA, profit_per_dollar = ne / (ex$capital * s))
    s1 <- ex$s1; dR <- (c$best_r - c$r_orig)[s1]
    dists[[en]] <- data.frame(engine = en, n_g1 = sum(s1), n_g2_sourced = ex$n_g2_sourced, n_g2_rescued = ex$n_g2_rescued, n_gate = ex$n_gate,
      mean_dR_g1 = mean(dR), share_hike_g1 = mean(dR > 0.25), share_term36_g1 = mean(c$best_term[s1] == 36), mean_frac_g1 = mean((c$best_l / c$l_orig)[s1]),
      mean_p_accept_g1 = mean(c$best_p[s1]), P1 = ex$P1, P2 = ex$P2, P3 = ex$P3, net_exp_ci_lo = ci[1], net_exp_ci_hi = ci[2],
      net_value_exp = (ex$gross_value - ex$cac) * s, gross_exp = ex$gross_exp * s, gross_real = ex$gross_real * s, cac_paid = ex$cac * s)
    if (en == "Multi-Stage Engine" && label == "base") { assign("contracts_base", c, envir = .GlobalEnv); assign("ex_base", ex, envir = .GlobalEnv) }
  }
  list(table = do.call(rbind, rows), dist = do.call(rbind, dists))
}
engine_net <- function(cf, em = emp) { c <- build_contracts(cf, em); ex <- expectation(c); list(net = (ex$gross_exp - ex$cac) * cf$scale, c = c, ex = ex) }
best_passive_net <- function(cf) { v <- c(); for (nm in c("pd_fico", "pd_grade", "pd_sub")) v[nm] <- max(passive_row(test[[nm]], cf, "npv")$net_exp, passive_row(test[[nm]], cf, "nam")$net_exp)
  p <- select_pd(cf); v["rf"] <- max(passive_row(p, cf, "npv")$net_exp, passive_row(p, cf, "nam")$net_exp, passive_row(p, cf, "nam_k")$net_exp); v }

# ============================ RUN ===========================================
options(width = 220); fmtM <- function(x) format(round(x), big.mark = ",")
cat("\n=== MASTER MATRIX v13 — BASE CASE (tau=12, LGD=0.8, calibrated RF, all fixes; rho_EC = Basel other-retail) ===\n")
extra <- list("Multi-Stage Engine (unconstrained)" = with_cfg(cfg, constrain = FALSE), "Multi-Stage Engine (no capital charge)" = with_cfg(cfg, capital_in_objective = FALSE))
res_base <- evaluate_all(cfg, TRUE, cfg$n_mc, extra, label = "base"); print(res_base$table, row.names = FALSE, digits = 4); print(res_base$dist, row.names = FALSE, digits = 4)
write.csv(res_base$table, "Master_Summary_Matrix_v13.csv", row.names = FALSE); write.csv(res_base$dist, "Engine_Contract_Distribution_v13.csv", row.names = FALSE)
stopifnot(!any(is.na(res_base$table$net_exp)), all(contracts_base$best_pd > 0 & contracts_base$best_pd < 1), nrow(test) == 2000,
          all(abs(contracts_base$P1 + contracts_base$P2 + contracts_base$P3 + contracts_base$Pw - 1) < 1e-9))
cat("Sanity assertions PASSED.\n")

cat("\n=== CONSERVATIVE BOUND (frozen loss spec tau=1, LGD=1; all other fixes on) ===\n")
res_cons <- evaluate_all(with_cfg(cfg, tau_bar = 1, lgd = 1.0), TRUE, cfg$n_mc / 2); print(res_cons$table, row.names = FALSE, digits = 4)
write.csv(res_cons$table, "Master_Summary_Matrix_v13_conservative.csv", row.names = FALSE)
cat("\n=== RHO = 0.15 STRESS (EC) ===\n")
res_rho <- evaluate_all(with_cfg(cfg, rho_ec_mode = "fixed"), TRUE, cfg$n_mc / 2); print(res_rho$table, row.names = FALSE, digits = 4)
write.csv(res_rho$table, "Master_Summary_Matrix_v13_rho015.csv", row.names = FALSE)
cat("\n=== FROZEN SPEC REPLICATION (sidur_12 mechanics inside v13 code; should match sidur_12 up to RF/split noise) ===\n")
res_frozen <- evaluate_all(do.call(with_cfg, c(list(cfg), frozen_overrides)), TRUE, 2000); print(res_frozen$table, row.names = FALSE, digits = 4)
write.csv(res_frozen$table, "Master_Summary_Matrix_v13_frozenspec.csv", row.names = FALSE)

cat("\n=== LOSS-PROFILE SCENARIOS (analytic) ===\n")
scen <- list(c(1, 1), c(6, 1), c(12, 1), c(12, 0.8), c(18, 0.7), c(24, 0.6)); srows <- list()
for (s in scen) { cf2 <- with_cfg(cfg, tau_bar = s[1], lgd = s[2]); e <- engine_net(cf2); bp <- best_passive_net(cf2)
  srows[[length(srows) + 1]] <- data.frame(tau_bar = s[1], lgd = s[2], engine_net_exp = e$net, engine_net_ps = (e$ex$gross_ps - e$ex$cac) * cfg$scale,
    engine_rate = e$ex$approved / nrow(test) * 100, best_passive = names(which.max(bp)), best_passive_net_exp = max(bp), ratio = e$net / max(bp)) }
scen_tab <- do.call(rbind, srows); print(scen_tab, row.names = FALSE, digits = 4); write.csv(scen_tab, "Scenario_Grid_v13.csv", row.names = FALSE)

cat("\n=== TORNADO (+/-25%) ===\n")
base_net <- engine_net(cfg)$net
tparams <- c("kappa", "lambda", "alpha", "phi", "gamma_cac", "eta0", "xi", "eta_cut", "hurdle", "r_disc", "tau_bar", "lgd", "gamma_emp", "delta_emp")
trows <- list()
for (p in tparams) {
  if (p %in% c("gamma_emp", "delta_emp")) { v <- sapply(c(0.75, 1.25), function(f) { em <- emp; em[[p]] <- em[[p]] * f; engine_net(cfg, em)$net }) }
  else { b <- cfg[[p]]; lo <- if (p == "tau_bar") max(1, round(b * 0.75)) else b * 0.75; hi <- if (p == "tau_bar") round(b * 1.25) else if (p == "lgd") min(1, b * 1.25) else b * 1.25
    v <- sapply(list(lo, hi), function(x) { cf2 <- cfg; cf2[[p]] <- x; engine_net(cf2)$net }) }
  trows[[length(trows) + 1]] <- data.frame(Parameter = p, Net_at_minus25 = v[1], Net_Base = base_net, Net_at_plus25 = v[2], Swing_USD = abs(v[2] - v[1]))
}
tornado <- do.call(rbind, trows); tornado <- tornado[order(-tornado$Swing_USD), ]; print(tornado, row.names = FALSE, digits = 4)
write.csv(tornado, "Engine_Sensitivity_Tornado_v13.csv", row.names = FALSE)

cat("\n=== BREAKEVEN DEMAND ELASTICITY (eta0) ===\n")
bp0 <- best_passive_net(cfg); brows <- list()
for (e0 in c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5)) { e <- engine_net(with_cfg(cfg, eta0 = e0))
  brows[[length(brows) + 1]] <- data.frame(eta0 = e0, price_elasticity_at_mean_rate = -e0 * R_bar, engine_net_exp = e$net, best_passive_net_exp = max(bp0),
    ratio = e$net / max(bp0), mean_dR_g1 = mean((e$c$best_r - e$c$r_orig)[e$ex$s1])) }
be <- do.call(rbind, brows); print(be, row.names = FALSE, digits = 4); write.csv(be, "Breakeven_Elasticity_v13.csv", row.names = FALSE)

cat("\n=== EXTRA SENSITIVITIES (attrition, xi, hurdle, max_hike) ===\n")
xrows <- list()
for (p in list(c("attrition", 0, 0.1, 0.25, 0.5), c("xi", 0, 0.25, 0.5, 1), c("hurdle", 0, 0.05, 0.1, 0.15, 0.2), c("max_hike", 2.5, 5, 7.5, 10)))
  for (v in as.numeric(p[-1])) { cf2 <- cfg; cf2[[p[1]]] <- v; e <- engine_net(cf2)
    xrows[[length(xrows) + 1]] <- data.frame(param = p[1], value = v, engine_net_exp = e$net, ratio = e$net / max(bp0), approved = e$ex$approved * cfg$scale,
      P1 = e$ex$P1, P2 = e$ex$P2, P3 = e$ex$P3, mean_dR_g1 = mean((e$c$best_r - e$c$r_orig)[e$ex$s1])) }
xs <- do.call(rbind, xrows); print(xs, row.names = FALSE, digits = 4); write.csv(xs, "Sensitivity_Extra_v13.csv", row.names = FALSE)

cat("\n=== AUC (OOS) — passives (PD) and engine profit score ===\n")
auc_of <- function(s) as.numeric(auc(roc(y_te, s, quiet = TRUE)))
auc_tab <- data.frame(Model = c("Integrated Logistic Rating", "Bank Sub-Grade Model", "Random Forest ML (Platt)", "Bank Grade Model", "FICO-Only Baseline", "Multi-Stage Engine (profit score)"),
  OOS_AUC = round(c(auc_of(test$pd_logit), auc_of(test$pd_sub), auc_of(test$pd_rf_platt), auc_of(test$pd_grade), auc_of(test$pd_fico), auc_of(contracts_base$best_npv_w)), 4))
print(auc_tab, row.names = FALSE); write.csv(auc_tab, "Model_AUC_Comparison_v13.csv", row.names = FALSE)

cat("\n=== APPENDIX: REPEATED RANDOM SPLITS (20 seeds; every model re-estimated on each train fold; analytic, no MC) ===\n")
sp <- list()
for (s_ in 0:19) { ok <- tryCatch({ fit_split(s_); TRUE }, error = function(e) { cat("  seed", s_, "skipped:", conditionMessage(e), "\n"); FALSE })
  if (!ok) next
  tt <- evaluate_all(cfg, do_mc = FALSE)$table; tt$seed <- s_
  pas <- tt[!grepl("^Multi", tt$regime), ]; e <- tt[tt$regime == "Multi-Stage Engine", ]
  tt$ratio_to_best_passive <- NA; tt$ratio_to_best_passive[tt$regime == "Multi-Stage Engine"] <- e$net_exp / max(pas$net_exp)
  tt$best_passive <- pas$regime[which.max(pas$net_exp)]
  sp[[length(sp) + 1]] <- tt[, c("seed", "regime", "approved", "rate", "capital", "net_exp", "net_real", "net_ps", "ratio_to_best_passive", "best_passive")]
  cat(sprintf("  seed %2d: engine %s | best passive %s (%s) | ratio %.3f\n", s_, fmtM(e$net_exp), fmtM(max(pas$net_exp)), tt$best_passive[1], e$net_exp / max(pas$net_exp)))
}
splits <- do.call(rbind, sp); write.csv(splits, "Repeated_Splits_v13.csv", row.names = FALSE)
summ <- do.call(rbind, lapply(split(splits, splits$regime), function(g) data.frame(regime = g$regime[1],
  net_exp_mean = mean(g$net_exp), net_exp_sd = stats::sd(g$net_exp), net_exp_q025 = quantile(g$net_exp, .025, names = FALSE),
  net_exp_q975 = quantile(g$net_exp, .975, names = FALSE), net_real_mean = mean(g$net_real), share_real_negative = mean(g$net_real < 0),
  net_ps_mean = mean(g$net_ps), ratio_mean = mean(g$ratio_to_best_passive), ratio_min = min(g$ratio_to_best_passive), ratio_max = max(g$ratio_to_best_passive))))
print(summ, row.names = FALSE, digits = 4); write.csv(summ, "Repeated_Splits_Summary_v13.csv", row.names = FALSE)
fit_split(cfg$seed_split)   # restore the canonical seed-42 split for anything run afterwards


# ============================ ADDENDUM (2026-09-25) =========================
# Inference, decomposition and identification. Append-only: nothing above changes, so every number
# printed before this line must reproduce the 2026-09-15 freeze to the dollar (checked below).
cat("\n=== ADDENDUM: freeze reproduction check ===\n")
stopifnot(abs(res_base$table$net_exp[res_base$table$regime == "Multi-Stage Engine"] - 12549535) < 5,
          abs(res_base$table$net_exp[res_base$table$regime == "FICO-Only Baseline +NAM"] - 8620023) < 5)
cat("Base-case engine and best-passive expected net reproduce the 2026-09-15 freeze.\n")

per_borrower_engine <- function(c) {            # expected net per test borrower (test scale), CAC charged if sourced
  g1 <- c$g1 == 1; gt <- c$gate == 1; s1 <- g1 & gt; s2 <- (!g1) & gt; ge <- numeric(nrow(c))
  ge[s1] <- (c$P1 * c$best_cond + c$P2 * c$ind_cond + c$P3 * c$orig_cond)[s1]
  ge[s2] <- (c$best_p * c$best_cond)[s2]
  ge - ifelse(gt, c$cac, 0)
}
per_borrower_passive <- function(pd, cf, rule) {
  L <- test$loan_amnt; R <- test$int_rate; r_m <- cf$r_disc / 12; cac <- cac_of(L, R, cf)
  cc <- pv_components(L, R, cf$term_orig, cf$tau_bar, cf$lgd, r_m); npv <- -L + (1 - pd) * cc$pv_nd + pd * cc$pv_d
  ap <- if (rule == "nam_k") (npv - cf$hurdle * basel_k(pd, cf$lgd, rho_of(pd, cf$rho_mode, cf$rho_fixed), cf$var_level) * L * balance_pv_years(L, R, cf$term_orig, r_m)) > cac
        else if (rule == "nam") npv > cac else npv > 0
  ifelse(ap, npv - cac, 0)
}
e_i <- per_borrower_engine(contracts_base)
stopifnot(abs(sum(e_i) * cfg$scale - res_base$table$net_exp[res_base$table$regime == "Multi-Stage Engine"]) < 1)

# ---- A1: zero acquisition cost (kappa = 0 for every regime) ----------------------------------
cat("\n=== A1: ZERO ACQUISITION COST (kappa = 0, all regimes; analytic) ===\n")
res_z <- evaluate_all(with_cfg(cfg, kappa = 0), do_mc = FALSE)
tz <- res_z$table[, c("regime", "approved", "rate", "capital", "net_exp", "net_real", "net_ps", "profit_per_dollar")]
pz <- tz[!grepl("^Multi", tz$regime), ]; ez <- tz$net_exp[tz$regime == "Multi-Stage Engine"]
tz$ratio_to_best_passive <- ifelse(tz$regime == "Multi-Stage Engine", ez / max(pz$net_exp), NA)
print(tz, row.names = FALSE, digits = 4); write.csv(tz, "ZeroCAC_v13.csv", row.names = FALSE)
cat(sprintf("Engine %s vs best passive %s (%s): ratio %.3f\n", fmtM(ez), fmtM(max(pz$net_exp)), pz$regime[which.max(pz$net_exp)], ez / max(pz$net_exp)))

# ---- A2: paired bootstrap over test borrowers (models and contracts held fixed) -----------------
cat("\n=== A2: PAIRED BOOTSTRAP OF THE ENGINE'S ADVANTAGE (B = 5,000 resamples of the 2,000 test borrowers) ===\n")
set.seed(2026); B <- 5000; n <- nrow(test); boot_idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B)
comp <- list("FICO-Only Baseline +NAM" = per_borrower_passive(test$pd_fico, cfg, "nam"),
             "Random Forest ML +NAM" = per_borrower_passive(select_pd(cfg), cfg, "nam"),
             "Random Forest ML +NAM+K" = per_borrower_passive(select_pd(cfg), cfg, "nam_k"),
             "Bank Sub-Grade Model +NAM" = per_borrower_passive(test$pd_sub, cfg, "nam"))
brows <- list(); eb <- apply(boot_idx, 1, function(ix) sum(e_i[ix])) * cfg$scale
for (nm in names(comp)) { p_i <- comp[[nm]]; pb <- apply(boot_idx, 1, function(ix) sum(p_i[ix])) * cfg$scale; d <- eb - pb; r <- eb / pb
  tt <- t.test(e_i, p_i, paired = TRUE)
  brows[[nm]] <- data.frame(comparator = nm, diff = (sum(e_i) - sum(p_i)) * cfg$scale, diff_lo = quantile(d, .025, names = FALSE), diff_hi = quantile(d, .975, names = FALSE),
    ratio = sum(e_i) / sum(p_i), ratio_lo = quantile(r, .025, names = FALSE), ratio_hi = quantile(r, .975, names = FALSE), p_boot_one_sided = mean(d <= 0),
    t_paired = unname(tt$statistic), p_t = tt$p.value, share_borrowers_engine_higher = mean(e_i > p_i)) }
boot_tab <- do.call(rbind, brows); print(boot_tab, row.names = FALSE, digits = 4); write.csv(boot_tab, "Bootstrap_Inference_v13.csv", row.names = FALSE)
cat(sprintf("Split-level sign test: engine beats best passive in %d of %d splits, one-sided p = %.2e\n",
            sum(splits$ratio_to_best_passive > 1, na.rm = TRUE), sum(!is.na(splits$ratio_to_best_passive)), 0.5^sum(!is.na(splits$ratio_to_best_passive))))

# ---- A3: value of each contract axis (grid restricted to a subset of rate / principal / term) ------
cat("\n=== A3: AXIS DECOMPOSITION (which contract dimensions create the value; analytic) ===\n")
bp0 <- best_passive_net(cfg); arows <- list()
for (m in c("none", "amount_term", "rate", "rate_term", "rate_amount", "all")) { e <- engine_net(with_cfg(cfg, grid_mask = m)); s1 <- e$ex$s1
  arows[[m]] <- data.frame(axes = m, engine_net_exp = e$net, engine_net_ps = (e$ex$gross_ps - e$ex$cac) * cfg$scale, ratio = e$net / max(bp0),
    approved = e$ex$approved * cfg$scale, mean_dR_g1 = mean((e$c$best_r - e$c$r_orig)[s1]), mean_frac_g1 = mean((e$c$best_l / e$c$l_orig)[s1]),
    share_term36_g1 = mean(e$c$best_term[s1] == 36)) }
axes_tab <- do.call(rbind, arows); print(axes_tab, row.names = FALSE, digits = 4); write.csv(axes_tab, "Axis_Decomposition_v13.csv", row.names = FALSE)

# ---- A4: robustness to the PD model feeding the engine ------------------------------------------
cat("\n=== A4: PD SOURCE ROBUSTNESS (engine and RF-based passives re-run with each PD; analytic) ===\n")
prow <- list()
for (src in c("rf_cal", "rf_raw", "logit")) { cf2 <- with_cfg(cfg, pd_source = src); e <- engine_net(cf2); b <- best_passive_net(cf2)
  prow[[src]] <- data.frame(pd_source = src, engine_net_exp = e$net, engine_net_ps = (e$ex$gross_ps - e$ex$cac) * cfg$scale, best_passive = names(which.max(b)),
    best_passive_net_exp = max(b), ratio = e$net / max(b), approved = e$ex$approved * cfg$scale) }
pd_tab <- do.call(rbind, prow); print(pd_tab, row.names = FALSE, digits = 4); write.csv(pd_tab, "PD_Source_Robustness_v13.csv", row.names = FALSE)

# ---- A5: identification of the rate channel — within-sub-grade (fixed-effects) estimate ----------
cat("\n=== A5: STIGLITZ-WEISS REGRESSION WITH SUB-GRADE FIXED EFFECTS (train) ===\n")
fit_fe <- glm(Default_numeric ~ int_rate + pmt_to_inc + fico + dti + revol_util + sub_grade, data = train, family = binomial)
cf_fe <- summary(fit_fe)$coefficients[c("int_rate", "pmt_to_inc", "fico", "dti", "revol_util"), ]; print(cf_fe)
r2_sg <- summary(lm(int_rate ~ sub_grade, data = train))$r.squared
sd_within <- mean(tapply(train$int_rate, train$sub_grade, sd), na.rm = TRUE)
cat(sprintf("Share of int_rate variance explained by sub_grade: %.4f; mean within-sub-grade SD %.3f pts (total SD %.3f)\n", r2_sg, sd_within, sd(train$int_rate)))
emp_fe <- list(gamma_emp = unname(coef(fit_fe)["int_rate"]), delta_emp = unname(coef(fit_fe)["pmt_to_inc"]))
e_fe <- engine_net(cfg, emp_fe); e_g0 <- engine_net(cfg, list(gamma_emp = 0, delta_emp = emp$delta_emp))
ci_fe <- emp_fe$gamma_emp + c(-1, 1) * 1.96 * cf_fe["int_rate", 2]
fe_tab <- data.frame(spec = c("pooled (base)", "sub-grade FE", "gamma = 0"), gamma_emp = c(emp$gamma_emp, emp_fe$gamma_emp, 0),
  gamma_se = c(summary(fit_emp)$coefficients["int_rate", 2], cf_fe["int_rate", 2], NA), delta_emp = c(emp$delta_emp, emp_fe$delta_emp, emp$delta_emp),
  engine_net_exp = c(base_net, e_fe$net, e_g0$net), ratio = c(base_net, e_fe$net, e_g0$net) / max(bp0),
  mean_dR_g1 = c(mean((contracts_base$best_r - contracts_base$r_orig)[ex_base$s1]), mean((e_fe$c$best_r - e_fe$c$r_orig)[e_fe$ex$s1]), mean((e_g0$c$best_r - e_g0$c$r_orig)[e_g0$ex$s1])))
print(fe_tab, row.names = FALSE, digits = 4); cat(sprintf("95%% CI for the FE rate coefficient: [%.4f, %.4f]\n", ci_fe[1], ci_fe[2]))
write.csv(cf_fe, "Empirical_StiglitzWeiss_FE_v13.csv"); write.csv(fe_tab, "Identification_Robustness_v13.csv", row.names = FALSE)

# ---- A6: who is re-priced — heterogeneity by grade and by income quartile ----------------------
cat("\n=== A6: HETEROGENEITY OF THE ENGINE'S TREATMENT (test book, base case) ===\n")
p_fn <- comp[["FICO-Only Baseline +NAM"]]; s1b <- ex_base$s1; s2b <- ex_base$s2; dRb <- contracts_base$best_r - contracts_base$r_orig
het_by <- function(g, lab) do.call(rbind, lapply(sort(unique(as.character(g))), function(v) { k <- as.character(g) == v; kk <- k & s1b
  data.frame(group = lab, level = v, n = sum(k), mean_pd = mean(contracts_base$orig_pd[k]), share_sourced = mean((s1b | s2b)[k]),
    mean_hike_g1 = if (any(kk)) mean(dRb[kk]) else NA, share_at_cap_g1 = if (any(kk)) mean(dRb[kk] >= cfg$max_hike - 1e-6) else NA,
    share_rate_mult_binding_g1 = if (any(kk)) mean(abs(contracts_base$best_r[kk] - cfg$max_rate_mult * contracts_base$r_orig[kk]) < 0.25) else NA,
    mean_accept_t1 = if (any(kk)) mean(contracts_base$best_p[kk]) else NA, engine_net = sum(e_i[k]) * cfg$scale, fico_nam_net = sum(p_fn[k]) * cfg$scale) }))
inc_q <- cut(test$annual_inc, quantile(test$annual_inc, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
grade_g <- as.character(test$grade); grade_g[grade_g %in% c("F", "G")] <- "F-G"
het_tab <- rbind(het_by(grade_g, "grade"), het_by(inc_q, "income quartile")); print(het_tab, row.names = FALSE, digits = 3)
write.csv(het_tab, "Heterogeneity_v13.csv", row.names = FALSE)

# ---- A7: are the passive models' AUC differences statistically meaningful? (DeLong 1988) ---------
cat("\n=== A7: DeLONG TESTS OF PAIRED AUC DIFFERENCES (test set) ===\n")
rocs <- list(FICO = roc(y_te, test$pd_fico, quiet = TRUE), Grade = roc(y_te, test$pd_grade, quiet = TRUE), SubGrade = roc(y_te, test$pd_sub, quiet = TRUE),
             RF_Platt = roc(y_te, test$pd_rf_platt, quiet = TRUE), Logit = roc(y_te, test$pd_logit, quiet = TRUE))
pairs <- list(c("SubGrade", "FICO"), c("Grade", "FICO"), c("RF_Platt", "FICO"), c("SubGrade", "Grade"), c("RF_Platt", "SubGrade"), c("Logit", "SubGrade"), c("Logit", "RF_Platt"))
dl <- do.call(rbind, lapply(pairs, function(p) { t <- roc.test(rocs[[p[1]]], rocs[[p[2]]], method = "delong", paired = TRUE)
  data.frame(model_a = p[1], model_b = p[2], auc_a = as.numeric(auc(rocs[[p[1]]])), auc_b = as.numeric(auc(rocs[[p[2]]])), z = unname(t$statistic), p_value = t$p.value) }))
print(dl, row.names = FALSE, digits = 4); write.csv(dl, "DeLong_AUC_Tests_v13.csv", row.names = FALSE)
# ---- A8 (2026-09-26): likelihood evidence and optimality of the grid solution --------------------
# A8a  which PD model fits the data best? out-of-sample log-likelihood, log-loss, Brier, pseudo-R2, in-sample AIC/BIC,
#      Vuong (1989) tests of paired per-borrower log-likelihood differences; set beside the profit each PD earns.
# A8b  likelihood-ratio tests of the Stiglitz-Weiss channel (rate and payment burden) in the default logit,
#      pooled and with sub-grade fixed effects; out-of-sample log-likelihood gain of the unrestricted model.
# A8c  optimality gap of the 1,178-contract grid: the engine re-solved on 2x and 5x finer grids.
cat("\n=== A8a: LIKELIHOOD COMPARISON OF PD MODELS (test set; probabilities clipped to [1e-4, 1-1e-4] as in the engine) ===\n")
ll_obs <- function(p, y = y_te) { p <- clip01(p); y * log(p) + (1 - y) * log(1 - p) }
pd_cols <- c(FICO = "pd_fico", Grade = "pd_grade", SubGrade = "pd_sub", Logit = "pd_logit", RF_raw = "pd_rf_raw", RF_Platt = "pd_rf_platt")
glm_of <- list(FICO = m_fico, Grade = m_grade, SubGrade = m_sub, Logit = m_logit)
eng_src <- c(Logit = "logit", RF_raw = "rf_raw", RF_Platt = "rf_cal")
LL0 <- sum(ll_obs(rep(mean(y_tr), length(y_te))))            # null model: train default rate
lik_tab <- do.call(rbind, lapply(names(pd_cols), function(k) { p <- test[[pd_cols[k]]]; l <- ll_obs(p); g <- glm_of[[k]]
  data.frame(model = k, oos_loglik = sum(l), log_loss = -mean(l), brier = mean((clip01(p) - y_te)^2), pseudo_r2_oos = 1 - sum(l) / LL0,
    auc = as.numeric(auc(roc(y_te, p, quiet = TRUE))), n_par = if (is.null(g)) NA else length(coef(g)),
    aic_train = if (is.null(g)) NA else AIC(g), bic_train = if (is.null(g)) NA else BIC(g),
    passive_nam_net = sum(per_borrower_passive(p, cfg, "nam")) * cfg$scale,
    engine_net_with_this_pd = if (k %in% names(eng_src)) pd_tab$engine_net_exp[pd_tab$pd_source == eng_src[[k]]] else NA) }))
lik_tab <- lik_tab[order(-lik_tab$oos_loglik), ]; print(lik_tab, row.names = FALSE, digits = 5)
write.csv(lik_tab, "Likelihood_PD_Models_v13.csv", row.names = FALSE)
vuong <- function(a, b) { d <- ll_obs(test[[pd_cols[a]]]) - ll_obs(test[[pd_cols[b]]]); z <- sqrt(length(d)) * mean(d) / sd(d)
  data.frame(model_a = a, model_b = b, loglik_a = sum(ll_obs(test[[pd_cols[a]]])), loglik_b = sum(ll_obs(test[[pd_cols[b]]])),
             mean_diff = mean(d), z = z, p_two_sided = 2 * pnorm(-abs(z))) }
vu_pairs <- list(c("RF_Platt", "FICO"), c("RF_Platt", "SubGrade"), c("RF_Platt", "Logit"), c("RF_Platt", "RF_raw"),
                 c("SubGrade", "FICO"), c("Logit", "SubGrade"), c("Logit", "FICO"))
vu_tab <- do.call(rbind, lapply(vu_pairs, function(p) vuong(p[1], p[2]))); print(vu_tab, row.names = FALSE, digits = 4)
write.csv(vu_tab, "Vuong_Tests_v13.csv", row.names = FALSE)
cat(sprintf("Rank correlation (Spearman) between OOS log-likelihood and passive +NAM net across the %d PD models: %.3f\n",
            nrow(lik_tab), cor(lik_tab$oos_loglik, lik_tab$passive_nam_net, method = "spearman")))

cat("\n=== A8b: LIKELIHOOD-RATIO TESTS OF THE STIGLITZ-WEISS CHANNEL (default logit, train; OOS gain on test) ===\n")
te2 <- test; te2$pmt_to_inc <- pmt_of(te2$loan_amnt, te2$int_rate, 60) / (te2$annual_inc / 12); te2$dti <- te2$revol_bal / te2$annual_inc * 100
sw_rhs <- c("int_rate", "pmt_to_inc", "fico", "dti", "revol_util")
restr <- list("drop int_rate + pmt_to_inc" = c("int_rate", "pmt_to_inc"), "drop int_rate" = "int_rate", "drop pmt_to_inc" = "pmt_to_inc")
lr_rows <- list()
for (spec in c("pooled", "sub-grade FE")) {
  fe_term <- if (spec == "pooled") character(0) else "sub_grade"
  f_full <- glm(reformulate(c(sw_rhs, fe_term), "Default_numeric"), data = train, family = binomial)
  for (nm in names(restr)) {
    f_r <- glm(reformulate(c(setdiff(sw_rhs, restr[[nm]]), fe_term), "Default_numeric"), data = train, family = binomial)
    LR <- 2 * (as.numeric(logLik(f_full)) - as.numeric(logLik(f_r))); dfr <- length(restr[[nm]])
    d <- ll_obs(predict(f_full, newdata = te2, type = "response")) - ll_obs(predict(f_r, newdata = te2, type = "response"))
    lr_rows[[length(lr_rows) + 1]] <- data.frame(spec = spec, restriction = nm, LR = LR, df = dfr, p_value = pchisq(LR, dfr, lower.tail = FALSE),
      dAIC_full_minus_restricted = AIC(f_full) - AIC(f_r), dBIC_full_minus_restricted = BIC(f_full) - BIC(f_r),
      oos_loglik_gain = sum(d), oos_vuong_z = sqrt(length(d)) * mean(d) / sd(d))
  }
}
lr_tab <- do.call(rbind, lr_rows); print(lr_tab, row.names = FALSE, digits = 4); write.csv(lr_tab, "LR_Tests_StiglitzWeiss_v13.csv", row.names = FALSE)

cat("\n=== A8c: OPTIMALITY GAP OF THE CONTRACT GRID (engine re-solved on finer grids; analytic) ===\n")
GRID_base <- GRID; grow <- list(); w0 <- contracts_base$best_npv_w; s1b <- ex_base$s1
for (g in list(c(0.05, 0.5), c(0.025, 0.25), c(0.01, 0.1))) { t0 <- Sys.time()
  GRID <- if (g[1] == 0.05) GRID_base else expand.grid(pct = round(seq(1.00, 0.10, by = -g[1]), 6), term = c(36, 60), shift = round(seq(-5, 10, by = g[2]), 6))
  e <- engine_net(cfg); w1 <- e$c$best_npv_w; gap <- ((w1 - w0) / abs(w0))[s1b]
  grow[[length(grow) + 1]] <- data.frame(pct_step = g[1], shift_step = g[2], n_contracts = nrow(GRID), engine_net_exp = e$net, gain_vs_base = e$net - base_net,
    rel_gain_vs_base = e$net / base_net - 1, ratio_to_best_passive = e$net / max(bp0), approved = e$ex$approved * cfg$scale,
    tier1_obj_gap_mean = mean(gap), tier1_obj_gap_median = median(gap), tier1_obj_gap_p95 = quantile(gap, .95, names = FALSE), share_g1_gap_above_1pct = mean(gap > 0.01),
    mean_dR_g1 = mean((e$c$best_r - e$c$r_orig)[e$ex$s1]), minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  cat(sprintf("  grid %5d contracts: engine %s (%+.3f%% vs base), ratio %.3f\n", nrow(GRID), fmtM(e$net), 100 * (e$net / base_net - 1), e$net / max(bp0)))
}
GRID <- GRID_base   # restore the official grid
gap_tab <- do.call(rbind, grow); print(gap_tab, row.names = FALSE, digits = 5); write.csv(gap_tab, "Optimality_Gap_v13.csv", row.names = FALSE)
stopifnot(abs(gap_tab$engine_net_exp[1] - base_net) < 1, nrow(GRID) == 1178)
cat("\n=== END OF ADDENDUM ===\n")


pdf("active_credit_engine_plots.pdf", width = 9, height = 6)
plot(roc(y_te, test$pd_rf_platt, quiet = TRUE), col = "darkred", lwd = 2, legacy.axes = TRUE, main = "ROC: calibrated RF vs engine profit score")
plot(roc(y_te, contracts_base$best_npv_w, quiet = TRUE), col = "darkgreen", lwd = 2, add = TRUE)
hist((contracts_base$best_r - contracts_base$r_orig)[ex_base$s1], breaks = seq(-5.25, 10.25, 0.5), col = "steelblue", main = "Tier-1 rate shift, Group 1 (base)", xlab = "pts")
barplot(tornado$Swing_USD / 1e6, names.arg = tornado$Parameter, las = 2, main = "Tornado swing ($M)")
plot(be$eta0, be$ratio, log = "x", type = "b", main = "Engine / best passive vs eta0", xlab = "eta0", ylab = "ratio"); abline(h = 1, lty = 2)
dev.off()
cat(sprintf("\nElapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins")))); print(sessionInfo()); sink()
