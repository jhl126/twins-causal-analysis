import numpy as np
from sklearn.linear_model import LogisticRegression

def doubly_robust_ate(X, T, Y):
    """Estimate the ATE using the AIPW doubly robust estimator.

    Both the propensity model and the outcome model are logistic regressions,
    appropriate for the binary treatment and binary outcome. Returns the
    point estimate of the average treatment effect.
    """
    X = np.asarray(X)
    T = np.asarray(T)
    Y = np.asarray(Y)

    # Propensity model: P(T=1 | X)
    ps_model = LogisticRegression(penalty=None, max_iter=1000)
    ps_model.fit(X, T)
    ps = ps_model.predict_proba(X)[:, 1]

    # Outcome model: P(Y=1 | X, T)
    # Concatenate T as a column onto X so treatment is a predictor
    XT = np.column_stack([X, T])
    outcome_model = LogisticRegression(penalty=None, max_iter=1000)
    outcome_model.fit(XT, Y)

    # Counterfactual predictions
    # mu1: predict with T set to 1 for everyone
    # mu0: predict with T set to 0 for everyone
    XT1 = np.column_stack([X, np.ones(len(T))])
    XT0 = np.column_stack([X, np.zeros(len(T))])
    mu1 = outcome_model.predict_proba(XT1)[:, 1]
    mu0 = outcome_model.predict_proba(XT0)[:, 1]

    # AIPW combination
    ate = (np.mean(T * (Y - mu1) / ps + mu1) - np.mean((1 - T) * (Y - mu0) / (1 - ps) + mu0))

    return ate

def cluster_bootstrap_ci(X, T, Y, pair_id, n_boot=500, seed=42):
    """Compute a 95% confidence interval for the DR ATE via cluster bootstrap.

    Resamples whole pairs (not individual twins) with replacement to preserve
    the within-pair correlation structure, recomputing the ATE on each
    resample. Returns the point estimate and the 2.5th/97.5th percentile CI.
    """
    rng = np.random.default_rng(seed)

    # Convert to arrays
    X = np.asarray(X)
    T = np.asarray(T)
    Y = np.asarray(Y)
    pair_id = np.asarray(pair_id)

    unique_pairs = np.unique(pair_id)
    n_pairs = len(unique_pairs)

    boot_estimates = []
    for b in range(n_boot):
        # Sample pairs with replacement
        sampled_pairs = rng.choice(unique_pairs, size=n_pairs, replace=True)

        # Build row indices for the sampled pairs
        rows = np.concatenate([np.where(pair_id == p)[0] for p in sampled_pairs])

        # Recompute the ATE on the resampled data
        ate_b = doubly_robust_ate(X[rows], T[rows], Y[rows])
        boot_estimates.append(ate_b)

    boot_estimates = np.array(boot_estimates)
    point = doubly_robust_ate(X, T, Y)
    ci_low, ci_high = np.percentile(boot_estimates, [2.5, 97.5])

    return point, ci_low, ci_high, boot_estimates

def regression_adjusted_ate(X, T, Y):
    """Estimate the ATE using outcome regression only (g-computation).

    Fits a single logistic outcome model on covariates and treatment, then
    averages the difference between predicted mortality under treatment and
    under control. Uses no propensity weighting.
    """
    X = np.asarray(X)
    T = np.asarray(T)
    Y = np.asarray(Y)

    # Fit outcome model: P(Y=1 | X, T)
    XT = np.column_stack([X, T])
    outcome_model = LogisticRegression(penalty=None, max_iter=1000)
    outcome_model.fit(XT, Y)

    # Counterfactual predictions: everyone treated vs everyone control
    XT1 = np.column_stack([X, np.ones(len(T))])
    XT0 = np.column_stack([X, np.zeros(len(T))])
    mu1 = outcome_model.predict_proba(XT1)[:, 1]
    mu0 = outcome_model.predict_proba(XT0)[:, 1]

    # ATE is the average difference in predicted outcomes
    return np.mean(mu1) - np.mean(mu0)

def naive_ate(T, Y):
    """Unadjusted difference in mean outcomes between treated and control."""
    T = np.asarray(T)
    Y = np.asarray(Y)
    return Y[T == 1].mean() - Y[T == 0].mean()

def t_learner(X, T, Y):
    """Estimate individual treatment effects via the T-Learner.

    Fits separate logistic outcome models on the treated and control arms,
    then predicts mortality under both models for every unit. Returns the
    estimated treatment effect tau(x) = mu1(x) - mu0(x) for each observation.
    """
    X = np.asarray(X)
    T = np.asarray(T)
    Y = np.asarray(Y)

    # Filter to treated and control rows (for training)
    X_treated, Y_treated = X[T == 1], Y[T == 1]
    X_control, Y_control = X[T == 0], Y[T == 0]

    # Train mu1 on treated arm, mu0 on control arm
    mu1_model = LogisticRegression(penalty=None, max_iter=1000)
    mu1_model.fit(X_treated, Y_treated)

    mu0_model = LogisticRegression(penalty=None, max_iter=1000)
    mu0_model.fit(X_control, Y_control)

    # Predict mortality under both models for the FULL dataset
    mu1 = mu1_model.predict_proba(X)[:, 1]
    mu0 = mu0_model.predict_proba(X)[:, 1]

    # calculate tau(x) = mu1(x) - mu0(x) for every twin
    tau = mu1 - mu0

    return tau

def t_learner_subgroup_bootstrap(X, T, Y, exposure, pair_id, n_boot=500, seed=42):
    """Bootstrap the subgroup CATEs and their difference for the T-Learner.

    Resamples whole pairs with replacement (preserving within-pair correlation),
    refits the T-Learner on each resample, and computes the exposed CATE, the
    unexposed CATE, and their difference. Returns point estimates and 95%
    percentile confidence intervals for all three quantities.
    """
    rng = np.random.default_rng(seed)

    X = np.asarray(X)
    T = np.asarray(T)
    Y = np.asarray(Y)
    exposure = np.asarray(exposure)
    pair_id = np.asarray(pair_id)

    unique_pairs = np.unique(pair_id)
    n_pairs = len(unique_pairs)

    # Pre-build a lookup from pair -> its row indices (faster than np.where in loop)
    pair_to_rows = {p: np.where(pair_id == p)[0] for p in unique_pairs}

    boot_exposed, boot_unexposed, boot_diff = [], [], []

    for b in range(n_boot):
        # Resample pairs with replacement
        sampled_pairs = rng.choice(unique_pairs, size=n_pairs, replace=True)
        rows = np.concatenate([pair_to_rows[p] for p in sampled_pairs])

        Xb, Tb, Yb, eb = X[rows], T[rows], Y[rows], exposure[rows]

        # Refit T-Learner on the resample
        tau_b = t_learner(Xb, Tb, Yb)

        # Subgroup CATEs on the resample
        ce = tau_b[eb == 1].mean()
        cu = tau_b[eb == 0].mean()

        boot_exposed.append(ce)
        boot_unexposed.append(cu)
        boot_diff.append(ce - cu)

    boot_exposed = np.array(boot_exposed)
    boot_unexposed = np.array(boot_unexposed)
    boot_diff = np.array(boot_diff)

    # Point estimates on the full data
    tau_full = t_learner(X, T, Y)
    point_exposed = tau_full[exposure == 1].mean()
    point_unexposed = tau_full[exposure == 0].mean()
    point_diff = point_exposed - point_unexposed

    results = {
        'exposed': (point_exposed, *np.percentile(boot_exposed, [2.5, 97.5])),
        'unexposed': (point_unexposed, *np.percentile(boot_unexposed, [2.5, 97.5])),
        'difference': (point_diff, *np.percentile(boot_diff, [2.5, 97.5])),
    }
    return results, boot_diff