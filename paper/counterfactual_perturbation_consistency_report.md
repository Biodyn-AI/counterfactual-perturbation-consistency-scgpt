# Counterfactual Perturbation Consistency Report

## Question
Do scGPT causal intervention effects align with real perturbation outcomes at matched gene-pair level?

## Evidence
- adamson (ablation): raw |rho|=-0.069, 95% CI [-0.416, 0.309], controlled |rho|=-0.045, sign agreement=0.484, sign 95% CI [0.323, 0.645], null p(|rho|)=0.5306, matched pairs=31
- dixit_13d (ablation): raw |rho|=0.051, 95% CI [-0.289, 0.379], controlled |rho|=0.282, sign agreement=0.222, sign 95% CI [0.053, 0.440], null p(|rho|)=0.3865, matched pairs=40
- dixit_7d (ablation): raw |rho|=0.230, 95% CI [-0.112, 0.520], controlled |rho|=0.187, sign agreement=0.500, sign 95% CI [0.267, 0.727], null p(|rho|)=0.2170, matched pairs=40
- shifrut (ablation): raw |rho|=0.064, 95% CI [-0.255, 0.359], controlled |rho|=0.093, sign agreement=0.250, sign 95% CI [0.000, 0.770], null p(|rho|)=0.5403, matched pairs=40

## Baseline Perturbation Benchmark
- adamson (ablation): AUPR=0.576, AUROC=0.483, n_pairs=31, n_pos=18
- dixit_13d (ablation): AUPR=0.839, AUROC=0.550, n_pairs=40, n_pos=33
- dixit_7d (ablation): AUPR=0.772, AUROC=0.589, n_pairs=40, n_pos=28
- shifrut (ablation): AUPR=0.794, AUROC=0.406, n_pairs=40, n_pos=36

## Timepoint Robustness
- dixit_13d_vs_7d (ablation): n_overlap=40, Spearman |delta|=0.387, Spearman signed delta=0.075, delta sign agreement=0.250

## Interpretation
Across evaluated rows, mean raw |rho|=0.069 and mean controlled |rho|=0.129.
Control adjustment improved |rho| in 3/4 evaluable rows.
This supports partial counterfactual consistency while indicating that source/target coverage and confounding still limit alignment.

## Limitations
- Causal scores are from low-compute pair subsets, not full perturbation universes.
- Composition control is coarse (group-level distribution shift), not cell-level causal adjustment.
- Shuffled-label nulls are source-permutation controls, not full generative null models.
