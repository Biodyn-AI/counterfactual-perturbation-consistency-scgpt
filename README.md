# Counterfactual Perturbation Consistency for scGPT

This repository contains code and result artifacts for evaluating whether scGPT causal intervention effects are consistent with perturbation-derived gene expression changes.

## Contents
- `analysis/run_counterfactual_consistency_study.py`: main analysis pipeline
- `analysis/test_counterfactual_consistency.py`: unit tests for analysis helpers
- `outputs/tables/`: generated metric tables and matched-pair exports
- `outputs/figures/`: generated figures
- `paper/`: project report draft

## Reproducibility
The original computations were run in a Python 3.11 environment with `numpy`, `pandas`, `scanpy`, `matplotlib`, `pyyaml`, and `torch` available.

To rerun from a compatible workspace containing the required perturbation data and scGPT assets:

```bash
python analysis/run_counterfactual_consistency_study.py --skip-causal --shuffle-n 400 --bootstrap-n 1000
```

If causal score files are not already present, remove `--skip-causal` and provide the expected `single_cell_mechinterp` dependency tree.

## Data access
The analysis expects processed perturbation datasets (Adamson, Dixit, Dixit 7-day, Shifrut) and scGPT assets aligned with the configuration paths used in the script.
