# Paper Quality Delivery: Counterfactual Perturbation Consistency

## A) Submission Readiness Verdict
**VERDICT: NOT READY**

### Central claim (one sentence)
Current low-compute scGPT intervention outputs show only weak and uncertain counterfactual alignment with perturbation-derived gene-level effects across Adamson, Dixit, and Shifrut settings.

### Main contributions (3-5)
1. A reproducible counterfactual consistency evaluation framework that directly matches intervention scores to perturbation deltas at source-target pair level.
2. Multi-dataset consistency quantification (Adamson, Dixit 13-day, Dixit 7-day, Shifrut) with rank, sign, and confound-adjusted metrics.
3. Falsification controls via source-shuffled null distributions and bootstrap uncertainty intervals.
4. Cross-timepoint robustness analysis for overlapping Dixit pairs.

### Evidence for each contribution
- Contribution 1: End-to-end pipeline and matched-pair exports are implemented and versioned in the external repository.
- Contribution 2: `counterfactual_consistency_metrics.tsv` reports ablation consistency ranging from -0.069 to 0.230 across datasets.
- Contribution 3: Null tests (`counterfactual_consistency_null.tsv`) show no ablation dataset significantly above shuffled controls.
- Contribution 4: `timepoint_consistency.tsv` reports moderate |delta| timepoint concordance (0.387) but weak signed concordance (0.075).

### Top 5 likely reviewer objections
1. Pair-level consistency is weak and mostly not significant versus null.
2. Pair coverage is too small to support strong mechanistic claims.
3. Confound control is source-level and does not fully model cell-level covariates.
4. Biological validation is hypothesis-level (few high-confidence concordant pairs) rather than confirmatory.
5. Results could reflect intervention sensitivity limits rather than true mechanistic mismatch.

### Research gaps (required because VERDICT=NOT READY)

#### Gap 1: Statistical power and effect stability
- Why it matters: without narrower uncertainty and stronger null separation, central faithfulness claims remain unsupported.
- Minimal additional analysis/experiment: increase pair and cell budgets for each dataset; rerun with multi-seed intervention estimates.
- Expected outcome patterns: if true alignment exists, confidence intervals should tighten around positive consistency and null p-values should decrease.
- Stop condition: at least two independent datasets show positive ablation consistency with null-rejecting empirical p-values and stable CIs.

#### Gap 2: Confound robustness at cell level
- Why it matters: source-level composition shifts can miss within-group technical and biological confounds.
- Minimal additional analysis/experiment: fit cell-level adjustment models (e.g., library size, batch, donor/cell state proxies) before pair aggregation.
- Expected outcome patterns: robust mechanistic signal should persist after richer adjustment; confounded signal should collapse.
- Stop condition: main directional conclusions hold under at least two distinct confound-adjustment strategies.

#### Gap 3: Biological external validation breadth
- Why it matters: isolated concordant pairs do not support broad biological mechanism claims.
- Minimal additional analysis/experiment: test enrichment of concordant pairs against external pathway/regulatory resources.
- Expected outcome patterns: biologically meaningful consistency should concentrate in coherent pathways/axes.
- Stop condition: concordant set enrichment remains significant under null-preserving permutation controls.

### Missing research executed in this pass
1. Added 1000-draw bootstrap confidence intervals for rank/sign consistency.
2. Added 400-draw source-shuffled null controls with empirical p-values.
3. Added timepoint robustness analysis (Dixit 13-day vs 7-day overlap).
4. Exported matched source-target tables and high-confidence concordant pair table for direct audit.

These additions improved integrity and uncertainty reporting but did not close the readiness gaps above.

## B) Top 5 Reviewer Objections and Responses
1. **Objection:** "Your main claim is not statistically supported."
   **Response:** We now explicitly report null-based empirical p-values and bootstrap CIs; claims were weakened to "partial/weak consistency" and submission readiness was downgraded to NOT READY.

2. **Objection:** "Confounding could explain the observed signal."
   **Response:** We added confound-adjusted residualized consistency metrics (baseline expression + composition shift) and report both raw and adjusted values side-by-side.

3. **Objection:** "Cross-condition robustness is missing."
   **Response:** We added an explicit Dixit 13-day vs 7-day overlap analysis showing magnitude-level but not direction-level robustness.

4. **Objection:** "This is an aggregate metric story without pair-level traceability."
   **Response:** We now export per-dataset matched pair tables and a high-confidence concordance table for pair-level inspection.

5. **Objection:** "Reproducibility is unclear."
   **Response:** We published a public repository with analysis code, outputs, and reproducibility-oriented README: https://github.com/Biodyn-AI/counterfactual-perturbation-consistency-scgpt

## C) Revised Paper in Full (Submission-Style Draft)
# Counterfactual Perturbation Consistency of scGPT Causal Interventions Across CRISPR Single-Cell Screens

## Abstract
Mechanistic interpretability for single-cell foundation models is often evaluated against curated network references, but this leaves a central question unresolved: do model interventions align with experimentally observed counterfactual perturbation effects at the gene level? We evaluated this question for scGPT by comparing intervention-derived source-target effects against perturbation-derived expression deltas in four settings: Adamson, Dixit (13-day), Dixit (7-day), and Shifrut CRISPR screens. We measured rank consistency, sign agreement, confound-adjusted consistency (baseline expression and composition shift), and source-shuffled null controls with bootstrap uncertainty. Across ablation interventions, raw consistency was modest (Spearman magnitude range: -0.069 to 0.230), confidence intervals were broad, and no dataset outperformed shuffled nulls at conventional significance. Confound adjustment increased consistency in three of four datasets, but the effect remained small and unstable. Cross-timepoint robustness in Dixit showed moderate agreement in effect magnitude ranking (Spearman on |delta| = 0.387) but weak signed agreement (Spearman = 0.075; sign agreement = 0.25). These results indicate partial but weak counterfactual alignment under current low-compute pair coverage, supporting cautious claims and motivating higher-coverage interventions with stronger confound controls.

## Introduction
Single-cell foundation models have made rapid progress in representation quality and transfer performance, but mechanistic evidence for their biological validity remains limited [1]. Most evaluation pipelines in this area rely on overlap with curated transcriptional edges, which can reward statistical alignment without demonstrating counterfactual faithfulness. For mechanistic claims to be biologically credible, intervention effects in model space should show directional or rank-level agreement with measured perturbation outcomes.

Perturb-seq style datasets provide exactly this opportunity by measuring expression changes after targeted perturbations [3-5]. The key challenge is that perturbation data are noisy, condition-specific, and confounded by composition and baseline expression effects [9,10]. A useful validation framework therefore needs both agreement metrics and controls that can reject trivial explanations.

This study asks: when scGPT gene-level interventions are applied to perturbation-matched source-target pairs, do model effects align with observed perturbation deltas strongly enough to support mechanistic claims?

## Related Work
Foundation-scale single-cell modeling has expanded rapidly, including scGPT and related transfer-learning frameworks [1]. In parallel, mechanistic interpretability research has emphasized causal interventions, circuit-level hypotheses, and falsification controls rather than purely descriptive attribution [6,7].

In single-cell biology, pooled CRISPR perturbation screens (including Adamson, Dixit, and Shifrut) are widely used to infer causal regulatory consequences at transcriptome scale [3-5]. Best-practice single-cell analysis also highlights the need to model technical and compositional confounders when drawing biological conclusions [9,10].

Our work sits at the intersection: we treat perturbation response as a counterfactual target for intervention-based mechanistic validation, and we explicitly include negative controls and confound-aware diagnostics.

## Methods
### Datasets
We used four processed perturbation datasets: Adamson [3], Dixit 13-day [4], Dixit 7-day [4], and Shifrut [5]. For each dataset, perturbation and control cells were identified from dataset-specific condition labels.

### Model and Intervention Outputs
We reused a pre-trained whole-human scGPT model and existing intervention pipeline [1]. Two intervention types were evaluated where available:
- `ablation`: value ablation at source-gene positions.
- `swap`: random source-value swap control intervention.

### Pair Matching and Counterfactual Targets
For each run, causal score tables were matched to perturbation-derived source-target deltas on exact gene-symbol overlap. Matched pair counts were:
- Adamson: 31 pairs,
- Dixit 13-day: 40 pairs,
- Dixit 7-day: 40 pairs,
- Shifrut: 40 pairs.

### Metrics
Primary metrics:
- Spearman correlation between `|effect_mean|` and `|delta|`.
- Spearman correlation between signed `effect_mean` and signed `delta`.
- Sign agreement between non-zero effect and non-zero delta signs.

Confound-aware metric:
- Residualized Spearman on absolute effects and deltas after regressing out log baseline expression and source-level composition shift (L1 distance from control composition).

Uncertainty and controls:
- 1000-sample bootstrap confidence intervals for primary metrics.
- 400 source-shuffled null draws to estimate empirical p-values and null z-scores.

Additional robustness:
- Dixit 13-day vs 7-day overlap consistency for matched ablation pairs.

## Results
### 1. Counterfactual consistency is modest and uncertain across datasets
Ablation-level consistency was weak-to-modest across datasets:
- Adamson: Spearman(|effect|, |delta|) = -0.069 (95% CI: -0.416, 0.309).
- Dixit 13-day: 0.051 (95% CI: -0.289, 0.379).
- Dixit 7-day: 0.230 (95% CI: -0.112, 0.520).
- Shifrut: 0.064 (95% CI: -0.255, 0.359).

No ablation dataset beat source-shuffled null at conventional significance (empirical p-values: 0.217 to 0.541).

Interpretation: under current pair coverage and intervention settings, we do not observe strong evidence that scGPT intervention magnitudes consistently track perturbation magnitudes.

### 2. Confound adjustment improves consistency in most datasets, but effect size remains limited
Confound-aware residualization increased absolute consistency in 3/4 ablation rows:
- Adamson: +0.025,
- Dixit 13-day: +0.230,
- Dixit 7-day: -0.042,
- Shifrut: +0.030.

Mean ablation consistency increased from 0.069 (raw) to 0.129 (controlled).

Interpretation: baseline expression and composition shift explain part of the raw discrepancy, but adjustment alone is insufficient to produce robust, high-confidence alignment.

### 3. Baseline perturbation ranking metrics are not sufficient evidence of mechanistic faithfulness
Ablation perturbation benchmark AUPR values were moderate to high:
- Adamson: 0.576,
- Dixit 13-day: 0.839,
- Dixit 7-day: 0.772,
- Shifrut: 0.794.

However, these benchmark scores coexisted with weak pair-level counterfactual consistency and broad uncertainty.

Interpretation: aggregate benchmark separability can overstate mechanistic faithfulness when pair-level directional alignment is weak.

### 4. Timepoint robustness in Dixit is asymmetric between magnitude and direction
For overlapping Dixit ablation pairs across 13-day and 7-day data:
- Spearman on |delta|: 0.387,
- Spearman on signed delta: 0.075,
- Sign agreement: 0.25.

Interpretation: relative perturbation magnitude is somewhat preserved across timepoints, but directionality is unstable, weakening causal sign-level claims.

## Biological Interpretation
The strongest cross-timepoint concordant pair under strict dual agreement was IRF1 -> BLVRB. This is hypothesis-generating rather than confirmatory evidence: it suggests a possible stress/immune-response axis where intervention effects and perturbation deltas agree at both timepoints, but the support size is too small for broad mechanistic claims.

More generally, concordant pairs were sparse, especially in Shifrut where many model effects were near zero. This pattern is more consistent with limited intervention sensitivity (or pair-selection mismatch) than with broad biological circuit recovery.

## Discussion
This study provides a conservative counterfactual validation layer for mechanistic interpretability claims in single-cell transformers. The central negative result is that current low-compute intervention runs do not provide strong pair-level counterfactual alignment against perturbation measurements.

What would change our conclusion:
1. Larger pair coverage with tighter confidence intervals.
2. Stable sign-level agreement across timepoints and datasets.
3. Null-rejecting consistency after stronger confound controls.

Limitations:
- Low-compute pair subsets likely under-sample relevant regulatory structure.
- Composition adjustment is source-level rather than full cell-level modeling.
- Source-shuffled null is informative but not a full generative null.

## Conclusion
Counterfactual consistency between scGPT interventions and perturbation outcomes is detectable in limited settings but currently too weak and uncertain to support strong submission-level mechanistic claims. The framework is useful, but evidence must be strengthened through higher-coverage interventions, stronger confound modeling, and more robust cross-condition validation.

## Data and Code Availability
All code, analysis scripts, generated tables, and figure-generation artifacts are publicly available at:

https://github.com/Biodyn-AI/counterfactual-perturbation-consistency-scgpt

The repository README provides instructions for reproducing the main figures and tables and documents required dataset access prerequisites.

## References
1. Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*. 2024;21:1470-1480. doi:10.1038/s41592-024-02201-0.
2. Tabula Sapiens Consortium, Jones RC, Karkanias J, et al. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science*. 2022;376(6594):eabl4896. doi:10.1126/science.abl4896.
3. Adamson B, Norman TM, Jost M, et al. A multiplexed single-cell CRISPR screening platform enables systematic dissection of the unfolded protein response. *Cell*. 2016;167(7):1867-1882.e21. doi:10.1016/j.cell.2016.11.048.
4. Dixit A, Parnas O, Li B, et al. Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. *Cell*. 2016;167(7):1853-1866.e17. doi:10.1016/j.cell.2016.11.038.
5. Shifrut E, Carnevale J, Tobin V, et al. Genome-wide CRISPR screens in primary human T cells reveal key regulators of immune function. *Cell*. 2018;175(7):1958-1971.e15. doi:10.1016/j.cell.2018.10.024.
6. Geva M, Schuster R, Berant J, Levy O. Transformer feed-forward layers are key-value memories. In: *Proceedings of EMNLP 2021*. doi:10.18653/v1/2021.emnlp-main.446.
7. Meng K, Sharma A, Andonian A, et al. Locating and editing factual associations in GPT. In: *Advances in Neural Information Processing Systems*. 2022. arXiv:2202.05262.
8. Sundararajan M, Taly A, Yan Q. Axiomatic attribution for deep networks. In: *Proceedings of ICML 2017*.
9. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*. 2019;15(6):e8746. doi:10.15252/msb.20188746.
10. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*. 2019;16:1289-1296. doi:10.1038/s41592-019-0619-0.
11. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature*. 2023;618:616-624. doi:10.1038/s41586-023-06139-9.
12. Rives A, Meier J, Sercu T, et al. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *Proceedings of the National Academy of Sciences*. 2021;118(15):e2016239118. doi:10.1073/pnas.2016239118.

## D) Figure/Table Audit
1. **Figure: consistency_vs_null.png**
   - Issues found: single-view summary was previously missing uncertainty context.
   - Fixes applied: uncertainty now provided in companion bootstrap-CI metrics table; figure retained as overview and interpreted conservatively in text.
   - Remaining improvement for submission: add panelized figure with CI whiskers and per-dataset null distributions.

2. **Table: counterfactual_consistency_metrics.tsv**
   - Issues found: initial version lacked confidence intervals and dataset-robustness fields.
   - Fixes applied: added bootstrap CIs for Spearman/sign agreement and retained null statistics.

3. **Table: timepoint_consistency.tsv**
   - Issues found: absent in initial draft.
   - Fixes applied: added overlap-based timepoint robustness metrics for Dixit 13d vs 7d.

4. **Table: perturbation_baseline_metrics.tsv**
   - Issues found: could be misread as mechanistic faithfulness evidence by itself.
   - Fixes applied: interpretation now explicitly separates baseline ranking quality from counterfactual pair-level consistency.

5. **Tables: matched_pairs_*.tsv and high_confidence_consistency_pairs.tsv**
   - Issues found: missing pair-level transparency previously.
   - Fixes applied: exported pair-level matched data and concordant-pair shortlist for direct audit and downstream biological validation.

## E) Claims Table
| Main claim | Support location | Strength |
|---|---|---|
| scGPT intervention-to-perturbation consistency is currently weak and uncertain. | Revised paper: Results §1; Table `counterfactual_consistency_metrics.tsv`; null table | **Strong** |
| Confound adjustment improves consistency in most datasets but not enough for strong claims. | Revised paper: Results §2; adjusted-vs-raw metrics columns | **Medium** |
| Baseline AUPR/AUROC alone can overstate mechanistic validity. | Revised paper: Results §3; baseline vs pair-level contrast | **Medium** |
| Dixit timepoint robustness is magnitude-stable but direction-unstable. | Revised paper: Results §4; `timepoint_consistency.tsv` | **Medium** |
| IRF1 -> BLVRB is a plausible hypothesis-level concordant pair, not confirmation. | Revised paper: Biological Interpretation; `high_confidence_consistency_pairs.tsv` | **Weak** |

## F) References (Complete List)
1. Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*. 2024;21:1470-1480. doi:10.1038/s41592-024-02201-0.
2. Tabula Sapiens Consortium, Jones RC, Karkanias J, et al. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science*. 2022;376(6594):eabl4896. doi:10.1126/science.abl4896.
3. Adamson B, Norman TM, Jost M, et al. A multiplexed single-cell CRISPR screening platform enables systematic dissection of the unfolded protein response. *Cell*. 2016;167(7):1867-1882.e21. doi:10.1016/j.cell.2016.11.048.
4. Dixit A, Parnas O, Li B, et al. Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. *Cell*. 2016;167(7):1853-1866.e17. doi:10.1016/j.cell.2016.11.038.
5. Shifrut E, Carnevale J, Tobin V, et al. Genome-wide CRISPR screens in primary human T cells reveal key regulators of immune function. *Cell*. 2018;175(7):1958-1971.e15. doi:10.1016/j.cell.2018.10.024.
6. Geva M, Schuster R, Berant J, Levy O. Transformer feed-forward layers are key-value memories. In: *Proceedings of EMNLP 2021*. doi:10.18653/v1/2021.emnlp-main.446.
7. Meng K, Sharma A, Andonian A, et al. Locating and editing factual associations in GPT. In: *Advances in Neural Information Processing Systems*. 2022. arXiv:2202.05262.
8. Sundararajan M, Taly A, Yan Q. Axiomatic attribution for deep networks. In: *Proceedings of ICML 2017*.
9. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*. 2019;15(6):e8746. doi:10.15252/msb.20188746.
10. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*. 2019;16:1289-1296. doi:10.1038/s41592-019-0619-0.
11. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature*. 2023;618:616-624. doi:10.1038/s41586-023-06139-9.
12. Rives A, Meier J, Sercu T, et al. Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *Proceedings of the National Academy of Sciences*. 2021;118(15):e2016239118. doi:10.1073/pnas.2016239118.
