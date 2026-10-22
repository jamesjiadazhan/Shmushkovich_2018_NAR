# sdrnapred

An independent R reimplementation of the linear regression-based positional
base-preference algorithm for predicting the efficacy of cholesterol-conjugated,
chemically modified, self-delivering siRNAs (sdRNAs), from:

> Shmushkovich T, Monopoli KR, Homsy D, Leyfer D, Betancur-Boissel M,
> Khvorova A, Wolfson AD (2018). Functional features defining the efficacy of
> cholesterol-conjugated, self-deliverable, chemically modified siRNAs.
> *Nucleic Acids Research* 46(20):10905–10916.
> [doi:10.1093/nar/gky745](https://doi.org/10.1093/nar/gky745)

**This is not the authors' code.** No software accompanied the paper. The
package was written from the published Methods and validated against the
performance figures the paper reports. Read the fidelity section below before
relying on it.

## Install

```r
# install.packages("devtools")
devtools::install_local("sdrnapred")
```

## Use

```r
library(sdrnapred)

# Build the weight matrix the paper selected: functional < 24% expression
# remaining, non-functional > 44%, positions 2-49.
wm <- sdrna_weight_matrix(functional_cutoff = 24, seed = 1)
wm
#> sdRNA positional base-preference weight matrix
#>   functional < 24% expression (n = 91)
#>   non-functional > 44% expression (n = 157)

# Score the training panel, or any 50-nt target regions of your own.
scores <- sdrna_score(sdrna356$sequence, wm)

# Higher score = predicted more potent = lower expected expression.
cor(scores, sdrna356$expression)

# Performance curves, as in Figure 2B.
curve <- sdrna_ppp_sensitivity(scores, sdrna356$expression, 24)
sdrna_ppp_at(curve, sensitivity = 40)

# Honest out-of-sample estimate.
sdrna_cv(functional_cutoff = 24, k = 10, seed = 1)$auc
```

Scoring your own candidates requires the **50-nucleotide mRNA target region**
in the paper's layout: the 20-nucleotide siRNA targeting site with 15
nucleotides of flanking sequence on each side (characters 16–35 are the
targeting site). The weights are specific to the paper's heavily modified,
asymmetric, cholesterol-conjugated chemical scaffold — the paper's central
result is precisely that algorithms do **not** transfer between chemical
scaffolds, so do not apply this matrix to unmodified siRNA.

## Data

`sdrna356` is the complete 356-compound training panel from Supplementary
Table 1 (IDs 1000–1355), redistributed under the paper's CC BY 4.0 licence. It
was extracted from the publisher's PDF and verified cell by cell against an
independent parser (1,424/1,424 cells exact). Two source quirks are preserved
rather than normalised: expression values `25.405` (ID 1095) and `29.2`
(ID 1107) depart from the table's two-decimal convention, and ID 1284 is
printed with no bold in the source, so its targeting site is unmarked there.

The bundled subset sizes reproduce the paper's stated bins exactly — 41, 91 and
138 functional at the <17%, <24% and <35% cutoffs, 157 non-functional at
>44% — which is the first check that the extraction is faithful.

## How closely does it replicate?

All figures below use `seed = 1`. The paper's own numbers are resubstitution
estimates (matrix built on all 356 compounds, then scored on those same
compounds), so the comparable rows are in-sample.

| Benchmark | Paper | This package |
|---|---|---|
| Pearson \|r\|, score vs expression (24%) | 0.55 | **0.565** |
| PPP at 40% sensitivity (24% cutoff) | ~85% | **84.1%** |
| AUC, functional <35% | 0.84 | **0.848** |
| AUC, functional <24% | 0.90 | 0.875 |
| AUC, functional <17% | 0.81 | 0.943 |
| Control algorithm AUC | 0.69 | 0.760 |

The two headline operating points the paper uses to justify its choice of
cutoff — r = 0.55 and ~85% accuracy at 40% sensitivity — reproduce closely, as
does the <35% AUC. Two rows do not:

- **The <17% AUC overshoots badly (0.943 vs 0.81).** That cutoff leaves only 41
  functional compounds to fit 125 non-zero weights, so this is where
  resubstitution optimism is largest and where any difference in how strictly
  cells are admitted shows up most. The likeliest cause is that this
  implementation's significance filter admits more cells than the authors'.
- **The control AUC is too high (0.760 vs 0.69).** The paper's control replaces
  the training set's per-position frequencies with "average RefSeq database
  frequencies adjusted for initial bias in the training set". Those adjusted
  frequencies are not published, so `sdrna_control_matrix()` defaults to a
  uniform 25% composition, which is a weaker control than theirs. Supply your
  own via the `refseq` argument to get closer.

### Out-of-sample performance (not reported in the paper)

10-fold cross-validation, rebuilding the matrix from scratch on each training
fold:

| Functional cutoff | In-sample AUC | Cross-validated AUC |
|---|---|---|
| <17% | 0.943 | **0.745** |
| <24% | 0.875 | **0.701** |
| <35% | 0.848 | **0.661** |

The algorithm retains real but much more modest predictive signal out of
sample. Anyone using this to prioritise candidates should plan against the
cross-validated numbers, not the published in-sample ones. Note also that folds
are drawn over compounds, and the panel targets only 17 genes, so compounds
against the same transcript share sequence context; a gene-grouped split would
be stricter still, but Supplementary Table 1 does not publish which gene each
compound targets.

## Where the published description is ambiguous

The Methods do not fully determine an implementation. The choices made here:

1. **`Md` in Equation 2** is read as the per-position base *count* in the
   functional or non-functional subset, so that `Md_F/q` and `Md_NF/z` are
   frequencies. This is the reading that makes the equation dimensionally
   coherent.
2. **σ in Equation 2** is the standard deviation, across simulated random sets,
   of the per-position base counts at the relevant subset size. The paper
   generates `N = 356` such sets per size; `n_replicates` controls this.
3. **Degrees of freedom** for the one-sample t-test are taken as
   `n_replicates - 1`.
4. **Bonferroni correction** is applied across all tested cells
   (4 bases × retained positions).
5. **The sign** of Equation 2's numerator is opposite to that of the weight
   (`freq_F - freq_NF`); since the test is two-sided this affects only the
   sign of the statistic, not the p-value.
6. **The control's RefSeq frequencies** are not published; see above.

Items 1–3 are the ones most likely to explain the residual gaps. All are
parameterised rather than hard-coded, so alternative readings can be tested
without editing the package.

## What this package does not do

- It does not reproduce the paper's AU-preference or thermodynamic-flexibility
  analyses (Figure 5), nor the position-14 2′-O-methyl experiments (Figure 6),
  which are wet-lab results.
- It does not include the two independent validation datasets (50 sequences by
  qPCR; 94 by QuantiGene), which the paper does not publish as sequences.
- It does not model chemical modification patterns. The scaffold is assumed
  fixed and uniform, as it is throughout the training panel.

## Testing

```r
devtools::test()
```

49 tests cover the extraction invariants, the weight-matrix construction
(including an independent recomputation of the weights), scoring, the
evaluation metrics against brute-force references, the control substitution,
and cross-validation.

## Licence

Code: MIT. Bundled data: CC BY 4.0, © the authors of
[doi:10.1093/nar/gky745](https://doi.org/10.1093/nar/gky745); please cite the
paper if you use it.
