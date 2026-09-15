# Chemistry-conditioned reliability of machine-learned polymer glass transition temperatures

A benchmark that asks a question point-accuracy tables cannot answer: **when a
Tg model meets a repeat-unit chemistry it was never trained on, does it know?**

Published Tg models routinely report R² above 0.85 on random splits of the same
public dataset. That number describes interpolation inside well-populated
chemistry. A screening campaign, by contrast, spends its time on repeat units
that are unlike anything in the training set — which is exactly where a random
split never looks. This project measures the size of that gap, separates its
causes, and tests whether uncertainty estimates remain trustworthy across it.

## What is here

| Stage | Question |
|---|---|
| Curation audit | Does the public record survive an independent re-derivation? |
| Split regimes | How much of the reported accuracy is interpolation? |
| Matched-pair design | Is the family-holdout penalty caused by *which* chemistry is missing, or just by *how much* data is missing? |
| Conformal comparison | Does an interval with a 90% guarantee actually cover 90% of the unfamiliar polymers? |
| Applicability domain | Is there a structural-similarity threshold beyond which predictions should not be trusted? |

## The three claims

1. **Chemistry-aware splits cost roughly twice the error of random splits**, and
   the penalty is family-specific rather than uniform.
2. **The penalty is causal, not a population artefact.** A matched design holds
   the test set and the training-set *size* fixed and varies only which family is
   removed, which no previous treatment of this dataset does. The confound is
   acknowledged in prior work and left open.
3. **Marginal conformal coverage hides the failure; conditional coverage exposes
   it — and a novelty-conditioned taxonomy repairs it.** Split conformal holds
   its nominal 90% overall while covering far less in the least-familiar
   structural band. Conditioning on polymer family cannot fix this (an unseen
   family has no calibration data, by construction). Conditioning on
   *nearest-training similarity* can, because that coordinate is defined for any
   repeat unit.

## Data

7,208 repeat-unit records (SMILES + reported Tg), the PoLyInfo-derived collection
redistributed with the POINT² benchmark, aggregated here to **7,174 unique
canonical structures**. Our independent curation reproduces a prior public audit
of this dataset exactly — 7,208 → 7,174 structures, 31 repeated canonical groups,
34 absorbed records, 28 groups with disagreeing Tg, 115 K maximum within-structure
range, 7,170/3/1 structures with 2/3/4 attachment points.

Polymer families are **derived, not inherited**: a backbone-aware SMARTS
classifier locates the chain path between attachment points and classifies the
linkages on it, so a poly(alkyl acrylate) is a polyacrylic rather than a
polyester. See `src/ptgbench/families.py`.

## Reproducing

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt

python scripts/prepare_features.py                 # cache descriptors + fingerprints
./scripts/run_all.sh random scaffold cluster family
N_REPEATS=3 ./scripts/run_all.sh matched
python scripts/analyze.py                          # tables + figures
python -m pytest tests -q
```

A cold run is dominated by descriptor computation (~2 min) and the matched-pair
stage (~1 h on 8 cores); everything else is minutes.

## Layout

```
src/ptgbench/
  data.py       curation, canonicalisation, structure-level aggregation
  families.py   backbone-aware polymer family assignment
  features.py   RDKit descriptors, Morgan fingerprints, leak-free column filter
  splits.py     random / scaffold / cluster / family / matched-pair partitions
  models.py     regressor pipelines (filtering fitted in-fold)
  conformal.py  split, normalised and Mondrian conformal predictors
  metrics.py    point, interval, subgroup-coverage and sparsification metrics
  pipeline.py   one split, end to end
  figures.py    validated colour-vision-safe figure style
scripts/        prepare_features.py, run_benchmark.py, run_all.sh, analyze.py
paper/          manuscript and generated figures
```

## Methodological choices worth knowing

- **Descriptor filtering is fitted inside each fold.** Dropping constant or
  degenerate columns using the full table — common in this literature, and an
  acknowledged limitation of the closest prior analysis — lets test structures
  influence the feature definition. Here it is a pipeline step.
- **Training pools are split three ways**, not two: the normalised conformal
  predictor's difficulty model is fitted on residuals the point model did not
  train on, so it cannot learn that the model is everywhere accurate.
- **Every regime is repeated.** Single-split model rankings on this dataset sit
  well inside run-to-run noise.
- **R² is never used alone to rank family difficulty.** It is normalised by the
  test set's own variance, so a narrow-Tg family can show a catastrophic R² at a
  modest absolute error.
