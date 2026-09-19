---
title: "Backbone Chemistry Governs the Reliability of Machine-Learned Glass Transition Temperatures: A Family-Resolved Benchmark with Novelty-Conditioned Conformal Uncertainty"
target_journal: "Polymers (MDPI)"
article_type: "Article"
keywords:
  - glass transition temperature
  - backbone chemistry
  - structure-property relationships
  - group-contribution methods
  - repeat-unit chemistry
  - polymer informatics
  - uncertainty quantification
  - conformal prediction
  - applicability domain
  - extrapolation
---

## Abstract

Which polymer backbones can a machine-learned glass transition temperature
(*T*~g~) model be trusted on? Random partitions, reporting *R*² above 0.85,
cannot answer: their test sets hold chemistry seen in training. On 7,174
repeat units in twenty backbone families plus an abstention class, histogram
gradient boosting gives 28.1 K mean absolute error (random), 33.8 K (scaffold),
44.5 K (fingerprint-cluster), 48.1 K (whole-family holdout), pooled *R*² falling
0.87 to 0.75; the two best models differ by 0.38 K. The failure is chemical:
withheld polyphosphazenes and polysiloxanes are overpredicted (bias −83.7,
−63.7 K), polydienes and polyolefins underpredicted (+48.0, +37.8 K);
deterioration tracks structural isolation, not family size. A matched-pair
design with fixed test structures attributes +14.97 K (95% CI [+11.21, +19.28]) to
removing the family, worse in all twenty-one families for both models, and under
1 K of unresolved sign to removing as many unrelated structures. Uncertainty
fails there: for extremely randomised trees, split conformal holds 0.904
marginal coverage (mean of splits) on random splits while covering 0.657 of the
least familiar band (pooled), and loses marginal validity under shift.
Family-conditioned calibration cannot repair this, a withheld family having no
calibration members; similarity-conditioned calibration holds 0.877–0.907 in all
four regimes (mean of splits).

## 1. Introduction

The glass transition temperature is the single most consulted thermal descriptor of
an amorphous polymer. It bounds the service window of a thermoplastic, sets
processing conditions, and governs the mechanical response of a matrix in a
composite. What fixes its value is the chain itself: *T*~g~ rises with the barrier
to segmental rotation about backbone bonds and with the strength of the
interactions between neighbouring chains, so a rigid, strongly interacting
heteroaromatic backbone vitrifies far above a flexible Si–O or P=N one. The
families at the two ends of the dataset used here show that spread directly:
aromatic polyimides have a median *T*~g~ of 249 °C, polysiloxanes and
polyphosphazenes of −8 °C (Table 2).

Measuring *T*~g~ by differential scanning calorimetry or dynamic mechanical
analysis is routine but slow, and it requires a synthesised, purified sample, so
predicting it from the structure of the repeat unit has been a target of
polymer theory since group-contribution methods were first formalised [1–3].
Those schemes write *T*~g~ as a composition-weighted sum of tabulated group
contributions over the repeat-unit
molar mass. The assumption is compositional additivity: no interaction terms,
and no dependence on which linkage a group sits in. That assumption is what lets
such a scheme be applied to a repeat unit its table was never fitted on, and it
is also what caps its accuracy. Refitted on the training data used here, an
additive scheme of this form reaches 43.8 K mean absolute error on random
partitions, where a descriptor ensemble reaches 28.1 K (Section 3.9).

Data-driven surrogates have largely displaced those schemes. Curated repositories
such as PoLyInfo [4] made supervised learning on experimental *T*~g~
practical, and a decade of work has explored the resulting design space: physically
interpretable descriptor sets with tree ensembles and kernel methods [5,6],
graph neural networks over the repeat-unit graph [7], and chemical language
models pretrained on millions of hypothetical polymers [8,9]. Recent community
infrastructure has consolidated the datasets and the evaluation code [10,11].
The accuracies reported are genuinely useful:
mean absolute errors in the region of 25–30 K against experimental values whose own
inter-source spread can exceed 100 K.

Those numbers describe interpolation. They are obtained on random partitions, in
which a test repeat unit typically has close structural analogues in training, and
they therefore estimate how well a model fills in gaps within chemistry it has
already seen. The motivating application is the reverse. A screening campaign
enumerates candidates precisely because they are not yet in the literature, and the
candidates that matter most are the ones furthest from the training distribution.
The mismatch between the validation regime and the deployment regime is well known
in cheminformatics, where scaffold and time-based splits were introduced for
exactly this reason [12,13], and it has begun to be examined for polymers
[14,15].

Three questions remain open, and this work addresses all three.

**Which backbone chemistries does a model fail on, and in which direction?** An
aggregate error over a heterogeneous collection is not actionable for a chemist
working on one class of polymer. Errors resolved by family, and reported with
their sign, say what a pooled MAE cannot: whether a backbone is placed
systematically too high or too low, and whether that displacement follows the
chemistry or merely follows the distance of the family's *T*~g~ from the rest of
the record. The two readings have to be separated before the error structure can
be called chemical.

**Is the family-holdout penalty caused by the missing chemistry?** Withholding an
entire polymer family and observing an error increase is the standard diagnostic,
but it changes two things at once. The test population changes — families differ in
how intrinsically hard they are and in how wide their *T*~g~ distribution is — and
the training set shrinks. A comparison against a random-split baseline therefore
measures the sum of a population effect, a data-volume effect and the effect of
interest. The closest prior analysis of this dataset states the confound
explicitly and leaves it unresolved [16]. Resolving it requires holding the test
set fixed and matching the training pools on size.

**Can a *T*~g~ prediction for a backbone the model has not seen be trusted?** A
point prediction without a trustworthy interval is not actionable for screening,
and split conformal
prediction [17–19] is an attractive answer because it
wraps any regressor in intervals with a finite-sample coverage guarantee under
exchangeability. Two difficulties arise here. Exchangeability is precisely what a
family holdout violates. And the guarantee is *marginal*: a predictor can achieve
90% coverage overall while systematically under-covering a subpopulation, and in
this problem that subpopulation is the unfamiliar chemistry the model was deployed
to explore. Ensemble uncertainty paired with an explicit applicability-domain
check has been applied to polymer property prediction [20], but without
a family-resolved conditional-coverage test. A recent analysis on a 410-sample
simulation-derived set reports this
pattern and characterises conformal intervals as conservative triage indicators
rather than fine-grained screening tools [15]; whether the failure is
repairable, and whether it persists at the scale and heterogeneity of the
experimental record, has not been tested.

### 1.1. Contributions

1. **An independently re-derived, chemically resolved dataset.** We re-curate
   7,208 experimental repeat-unit records to 7,174 canonical structures,
   reproducing a prior public analysis of this collection [16] step for step, and assign
   each structure to one of twenty backbone families or to an explicit
   abstention class using a *backbone-aware* substructure classifier.
   Classifying on the chain path rather than on whole-molecule matches keeps
   poly(alkyl acrylate)s out of the polyester class, which matters because every
   family-holdout conclusion inherits the labels.
2. **A family-resolved account of where prediction fails, and in which
   direction.** Holdout deterioration spans a factor of nearly three across the
   twenty-one groups, and the signed error is systematic: withheld
   polyphosphazenes and polysiloxanes are overpredicted, by 83.7 K and 63.7 K on
   average, and withheld polydienes and polyolefins underpredicted, by 48.0 K and
   37.8 K. Regression toward the training mean accounts for part of the first
   pair and none of the second, so the chemical reading rests on the
   underpredicted hydrocarbon families rather than on the inorganic ones
   (Section 3.3). Deterioration tracks how structurally isolated a family is
   (Pearson *r* between −0.62 and −0.69 across four models) while a family-size
   effect is not resolved at twenty-one families.
3. **A matched-pair design that isolates the family effect.** Three training
   pools — one retaining the family, one with it removed, and one size-matched
   control from which an equal number of unrelated structures was removed —
   predict an identical held-out test set, separating the cost of missing
   chemistry from the cost of missing data.
4. **A test of whether a *T*~g~ prediction for an unseen backbone can be
   trusted**, by conditional-coverage evaluation of conformal uncertainty, with a
   novelty-conditioned Mondrian taxonomy that substantially narrows the subgroup
   coverage gap where a family-conditioned taxonomy structurally cannot.
5. **An applicability-domain diagnostic in units a practitioner can act on.**
   Using conformal prediction itself as the applicability-domain statement, in
   the tradition established by Norinder et al. [21], coverage, interval width
   and absolute error are resolved by nearest-neighbour Tanimoto similarity, and
   the similarity itself is reported alongside each prediction so that a reader
   can see which regime a prediction sits in. We do
   not offer a single similarity below which an interval should be disbelieved,
   because under family holdout the error is not monotone in similarity: for
   extremely randomised trees it falls from 55.5 K below Tanimoto 0.4 to 30.3 K
   in the 0.6–0.7 band and then rises again to 36.7 K above 0.8 (Section 3.5), so
   no cut-point separates trustworthy from untrustworthy predictions.
6. **A Van Krevelen group-contribution reference point and a taxonomy
   sensitivity check**, so that the machine-learned result is placed against the
   additive schemes it proposes to replace and against a perturbation of its own
   labels.
7. **A reproducible open implementation**, with leak-free in-fold preprocessing
   and a test suite pinning both the chemistry of the family assignment and the
   finite-sample behaviour of the conformal predictors.

## 2. Materials and Methods

### 2.1. Dataset and curation

We use the experimental homopolymer *T*~g~ collection redistributed with the
POINT² polymer informatics benchmark [11], which derives from the PoLyInfo
record [4]. Each entry pairs a repeat unit written as a PSMILES string —
SMILES with two ``*`` atoms marking the polymerisation attachment points — with a
reported *T*~g~.

The raw export is record-level: one repeat unit may appear several times carrying
values from different primary sources. Modelling on records would allow identical
inputs to fall on both sides of a partition and would weight frequently reported
polymers more heavily, so all analysis is performed on a structure-level table.
Repeat units were canonicalised with RDKit [22] and grouped; where a canonical structure
carried several reported values, the median was taken, and the number of
contributing records and the full reported range were retained as provenance.

### 2.2. Backbone-aware family assignment

Polymer family labels were derived rather than inherited, so that the taxonomy is
reproducible and its chemistry is inspectable. For each repeat unit the *backbone*
is taken as the union of shortest paths between attachment points, extended to
include any ring the path enters. Characteristic linkages are then matched against
SMARTS patterns and accepted only when their core atoms lie on that backbone.

Because a repeat unit may be written cut at its own characteristic linkage —
nylon-6 as ``*NCCCCCC(=O)*`` divides its amide between the two ends — the
classifier operates on a head-to-tail dimer rather than on the single unit.
Joining two copies restores the linkage the cut destroyed while keeping the
junction a bond that exists in the real chain.

Backbone restriction is what separates a poly(alkyl acrylate) from a polyester:
both contain an ester group, but in the acrylate it decorates an all-carbon chain.
Patterns are tested in order of specificity, because the groups are nested — every
cyclic imide contains two amide-like motifs, and a urethane contains both
ester-like and amide-like fragments. Repeat units whose backbone carries no
heteroatom linkage are classified by the groups pendant to the chain
(aromatic-dominated, acrylic, styrenic, halogenated, other heteroatom-substituted
vinyl, unsaturated, or saturated hydrocarbon). A unit whose backbone does carry a
heteroatom but matches no linkage pattern is assigned to an explicit abstention
class, **Other backbone**, rather than forced into the nearest hydrocarbon
category. The twenty-five chemistry cases in Table S1 are pinned by unit tests.

### 2.3. Representations

Two representations were computed. The full RDKit descriptor block (217 columns)
provides physically interpretable features; the topological ``Ipc`` index was
replaced by its base-10 logarithm because it overflows for larger repeat units.
Binary Morgan fingerprints (radius 2, 2,048 bits) [23] provide the
similarity space in which structural novelty is measured. Attachment points were
retained in both, since they identify where the chain continues and the environments around them
encode exactly the linkage motifs that set *T*~g~.

Descriptor column filtering — removal of non-finite, constant and exactly
duplicated columns — is implemented as the first step of a scikit-learn pipeline
[24] and is therefore **fitted on the training fold only**. Performing this filtering
once over the full table, as is common, lets test structures influence the feature
definition; the closest prior analysis of this dataset notes this as an unquantified
limitation of its own results [16].

### 2.4. Partitioning regimes

Four regimes were used, in increasing order of enforced novelty: **random**
partitions (repeated, 80:20); **scaffold** partitions, in which generic
Bemis–Murcko scaffolds [25] are assigned whole; **cluster** partitions,
using average-linkage agglomerative clusters in Jaccard space over Morgan bits;
and **family** holdouts, in which every member of one backbone family is withheld.
All non-family regimes were repeated ten times with different seeds, because
single-split rankings among tree ensembles on this dataset fall inside run-to-run
noise. The family regime contributes twenty-one splits, one per family.

### 2.5. Matched-pair design

A family is eligible if it holds at least forty structures, which all twenty-one
groups in this dataset do. For each eligible family *F*, a test set *T* ⊂ *F* is
drawn (30% of *F*) and three training pools are constructed: **informed**
(everything but *T*, so the remainder of *F* is available), **naive** (everything
but *F* entirely), and **control** (**informed** minus a randomly chosen block of
structures from other families, equal in number to *F* \ *T*). All three predict
the same *T*. The draw of *T* and of the control block is repeated three times
per family with different seeds, giving sixty-three matched triples per model.

The contrast **naive** − **control** is the quantity of interest. Both pools are
identical in size and both exclude *T*; they differ only in *which* structures
were removed. The contrast **control** − **informed** measures the pure
data-volume effect over the same test set.

Two properties of the design bound what it can establish. First, the removed
block is small: it is 70% of one family, a median of 1.85% of the informed pool
across the twenty-one families, ranging from 0.43% for polyureas to 17.94% for
polyimides, with sixteen of twenty-one families under 5%. Second, the two arms
remove the same *quantity* but not the same *kind* of structure: the naive arm
removes a coherent chemical neighbourhood, the control arm a diffuse random
sample spread over the remaining families. The contrast therefore measures the
effect of removal *locality* at fixed volume. That chemistry is what locality
means here is an interpretation, argued in Section 4.2, not a separate
measurement.

### 2.6. Models

Five regressors were evaluated on the descriptor representation: a training-median
baseline, support vector regression, random forest, histogram gradient boosting and
extremely randomised trees. Every model is wrapped in the same pipeline so that all
data-dependent preprocessing is fitted in-fold.

A sixth model, a group-contribution baseline, was run through the same four
regimes. It implements the classical Van Krevelen additive form
*T*~g~ = *Y*~g~/*M* with *Y*~g~ = Σ*n*~i~*Y*~i~, where *M* is the repeat-unit
molar mass, *n*~i~ the count of structural group *i* and *Y*~i~ its molar
glass-transition contribution [1–3]. Groups
are RDKit fragment counts supplemented by element, ring, rotatable-bond,
hydrogen-bonding and backbone-length counts, which carry the chain skeleton that
fragment counts alone would miss. The contributions *Y*~i~ are fitted by ridge
regression on the training fold rather than taken from a published table, so the
additive scheme receives the same data the ensembles receive and any difference
between them reflects the functional form rather than the calibration set.
Fitting is performed in *Y*~g~ = *T*~g~·*M* space and in absolute temperature,
since an additive law expressed in Celsius would make the contributions depend on
an arbitrary zero.

### 2.7. Conformal prediction and coverage diagnostics

The screening question this machinery is built to answer is a narrow one: given a
predicted *T*~g~ for a repeat unit whose backbone chemistry the model has not
seen, how wide an interval has to be quoted before the number can be acted on?
The construction below is standard; what it is asked to do here is chemical.

Each training pool is divided three ways: a **fit** partition (60%) trains the
point regressor, a **difficulty** partition (20%) supplies out-of-sample
residuals, and a **calibration** partition (20%) supplies the conformal quantile.
Fitting the difficulty model on the point model's own training residuals would
teach it that the model is everywhere accurate, so that partition must be held
out. All conformal variants share one point predictor and one calibration set, so
differences between them reflect only how residuals are converted into intervals.

The three-way division costs accuracy. On a random 80:20 partition the point
regressor is fitted on about 48% of the dataset, three-fifths of what an
80%-trained literature model sees, so the random-split errors reported below are
mildly pessimistic relative to published figures obtained without a calibration
carve-out. Data-reusing constructions such as the jackknife+ avoid the carve-out
at the cost of refitting [26]; the split construction is kept here so
that all four interval methods share one fitted model and one calibration set and
differ in nothing else.

Given calibration scores *s*~1~…*s*~n~ and a miscoverage level α, the conformal
quantile is the ⌈(*n*+1)(1−α)⌉-th order statistic; when that index exceeds *n* the
honest interval is unbounded. We report α = 0.1 throughout. Three predictors are
compared:

* **SCP** — one global quantile of absolute residuals.
* **Normalised SCP** — residuals scaled by a fitted difficulty estimate, so
  intervals widen where the model expects to struggle rather than everywhere
  [27].
* **Mondrian SCP** — a separate quantile per category [28,29],
  under two taxonomies: *polymer family*, and *nearest-training-similarity band*.

The asymmetry between the two taxonomies is the point. A family-conditioned
taxonomy is undefined for a family with no calibration members, which is exactly
the family-holdout case; a novelty-conditioned taxonomy is defined for any repeat
unit, seen family or not. A third route to validity under shift is to reweight
the calibration scores by the likelihood ratio between the test and training
covariate distributions [30]; Section 4.3 says why that route is not
taken here.

A fourth, an **oracle inflation control**, is used only to test whether the
novelty-conditioned predictor's advantage reduces to its intervals being wider.
It multiplies every split conformal interval about its centre by a single
constant *k*. Two choices of *k* are reported. The **coverage-matched** constant
reproduces the novelty-conditioned predictor's marginal coverage on the same test
set, and is the primary control, because it fixes the coverage budget and leaves
only its allocation across the similarity range to differ. The **width-matched**
constant equalises the mean interval width instead; it answers a different
question, and it hands the inflated predictor more marginal coverage than the
predictor it is compared against. Both are oracles and neither is a deployable
method, because *k* is read off the test set it is then compared against; their
role is to bound what *any* global inflation of split conformal could achieve.

Alongside marginal coverage and mean interval width we report **conditional**
coverage within each similarity band and each family, the **maximum coverage
gap** (the largest shortfall below nominal across subgroups with at least ten
members), and the coverage–width criterion. Ensemble spread is additionally
assessed by sparsification error.

Coverage can be aggregated in four ways that differ here by several points, so
every coverage figure below is labelled with the one used.

* **Pooled** marginal coverage is the fraction of all held-out predictions in a
  regime that fall inside their interval, taking every split together. It is the
  quantity a single large test set would give.
* **Mean-of-splits** marginal coverage is the per-split coverage averaged over
  splits, so a 44-structure family weighs as much as a 1,707-structure one. Under
  family holdout split conformal gives 0.701 this way against 0.743 pooled.
* **Band-resolved** coverage is pooled within a similarity band across splits,
  since a single split may contribute few or no structures to a band.
* **Split-averaged worst-band** coverage takes the lowest band coverage within
  each split, over bands holding at least ten structures, and averages those
  minima. Taking the minimum before the average makes it lower than the
  band-resolved coverage of the least-similar band. The size restriction is
  material: dropping it lowers the split conformal figure of Table 13 from 0.564
  to 0.539, because a band of two or three structures can read 0.0 on binomial
  noise alone.

The convention adopted is mean-of-splits for marginal coverage and band-resolved
for coverage within a similarity band. Two places depart from it and say so: the
split-averaged worst band is used in Figure 4b and in the calibration-construction
comparison of Section 3.8, and the oracle controls of Section 3.4 are quoted
pooled throughout, because they compare two predictors structure by structure on
the same held-out set and a split-level average would discard that pairing.

The coverage–width criterion (CWC) is the mean interval width multiplied by
exp(η(1 − α − coverage)) when coverage falls below nominal and left unpenalised
otherwise, with η = 30. Because the penalty is exponential in the shortfall, CWC
differences between methods are dominated by that exponent and not by width: the
criterion is close to a recoding of the coverage gap already reported, and it is
quoted only for comparability with the conformal literature.

### 2.8. Statistical treatment

Confidence intervals on the matched-pair effects are percentile bootstrap
intervals over the twenty-one family-level means, each of which averages the
three repeats. Those twenty-one values are neither independent nor a sample. The
naive and control pools of a given family share at least 78% of their training
data by construction — the minimum is polyimides, the largest family, and the
figure exceeds 95% in sixteen of the twenty-one families; the families are chemically
nested, so that polycarbonates sit inside the same carbonyl chemistry as
polyesters; and the twenty-one are the entire population of families in this
dataset rather than a draw from a larger one. The intervals should therefore be
read as descriptive dispersion across families, not as inference about a
population of families.

The Wilcoxon signed-rank test on twenty-one paired values has a smallest
attainable two-sided *p* of 2⁻²⁰ = 9.5 × 10⁻⁷, reached whenever all twenty-one
signs agree. Where that value appears below it is the floor of the test, not a
measured quantity, and the informative statement is the sign count.

*p*-values are reported unadjusted for multiplicity. The one place where this
matters is the driver analysis of Section 3.3, which runs sixteen correlations
(four candidate drivers × four models). Under a Bonferroni correction at
α = 0.05/16 = 0.0031, the nearest-training-similarity correlation survives for
all four models and no other driver survives for any.

### 2.9. Reproducibility

All code, the curated structure-level table and every result file are available at
https://github.com/Rugvedhs/spectrasense, in the directory
``polymer-tg-benchmark``. Seeds are fixed, and every figure and every result CSV
is written by a script from the stage outputs, the oracle controls of
Section 3.4 and the shrinkage regression of Section 3.3 included: eleven of the
fourteen table sources by ``scripts/analyze.py``, the curation and composition
tables by ``scripts/prepare_features.py`` and the calibration study by
``scripts/run_calibration_study.py``. A cold run reproduces every number in this
paper. The fourteen tables typeset in the body are transcriptions of the CSV
named in each caption.
``scripts/check_manuscript_tables.py`` re-parses those tables out of the
manuscript source and asserts every cell against its named CSV at the rounding
quoted; it runs as part of the test suite.

## 3. Results

### 3.1. The public record survives re-derivation

Independent re-curation of the 7,208 raw records reproduces a prior public
analysis of this collection [16] step for step, and Table 1 sets the two side by
side so that the correspondence can be read rather than asserted: 7,174 unique
raw SMILES, no RDKit parse failures, 7,174 canonical structures, 31 canonical
groups containing more than one record, 34 records absorbed by aggregation, 28 of
those groups carrying disagreeing *T*~g~ values, and a maximum within-structure
range of 115 K. Attachment-point counts likewise agree once the level is matched:
that analysis counts 7,204 records with two attachment points, three with three
and one with four, and we count 7,170 structures with two, three with three and
one with four, the difference being the 34 duplicate records aggregated away. The
dataset is therefore a stable object to build on, and the numbers below are not
sensitive to curation choices.

**Table 1.** Curation audit of the POINT² experimental *T*~g~ export
[11]: reported by [16] | this work. Counts describe the whole collection,
not a model run; the within-structure range is the largest spread between
independent reports of one canonical repeat unit. The attachment-point rows are
record-level in the prior analysis [16] and structure-level here, which is the
whole of their disagreement. The values in the prior-analysis column are
transcribed from that work and are not an output of this pipeline; our column is
`results/table_curation.csv` and `results/table_attachment_points.csv`.

| Step | Reported by [16] | This work |
|---|---|---|
| Raw records | 7,208 | 7,208 |
| Unique raw SMILES | 7,174 | 7,174 |
| RDKit parse failures | 0 | 0 |
| Canonical structures retained | 7,174 | 7,174 |
| Structures aggregated from >1 record | 31 | 31 |
| Records absorbed by aggregation | 34 | 34 |
| Aggregated groups with disagreeing *T*~g~ | 28 | 28 |
| Max within-structure *T*~g~ range (K) | 115 | 115 |
| Units with 2 / 3 / 4 attachment points | 7,204 / 3 / 1 (records) | 7,170 / 3 / 1 (structures) |

That 115 K figure deserves emphasis before any model error is quoted. It is the
spread between independent literature reports of the *same* repeat unit, and it
places a floor under what structure-based prediction can achieve: a repeat-unit
representation cannot resolve differences in molecular weight, tacticity,
crosslinking, thermal history or measurement protocol, and those differences are
folded into the target.

The derived taxonomy places all twenty-five reference polymers of Table S1 in the
family a polymer chemist would assign, and yields the family composition of
Table 2. The resulting family *T*~g~ ordering
(Figure 1) is likewise the expected one, which is a further check that the labels
are chemically meaningful rather than merely self-consistent: among the named
families, polyphosphazenes (median −8 °C) and polysiloxanes (−8 °C) sit lowest,
reflecting the low rotational barriers of P=N and Si–O backbones, and aromatic
polyimides sit highest (249 °C), reflecting rigid, strongly interacting
heteroaromatic chains. The only group below them is the abstention class
Other backbone (median −19.6 °C, 58 structures), of which 57 contain silicon;
its position is consistent with the silane and carbosilane chemistry it holds.

**Table 2.** Composition of the twenty backbone families and the *Other
backbone* abstention class, over the 7,174 structure-level records; *T*~g~ in °C
as reported, ordered by family size. Means and standard deviations are in
`results/table_families.csv`.

| Family | *n* | Median *T*~g~ (°C) | Min (°C) | Max (°C) |
|---|---|---|---|---|
| Polyimides | 1707 | 249.0 | -9.0 | 490.0 |
| Polyamides | 999 | 154.0 | -62.0 | 395.0 |
| Polyesters | 761 | 75.0 | -78.6 | 400.0 |
| Polyethers | 713 | 167.0 | -116.0 | 437.0 |
| Polyacrylics | 661 | 60.0 | -100.0 | 330.0 |
| Polyphenylenes | 423 | 133.0 | -71.4 | 495.0 |
| Polystyrenes | 250 | 95.0 | -75.0 | 330.0 |
| Polysulfones | 232 | 192.0 | 23.0 | 422.0 |
| Polysiloxanes | 209 | -8.0 | -139.0 | 177.5 |
| Polyvinyls | 197 | 57.0 | -79.5 | 174.6 |
| Polyurethanes | 189 | 62.0 | -60.0 | 225.0 |
| Polycarbonates | 189 | 113.0 | -50.0 | 355.0 |
| Polyolefins | 126 | 60.0 | -88.0 | 350.0 |
| Polyphosphazenes | 118 | -8.0 | -105.0 | 291.0 |
| Polydienes | 77 | 101.0 | -99.6 | 420.0 |
| Polysulfides | 60 | 130.5 | -118.0 | 386.0 |
| Other backbone | 58 | -19.6 | -90.0 | 152.0 |
| Polyhalo-olefins | 56 | 33.2 | -108.0 | 232.5 |
| Polyanhydrides | 53 | 98.0 | 20.0 | 260.0 |
| Polyimines | 52 | 159.0 | 10.0 | 327.0 |
| Polyureas | 44 | 169.0 | -50.0 | 292.0 |

### 3.2. Accuracy is governed by the split, not by the model

On random partitions the best model, histogram gradient boosting, reaches
28.1 ± 1.0 K MAE (*R*² = 0.87), consistent with published values for
descriptor-based ensembles on this dataset once the calibration carve-out of
Section 2.7 is allowed for. That ranking is fragile. Across ten repeated splits,
histogram gradient boosting beats extremely randomised trees by 0.38 K,
consistently in sign (9 of 10 splits, Wilcoxon *p* = 0.049) but negligible in
magnitude; a paired *t*-test does not resolve it (*p* = 0.069) because one repeat
reverses. The magnitude, not the significance, is what matters: 0.38 K against a
20.0 K protocol effect. Random forests and support vector regression are worse by
2.17 K and 2.74 K (*p* < 0.0001 each). Reporting a single split, as is common,
would allow either of the leading models to be declared the winner.

Against that sensitivity, the choice of partition moves the error by far more
(Figure 2, Table 3). Holding the model at histogram gradient boosting, MAE rises
from 28.1 K on random splits to 33.8 K on scaffold splits, 44.5 K on
fingerprint-cluster splits and 48.1 K under whole-family holdout. The
random-to-family contrast is therefore 20.0 K: about fifty times the 0.38 K that
separates the two best models, and about seven times the 2.74 K that spans all
four non-trivial regressors. Only the training-median baseline, at +65.7 K, is
moved more by the choice of model than by the choice of protocol. Whatever a
leaderboard on this dataset is ranking, it is not the property that changes most.

**Table 3.** Test mean absolute error in K by model and split regime. Each cell
is the mean over splits, with the standard deviation across splits after ±.
Random, scaffold and cluster are ten repeated splits; family holdout contributes
twenty-one splits, one per family, so its spread is dispersion between families
rather than run-to-run noise. The last row is the mean nearest-training Tanimoto
similarity of the held-out structures. RMSE and *R*² per cell are in
`results/table_point_by_regime.csv`.

| Model | Random (10) | Scaffold (10) | Cluster (10) | Family (21) |
|---|---|---|---|---|
| Histogram gradient boosting | 28.13 ± 1.01 | 33.81 ± 7.76 | 44.47 ± 10.37 | 48.11 ± 15.39 |
| Extremely randomised trees | 28.51 ± 0.93 | 35.50 ± 8.48 | 48.71 ± 11.27 | 49.56 ± 16.19 |
| Random forest | 30.30 ± 1.00 | 36.95 ± 8.37 | 48.79 ± 11.81 | 50.30 ± 15.66 |
| Support vector regression | 30.88 ± 0.90 | 35.34 ± 6.15 | 53.94 ± 17.06 | 53.06 ± 18.28 |
| Training-median baseline | 93.79 ± 1.65 | 95.71 ± 13.67 | 121.61 ± 24.23 | 95.46 ± 28.61 |
| Mean nearest-training similarity | 0.762 | 0.671 | 0.403 | 0.454 |

Mean nearest-training Tanimoto similarity falls alongside the error, from 0.76 on
random splits to 0.67, 0.40 and 0.45. That ordering is not strictly monotone:
family holdout leaves test structures slightly *less* isolated on average (0.45)
than cluster holdout (0.40), because a withheld family may still have near
analogues in a chemically adjacent one, yet family holdout produces the higher
error. Structural distance alone therefore does not determine difficulty;
*which* chemistry is missing also matters, a point Section 3.3 takes up directly.

*R*² must be quoted with its aggregation stated, because the two natural choices
disagree by more than any effect in this paper (Table 4). Pooling residuals over
all held-out structures in a regime, which is the quantity comparable to a
random-split *R*², gives 0.87 on random splits, 0.81 on scaffold, 0.68 on cluster
and 0.75 under family holdout. Averaging *R*² over the twenty-one individual
family holdouts instead gives 0.26, because each family is scored against its own
*T*~g~ variance and a family with a narrow *T*~g~ spread contributes a large
negative value at a modest absolute error. That 0.26 is a real number but it is
not commensurable with 0.87, and the fall from 0.87 to 0.26 that a naive
comparison suggests is an aggregation artefact. The comparable fall is 0.87 to
0.75.

**Table 4.** The two aggregations of *R*², for histogram gradient boosting.
Pooled *R*² is computed once over all held-out predictions in a regime; the
mean-of-splits value averages the per-split *R*², each scored against its own
split's *T*~g~ variance. The two are not interchangeable, and the family-holdout
row is where they diverge. All five models are in `results/table_pooled_r2.csv`.

| Regime | Predictions pooled | Pooled *R*² | Mean of per-split *R*² |
|---|---|---|---|
| Random (10 splits) | 14,350 | 0.866 | 0.865 |
| Scaffold (10 splits) | 14,863 | 0.808 | 0.792 |
| Cluster (10 splits) | 24,985 | 0.679 | 0.538 |
| Family (21 splits) | 7,174 | 0.746 | 0.256 |

The dispersion is as informative as the mean. Random-split MAE varies by ±1.0 K
across ten repeats, and scaffold and cluster splits by ±7.8 K and ±10.4 K; these
are repeat-to-repeat quantities. The ±15.4 K quoted for family holdout is not:
each family is held out exactly once, so it is the standard deviation *between
families*, and it describes how unevenly the penalty is distributed over
chemistry rather than how reproducible any one holdout is. Either way, a point
estimate from one chemistry-aware split should not be trusted.

### 3.3. Transfer failure is a property of backbone chemistry

Family-holdout deterioration, defined as the ratio of family-holdout MAE to
random-split MAE for the same model, spans a factor of nearly three across the
twenty-one families (Table 5; histogram gradient boosting throughout this
section). At the difficult end sit polyphosphazenes (3.02×, 84.8 K MAE),
polydienes (2.83×, 79.6 K) and polysiloxanes (2.48×, 69.7 K). At the easy end,
polyvinyls deteriorate by 1.07×, polyethers by 1.19× and polycarbonates by
1.20×; no family is predicted better when withheld than the average random-split
structure, though polyvinyls come within 8%.

**Table 5.** Family-holdout error for histogram gradient boosting: the three
most and three least deteriorated of the twenty-one families. Each family is
withheld exactly once, so *n* is the number of held-out structures and there is
no repeat-to-repeat spread. Deterioration is the family's holdout MAE divided by
the same model's mean random-split MAE of 28.13 K. All twenty-one families are
listed in Appendix A and in `results/table_family_holdout.csv`.

| Family | *n* | MAE (K) | Deterioration | Mean NN similarity |
|---|---|---|---|---|
| Polyphosphazenes | 118 | 84.8 | 3.02× | 0.304 |
| Polydienes | 77 | 79.6 | 2.83× | 0.367 |
| Polysiloxanes | 209 | 69.7 | 2.48× | 0.388 |
| … | | | | |
| Polycarbonates | 189 | 33.7 | 1.20× | 0.558 |
| Polyethers | 713 | 33.4 | 1.19× | 0.553 |
| Polyvinyls | 197 | 30.2 | 1.07× | 0.430 |

The signed error carries the mechanism, but it has to be separated from a purely
statistical effect first. Defining bias as mean(*T*~g~ observed − *T*~g~
predicted), so that a negative value is overprediction, the withheld inorganic
backbones are displaced by wide margins: polyphosphazenes by −83.7 K,
polysiloxanes by −63.7 K and the Other backbone silane class by −45.3 K, against
MAEs of 84.8, 69.7 and 52.5 K. Across the twenty-one families the median ratio of
|bias| to MAE is 0.54, so about half of a typical family-holdout error is
systematic displacement rather than scatter, and its magnitude tracks structural
isolation, the absolute bias correlating with mean nearest-training similarity at
*r* = −0.53 (*p* = 0.014, unadjusted).

A model fitted almost entirely on other families will pull any withheld family
toward the bulk of the training distribution, so a family whose *T*~g~ sits far
from that bulk is displaced toward it for reasons that have nothing to do with
its backbone. Regressing bias on the family's *T*~g~ offset from the rest of the
dataset measures how much of the signal that accounts for (Table 6). The two are
positively correlated, at *r* = +0.425 (*p* = 0.055, twenty-one families):
marginal, but large enough to matter here, because it is the three families just
quoted that the line fits. Polyphosphazenes sit 145 K below the rest of the
dataset, polysiloxanes 148 K below and the silane class 155 K below, and
shrinkage alone accounts for −23.4 K, −24.0 K and −25.3 K of their bias, which is
28%, 38% and 56% of what is observed. They remain the most displaced families
in absolute terms, but they cannot carry the chemical argument by themselves.

The families that can carry it are the ones moving against the shrinkage line.
Polydienes (bias +48.0 K, residual +49.9 K), polyolefins (+37.8 K, +47.7 K),
polyhalo-olefins (+28.4 K, +43.5 K) and polystyrenes (+19.9 K, +23.6 K) are all
*under*predicted even though every one of them sits below the dataset median,
where regression toward the training mean would push the prediction up rather
than down. A model deprived of its flexible unsaturated and saturated hydrocarbon
chains places them too low, having learned its temperature scale from stiffer
condensation backbones, and no shrinkage toward the training mean produces that
sign. The chemical reading of the signed error therefore rests on these four
families rather than on the inorganic ones, and this line of evidence is partial:
for the families with the most extreme *T*~g~ the two explanations are not
separated by the data available here.

**Table 6.** Family-holdout bias against regression toward the training mean,
histogram gradient boosting, one row per withheld family. Offset is the family's
median *T*~g~ minus the median of the rest of the dataset; bias is
mean(*T*~g~ observed − *T*~g~ predicted) over the withheld family; both in K.
Expected bias is the least-squares fit of bias on offset across all twenty-one
families (Pearson *r* = +0.425, *p* = 0.055), and residual is bias minus that
fit. Shown are the four families with the largest positive residuals, the four
with the largest negative, and the abstention class. All twenty-one are in
`results/table_bias_shrinkage.csv`.

| Family | Offset (K) | Bias (K) | Expected from offset (K) | Residual (K) |
|---|---|---|---|---|
| Polydienes | −34.0 | +48.0 | −1.9 | +49.9 |
| Polyolefins | −75.0 | +37.8 | −9.9 | +47.7 |
| Polyhalo-olefins | −101.8 | +28.4 | −15.1 | +43.5 |
| Polystyrenes | −43.0 | +19.9 | −3.7 | +23.6 |
| Other backbone | −154.6 | −45.3 | −25.3 | −20.0 |
| Polyimides | +150.5 | +3.9 | +33.8 | −29.9 |
| Polysiloxanes | −148.0 | −63.7 | −24.0 | −39.6 |
| Polyureas | +35.0 | −43.1 | +11.4 | −54.5 |
| Polyphosphazenes | −145.0 | −83.7 | −23.4 | −60.3 |

Polycarbonates illustrate the benign case from the same angle: they are
surrounded in structure space by polyesters and polyethers, which remain in
training and carry nearly the same carbonyl and ether linkages, so their
withdrawal removes little the model cannot recover, and their bias is a modest
+15.8 K.

Four candidate explanations were tested against deterioration across the
twenty-one families, and only structural isolation survives. Mean
nearest-training similarity correlates strongly and negatively for every model
(Pearson *r* between −0.62 and −0.69, *p* ≤ 0.0027, significant for all four
models after Bonferroni correction over the sixteen tests). The shift between the
family's median *T*~g~ and the rest of the dataset gives Pearson *r* between 0.32
and 0.39 (*p* ≥ 0.079), and the ratio of the family's *T*~g~ interquartile range
to the dataset's gives *r* between 0.15 and 0.32 (*p* ≥ 0.16); neither is
resolved here. Family size is likewise not resolved, at *r* between −0.05 and
−0.35 (*p* ≥ 0.115; for extremely randomised trees, *r* = −0.27, 95% CI
[−0.63, +0.18]). That interval is wide enough to contain a moderate effect in
either direction, so the honest statement is that a size effect is undetectable
at twenty-one families, not that it is absent.

One feature of the design protects the headline contrast. Mean nearest-training
similarity and family size are nearly orthogonal across the twenty-one families
(*r* = 0.10, *p* = 0.67), so the similarity result is not a relabelled size
result and the two candidate drivers are not competing for the same variance.

### 3.4. Conformal intervals lose validity exactly where they are needed

All conformal results below use extremely randomised trees as the point
predictor, not the histogram gradient boosting of Section 3.2. Of the two leading
models it is the only one that exposes a per-tree ensemble spread, which the
sparsification diagnostic of Section 2.7 needs, and the accuracy it gives up is
0.38 K of random-split MAE, consistent in sign across nine of the ten repeats but
small against every effect this paper measures (Section 3.2). The paper's two
headline results therefore come from different
models: the 28.1 K to 48.1 K accuracy gap from histogram gradient boosting, every
coverage figure from extremely randomised trees.

Split conformal behaves as advertised on random partitions: 0.904 empirical
coverage against a nominal 0.90, as a mean over the ten splits. Under shift the
guarantee degrades with structural novelty — 0.834 on scaffold splits, 0.730 on
cluster splits and 0.701 under family holdout, all mean-of-splits (Figure 4a;
Table 7). This is not a defect of the method but of its premise: conformal
validity requires calibration and test data to be exchangeable, and a
chemistry-aware holdout is constructed precisely to break that.

Marginal coverage also conceals a failure that is present even when it holds. On
random splits, where overall coverage is a healthy 0.904, coverage among the
structures whose nearest training analogue lies below Tanimoto 0.4 is 0.657
pooled over the ten splits (Figure 5); the split-averaged worst band, which takes
each split's weakest band before averaging, is 0.650 (Figure 4b). A nominally 90%
interval covers two-thirds of the repeat units that a screening campaign would
actually be evaluating. A model can pass every validity check in current practice
and still be wrong about a third of the unfamiliar polymers it is deployed on.

Neither obvious remedy works. The normalised predictor, which scales residuals by
a fitted difficulty estimate, does lift the least-similar band: from 0.657 to
0.759 on random splits and from 0.552 to 0.629 on scaffold splits, both
band-resolved. Its marginal coverage still falls away under shift, to 0.825 on
scaffold, 0.711 on cluster and 0.709 under family holdout, mean-of-splits, and
under cluster holdout that is the worst marginal coverage of the four methods.
Conditioning the quantile on polymer family fails for a structural reason rather
than an empirical one: under family holdout the withheld family has no
calibration members by construction, so the predictor falls back to the global
quantile for 100% of test structures and reproduces split conformal exactly
(0.701 mean-of-splits, identical interval widths).

Conditioning on nearest-training similarity is the variant that survives, because
that coordinate is defined for every repeat unit whether or not its family was
seen. It holds 0.907, 0.895, 0.886 and 0.877 mean-of-splits across the four
regimes, never requiring a fallback under family holdout. On random splits it
achieves this while being *narrower* than split conformal (129.1 K against
132.3 K), and it widens only where the chemistry demands it, to 155 K under
scaffold, 232 K under cluster and 208 K under family holdout (Table 7). The
coverage–width criterion registers the same thing, at 411 under family holdout
against 48,020 for split conformal, but that ratio is almost entirely the
exponential penalty of Section 2.7 acting on a 0.20 coverage shortfall rather
than a statement about width, and it should be read as a restatement of the
coverage gap.

**Table 7.** Conformal performance by split regime and method, extremely
randomised trees, nominal coverage 0.90. Marginal coverage is the mean over
splits (ten for random, scaffold and cluster; twenty-one for family). The
least-similar band is Tanimoto < 0.4, pooled over splits. Width is in K. Fallback
is the mean fraction of test structures whose category had no calibration members
and which therefore received the global quantile. CWC is the coverage–width
criterion with η = 30. Source: `results/table_conformal.csv` and
`results/table_applicability_domain.csv`.

| Regime | Method | Coverage (mean of splits) | Band < 0.4 (pooled) | Width (K) | Fallback | CWC |
|---|---|---|---|---|---|---|
| Random | SCP | 0.904 | 0.657 | 132.3 | 0.00 | 132 |
| Random | Normalised SCP | 0.901 | 0.759 | 127.0 | 0.00 | 127 |
| Random | Mondrian, family | 0.899 | 0.696 | 133.1 | 0.08 | 136 |
| Random | Mondrian, similarity | 0.907 | 0.888 | 129.1 | 0.00 | 129 |
| Scaffold | SCP | 0.834 | 0.552 | 126.5 | 0.00 | 923 |
| Scaffold | Normalised SCP | 0.825 | 0.629 | 124.1 | 0.00 | 1,165 |
| Scaffold | Mondrian, family | 0.832 | 0.576 | 127.2 | 0.10 | 990 |
| Scaffold | Mondrian, similarity | 0.895 | 0.938 | 154.8 | 0.00 | 179 |
| Cluster | SCP | 0.730 | 0.655 | 126.4 | 0.00 | 20,920 |
| Cluster | Normalised SCP | 0.711 | 0.667 | 127.7 | 0.00 | 37,007 |
| Cluster | Mondrian, family | 0.758 | 0.694 | 141.8 | 0.21 | 10,101 |
| Cluster | Mondrian, similarity | 0.886 | 0.825 | 232.2 | 0.13 | 356 |
| Family | SCP | 0.701 | 0.636 | 121.7 | 0.00 | 48,020 |
| Family | Normalised SCP | 0.709 | 0.650 | 129.5 | 0.00 | 39,813 |
| Family | Mondrian, family | 0.701 | 0.636 | 121.7 | 1.00 | 48,020 |
| Family | Mondrian, similarity | 0.877 | 0.890 | 208.3 | 0.00 | 411 |

Whether a global widening could have bought the same result is worth testing
directly, because if it could the taxonomy would be doing no work. The oracle
control of Section 2.7 multiplies every split conformal interval about its centre
by a single constant *k* (Table 8, Figure 6). The constant can be chosen two
ways, and the choice decides what the comparison means. Matching on mean *width*
leaves the
inflated predictor holding more marginal coverage than the novelty-conditioned
one, 0.908 against 0.898 under family holdout and 0.905 against 0.870 under
cluster, so a band comparison at that setting is not like for like: part of any
advantage it shows is simply the extra coverage it was handed. Matching on
marginal *coverage* fixes the coverage budget and asks only how each predictor
allocates it across the similarity range, which is the question the taxonomy
exists to answer. The coverage-matched comparison is the primary one here.

At matched marginal coverage the novelty-conditioned predictor is not beaten in
either regime. Under family holdout the constant that reproduces its 0.898 pooled
coverage is *k* = 1.599, and it covers 0.835 of the least-similar band
(Tanimoto < 0.4, 2,614 pooled predictions) against 0.890. Those are paired
predictions on identical structures, so the discordant pairs settle the
comparison: the novelty-conditioned predictor covers 164 structures the inflated
one misses and misses 21 that it covers, McNemar *p* = 1.1 × 10⁻²⁸. Under cluster
holdout the matched constant is *k* = 1.568 and the two are indistinguishable in
that band, 0.822 against 0.825, on 685 discordant pairs one way and 644 the other
(McNemar *p* = 0.27). The taxonomy therefore wins decisively under family
holdout, where it never falls back to the global quantile, and ties under cluster
holdout, where 12.8% of test structures fall in a similarity band with no
calibration members and receive the global quantile anyway (Table 7;
Section 4.3).

What it wins is allocation and not efficiency, and that should be said plainly.
At matched coverage the novelty-conditioned intervals are the *wider* of the two
in both regimes: 202.5 K against 194.1 K under family holdout and 226.3 K against
197.9 K under cluster. A uniform inflation buys the same marginal coverage for
less total width. What it cannot buy is the same distribution of that width over
the similarity range.

A reader may reasonably ask the width-matched question as well, since width is
the resource an interval spends, and the answer points the same way under family
holdout. At equal mean width there, *k* = 1.669, the inflated predictor still
covers only 0.849 of the least-similar band against 0.890 (134 discordant pairs
against 27, McNemar *p* = 3.0 × 10⁻¹⁸), and it does so while over-covering
marginally at 0.908. Matching the novelty-conditioned band coverage by inflation
alone takes *k* = 1.93, which carries the mean interval from 203 K to 234 K;
reaching nominal 0.90 in that band takes *k* = 2.0, at 243 K and 0.945 marginal.
Under cluster holdout the equal-width constant *k* = 1.792 does reach 0.864 in
that band against 0.825, but at 0.905 marginal coverage against 0.870, which is
where the difference comes from: the coverage-matched row above removes it.

Both controls are oracles. The constant is computed from the test set it is then
compared against, so none of 1.568, 1.599, 1.669, 1.792, 1.93 or 2.0 is knowable
at prediction time, and it is regime-dependent, so a constant chosen on one
holdout has no claim on another. Their role is to bound what any global inflation
of split conformal could achieve, not to be deployed.

**Table 8.** Oracle inflation controls, extremely randomised trees, nominal
coverage 0.90. Split conformal intervals are multiplied about their centre by the
constant *k*. The coverage-matched row in each block is the *k* that reproduces
the novelty-conditioned predictor's marginal coverage on the same test set and is
the primary comparison; the width-matched row is the *k* that equalises mean
width instead, and is handed more marginal coverage in the process. Neither is
available at prediction time. Coverage is pooled over all held-out predictions in
the regime, 7,174 under family holdout and 24,985 under cluster holdout; the
least-similar band is Tanimoto < 0.4 and holds 2,614 and 12,846 of them. Width in
K. The full sweep over *k* is in `results/table_width_matched.csv`.

| Regime | Predictor | *k* | Coverage (pooled) | Band < 0.4 (pooled) | Width (K) |
|---|---|---|---|---|---|
| Family | Split conformal | 1.000 | 0.743 | 0.636 | 121.4 |
| Family | Coverage-matched oracle | 1.599 | 0.898 | 0.835 | 194.1 |
| Family | Mondrian, similarity | — | 0.898 | 0.890 | 202.5 |
| Family | Width-matched oracle | 1.669 | 0.908 | 0.849 | 202.5 |
| Family | Inflated to match band | 1.930 | 0.939 | 0.890 | 234.2 |
| Family | Inflated to nominal band | 2.000 | 0.945 | 0.900 | 242.7 |
| Cluster | Split conformal | 1.000 | 0.716 | 0.655 | 126.3 |
| Cluster | Coverage-matched oracle | 1.568 | 0.870 | 0.822 | 197.9 |
| Cluster | Mondrian, similarity | — | 0.870 | 0.825 | 226.3 |
| Cluster | Width-matched oracle | 1.792 | 0.905 | 0.864 | 226.3 |

### 3.5. An applicability domain, and an honest limit on the repair

Resolving coverage by similarity band under family holdout (Table 9) turns the
result into an operational rule, and into an applicability-domain statement of
the kind conformal prediction has been used to make since Norinder et al. [21]:
the interval, rather than a separate distance threshold, carries the domain
information. All coverage in this section is band-resolved:
pooled within a band over the twenty-one splits, which together hold out each of
the 7,174 structures exactly once. Split conformal covers 0.636 of structures
whose nearest training analogue lies below Tanimoto 0.4, and 0.771, 0.809 and
0.874 in the next three bands. Its intervals are essentially constant in width
(120.7–122.9 K) across all six bands, while the absolute error across those bands
varies by a factor of 1.83, from 55.5 K in the least-similar band to 30.3 K in
the 0.6–0.7 band, its lowest (Figure 7). Applying one width to a 1.8-fold spread
in error is precisely the failure mode. Novelty-conditioned calibration instead
holds 0.890, 0.919, 0.907 and 0.937 in the four least-familiar bands, at widths
that scale from 246 K down to 149 K.

The repair has a cost that should be stated plainly. In the most-similar band
(Tanimoto ≥ 0.8) the novelty-conditioned predictor *under*-covers, at 0.673
against split conformal's 0.788, because it narrows those intervals to 88 K. Under
a family holdout a test structure that closely resembles training data is unusual,
and the calibration points populating that band come from families where high
similarity genuinely did signal an easy prediction. The trade is nonetheless
favourable in practice: that band holds 156 of 7,174 held-out predictions
(2.2%), while the two least-similar bands hold 4,647 (64.8%). The recommendation
is therefore to report novelty-conditioned intervals together with the
nearest-training similarity itself, so that a reader can see which regime a given
prediction sits in rather than trusting a single interval uniformly.

**Table 9.** Coverage, interval width and absolute error resolved by
nearest-training Tanimoto similarity band under family holdout, extremely
randomised trees, nominal coverage 0.90. The twenty-one splits hold out each of
the 7,174 structures exactly once, so *n* sums to 7,174 and all quantities are
pooled within a band across splits. Widths and errors in K. Error is not monotone
in similarity: it falls to 30.3 K in the 0.6–0.7 band and rises again above 0.8.
Sources: `results/table_applicability_domain.csv`, `results/table_band_mae.csv`.

| Band | *n* | MAE (K) | SCP coverage | SCP width (K) | Mondrian coverage | Mondrian width (K) |
|---|---|---|---|---|---|---|
| < 0.4 | 2,614 | 55.5 | 0.636 | 121.2 | 0.890 | 245.5 |
| 0.4–0.5 | 2,033 | 42.1 | 0.771 | 121.3 | 0.919 | 205.5 |
| 0.5–0.6 | 1,403 | 37.8 | 0.809 | 121.3 | 0.907 | 176.7 |
| 0.6–0.7 | 605 | 30.3 | 0.874 | 121.5 | 0.937 | 148.7 |
| 0.7–0.8 | 363 | 31.2 | 0.854 | 122.9 | 0.832 | 115.0 |
| ≥ 0.8 | 156 | 36.7 | 0.788 | 120.7 | 0.673 | 87.9 |

### 3.6. The penalty follows which structures were removed, not how many

The matched design separates the two effects that a conventional family holdout
confounds (Figure 3, Table 10). Across twenty-one families and three repeats, with
the test set held identical and the two comparison pools held identical in size:

* Removing an equal quantity of *unrelated* structures — the **control** arm —
  changes MAE by +0.29 K (95% CI [−0.63, +1.24], Wilcoxon *p* = 0.71, worse in
  9 of 21 families) for extremely randomised trees, and by −0.58 K
  (95% CI [−1.81, +0.64], *p* = 0.45, worse in 10 of 21) for histogram gradient
  boosting. The two models do not agree even on the sign.
* Removing *the family itself* — the **naive** arm — costs +14.97 K
  (95% CI [+11.21, +19.28]) for extremely randomised trees and +15.15 K
  (95% CI [+11.02, +20.03]) for histogram gradient boosting, and **is worse in
  all twenty-one families for both models**. With forty-two of forty-two signs
  positive the Wilcoxon test returns its floor, *p* = 9.5 × 10⁻⁷ per model; the
  sign count is the result, not the *p*-value.

**Table 10.** Matched-pair effects in K, over twenty-one families and three
repeats per family, with the test set held identical within each triple. Family
effect is naive − control, size effect is control − informed, total effect is
naive − informed; positive means worse. Each family contributes the mean of its
three repeats, the mean column averages those twenty-one family means, and the
95% interval is a percentile bootstrap over them, to be read as dispersion across
families rather than inference about a population (Section 2.8). The Wilcoxon
*p* of 9.5 × 10⁻⁷ is the floor of the test at twenty-one paired values, reached
whenever all signs agree. Source: `results/table_matched_statistics.csv`.

| Model | Effect | Mean (K) | 95% CI (K) | Median (K) | Wilcoxon *p* | Positive |
|---|---|---|---|---|---|---|
| Extremely randomised trees | Family | +14.97 | [+11.21, +19.28] | +12.86 | 9.5 × 10⁻⁷ | 21 / 21 |
| Extremely randomised trees | Size | +0.29 | [−0.63, +1.24] | −0.15 | 0.71 | 9 / 21 |
| Extremely randomised trees | Total | +15.26 | [+11.35, +19.83] | +12.55 | 9.5 × 10⁻⁷ | 21 / 21 |
| Histogram gradient boosting | Family | +15.15 | [+11.02, +20.03] | +11.28 | 9.5 × 10⁻⁷ | 21 / 21 |
| Histogram gradient boosting | Size | −0.58 | [−1.81, +0.64] | −0.02 | 0.45 | 10 / 21 |
| Histogram gradient boosting | Total | +14.57 | [+10.63, +19.35] | +12.11 | 9.5 × 10⁻⁷ | 21 / 21 |

The two effects are therefore not merely different in size but different in
character. The family effect is roughly 15 K and unanimous in sign across every
family and both models. The size effect is under 1 K in magnitude, has opposite
signs in the two models, and scatters around zero per family, from −3.6 K to
+4.2 K for extremely randomised trees and from −8.9 K to +7.3 K for histogram
gradient boosting. We deliberately do not quote a ratio between them: with the
denominator's sign unresolved, a ratio is not a quantity. What the design
establishes is that the total effect, +15.26 K and +14.57 K for the two models,
is accounted for by the family term alone, and it does so on matched test
structures, so it is not an artefact of families differing in intrinsic
difficulty or in *T*~g~ spread.

What the two arms actually differ in is the *locality* of the removal, not its
volume: the naive arm deletes a connected region of structure space, the control
arm a diffuse sample of the same size scattered across the other twenty families.
The chemical reading of that contrast is developed in Section 4.2; the
measurement is the locality contrast.

Family-level magnitudes rank as the chemistry predicts. For extremely randomised
trees, polyphosphazenes suffer most (+43.7 K), followed by polysiloxanes
(+30.6 K) and polyamides (+25.2 K); at the other end polyimines (+1.4 K),
polyvinyls (+4.1 K) and polyethers (+5.8 K) are nearly unaffected, because the
linkages that set their backbone mobility survive elsewhere in the training set.
Note that polyamides rank third here despite being the second-largest family in
the dataset (999 structures) — further evidence that abundance does not confer
transferability when the chemistry is distinctive.

A near-duplicate objection has to be answered separately. The informed arm keeps
the rest of the family, so for some test structures it holds a very close
analogue that the naive arm cannot have, and the family effect might be no more
than the loss of those analogues. Splitting the effect by how close the informed
arm's nearest training neighbour actually was tests that directly (Table 11). The
effect does not vanish where no near analogue existed. For extremely randomised
trees it is +11.4 K over the 1,426 test structures whose nearest informed-arm
neighbour lay below Tanimoto 0.6, against +15.6 K over the 1,655 whose neighbour
lay above 0.9; for histogram gradient boosting the same two bands give +11.6 K
and +13.7 K. The effect is larger where a near-duplicate was available, which is
what one would expect, but most of it is present where none was.

**Table 11.** Matched-pair family effect resolved by how close the informed arm's
nearest training neighbour was, in K, over twenty-one families and three repeats.
Bands are the informed arm's nearest-training Tanimoto similarity for the same
test structure; *n* is the number of test structures in the band, each appearing
once per repeat, paired across arms. Arm means are MAE in K over those
structures, and the family effect is the paired mean of naive − control. Source:
`results/table_matched_by_analogue.csv`.

| Model | Informed-arm band | *n* | Informed | Control | Naive | Family effect |
|---|---|---|---|---|---|---|
| Extremely randomised trees | < 0.6 | 1,426 | 45.48 | 43.63 | 55.05 | +11.42 |
| Extremely randomised trees | 0.6–0.8 | 2,307 | 29.89 | 30.44 | 46.18 | +15.74 |
| Extremely randomised trees | 0.8–0.9 | 1,071 | 23.47 | 23.87 | 38.52 | +14.65 |
| Extremely randomised trees | ≥ 0.9 | 1,655 | 17.83 | 20.48 | 36.09 | +15.61 |
| Histogram gradient boosting | < 0.6 | 1,426 | 43.95 | 41.63 | 53.27 | +11.63 |
| Histogram gradient boosting | 0.6–0.8 | 2,307 | 29.35 | 29.28 | 42.61 | +13.32 |
| Histogram gradient boosting | 0.8–0.9 | 1,071 | 23.88 | 24.03 | 36.10 | +12.07 |
| Histogram gradient boosting | ≥ 0.9 | 1,655 | 18.46 | 20.12 | 33.83 | +13.71 |

### 3.7. The gap is not an artefact of the descriptor set

A natural objection is that the generalisation gap reflects a limitation of RDKit
descriptors rather than of the data, and that a representation better suited to
novel chemistry would close it. Repeating the random and family-holdout regimes
with binary Morgan fingerprints (radius 2, 2,048 bits) in place of the descriptor
block tests this directly, since the two representations share almost nothing:
one is a set of physically motivated aggregate quantities, the other a sparse
record of local substructural environments. Both use extremely randomised trees.
The two family arms are the full twenty-one holdouts, but the random arms differ
in repeat count: five splits for Morgan against ten for descriptors.

The gap survives the substitution and is not narrower (Table 12). Fingerprints are
worse in absolute terms in both regimes — 35.8 K against 28.5 K on random splits
and 62.7 K against 49.6 K under family holdout — and their family-holdout penalty
is larger in absolute terms (+26.9 K against +21.1 K) while being similar as a
ratio (1.752 against 1.738). We do not attach an interval to the difference
between those two ratios, because the random arms behind them rest on different
numbers of repeats; the statement supported is that substructural environments
offer no protection here in aggregate.

**Table 12.** Representation ablation, extremely randomised trees. MAE in K as the
mean over splits, with the standard deviation across splits after ± and the
number of splits in parentheses. Penalty is family MAE minus random MAE and the
ratio is their quotient. Sources:
`results/table_representation_summary.csv`,
`results/table_representation_ablation.csv`.

| Representation | Random MAE (K) | Family MAE (K) | Penalty (K) | Penalty ratio |
|---|---|---|---|---|
| RDKit descriptors, 217 columns | 28.51 ± 0.93 (10) | 49.56 ± 16.19 (21) | +21.05 | 1.738 |
| Morgan, radius 2, 2,048 bits | 35.81 ± 0.92 (5) | 62.73 ± 15.38 (21) | +26.92 | 1.752 |

The per-family picture is not uniform, and one plausible mechanism is ruled out
by it. Fingerprints are worse than descriptors on eighteen of the twenty-one
withheld families, but they are *better* on three, and the largest single
improvement is on polyphosphazenes (62.7 K against 87.1 K), the family whose
absent P=N bits would be the obvious reason for a fingerprint model to fail. Nor
are they better on polysiloxanes (77.6 K against 70.4 K). An explanation in terms
of missing backbone bits therefore does not survive contact with the per-family
data, and we do not offer one. What the ablation shows is narrower and still
useful: the gap is not specific to the descriptor block, so it cannot be
dismissed as an artefact of that feature set. Whether some other representation —
a model pretrained on millions of hypothetical polymers, for instance — could
narrow it is untested here and is the natural next step (Section 4.5).

### 3.8. The calibration set is a second lever, with a caveat

Sections 3.4 and 3.5 repair conditional coverage by changing how the conformal
*quantile* is computed. A second possibility is to leave the quantile alone and
change what the calibration set contains. Standard split conformal draws
calibration points at random from the training pool, so every calibration residual
describes an interpolation — the model had relatives of that repeat unit in
training. Under family holdout the test residuals describe extrapolation, so the
quantile is calibrated against the wrong population. Rebuilding the calibration
set by holding families out *within* the training pool produces residuals of the
kind the test set will actually demand.

Both levers work, and to a similar degree (Table 13).

**Table 13.** Calibration-set construction under family holdout, extremely
randomised trees, nominal coverage 0.90, with the point model identical in all
four arms. Coverage is the mean over the twenty-one splits. Worst band is the
split-averaged worst band defined in Section 2.7: within each split the weakest
similarity band holding at least ten structures is taken and those minima are
then averaged, which is why it sits below the band-resolved coverage of the
least-similar band in Table 9. Dropping the ten-structure restriction moves the
first row from 0.564 to 0.539. Width in K, as a mean over splits; *n* calibration
is the mean calibration-set size per split. Source:
`results/table_calibration_study.csv`.

| Calibration set | Quantile | Coverage (mean of splits) | Worst band (split-averaged) | Width (K) | *n* calibration |
|---|---|---|---|---|---|
| Random | Global (SCP) | 0.701 | 0.564 | 122 | 1,366 |
| Random | Novelty-conditioned | 0.877 | 0.753 | 208 | 1,366 |
| Family-out | Global (SCP) | 0.855 | 0.753 | 187 | 3,349 |
| Family-out | Novelty-conditioned | 0.863 | 0.778 | 189 | 3,349 |

Rebuilding the calibration set lifts plain split conformal from 0.701 to 0.855
mean-of-splits without touching the quantile rule, the model or the representation. It does so by
widening intervals from 122 K to 187 K, which is the appropriate response: the
residuals it now calibrates against are genuinely larger. The construction is not,
however, the only thing that changes in that swap. The family-out construction
also enlarges the calibration set from a mean of 1,366 structures to 3,349,
because residuals are pooled over the held-out training families rather than
drawn once. A larger calibration set tightens the quantile's own sampling error,
so the 0.701 → 0.855 improvement is attributable to the two changes jointly and
this experiment cannot separate them.

The two repairs are largely substitutes rather than complements. Applying both
gives 0.863 marginal coverage, no better than the novelty-conditioned quantile
alone at 0.877, because they address the same root cause from opposite ends —
that calibration residuals understate extrapolation error. The combination has
the best split-averaged worst band of the four arms, 0.778 against 0.753 and
0.753, but that margin is not resolved: across the twenty-one families the combination beats
the novelty-conditioned quantile alone on worst-band coverage in only 9, with a
mean difference of +0.025 and a Wilcoxon *p* of 0.40. It should be read as no
worse rather than as better.

That comparison depends on how band coverage is aggregated, and the dependence
should be disclosed. Under the band-resolved convention used in Table 9, the
least-similar band gives 0.890 for the novelty-conditioned quantile alone, 0.824
for family-out split conformal and 0.876 for the two combined, against 0.636 for
plain split conformal. Read that way the novelty-conditioned quantile is ahead of
both family-out arms rather than level with them, and the statement that both
levers work to a similar degree holds only under the split-averaged convention.
The split-averaged figures are quoted here because they are the more conservative
of the two against this section's own recommendation, which is to use the
novelty-conditioned quantile.

The practical guidance follows the cost. Novelty-conditioned calibration is free:
it reuses the same fitted model and calibration set, requiring only that
similarity be computed. Family-out calibration requires refitting the model once
per held-out training family — here eight additional fits per split. For most
purposes the novelty-conditioned quantile alone is the better trade.

### 3.9. A Van Krevelen baseline deteriorates less and predicts worse

The additive schemes that machine-learned models propose to replace make a
different inductive assumption: that *T*~g~ is a composition-weighted sum of
group contributions, with no interaction terms and no dependence on which
linkages the groups sit in. If the family-holdout penalty is really about absent
linkage chemistry, a model that never learns linkage-specific behaviour in the
first place should lose less when a linkage type is withdrawn. It does, and it
pays for that in accuracy everywhere (Table 14).

**Table 14.** Van Krevelen group-contribution baseline against the ensemble, MAE
in K as the mean over splits (ten for random, scaffold and cluster; twenty-one
for family). Penalty is family MAE minus random MAE and the ratio is their
quotient. The additive contributions *Y*~i~ are refitted on each training fold,
so the baseline receives the same data as the ensemble. Source:
`results/table_baseline_comparison.csv`.

| Model | Random | Scaffold | Cluster | Family | Penalty (K) | Ratio |
|---|---|---|---|---|---|---|
| Group contribution | 43.84 | 43.80 | 54.86 | 59.80 | +15.96 | 1.364 |
| Histogram gradient boosting | 28.13 | 33.81 | 44.47 | 48.11 | +19.98 | 1.710 |
| Training-median baseline | 93.79 | 95.71 | 121.61 | 95.46 | +1.67 | 1.018 |

The group-contribution baseline of Section 2.6 gives 43.8 K MAE on random splits,
43.8 K under scaffold holdout, 54.9 K under cluster holdout and 59.8 K under
family holdout. The absolute penalty, random to family, is +16.0 K against the
ensemble's +20.0 K, and the ratio is 1.36 against 1.71. The flat random-to-scaffold
step is the sharper observation: 43.84 K to 43.80 K, no change at all, where the
ensemble loses 5.7 K. A model whose only inputs are group counts is indifferent
to whether it has seen the scaffold before.

The honest way to read this is absolute first. The baseline is worse than the
ensemble in every regime, by 15.7 K on random splits and by 11.7 K under family
holdout, and its smaller penalty ratio partly reflects a worse starting point:
it has less accuracy to lose. Resolved by family, it beats the ensemble on only
two of the twenty-one withheld families, polysiloxanes (68.2 K against 69.7 K)
and polysulfones (35.3 K against 36.7 K), both narrowly. It is not an
alternative to the ensemble for prediction.

It is nonetheless independent support for the mechanism, from a model class with
a different inductive bias and a different failure mode. Compositional additivity
buys transferability at the cost of accuracy; the ensemble makes the opposite
trade, and the size of the trade is what Sections 3.6 and 4.2 are about.

### 3.10. Sensitivity to the family taxonomy

Every family-dependent result inherits the labels, so the taxonomy is a
legitimate target for a referee. During this work a boundary-split bug was found
and fixed in the classifier: because the backbone was taken as the path between
the two attachment points, a repeat unit written cut at its own characteristic
linkage lost that linkage, so nylon-6 written ``*NCCCCCC(=O)*`` read as a plain
hydrocarbon chain and was labelled a polyvinyl. The dimer construction of
Section 2.2 is the fix, and the same revision replaced a catch-all hydrocarbon
assignment with the Other backbone abstention class. Together these reassigned
204 of the 7,174 structures, 2.8% of the dataset, the largest movements being 88
polyesters that are really polycarbonates, 46 polyolefins carrying a backbone
heteroatom, 23 polyethers that are polyesters, and 15 polyvinyls that are
polyamides.

Rerunning the whole pipeline against both label sets gives a direct sensitivity
measurement. The random, scaffold and cluster results are unchanged to three
decimal places for every model, as they must be, since those regimes never
consult a family label; pooled *R*² in those three regimes is likewise identical
to three decimals. Under family holdout, where the labels do enter, pooled *R*²
moves by at most 0.006 across the five models — the largest is support vector
regression at 0.0052 — and is unchanged at 0.746 for histogram gradient boosting.
The remaining family-dependent results move by about a kelvin and no conclusion
changes sign: family-holdout MAE goes from 46.8 K to 48.1 K for histogram
gradient boosting and from 48.9 K to 49.6 K for extremely randomised trees, and
the matched-design family effect from +14.61 K to +14.97 K and from +13.28 K to
+15.15 K. The 1.29 K move in family-holdout MAE has two parts, because the
revision also created a twenty-first holdout group. Averaged over the same twenty
groups the old labels supplied, the corrected labels give 47.89 K, so relabelling
accounts for +1.07 K of the move and the new Other backbone group for the
remaining +0.22 K. The one qualitative change
is in the paper's favour. Under the old labels the family effect was positive in
20 of 20 families for one model but only 18 of 20 for the other, the two
exceptions being polycarbonates and polysulfides; polycarbonates was also the
family most altered by the correction, nearly doubling from 98 to 189 structures
as the misassigned polyesters returned to it. Under the corrected labels the
effect is positive in 21 of 21 families for both models. The earlier anomaly was
a labelling error, not a counterexample.

This is a sensitivity check against one specific perturbation, not a proof of
robustness to any taxonomy. A chemically different but defensible scheme — one
that split polyimides by dianhydride class, say — would move the family-level
numbers again. What the check establishes is the order of magnitude of that
movement for a 2.8% relabelling: about a kelvin on aggregate errors, with the
sign structure of the matched design intact.

## 4. Discussion

The result a polymer chemist can act on is the family-resolved one. A *T*~g~
model does not have a single accuracy: histogram gradient boosting predicts
withheld polyvinyls to 30.2 K and withheld polyphosphazenes to 84.8 K, and it
displaces whole backbone classes in a consistent direction rather than
scattering about them. What follows takes that structure apart: what an aggregate leaderboard
number conceals (Section 4.1), why the penalty follows which chemistry is missing
rather than how much data is missing (Section 4.2), why an interval conditioned
on structural novelty survives where one conditioned on polymer family cannot
(Section 4.3), what to do about it (Section 4.4), and what the evidence here does
not support (Section 4.5).

### 4.1. What a *T*~g~ leaderboard is actually measuring

The most consequential number in this study is not an error but a contrast.
Changing the model changes MAE by 0.38 K between the two best and by 2.74 K
across all four non-trivial regressors; changing the partitioning protocol from
random to family holdout changes it by 20.0 K. Any comparison of architectures on
this dataset that does not fix and report the split regime is dominated by a
variable it does not control.

This reframes a decade of incremental gains. Reported improvements from kernel
methods to tree ensembles to graph networks to pretrained chemical language models
are real, but they are improvements in interpolation, measured where structural
analogues are plentiful. None of the comparisons establishes that the newer model
extrapolates better, because none was evaluated under a regime that requires
extrapolation. That is not a criticism of the models; it is a criticism of the
protocol they were scored under, and it is cheap to fix.

### 4.2. Why chemistry and not volume

The matched design measures a locality effect: at fixed training-set size,
deleting a coherent chemical neighbourhood costs about 15 K on that
neighbourhood, deleting a diffuse sample of the same size costs nothing that
survives a change of model. The chemical reading of that asymmetry is the
straightforward one. Descriptor-based models of *T*~g~ learn how particular
backbone linkages translate into chain stiffness and interchain interaction, and
those mappings are local to a linkage type. Remove every phosphazene, and no
amount of additional polyester data teaches the model what a P=N backbone does to
segmental mobility; it extrapolates organic intuition into inorganic chemistry
and overpredicts those chains by 84 K on average (Section 3.3). Four lines of
evidence support that reading, one of them only in part. Deterioration correlates
with structural isolation while no size effect is resolved (Section 3.3). The
family effect is not the loss of near-duplicates: it is +11.4 K even among test
structures for which the informed arm itself held no training neighbour above
Tanimoto 0.6 (Section 3.6). An additive model that never learns linkage-specific
behaviour deteriorates proportionally less when a linkage type is withdrawn
(Section 3.9). The fourth line, that the withheld-family errors are systematic
displacements of a consistent sign, holds only partly, because for the families
with the most extreme *T*~g~ regression toward the training mean predicts the
same sign. It rests on the families that move against that line: polydienes,
polyolefins, polyhalo-olefins and polystyrenes, all underpredicted while sitting
below the dataset median (Section 3.3).

The consequence for data-collection strategy needs stating carefully, because the
experiment removes data rather than adding it. What was tested is that
withdrawing a median of 1.85% of the training pool at random, from chemistry
already covered, does not measurably change accuracy on a held-out family, while
withdrawing the same quantity as a coherent family does. That says a marginal
fraction of the existing pool carries little information about a new backbone. It
does not establish that *adding* data within covered chemistry would be useless,
which is a different experiment at a different scale and one these results cannot
settle. The defensible version of the recommendation is the comparative one: at
equal cost, structures from an uncovered family are worth more than structures
from a covered one, and a screening programme that plans its data collection
around coverage of the structure space rather than around volume is acting on
what was measured.

### 4.3. Why novelty-conditioning works and family-conditioning cannot

The contrast between the two Mondrian taxonomies is a structural fact rather than
an empirical accident, and it generalises beyond this dataset. A conditional
conformal predictor can only offer category-specific validity for categories it
has calibration data for. Polymer family fails this test in exactly the
deployment scenario of interest: the family is new, so it has no calibration
members, and the predictor silently degrades to the unconditional one, here for
100% of test structures. Nearest-training similarity is computable for any repeat
unit against any training set, its value is available before the label is, and
its bands are populated by calibration structures drawn from every family, so the
taxonomy remains well defined precisely where the family-based one collapses. The
same argument recommends similarity — or any other deployment-computable novelty
coordinate — over class-based taxonomies for conditional conformal prediction in
materials problems generally, wherever the classes of interest are the ones the
model has not seen.

The similarity taxonomy is not, however, free of the same defect in kind. A
similarity band can also be empty of calibration members, and under cluster
holdout it is: 12.8% of test structures fall back to the global quantile there,
against 21.5% for the family taxonomy. The difference is one of degree and of
prognosis. A family holdout empties the relevant family category by construction,
every time, and no amount of additional calibration data fixes it; a similarity
band is emptied only when the calibration pool happens to contain nothing that
novel, which calibration data drawn more widely can in principle remedy. The
advantage is
real but it should be stated as a smaller and repairable blind spot rather than
as none.

Section 3.4 tested whether a global widening explains the advantage, and at
matched marginal coverage it does not. An oracle inflation of split conformal set
to the same marginal coverage covers 0.835 of the least-similar band under family
holdout against 0.890, losing the paired comparison decisively (McNemar
*p* = 1.1 × 10⁻²⁸), and 0.822 against 0.825 under cluster holdout, which is a tie
(*p* = 0.27). The taxonomy is ahead where it never falls back and level where it
falls back for 12.8% of structures. What it buys is allocation rather than
efficiency: at matched coverage its intervals are the wider of the two, 202.5 K
against 194.1 K under family holdout, so the claim is about where the width goes
and not about spending less of it.

A weighted conformal correction is the other principled response to a covariate
shift of this kind [30], and it is not used here for a practical
reason. The weights it needs are the likelihood ratio between the test and
training covariate distributions, estimated from unlabelled test data. With 44 to
1,707 structures per family in a 217-dimensional descriptor space, and with the
shift being the wholesale absence of a chemistry rather than a reweighting of a
shared support, that ratio is not estimable at the sample sizes available. A
novelty-conditioned taxonomy asks less of the data: it needs only a coordinate
that is computable at prediction time and calibration members in each band.

### 4.4. Practical recommendations

For anyone building or reporting a *T*~g~ model intended for screening:

1. **Report the split regime as a primary result, not a methods detail**, and
   report at least one chemistry-aware regime alongside any random split. The gap
   between them is the honest estimate of screening performance.
2. **Repeat every regime that can be repeated, and say what the spread means.**
   Scaffold and cluster splits vary by ±8–10 K across repeats. A family holdout
   is run once per family, so its ±15 K is dispersion between families, not
   run-to-run noise; both are worth reporting, and conflating them is not.
3. **Report errors resolved by polymer family, with the sign of the error.** An
   aggregate MAE of 48 K conceals a range from 30 K to 85 K, and a user working
   on polysiloxanes is served by the wrong number; the −64 K bias on that family
   tells them more than either.
4. **Ship intervals, and condition them on a deployment-computable novelty
   coordinate.** Report the nearest-training similarity next to each prediction so
   a reader can see which regime it belongs to.
5. **Fit every data-dependent preprocessing step inside the fold**, including
   descriptor filtering. It is a one-line change to a pipeline and removes an
   unquantified optimism from the reported numbers.

### 4.5. Limitations

The target itself is noisy in a way no structural model can address. Independent
literature reports of the same repeat unit differ by up to 115 K in this dataset,
because a repeat-unit graph cannot express molecular weight, tacticity,
crosslink density, thermal history or measurement protocol. A 28 K random-split
MAE should be read against that floor, and part of the residual is irreducible.

The family taxonomy is our own. Section 3.10 measures its sensitivity to one
2.8% relabelling and finds about a kelvin of movement with every conclusion
intact, but that is one perturbation, not a survey of chemically defensible
alternatives. The validation in Table S1 should also be described accurately: it
is a pinned regression suite, and its twenty-five cases are the ones the
classifier's patterns were written to handle. It protects against future
regressions and it demonstrates that the rules encode the chemistry we claim they
encode; it is not an independent held-out validation set, and it cannot be read
as an accuracy estimate for the taxonomy on unseen chemistry.

The matched design removes data rather than adding it, and removes a median of
1.85% of the pool in the control arm (Section 2.5). Its conclusions are about
where a fixed data budget should be spent, not about the returns to enlarging it.

Our models are descriptor-based ensembles and one additive baseline. We show in
Section 3.7 that the generalisation gap survives a change of representation to
Morgan fingerprints, but we have not tested a pretrained chemical language model
or a graph neural network. Such a model, having seen millions of hypothetical
polymers during pretraining, might have effective coverage of families absent
from this dataset's *labelled* portion. Testing whether pretraining narrows the
family effect is the natural next experiment, and the matched design transfers to
it unchanged.

The signed-error evidence for the chemical mechanism is partial. For families
whose *T*~g~ sits far from the rest of the dataset, regression toward the
training mean and an absent linkage chemistry predict displacement of the same
sign, and the twenty-one family means available here do not separate them
(Section 3.3). The argument rests on the families that move against the shrinkage
line; a design that varied *T*~g~ offset and backbone novelty independently would
settle it, and this one does not.

Finally, the novelty-conditioned predictor under-covers the most-similar band
under family holdout (Section 3.5) and falls back to the global quantile for
12.8% of structures under cluster holdout, where it holds no advantage over an
oracle inflation of split conformal at matched marginal coverage (Section 3.4).
It is also the wider of the two at matched coverage in both regimes, so what it
offers is a better allocation of width rather than less of it. The affected
groups are small and the direction of the trade favours screening, but the method
is an improvement rather than a solution, and conditional coverage should
continue to be reported rather than assumed.

## 5. Conclusions

Machine-learned glass transition temperatures are far more reliable inside
familiar backbone chemistry than outside it, and standard practice measures only
the former. On 7,174 experimental repeat units, moving from random to
chemistry-aware partitions raises the error of the best model from 28.1 K to as
much as 48.1 K and lowers pooled *R*² from 0.87 to 0.75. That 20.0 K shift is
roughly fifty times the difference between the best and second-best model and
seven times the spread across all four non-trivial regressors.

A matched-pair design shows that the penalty follows the locality of what is
removed rather than its quantity: deleting an equal number of unrelated
structures moves MAE by under 1 K, with the two models disagreeing on the sign,
while deleting the family costs about +15 K and does so in all twenty-one
families for both models. Transferability tracks how structurally isolated a
family is; family size is not resolved at this sample size, and the two are
nearly orthogonal, so the first result is not the second in disguise. An additive
group-contribution baseline, which never learns linkage-specific behaviour,
loses proportionally less across the same regimes and predicts worse everywhere,
which is the same mechanism seen from a model class with a different inductive
bias.

The uncertainty attached to these predictions fails in the same place. Split
conformal intervals hold their nominal 90% marginally on random splits while
covering 66% of the least familiar structures, and lose marginal validity
under shift. Conditioning the calibration on polymer family cannot repair this,
because a withheld family has no calibration data by construction. Conditioning
on nearest-training structural similarity does, holding 0.88–0.91 as a mean over
splits in every regime tested and allocating width to the bands where the error
actually is: at matched marginal coverage it covers more of the least familiar
band than a uniform inflation of split conformal under family holdout, 0.890
against 0.835, and as much under cluster holdout, though it spends more total
width to do so. Rebuilding the calibration set from within-training family
holdouts is a second lever that reaches 0.855 with the plain global quantile,
though in this experiment it also enlarges the calibration set from 1,366 to
3,349 structures, so construction and size cannot be separated; combining the two
levers is no better than the novelty-conditioned quantile alone.

The practical implication is a short list: report a chemistry-aware split, repeat
what can be repeated, resolve errors by family and by sign, and attach
novelty-conditioned intervals together with the similarity value itself. None of
these is expensive, and together they turn a *T*~g~ model from something that
scores well into something a screening campaign can act on.

## Figure Captions

These are the captions carried by the seven main figures. They are held here, in
the manuscript source, and the typesetting script reads them from this section
rather than holding a second copy.

**Figure 1.** Reported glass transition temperature by repeat-unit family,
ordered by median; boxes span the interquartile range, whiskers the 5th to 95th
percentiles, and family sizes are listed at the right. The ordering is the one
polymer chemistry predicts, which is evidence that the derived taxonomy is
meaningful rather than merely self-consistent.

**Figure 2.** Test MAE in K by model and split regime; error bars are the
standard deviation across splits. For random, scaffold and cluster that is spread
across ten repeats of the same protocol; for family holdout each family is
withheld exactly once, so the bar is dispersion between families and not
run-to-run noise. The spread between regimes dwarfs the spread between models.

**Figure 3.** Matched-pair design, extremely randomised trees, twenty-one
families with three repeats of the test draw and the control block in each.
(a) The three training arms predicting an identical set of held-out structures,
one row per family, each point the mean of that family's three repeats. Informed
and control are coincident to within a few kelvin and of inconsistent sign
(polydienes 63.1 against 67.4 K, polyhalo-olefins 55.1 against 51.5 K,
polyolefins 44.9 against 41.7 K, largest separation 4.3 K), while naive is
displaced to the right in every family. (b) The family effect, naive minus
control, with training-set size held constant.

**Figure 4.** Coverage of nominally 90% conformal intervals, extremely randomised
trees. (a) Marginal coverage, as a mean over the splits of each regime. (b) The
split-averaged worst band: within each split the weakest similarity band holding
at least ten structures is taken and those minima are averaged. Line style
encodes method alongside colour, so the exact coincidence of split conformal and
the family-conditioned quantile under family holdout remains visible.

**Figure 5.** Conditional coverage against nearest-training Tanimoto similarity,
one panel per split regime, pooled within each band across the splits of that
regime. Split conformal degrades monotonically with structural novelty; the
novelty-conditioned quantile stays near nominal.

**Figure 6.** Coverage in the least-similar band (Tanimoto < 0.4), pooled over
splits, against mean interval width in K, under (a) family and (b) cluster
holdout. The grey line is split conformal with every interval multiplied about
its centre by a constant *k*, swept from 1.0 to 2.2. The two open markers are
oracle choices of *k*: the square reproduces the novelty-conditioned predictor's
marginal coverage, the diamond equalises its mean width. Neither is knowable at
prediction time. The filled point is the novelty-conditioned predictor itself.
The square is the like-for-like comparison; reading the figure along the width
axis alone compares predictors that do not hold the same marginal coverage, which
is the confound Section 3.4 sets out. Source:
`results/table_width_matched.csv`.

**Figure 7.** Mean absolute error in K against nearest-training similarity, for
random and family-holdout splits, extremely randomised trees, pooled within each
band. Error bars are the standard error of the mean of within-band absolute
error. The random-split bars are anticonservative, because the ten repeats draw
on the same 7,174 structures and a structure recurs across them, so the
within-band observations are not independent. Band sizes are not annotated on the
figure; they are in Table 9 for the family holdout and in
`results/table_band_mae.csv` for both regimes.

## Appendix A. Full per-family listing

**Table A1.** Family-holdout error and deterioration for histogram gradient
boosting, all twenty-one families, in full. MAE and RMSE in K; deterioration is
the family's holdout MAE over the model's mean random-split MAE. Table 5 quotes
the three most and three least deteriorated rows of this listing. Source:
`results/table_family_holdout.csv`.

## Supplementary Materials

The following are available in the repository named in the Data Availability
Statement: Table S1, the twenty-five reference polymers with expected and
assigned family; Table S2, family-holdout deterioration against structural
isolation, *T*~g~ shift, *T*~g~ interquartile ratio and family size, with Pearson
and Spearman correlations; Table S3, attachment-point counts per structure;
Table S4, per-band coverage under each calibration construction; Table S5,
absolute error resolved by nearest-training similarity band with the
worst-to-best band ratio; Table S6, per-family group-contribution error under
family holdout.

## Author Contributions

Conceptualisation, [AUTHOR]; methodology, [AUTHOR]; software, [AUTHOR];
validation, [AUTHOR]; formal analysis, [AUTHOR]; investigation, [AUTHOR]; data
curation, [AUTHOR]; writing—original draft preparation, [AUTHOR];
writing—review and editing, [AUTHOR]; visualisation, [AUTHOR]. All authors have
read and agreed to the published version of the manuscript.

## Funding

This research received no external funding.

## Institutional Review Board Statement

Not applicable.

## Informed Consent Statement

Not applicable.

## Data Availability Statement

All code, the curated structure-level table and every result file underlying this
work are available at https://github.com/Rugvedhs/spectrasense, in the directory
``polymer-tg-benchmark``. The underlying experimental *T*~g~ records are
redistributed with the POINT² benchmark [11] and derive from PoLyInfo
[4].

## Acknowledgments

The authors thank the maintainers of the POINT² benchmark and of PoLyInfo for
making the underlying experimental records available, and the RDKit and
scikit-learn developers for the software this work is built on.

## Conflicts of Interest

The authors declare no conflict of interest.

## Abbreviations

The following abbreviations are used in this manuscript:

| | |
|---|---|
| *T*~g~ | glass transition temperature |
| MAE | mean absolute error |
| RMSE | root mean squared error |
| CI | confidence interval |
| SCP | split conformal prediction |
| CWC | coverage–width criterion |
| SMILES | simplified molecular-input line-entry system |
| PSMILES | polymer SMILES, with attachment points marked |
| SMARTS | SMILES arbitrary target specification |
| DSC | differential scanning calorimetry |
| DMA | dynamic mechanical analysis |
