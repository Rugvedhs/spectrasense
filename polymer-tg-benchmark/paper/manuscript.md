---
title: "Backbone Chemistry Governs the Reliability of Machine-Learned Glass Transition Temperatures: A Family-Resolved Benchmark with Novelty-Conditioned Conformal Uncertainty"
target_journal: "Polymers (MDPI)"
article_type: "Article"
keywords:
  - glass transition temperature
  - polymer informatics
  - structure-property relationships
  - repeat-unit chemistry
  - uncertainty quantification
  - conformal prediction
  - applicability domain
  - extrapolation
---

## Abstract

Machine-learned models of the polymer glass transition temperature (*T*~g~)
routinely report coefficients of determination above 0.85 and are proposed as
screening tools, but those figures come from random partitions whose test sets
are filled with chemistry already well represented in training. Using 7,174
canonical repeat units in twenty-one backbone families, we measure what
deployment outside it costs. For histogram gradient boosting, mean
absolute error rises from 28.1 K on random partitions to 33.8 K under scaffold
holdout, 44.5 K under fingerprint-cluster holdout and 48.1 K under whole-family
holdout, and pooled *R*² falls from 0.87 to 0.75, while the two best models
differ by 0.38 K. A matched-pair design, fixing the test set and varying only
which structures the training pool loses, attributes +14.97 K (95% CI [+11.21,
+19.28]) to removing the family and +0.29 K to removing an equal number of
unrelated structures, the family arm being worse in all twenty-one families for
both models. Uncertainty fails in the same place: split conformal holds 0.904
marginal coverage on random splits while covering 0.650 of the least familiar
structural band, and loses marginal validity under shift. Family-conditioned
calibration cannot repair this, since a withheld family has no calibration
members; similarity-conditioned calibration holds 0.877–0.907 across all four
regimes.

## 1. Introduction

The glass transition temperature is the single most consulted thermal descriptor of
an amorphous polymer. It bounds the service window of a thermoplastic, sets
processing conditions, and governs the mechanical response of a matrix in a
composite. Measuring it by differential scanning calorimetry or dynamic mechanical
analysis is routine but slow, and it requires a synthesised, purified sample, so
predicting *T*~g~ from the structure of the repeat unit has been a target of
polymer theory since group-contribution methods were first formalised
[VanKrevelen2009, Bicerano, Askadskii2003].

Data-driven surrogates have largely displaced those schemes. Curated repositories
such as PoLyInfo [Otsuka] made supervised learning on experimental *T*~g~
practical, and a decade of work has explored the resulting design space: physically
interpretable descriptor sets with tree ensembles and kernel methods [Tao2021],
graph neural networks over the repeat-unit graph [Qiu2023], and chemical language
models pretrained on millions of hypothetical polymers [Kuenneth2023, Xu2023].
Recent community infrastructure has consolidated the datasets and the evaluation
code [Jablonka2025, Xu2025]. The accuracies reported are genuinely useful:
mean absolute errors in the region of 25–30 K against experimental values whose own
inter-source spread can exceed 100 K.

Those numbers describe interpolation. They are obtained on random partitions, in
which a test repeat unit typically has close structural analogues in training, and
they therefore estimate how well a model fills in gaps within chemistry it has
already seen. The motivating application is the reverse. A screening campaign
enumerates candidates precisely because they are not yet in the literature, and the
candidates that matter most are the ones furthest from the training distribution.
The mismatch between the validation regime and the deployment regime is well known
in cheminformatics, where scaffold and time-based splits were introduced for exactly
this reason [Sheridan2013, Wu2018], and it has begun to be examined for polymers
[Babbar2023, Akdogan2026].

Two questions remain open, and this work addresses both.

**Is the family-holdout penalty caused by the missing chemistry?** Withholding an
entire polymer family and observing an error increase is the standard diagnostic,
but it changes two things at once. The test population changes — families differ in
how intrinsically hard they are and in how wide their *T*~g~ distribution is — and
the training set shrinks. A comparison against a random-split baseline therefore
measures the sum of a population effect, a data-volume effect and the effect of
interest. Published treatments of this dataset state the confound explicitly and
leave it unresolved. Resolving it requires holding the test set fixed and matching
the training pools on size.

**Does the uncertainty survive the shift?** A point prediction without a
trustworthy interval is not actionable for screening, and split conformal
prediction [Papadopoulos2002, Vovk2005, Lei2018] is an attractive answer because it
wraps any regressor in intervals with a finite-sample coverage guarantee under
exchangeability. Two difficulties arise here. Exchangeability is precisely what a
family holdout violates. And the guarantee is *marginal*: a predictor can achieve
90% coverage overall while systematically under-covering a subpopulation, and in
this problem that subpopulation is the unfamiliar chemistry the model was deployed
to explore. A recent analysis on a 410-sample simulation-derived set reports this
pattern and characterises conformal intervals as conservative triage indicators
rather than fine-grained screening tools [Akdogan2026]; whether the failure is
repairable, and whether it persists at the scale and heterogeneity of the
experimental record, has not been tested.

### 1.1. Contributions

1. **An independently re-derived, chemically resolved dataset.** We re-curate
   7,208 experimental repeat-unit records to 7,174 canonical structures,
   reproducing the published audit of this collection step for step, and assign
   each structure to one of twenty backbone families or to an explicit
   abstention class using a *backbone-aware* substructure classifier.
   Classifying on the chain path rather than on whole-molecule matches keeps
   poly(alkyl acrylate)s out of the polyester class, which matters because every
   family-holdout conclusion inherits the labels.
2. **A matched-pair design that isolates the family effect.** Three training
   pools — one retaining the family, one with it removed, and one size-matched
   control from which an equal number of unrelated structures was removed —
   predict an identical held-out test set, separating the cost of missing
   chemistry from the cost of missing data.
3. **A conditional-coverage evaluation of conformal uncertainty**, and a
   novelty-conditioned Mondrian taxonomy that restores subgroup validity where a
   family-conditioned taxonomy structurally cannot.
4. **A validated applicability-domain criterion** expressed in units a
   practitioner can act on: the nearest-neighbour Tanimoto similarity below which
   the stated interval should not be believed.
5. **A group-contribution reference point and a taxonomy sensitivity check**,
   so that the machine-learned result is placed against the additive schemes it
   proposes to replace and against a perturbation of its own labels.
6. **A reproducible open implementation**, with leak-free in-fold preprocessing
   and a test suite pinning both the chemistry of the family assignment and the
   finite-sample behaviour of the conformal predictors.

## 2. Materials and Methods

### 2.1. Dataset and curation

We use the experimental homopolymer *T*~g~ collection redistributed with the
POINT² polymer informatics benchmark [Xu2025], which derives from the PoLyInfo
record [Otsuka]. Each entry pairs a repeat unit written as a PSMILES string —
SMILES with two ``*`` atoms marking the polymerisation attachment points — with a
reported *T*~g~.

The raw export is record-level: one repeat unit may appear several times carrying
values from different primary sources. Modelling on records would allow identical
inputs to fall on both sides of a partition and would weight frequently reported
polymers more heavily, so all analysis is performed on a structure-level table.
Repeat units were canonicalised with RDKit and grouped; where a canonical structure
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
Binary Morgan fingerprints (radius 2, 2,048 bits) provide the similarity space in
which structural novelty is measured. Attachment points were retained in both,
since they identify where the chain continues and the environments around them
encode exactly the linkage motifs that set *T*~g~.

Descriptor column filtering — removal of non-finite, constant and exactly
duplicated columns — is implemented as the first step of a scikit-learn pipeline
and is therefore **fitted on the training fold only**. Performing this filtering
once over the full table, as is common, lets test structures influence the feature
definition; the closest prior analysis of this dataset notes this as an unquantified
limitation of its own results.

### 2.4. Partitioning regimes

Four regimes were used, in increasing order of enforced novelty: **random**
partitions (repeated, 80:20); **scaffold** partitions, in which generic
Bemis–Murcko scaffolds [Bemis1996] are assigned whole; **cluster** partitions,
using average-linkage agglomerative clusters in Jaccard space over Morgan bits;
and **family** holdouts, in which every member of one backbone family is withheld.
All non-family regimes were repeated ten times with different seeds, because
single-split rankings among tree ensembles on this dataset fall inside run-to-run
noise. The family regime contributes twenty-one splits, one per family.

### 2.5. Matched-pair design

For each eligible family *F*, a test set *T* ⊂ *F* is drawn (30% of *F*) and three
training pools are constructed: **informed** (everything but *T*, so the remainder
of *F* is available), **naive** (everything but *F* entirely), and **control**
(**informed** minus a randomly chosen block of structures from other families,
equal in number to *F* \ *T*). All three predict the same *T*.

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
glass-transition contribution [VanKrevelen2009, Bicerano, Askadskii2003]. Groups
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
carve-out.

Given calibration scores *s*~1~…*s*~n~ and a miscoverage level α, the conformal
quantile is the ⌈(*n*+1)(1−α)⌉-th order statistic; when that index exceeds *n* the
honest interval is unbounded. We report α = 0.1 throughout. Three predictors are
compared:

* **SCP** — one global quantile of absolute residuals.
* **Normalised SCP** — residuals scaled by a fitted difficulty estimate, so
  intervals widen where the model expects to struggle rather than everywhere.
* **Mondrian SCP** — a separate quantile per category [Vovk2012], under two
  taxonomies: *polymer family*, and *nearest-training-similarity band*.

The asymmetry between the two taxonomies is the point. A family-conditioned
taxonomy is undefined for a family with no calibration members, which is exactly
the family-holdout case; a novelty-conditioned taxonomy is defined for any repeat
unit, seen family or not.

A fourth, **width-matched control**, is used only to test whether the
novelty-conditioned predictor's advantage reduces to its being wider. It
multiplies every split conformal interval about its centre by a single constant
*k*, chosen so that the mean width equals the novelty-conditioned mean width on
the same test set. It is an oracle and not a deployable method, because *k* is
read off the test-set widths it is then compared against; its role is to bound
what *any* global inflation of split conformal could achieve.

Alongside marginal coverage and mean interval width we report **conditional**
coverage within each similarity band and each family, the **maximum coverage
gap** (the largest shortfall below nominal across subgroups with at least ten
members), and the coverage–width criterion. Ensemble spread is additionally
assessed by sparsification error.

### 2.8. Statistical treatment

Confidence intervals on the matched-pair effects are percentile bootstrap
intervals over the twenty-one family-level means. Those twenty-one values are
neither independent nor a sample. The arms of a given family share at least 82%
of their training data by construction, and more than 95% in sixteen of the
twenty-one families; the families are chemically
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
``polymer-tg-benchmark``. Seeds are fixed; a cold run reproduces every number in
this paper.

## 3. Results

### 3.1. The public record survives re-derivation

Independent re-curation of the 7,208 raw records reproduces the published audit of
this collection step for step: 7,174 unique raw SMILES, no RDKit parse failures,
7,174 canonical structures, 31 canonical groups containing more than one record,
34 records absorbed by aggregation, 28 of those groups carrying disagreeing
*T*~g~ values, and a maximum within-structure range of 115 K (Table 1).
Attachment-point
counts likewise agree, at 7,170 structures with two, three with three and one
with four. The dataset is therefore a stable object to build on, and the numbers
below are not sensitive to curation choices.

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

### 3.2. Accuracy is governed by the split, not by the model

On random partitions the best model, histogram gradient boosting, reaches
28.1 ± 1.0 K MAE (*R*² = 0.87), consistent with published values for
descriptor-based ensembles on this dataset once the calibration carve-out of
Section 2.7 is allowed for. That ranking is fragile. Across ten repeated splits,
histogram gradient boosting beats extremely randomised trees by 0.38 K, which a
paired *t*-test does not resolve (*p* = 0.069); random forests and support vector
regression are worse by 2.17 K and 2.74 K (*p* < 0.0001 each). Reporting a single
split, as is common, would allow either of the leading models to be declared the
winner.

Against that sensitivity, the choice of partition moves the error by far more
(Figure 2, Table 3). Holding the model at histogram gradient boosting, MAE rises
from 28.1 K on random splits to 33.8 K on scaffold splits, 44.5 K on
fingerprint-cluster splits and 48.1 K under whole-family holdout. The
random-to-family contrast is therefore 20.0 K: about fifty times the 0.38 K that
separates the two best models, and about seven times the 2.74 K that spans all
four non-trivial regressors. Only the training-median baseline, at +65.7 K, is
moved more by the choice of model than by the choice of protocol. Whatever a
leaderboard on this dataset is ranking, it is not the property that changes most.

Mean nearest-training Tanimoto similarity falls alongside the error, from 0.76 on
random splits to 0.67, 0.40 and 0.45. That ordering is not strictly monotone:
family holdout leaves test structures slightly *less* isolated on average (0.45)
than cluster holdout (0.40), because a withheld family may still have near
analogues in a chemically adjacent one, yet family holdout produces the higher
error. Structural distance alone therefore does not determine difficulty;
*which* chemistry is missing also matters, a point Section 3.3 takes up directly.

*R*² must be quoted with its aggregation stated, because the two natural choices
disagree by more than any effect in this paper (Table 10). Pooling residuals over
all held-out structures in a regime, which is the quantity comparable to a
random-split *R*², gives 0.87 on random splits, 0.81 on scaffold, 0.68 on cluster
and 0.75 under family holdout. Averaging *R*² over the twenty-one individual
family holdouts instead gives 0.26, because each family is scored against its own
*T*~g~ variance and a family with a narrow *T*~g~ spread contributes a large
negative value at a modest absolute error. That 0.26 is a real number but it is
not commensurable with 0.87, and the fall from 0.87 to 0.26 that a naive
comparison suggests is an aggregation artefact. The comparable fall is 0.87 to
0.75.

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
twenty-one families (Table 4; histogram gradient boosting throughout this
section). At the difficult end sit polyphosphazenes (3.02×, 84.8 K MAE),
polydienes (2.83×, 79.6 K) and polysiloxanes (2.48×, 69.7 K). At the easy end,
polyvinyls deteriorate by 1.07×, polyethers by 1.19× and polycarbonates by
1.20×; no family is predicted better when withheld than the average random-split
structure, though polyvinyls come within 8%.

The signed error makes the mechanism visible rather than merely plausible.
Defining bias as mean(*T*~g~ observed − *T*~g~ predicted), so that a negative
value is overprediction, the three families a model has least chemical warrant
for are overpredicted by a wide margin when withheld: polyphosphazenes by
−83.7 K, polysiloxanes by −63.7 K and the Other backbone silane class by
−45.3 K. Set against their MAEs of 84.8, 69.7 and 52.5 K, the error in these
families is dominated by displacement rather than by scatter. This is the
expected failure for an inorganic backbone
withdrawn from a training set dominated by organic condensation polymers: with
nothing to indicate that P=N or Si–O rotation is nearly free, the model assigns
these chains the stiffness of the organic chemistry it does know. The mirror
case exists and points the other way. Polydienes are *under*predicted by
+48.0 K and polyolefins by +37.8 K; a model deprived of its flexible unsaturated
and saturated hydrocarbon chains places them too low, having learned its
temperature scale from stiffer condensation backbones. Bias of either sign is
what an absent linkage type produces, and it is not confined to the extremes:
across the twenty-one families the median ratio of |bias| to MAE is 0.54, so
about half of a typical family-holdout error is systematic displacement.
Its magnitude tracks structural isolation, the absolute bias correlating with
mean nearest-training similarity at *r* = −0.53 (*p* = 0.014, unadjusted).

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
predictor. Split conformal behaves as advertised on random partitions: 0.904
empirical coverage against a nominal 0.90. Under shift the guarantee degrades
with structural novelty — 0.834 on scaffold splits, 0.730 on cluster splits and
0.701 under family holdout (Figure 4, left; Table 8). This is not a defect of the
method
but of its premise: conformal validity requires calibration and test data to be
exchangeable, and a chemistry-aware holdout is constructed precisely to break
that.

Marginal coverage also conceals a failure that is present even when it holds. On
random splits, where overall coverage is a healthy 0.904, coverage in the
least-similar structural band is 0.650 (Figure 4, right; Figure 5). A nominally
90% interval covers two-thirds of the repeat units that a screening campaign would
actually be evaluating. A model can pass every validity check in current practice
and still be wrong about a third of the unfamiliar polymers it is deployed on.

Neither obvious remedy works. The normalised predictor, which scales residuals by
a fitted difficulty estimate, improves the worst band on random splits (0.741 vs
0.650) but fails with everything else under shift (0.825 scaffold, 0.711 cluster,
0.709 family) and is the *worst* of the four under cluster holdout. Conditioning
the quantile on polymer family fails for a structural reason rather than an
empirical one: under family holdout the withheld family has no calibration members
by construction, so the predictor falls back to the global quantile for 100% of
test structures and reproduces split conformal exactly (0.701, identical interval
widths).

Conditioning on nearest-training similarity is the variant that survives, because
that coordinate is defined for every repeat unit whether or not its family was
seen. It holds 0.907, 0.895, 0.886 and 0.877 across the four regimes, never
requiring a fallback under family holdout. On random splits it achieves this
while being *narrower* than split conformal (129.1 K vs 132.3 K), and it widens
only where the chemistry demands it — to 155 K under scaffold, 232 K under
cluster and 208 K under family holdout. That adaptivity is what the
coverage–width criterion rewards (Figure 6): under family holdout, 411 against
48,020 for split conformal.

Whether width alone could have bought the same result is worth testing directly,
because if it could the taxonomy would be doing no work. Under family holdout the
width-matched control of Section 2.7 — split conformal inflated by the single
constant *k* = 1.67 that equalises mean width with the novelty-conditioned
predictor — reaches 0.908 pooled marginal coverage against the
novelty-conditioned predictor's 0.898. Marginally, therefore, a global inflation
does match and slightly exceed it. Conditionally it does not: in the least
similar band (Tanimoto < 0.4) the width-matched control covers 0.849 against
0.890. Recovering that band by global inflation alone requires *k* = 1.93, which
carries the mean interval to 234 K against the novelty-conditioned 203 K, over-
covering the bands that are already adequate in order to reach the one that is
not. Two caveats attach to this control and both weaken the inflation strategy
further. It is an oracle: *k* is computed from the test-set widths it is being
compared against, so neither 1.67 nor 1.93 is knowable at prediction time. And
the factor is regime-dependent, so a constant chosen on one holdout has no claim
on another. The defensible statement is the conditional one: no single global
inflation of split conformal reaches nominal coverage in the least-familiar band
without over-covering everywhere else, and the factor that would be required
cannot be computed when the prediction is made.

### 3.5. An applicability domain, and an honest limit on the repair

Resolving coverage by similarity band under family holdout (Table 5) turns the
result into an operational rule. Split conformal covers 0.636 of structures whose
nearest training analogue lies below Tanimoto 0.4, and 0.771, 0.809 and 0.874 in
the next three bands. Its intervals are essentially constant in
width (121–123 K) across all six bands, while the absolute error across those
bands varies by a factor of 1.83, from 55.5 K in the least-similar band to
30.3 K in the best (Figure 7). Applying one width to a 1.8-fold spread in error is precisely
the failure mode. Novelty-conditioned calibration instead holds 0.890, 0.919,
0.907 and 0.937 in the four least-familiar bands, at widths that scale from
246 K down to 149 K.

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

### 3.6. The penalty follows which structures were removed, not how many

The matched design separates the two effects that a conventional family holdout
confounds (Figure 3, Table 6). Across twenty-one families and three repeats, with
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

The practical reading follows from the locality contrast. At a fixed
training-set size, a pool that excludes the target family predicts it about 15 K
worse than a pool of the same size that retains the family and excludes an
equivalent quantity of unrelated structures instead.
Either the family must be represented in training, or the prediction must
carry an interval that widens honestly when it is not — which is the case for
novelty-conditioned calibration made in Section 3.4.

### 3.7. The gap is not an artefact of the descriptor set

A natural objection is that the generalisation gap reflects a limitation of RDKit
descriptors rather than of the data, and that a representation better suited to
novel chemistry would close it. Repeating the random and family-holdout regimes
with binary Morgan fingerprints (radius 2, 2,048 bits) in place of the descriptor
block tests this directly, since the two representations share almost nothing:
one is a set of physically motivated aggregate quantities, the other a sparse
record of local substructural environments. Both use extremely randomised trees.

The gap survives the substitution and is not narrower (Table 7). Fingerprints are
worse in absolute terms in both regimes — 35.8 K against 28.5 K on random splits
and 62.7 K against 49.6 K under family holdout — and their family-holdout penalty
is larger in absolute terms (+26.9 K against +21.1 K) while being
indistinguishable as a ratio (1.752 against 1.738). Substructural environments
offer no protection here in aggregate.

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

Both levers work, and to a similar degree (Table 9). Against a nominal 0.90 under
family holdout, with the point model identical in all four arms:

| Calibration set | Quantile | Coverage | Worst band | Width | *n* calibration |
|---|---|---|---|---|---|
| Random | Global (SCP) | 0.701 | 0.564 | 122 K | 1,366 |
| Random | Novelty-conditioned | 0.877 | 0.753 | 208 K | 1,366 |
| Family-out | Global (SCP) | 0.855 | 0.753 | 187 K | 3,349 |
| Family-out | Novelty-conditioned | 0.863 | 0.778 | 189 K | 3,349 |

Rebuilding the calibration set lifts plain split conformal from 0.701 to 0.855
without touching the quantile rule, the model or the representation. It does so by
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
the best worst-band coverage of the four arms, 0.778 against 0.753 and 0.753, but
that margin is not resolved: across the twenty-one families the combination beats
the novelty-conditioned quantile alone on worst-band coverage in only 9, with a
mean difference of +0.025 and a Wilcoxon *p* of 0.40. It should be read as no
worse rather than as better.

The practical guidance follows the cost. Novelty-conditioned calibration is free:
it reuses the same fitted model and calibration set, requiring only that
similarity be computed. Family-out calibration requires refitting the model once
per held-out training family — here eight additional fits per split. For most
purposes the novelty-conditioned quantile alone is the better trade.

### 3.9. An additive baseline deteriorates less and predicts worse

The additive schemes that machine-learned models propose to replace make a
different inductive assumption: that *T*~g~ is a composition-weighted sum of
group contributions, with no interaction terms and no dependence on which
linkages the groups sit in. If the family-holdout penalty is really about absent
linkage chemistry, a model that never learns linkage-specific behaviour in the
first place should lose less when a linkage type is withdrawn. It does, and it
pays for that in accuracy everywhere (Table 11).

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
moves by at most 0.005 across the five models and is unchanged at 0.746 for
histogram gradient boosting. The remaining family-dependent results move by about
a kelvin and no conclusion changes sign: family-holdout MAE goes
from 46.8 K to 48.1 K for histogram gradient boosting and from 48.9 K to 49.6 K
for extremely randomised trees, and the matched-design family effect from
+14.61 K to +14.97 K and from +13.28 K to +15.15 K. The one qualitative change
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
and overpredicts those chains by 84 K on average (Section 3.3). Three lines of
evidence support that reading: deterioration correlates with structural isolation
while no size effect is resolved (Section 3.3); the withheld-family errors are
systematic displacements of a consistent sign, about half the error on a typical
family (Section 3.3); and an additive model that never learns linkage-specific
behaviour deteriorates proportionally less when a linkage type is withdrawn
(Section 3.9).

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

Section 3.4 also tested whether width alone explains the advantage. Marginally it
can: an oracle inflation of split conformal matched to the same mean width
reaches 0.908 against 0.898. Conditionally it cannot, missing the least-similar
band at 0.849 against 0.890, and the global factor that would close that band
costs 234 K of mean width against 203 K and is not knowable at prediction time.
The claim we make is therefore about allocation, not about total width: the
novelty-conditioned predictor puts its width where the error is, and no constant
multiple of split conformal does that.

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

Finally, the novelty-conditioned predictor under-covers the most-similar band
under family holdout (Section 3.5) and falls back to the global quantile for
12.8% of structures under cluster holdout (Section 4.3). The affected groups are
small and the direction of the trade favours screening, but the method is an
improvement rather than a solution, and conditional coverage should continue to
be reported rather than assumed.

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
covering only 65% of the least familiar structures, and lose marginal validity
under shift. Conditioning the calibration on polymer family cannot repair this,
because a withheld family has no calibration data by construction. Conditioning
on nearest-training structural similarity does, holding 0.88–0.91 across every
regime tested and allocating width to the bands where the error actually is.
Rebuilding the calibration set from within-training family holdouts is a second
lever that reaches 0.855 with the plain global quantile, though in this
experiment it also enlarges the calibration set from 1,366 to 3,349 structures,
so construction and size cannot be separated; combining the two levers is no better than the novelty-conditioned
quantile alone.

The practical implication is a short list: report a chemistry-aware split, repeat
what can be repeated, resolve errors by family and by sign, and attach
novelty-conditioned intervals together with the similarity value itself. None of
these is expensive, and together they turn a *T*~g~ model from something that
scores well into something a screening campaign can act on.

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
redistributed with the POINT² benchmark [Xu2025] and derive from PoLyInfo
[Otsuka].

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
