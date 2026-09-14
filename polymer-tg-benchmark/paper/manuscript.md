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

Machine-learned models of the glass transition temperature (*T*~g~) of homopolymers
routinely report coefficients of determination above 0.85, and are increasingly
proposed as screening tools for polymer discovery. Those figures are obtained on
random partitions of the available data, which populate the test set with repeat
units whose backbone chemistry is already well represented in training. Screening
campaigns operate in the opposite regime. Here we quantify what that difference
costs, attribute it, and ask whether the uncertainty attached to a prediction
survives it. Using 7,174 canonical repeat units spanning twenty backbone families,
we compare random, Bemis–Murcko scaffold, fingerprint-cluster and whole-family
partitions, and introduce a matched-pair design in which three models predict an
identical set of held-out structures from training pools that differ only in what
was removed. Because a size-matched control arm is included, the penalty for
removing a family is separated from the penalty for training on less data — a
confound that previous family-holdout analyses of this dataset acknowledge but do
not resolve. We then evaluate split conformal prediction and two conditional
variants. The central finding is that reliability in this problem is a property of
backbone chemistry rather than of the model: siloxane, phosphazene and imide
backbones behave very differently from one another under identical treatment, and
conventional conformal intervals retain their nominal marginal coverage while
covering substantially less than nominal for the least familiar repeat units.
Conditioning the conformal calibration on polymer family cannot repair this,
because an unseen family has no calibration data by construction; conditioning on
nearest-neighbour structural similarity can, because that coordinate is defined
for every repeat unit. We recommend reporting family-resolved errors and
novelty-conditioned intervals alongside any *T*~g~ model intended for screening,
and provide an open, fully reproducible implementation.

## 1. Introduction

The glass transition temperature is the single most consulted thermal descriptor of
an amorphous polymer. It bounds the service window of a thermoplastic, sets
processing conditions, and governs the mechanical response of a matrix in a
composite. Measuring it by differential scanning calorimetry or dynamic mechanical
analysis is routine but slow, and it requires a synthesised, purified sample, so
predicting *T*~g~ from the structure of the repeat unit has been a target of
polymer theory since group-contribution methods were first formalised [Bicerano].

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
   each structure to one of twenty backbone families using a *backbone-aware*
   substructure classifier. Classifying on the chain path rather than on
   whole-molecule matches keeps poly(alkyl acrylate)s out of the polyester class,
   which matters because every family-holdout conclusion inherits the labels.
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
5. **A reproducible open implementation**, with leak-free in-fold preprocessing
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

Backbone restriction is what separates a poly(alkyl acrylate) from a polyester:
both contain an ester group, but in the acrylate it decorates an all-carbon chain.
Patterns are tested in order of specificity, because the groups are nested — every
cyclic imide contains two amide-like motifs, and a urethane contains both
ester-like and amide-like fragments. Repeat units whose backbone carries no
heteroatom linkage are classified by the groups pendant to the chain
(aromatic-dominated, acrylic, styrenic, halogenated, other heteroatom-substituted
vinyl, unsaturated, or saturated hydrocarbon). The twenty-two chemistry cases in
Table S1 are pinned by unit tests.

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
noise.

### 2.5. Matched-pair design

For each eligible family *F*, a test set *T* ⊂ *F* is drawn and three training
pools are constructed: **informed** (everything but *T*, so the remainder of *F*
is available), **naive** (everything but *F* entirely), and **control**
(**informed** minus a randomly chosen block of structures from other families,
equal in number to *F* \ *T*). All three predict the same *T*.

The contrast **naive** − **control** is the quantity of interest. Both pools are
identical in size and both exclude *T*; they differ only in *which* structures
were removed. The contrast **control** − **informed** measures the pure
data-volume effect over the same test set.

### 2.6. Models

Five regressors were evaluated on the descriptor representation: a training-median
baseline, support vector regression, random forest, histogram gradient boosting and
extremely randomised trees. Every model is wrapped in the same pipeline so that all
data-dependent preprocessing is fitted in-fold.

### 2.7. Conformal prediction and coverage diagnostics

Each training pool is divided three ways: a **fit** partition trains the point
regressor, a **difficulty** partition supplies out-of-sample residuals, and a
**calibration** partition supplies the conformal quantile. Fitting the difficulty
model on the point model's own training residuals would teach it that the model is
everywhere accurate, so that partition must be held out. All conformal variants
share one point predictor and one calibration set, so differences between them
reflect only how residuals are converted into intervals.

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

Alongside marginal coverage and mean interval width we report **conditional**
coverage within each similarity band and each family, the **maximum coverage
gap** (the largest shortfall below nominal across subgroups with at least ten
members), and the coverage–width criterion. Ensemble spread is additionally
assessed by sparsification error.

### 2.8. Reproducibility

All code, the curated structure-level table and every result file are available at
[repository]. Seeds are fixed; a cold run reproduces every number in this paper.

## 3. Results

### 3.1. The public record survives re-derivation

Independent re-curation of the 7,208 raw records reproduces the published audit of
this collection step for step: 7,174 unique raw SMILES, no RDKit parse failures,
7,174 canonical structures, 31 canonical groups containing more than one record,
34 records absorbed by aggregation, 28 of those groups carrying disagreeing
*T*~g~ values, and a maximum within-structure range of 115 K. Attachment-point
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
family a polymer chemist would assign. The resulting family *T*~g~ ordering
(Figure 1) is likewise the expected one, which is a further check that the labels
are chemically meaningful rather than merely self-consistent: polyphosphazenes
(median −8 °C) and polysiloxanes (−8 °C) sit lowest, reflecting the low rotational
barriers of P=N and Si–O backbones, and aromatic polyimides sit highest (249 °C),
reflecting rigid, strongly interacting heteroaromatic chains.

### 3.2. Accuracy is governed by the split, not by the model

On random partitions the best model reaches 28.1 ± 1.0 K MAE
(*R*² = 0.87), consistent with published values for descriptor-based ensembles on
this dataset. That ranking is fragile. Across ten repeated splits, histogram
gradient boosting beats extremely randomised trees by 0.38 K, which a paired
*t*-test does not resolve (*p* = 0.069); random forests and support vector
regression are worse by 2.17 K and 2.74 K (*p* < 0.0001 each). Reporting a single
split, as is common, would allow either of the leading models to be declared the
winner.

Against that 0.38 K of model-choice sensitivity, the choice of partition moves the
error by an order of magnitude more (Figure 2, Table 3). Holding the model fixed,
MAE rises from 28.1 K on random splits to 33.8 K on scaffold splits, 44.5 K on
fingerprint-cluster splits and 46.8 K under whole-family holdout; *R*² falls from
0.87 to 0.79, 0.54 and 0.32. Mean nearest-training Tanimoto similarity falls
alongside it, from 0.76 on random splits to 0.67, 0.40 and 0.46. That ordering is
not strictly monotone: family holdout leaves test structures slightly *less*
isolated on average (0.46) than cluster holdout (0.40), because a withheld family
may still have near analogues in a chemically adjacent one, yet family holdout
produces the higher error. Structural distance alone therefore does not determine
difficulty; *which* chemistry is missing also matters, a point Section 3.3 takes
up directly. **The evaluation protocol is worth about forty times more than the
model choice** (16.3 K against 0.38 K), which reframes what a leaderboard on this
dataset is measuring.

The dispersion is as informative as the mean. Random-split MAE varies by ±1.0 K
across repeats; scaffold and cluster splits vary by ±7.8 K and ±10.4 K, and
family holdouts by ±15.6 K. Performance under novelty is not a single number, and
a point estimate from one chemistry-aware split should not be trusted.

### 3.3. Transfer failure is a property of backbone chemistry

Family-holdout deterioration, defined as the ratio of family-holdout MAE to
random-split MAE for the same model, spans a factor of nearly four across the
twenty families (Table 4). At the difficult end sit polyphosphazenes (3.02×,
84.8 K MAE), polydienes (2.58×) and polysiloxanes (2.48×). At the easy end,
polycarbonates deteriorate by 0.80× — that is, they are predicted *better* when
their own family is withheld than the average random-split structure is.

The pattern is chemically coherent. Polyphosphazene and polysiloxane backbones are
inorganic, and nothing else in a dataset dominated by organic condensation
polymers constrains the relationship between their descriptors and chain
flexibility; the models fall back on organic intuition and overpredict badly.
Polycarbonates, by contrast, are surrounded in structure space by polyesters and
polyethers, which remain in training and carry nearly the same carbonyl and ether
linkages, so their withdrawal removes little the model cannot recover.

Testing three candidate explanations against deterioration across the twenty
families, only structural isolation survives. Mean nearest-training similarity
correlates strongly and significantly for every model (Pearson *r* between −0.63
and −0.68, *p* < 0.005). The shift between the family's median *T*~g~ and the rest
of the dataset is weaker and marginal (*r* ≈ 0.40–0.46, *p* ≈ 0.04–0.08), and
family size does not predict deterioration at all (*p* > 0.2). Larger families are
not easier to extrapolate to; isolated ones are harder.

This also settles a practical question. Because the effect tracks structural
isolation rather than family size, collecting more data *within* well-represented
families will not close the gap. Coverage of the structure space, not volume, is
the binding constraint.

### 3.4. Conformal intervals lose validity exactly where they are needed

Split conformal behaves as advertised on random partitions: 0.904 empirical
coverage against a nominal 0.90. Under shift the guarantee degrades monotonically
with structural novelty — 0.834 on scaffold splits, 0.730 on cluster splits and
0.711 under family holdout (Figure 4, left). This is not a defect of the method
but of its premise: conformal validity requires calibration and test data to be
exchangeable, and a chemistry-aware holdout is constructed precisely to break
that.

Marginal coverage also conceals a failure that is present even when it holds. On
random splits, where overall coverage is a healthy 0.904, coverage in the
least-similar structural band is 0.650 (Figure 4, right; Figure 5). A nominally
90% interval covers two-thirds of the repeat units that a screening campaign would
actually be evaluating. **A model can pass every validity check in current
practice and still be wrong about a third of the unfamiliar polymers it is
deployed on.**

Neither obvious remedy works. The normalised predictor, which scales residuals by
a fitted difficulty estimate, improves the worst band on random splits (0.741 vs
0.650) but fails with everything else under shift (0.826 scaffold, 0.711 cluster,
0.712 family) and is the *worst* of the four under cluster holdout. Conditioning
the quantile on polymer family fails for a structural reason rather than an
empirical one: under family holdout the withheld family has no calibration members
by construction, so the predictor falls back to the global quantile for 100% of
test structures and reproduces split conformal exactly (0.711, identical interval
widths).

Conditioning on nearest-training similarity is the variant that survives, because
that coordinate is defined for every repeat unit whether or not its family was
seen. It holds 0.907, 0.895, 0.886 and 0.883 across the four regimes, never
requiring a fallback under family holdout. Crucially it does not buy this with
width: on random splits it is *narrower* than split conformal (129.1 K vs
132.3 K) while covering far more of the difficult band, and it widens only where
the chemistry demands it — to 155 K under scaffold, 232 K under cluster and 208 K
under family holdout. That adaptivity is the behaviour a screening workflow needs,
and it is what the coverage–width criterion rewards (Figure 6): under family
holdout, 342 against 35,556 for split conformal.

### 3.5. An applicability domain, and an honest limit on the repair

Resolving coverage by similarity band under family holdout (Table 5) turns the
result into an operational rule. Split conformal covers 0.650 of structures whose
nearest training analogue lies below Tanimoto 0.4, rising through 0.778, 0.828 and
0.868 as similarity increases. Its intervals are essentially constant in width
(123–127 K) across all six bands — it applies one width to structures whose error
differs by a factor of three, which is precisely the failure mode.
Novelty-conditioned calibration instead holds 0.900, 0.931, 0.903 and 0.906 in the
four least-familiar bands, at widths that scale from 250 K down to 158 K.

The repair has a cost that should be stated plainly. In the most-similar band
(Tanimoto ≥ 0.8) the novelty-conditioned predictor *under*-covers, at 0.622
against split conformal's 0.803, because it narrows those intervals to 86 K. Under
a family holdout a test structure that closely resembles training data is unusual,
and the calibration points populating that band come from families where high
similarity genuinely did signal an easy prediction. The trade is nonetheless
strongly favourable in practice: that band holds 127 of 7,174 held-out predictions
(1.8%), while the two least-similar bands hold 4,759 (66%). The recommendation is
therefore to report novelty-conditioned intervals together with the
nearest-training similarity itself, so that a reader can see which regime a given
prediction sits in rather than trusting a single interval uniformly.

### 3.6. The penalty is caused by the missing chemistry, not the missing data

The matched design separates the two effects that a conventional family holdout
confounds (Figure 3, Table 6). Across twenty families and three repeats, with the
test set held identical and the two comparison pools held identical in size:

* Removing an equal quantity of *unrelated* structures — the **control** arm —
  changes MAE by +0.32 K (95% CI [−0.63, +1.23], Wilcoxon *p* = 0.47, worse in
  12 of 20 families) for extremely randomised trees, and by +0.21 K
  (*p* = 0.93, 9 of 20) for histogram gradient boosting. Per family the effect
  ranges from −4.0 K to +5.1 K and scatters around zero.
* Removing *the family itself* — the **naive** arm — costs +14.61 K
  (95% CI [+10.53, +19.41], *p* = 2 × 10⁻⁶) and **is worse in all twenty of
  twenty families**, and +13.28 K (*p* = 1 × 10⁻⁵, 18 of 20) for the second model.

The family effect is roughly forty-six times the size of the data-volume effect
and, unlike it, is consistent in sign across every family tested. Essentially the
entire family-holdout penalty is attributable to the absent chemistry: the total
effect of +14.93 K decomposes into +14.61 K of chemistry and +0.32 K of data
volume. This resolves the confound that motivated the design, and it does so on
matched test structures, so it is not an artefact of families differing in
intrinsic difficulty or in *T*~g~ spread.

Family-level magnitudes rank as the chemistry predicts. Polyphosphazenes suffer
most (+45.4 K), followed by polyamides (+30.0 K) and polysiloxanes (+28.7 K); at
the other end polysulfides (+1.5 K), polyvinyls (+3.7 K) and polyimines (+4.1 K)
are nearly unaffected, because the linkages that set their backbone mobility
survive elsewhere in the training set. Note that polyamides rank second here
despite being the second-largest family in the dataset (982 structures) — further
evidence that abundance does not confer transferability when the chemistry is
distinctive.

The practical reading is uncomfortable but clear. Because the penalty comes from
absent chemistry rather than absent volume, a *T*~g~ model cannot be made
trustworthy on a new backbone family by training it on more of what it already
has. Either the family must be represented in training, or the prediction must
carry an interval that widens honestly when it is not — which is the case for
novelty-conditioned calibration made in Section 3.4.
