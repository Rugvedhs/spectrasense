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
