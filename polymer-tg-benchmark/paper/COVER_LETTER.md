# Cover letter

To the Editorial Office, *Polymers* (MDPI)

Dear Editor,

We submit "Backbone Chemistry Governs the Reliability of Machine-Learned Glass
Transition Temperatures: A Family-Resolved Benchmark with Novelty-Conditioned
Conformal Uncertainty" for consideration as an Article.

The question the paper asks is a polymer-science one: for which backbone
chemistries can a structure–property model of *T*~g~ be trusted, and for which
not. Working from 7,174 experimental repeat units assigned by a backbone-aware
classifier to twenty families and an explicit abstention class, we withhold each
family in turn and report the error resolved by chemistry and by sign. The
structure of the failure is chemical rather than statistical. Withheld
polyphosphazenes and polysiloxanes are overpredicted, by 83.7 K and 63.7 K on
average, because a model trained on stiffer organic backbones has no way to learn
what the low rotational barriers of P=N and Si–O chains do to segmental mobility;
withheld polydienes and polyolefins are underpredicted, by 48.0 K and 37.8 K.
Deterioration tracks how structurally isolated a family is, not how large it is.
A matched-pair design, in which the test set and the training-set size are held
fixed and only the identity of the removed structures changes, attributes about
15 K to the missing chemistry and under 1 K, of unresolved sign, to the missing
data. A Van Krevelen group-contribution baseline refitted on the same folds
deteriorates proportionally less and predicts worse everywhere (43.8 K against
28.1 K on random partitions), which is the same mechanism seen from a model class
with the opposite inductive bias. We also show that a prediction interval
conditioned on structural novelty remains usable for a backbone the model has
never seen, where one conditioned on polymer family cannot be.

We believe this fits the journal's scope in polymer structure–property
relationships and computational polymer science: the subject is the repeat-unit
chemistry, and the machine learning is the instrument.

The underlying experimental *T*~g~ records are public, redistributed with the
POINT² benchmark and derived from PoLyInfo. All code, the curated
structure-level table and every result file are released openly, and a cold run
reproduces every number in the paper.

No preprint of this work has been posted. The manuscript is not under
consideration elsewhere, all authors have approved the submission, and we declare
no conflict of interest.

Yours sincerely,

[AUTHOR], on behalf of [AUTHOR(S)]
[AFFILIATION]
[CORRESPONDING EMAIL]

---

## Suggested reviewers

**The author must verify every name, affiliation and email below before
submission.** They are drawn from the reference list of this manuscript and the
affiliations are stated as the author's best recollection, not as a checked fact;
outbound access to Crossref and to publisher pages was blocked in the environment
this letter was drafted in. Confirm each from the cited paper's title page and
supply a current institutional email.

| Suggested reviewer | Basis in our reference list | Affiliation to confirm |
|---|---|---|
| Rampi Ramprasad | Ref. 8, polyBERT, chemical language models for polymer informatics | Georgia Institute of Technology (confirm) |
| Kevin M. Jablonka | Ref. 10, PolyMetriX, polymer informatics infrastructure and evaluation | Friedrich Schiller University Jena (confirm) |
| Amir Barati Farimani | Ref. 9, TransPolymer, transformer models for polymer properties | Carnegie Mellon University (confirm) |
| Ying Li | Ref. 5, benchmarking machine learning models for polymer *T*~g~ | University of Wisconsin–Madison (confirm) |
| M. Akdoğan | Ref. 15, conformal uncertainty for polymer *T*~g~; the closest prior work | Confirm from the *ACS Omega* paper |

None of the suggested reviewers is a collaborator or co-author of ours. If the
editorial office prefers reviewers without a stake in the methods being
criticised, note that Refs. 5, 8, 9 and 15 are among the works whose evaluation
protocol this paper argues is incomplete.
