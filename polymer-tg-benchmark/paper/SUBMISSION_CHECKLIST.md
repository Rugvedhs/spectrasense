# Pre-submission checklist — *Polymers* (MDPI)

Author-facing notes. This file is not part of the manuscript and must not be
submitted with it. `paper/references.md` is the bibliography and is kept clean of
these annotations.

Reference numbers below are the MDPI numeric ones produced by
`scripts/renumber_citations.py`, that is, the order of first appearance in the
text. They changed in this revision; the old author-year keys no longer exist.

---

## A. Blocking items — only the author can supply these

In priority order. Nothing below can be inferred from the repository.

1. **Author names**, in the intended order, with the spelling to appear in
   print.
2. **Affiliations** for every author, in MDPI's numbered form (department,
   institution, city, postcode, country).
3. **ORCID iDs** for every author who has one. MDPI displays them and the
   submission system asks for them at author-registration time.
4. **Corresponding author**: which author, with an institutional email address
   that will remain valid after publication, and a postal address.
5. **CRediT contribution statement.** The Author Contributions section of the
   manuscript currently carries ten `[AUTHOR]` placeholders, one per CRediT
   role. Replace each with initials, and delete any role that no author
   performed rather than assigning it to someone who did not do it. Roles used:
   conceptualisation, methodology, software, validation, formal analysis,
   investigation, data curation, writing—original draft, writing—review and
   editing, visualisation.
6. **Funding status.** The manuscript states "This research received no external
   funding." Confirm that this is true for every author, including institutional
   studentships, fellowships and equipment grants, and replace it with the
   funder name and grant number if it is not.
7. **Article processing charge.** The APC for *Polymers* is CHF 2,700, payable
   on acceptance. Confirm the current rate on the journal's own page before
   submitting, as it is revised periodically.
   - **A waiver or discount must be requested before submission and obtained in
     writing.** MDPI does not reliably grant waivers after a manuscript is
     accepted, and acceptance without a confirmed waiver leaves the full charge
     due. Keep the emailed confirmation, and quote its reference in the
     submission comments.
8. **Dataset redistribution licence.** The experimental *T*~g~ records come from
   PoLyInfo [4] by way of the POINT² benchmark [11], and this work redistributes
   a curated structure-level derivative of them in the repository named in the
   Data Availability Statement. Confirm that POINT²'s licence permits that
   redistribution and that PoLyInfo's terms permit POINT²'s. If either does not,
   the repository must ship the derivation code and the record identifiers
   rather than the values, and the Data Availability Statement must be reworded
   to match.
9. **Preprint position.** The cover letter states that no preprint has been
   posted. If one is posted before or during review, tell the editorial office
   and give the DOI.

---

## B. Bibliographic entries still to be verified

The following entries were cited from established knowledge and were **not**
checked against the publisher or an indexing service; outbound access to
Crossref and to GitHub was blocked in the environment the revision was prepared
in. Confirm volume, issue, page range, year and DOI for each before submission.

| # | Entry | What to confirm |
|---|---|---|
| 1 | Van Krevelen & te Nijenhuis, *Properties of Polymers*, 4th ed. | Publisher and year (Elsevier, 2009). |
| 2 | Bicerano, *Prediction of Polymer Properties*, 3rd ed. | Publisher and year of the edition cited (Marcel Dekker, 2002). |
| 3 | Askadskii, *Computational Materials Science of Polymers* | Publisher and year (Cambridge International Science Publishing, 2003). |
| 4 | Otsuka et al., PoLyInfo | Conference name, publisher and page range (pp 22–29). |
| 12 | Sheridan (2013) | Volume 53, issue 4, pages 783–790. |
| 13 | Wu et al. (2018) | Volume 9, pages 513–530, and the full author list. |
| 16 | Teh | See Section C. |
| 17 | Papadopoulos et al., ECML 2002 | LNCS volume and page range (2430, pp 345–356). |
| 18 | Vovk, Gammerman & Shafer (2005) | Publisher and year. |
| 19 | Lei et al. (2018) | Volume 113, issue 523, pages 1094–1111. |
| 21 | Norinder, Carlsson, Boyer & Eklund (2014) | Volume 54, issue 6, pages 1596–1603, DOI 10.1021/ci5001168. Transcribed, not verified. |
| 22 | RDKit | Version used and a citable release. |
| 23 | Rogers & Hahn (2010) | Volume 50, issue 5, pages 742–754. |
| 24 | Pedregosa et al. (2011) | Volume 12, pages 2825–2830. |
| 25 | Bemis & Murcko (1996) | Volume 39, issue 15, pages 2887–2893. |
| 26 | Barber et al. (2021) | Volume 49, issue 1, pages 486–507. |
| 27 | Papadopoulos & Haralambous (2011) | Volume 24, issue 8, pages 842–851. |
| 28 | Vovk (2013) | Volume 92, pages 349–376. |
| 29 | Boström & Johansson (2020) | PMLR volume 128, pages 114–133. |
| 30 | Tibshirani, Barber, Candès & Ramdas (2019) | NeurIPS 32 proceedings pagination and the arXiv identifier 1904.06019. Transcribed, not verified. |

Entries 5, 6, 7, 8, 9, 10, 11, 14, 15 and 20 were checked against the publisher
or an indexing service during earlier work on this manuscript.

---

## C. The prior-analysis entry, reference [16]

1. **Archive it.** A GitHub repository is not a published source and its state
   can change. Deposit a snapshot with Zenodo (DOI) or record a Software
   Heritage identifier, and cite that rather than the bare repository URL.
2. **Verify the author name and spelling** from the repository itself or from
   the author directly. The name in the entry was taken from the repository
   handle and has not been confirmed.
3. **Verify the reported values transcribed into Table 1.** The middle column of
   Table 1 (7,208 records; 7,174 unique SMILES; 0 parse failures; 7,174
   canonical structures; 31 repeated canonical groups; 34 additional occurrences
   consolidated; 28 groups with differing *T*~g~; 115 K maximum within-group
   range; record-level attachment-point counts 7,204 / 3 / 1) is transcribed
   from that analysis. It is **not** an output of this pipeline and is not
   checked by `scripts/check_manuscript_tables.py`. Re-read it off the archived
   snapshot before submission.
4. **Confirm the accessed-on date** once the archived version is cited.

---

## D. Suggested reviewers

`paper/COVER_LETTER.md` proposes five reviewers drawn from the reference list.
Every name, affiliation and email there is unverified and must be confirmed from
the cited paper's title page before the cover letter is sent.

---

## E. Mechanical checks that are already automated

These pass in the repository and do not need author input, but should be re-run
after any edit to the manuscript:

- `scripts/check_manuscript_tables.py` — all fourteen inline tables against
  their source CSVs.
- `scripts/renumber_citations.py --check` — every citation resolves and every
  reference is cited.
- `pytest tests -q` — the above plus the chemistry and conformal unit tests.
- `scripts/build_pdf.py` — typesets the current manuscript source.


## Why each annotated reference is cited

Moved out of `references.md`, which a submitted manuscript needs clean.

- **[1]** the additive form *T*~g~ = *Y*~g~/*M* used for the group-contribution baseline in Section 2.6.
- **[3]** an alternative additive scheme based on repeat-unit volume increments.
- **[11]** source of the *T*~g~ collection used here.
- **[13]** scaffold splits.
- **[15]** the closest prior work; 410 simulation-derived samples, marginal split conformal only.
- **[16]** a prior public analysis of this dataset: it reports the 7,208 → 7,174 curation audit, the 31/34/28 aggregation counts and the 115 K maximum within-structure range reproduced in Section 3.1, and states the family-holdout confound that Section 2.5 sets out to resolve.
- **[21]** conformal prediction used as the applicability-domain statement itself.
- **[26]** the data-reusing alternative to the split construction, noted in Section 2.7.
- **[27]** normalised nonconformity.
- **[28]** Mondrian conformal prediction.
- **[30]** weighted conformal prediction; discussed in Sections 2.7 and 4.3.
