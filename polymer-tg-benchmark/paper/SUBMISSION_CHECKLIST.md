# Pre-submission checklist

Author-facing notes. This file is not part of the manuscript and must not be
submitted with it. `paper/references.md` is the bibliography and is kept clean of
these annotations.

## Bibliographic entries still to be verified

The following entries were cited from established knowledge and were **not**
checked against the publisher or an indexing service; outbound access to
Crossref and to GitHub was blocked in the environment the revision was prepared
in. Confirm volume, issue, page range, year and DOI for each before submission.

| # | Entry | What to confirm |
|---|---|---|
| 9 | Bicerano, *Prediction of Polymer Properties*, 3rd ed. | Publisher and year of the edition cited (Marcel Dekker, 2002). |
| 10 | Van Krevelen & te Nijenhuis, *Properties of Polymers*, 4th ed. | Publisher and year (Elsevier, 2009). |
| 11 | Askadskii, *Computational Materials Science of Polymers* | Publisher and year (Cambridge International Science Publishing, 2003). |
| 12 | Otsuka et al., PoLyInfo | Conference name, publisher and page range (pp 22–29). |
| 15 | Norinder, Carlsson, Boyer, Eklund (2014) | Volume 54, issue 6, pages 1596–1603, DOI 10.1021/ci5001168. Added in this revision; transcribed, not verified. |
| 16 | Papadopoulos et al., ECML 2002 | LNCS volume and page range (2430, pp 345–356). |
| 17 | Vovk, Gammerman & Shafer (2005) | Publisher and year. |
| 18 | Vovk (2013) | Volume 92, pages 349–376. |
| 19 | Lei et al. (2018) | Volume 113, issue 523, pages 1094–1111. |
| 20 | Papadopoulos & Haralambous (2011) | Volume 24, issue 8, pages 842–851. |
| 21 | Boström & Johansson (2020) | PMLR volume 128, pages 114–133. |
| 22 | Barber et al. (2021) | Volume 49, issue 1, pages 486–507. |
| 23 | Tibshirani, Barber, Candès & Ramdas (2019) | NeurIPS 32 proceedings pagination and the arXiv identifier 1904.06019. Added in this revision; transcribed, not verified. |
| 24 | Bemis & Murcko (1996) | Volume 39, issue 15, pages 2887–2893. |
| 25 | Rogers & Hahn (2010) | Volume 50, issue 5, pages 742–754. |
| 26 | Sheridan (2013) | Volume 53, issue 4, pages 783–790. |
| 27 | Wu et al. (2018) | Volume 9, pages 513–530, and the full author list. |
| 28 | Pedregosa et al. (2011) | Volume 12, pages 2825–2830. |
| 29 | RDKit | Version used and a citable release. |
| 30 | Teh | See below. |

Entries 1–8, 13 and 14 were checked against the publisher or an indexing service
during earlier work on this manuscript.

## The [Teh] entry

1. **Archive it.** A GitHub repository is not a published source and its state
   can change. Deposit a snapshot with Zenodo (DOI) or record a Software
   Heritage identifier, and cite that rather than the bare repository URL.
2. **Verify the author name and spelling** from the repository itself or from
   the author directly. The name in the entry was taken from the repository
   handle and has not been confirmed.
3. **Verify the reported values transcribed into Table 1.** The left-hand column
   of Table 1 (7,208 records; 7,174 unique SMILES; 0 parse failures; 7,174
   canonical structures; 31 repeated canonical groups; 34 additional occurrences
   consolidated; 28 groups with differing *T*~g~; 115 °C maximum within-group
   range; record-level attachment-point counts 7,204 / 3 / 1) is transcribed
   from that analysis. It is **not** an output of this pipeline and is not
   checked by `scripts/check_manuscript_tables.py`. Re-read it off the archived
   snapshot before submission.

## Other pre-submission items

- Fill in the `[AUTHOR]` placeholders in the Author Contributions section.
- Confirm the accessed-on date for entry 30 once the archived version is cited.
