# References

`[verified]` marks entries whose bibliographic details were checked against the
publisher or an indexing service during this work. Unmarked entries are canonical
sources cited from established knowledge; **check volume, page and year before
submission** — they were not re-verified here.

## Polymer Tg prediction and polymer informatics

1. **[verified]** Tao, L.; Varshney, V.; Li, Y. Benchmarking Machine Learning
   Models for Polymer Informatics: An Example of Glass Transition Temperature.
   *J. Chem. Inf. Model.* **2021**, *61* (11), 5395–5413.
   DOI: 10.1021/acs.jcim.1c01031
2. **[verified]** Casanola-Martin, G. M.; Karuth, A.; Pham-The, H.; et al. Machine
   learning analysis of a large set of homopolymers to predict glass transition
   temperatures. *Commun. Chem.* **2024**, *7*, 226.
   DOI: 10.1038/s42004-024-01305-0
3. **[verified]** Kuenneth, C.; Ramprasad, R. polyBERT: a chemical language model
   to enable fully machine-driven ultrafast polymer informatics.
   *Nat. Commun.* **2023**, *14*, 4099. DOI: 10.1038/s41467-023-39868-6
4. **[verified]** Xu, C.; Wang, Y.; Barati Farimani, A. TransPolymer: a
   Transformer-based language model for polymer property predictions.
   *npj Comput. Mater.* **2023**, *9*, 64. DOI: 10.1038/s41524-023-01016-5
5. **[verified]** Kunchapu, S.; Jablonka, K. M. PolyMetriX: an ecosystem for
   digital polymer chemistry. *npj Comput. Mater.* **2025**, *11*, 312.
6. **[verified]** Xu, J.; et al. POINT²: A Polymer Informatics Training and
   Testing Database. *arXiv* **2025**, arXiv:2503.23491.
   — source of the *T*~g~ collection used here.
7. **[verified]** Babbar, A.; et al. Explainability and Transferability of Machine
   Learning Models for Predicting the Glass Transition Temperature of Polymers.
   *arXiv* **2023**, arXiv:2308.09898.
8. **[verified]** Qiu, H.; et al. Design of Polyimides with Targeted Glass
   Transition Temperature Using a Graph Neural Network.
   *J. Mater. Chem. C* **2023**. DOI: 10.1039/D2TC05174E
9. Bicerano, J. *Prediction of Polymer Properties*, 3rd ed.; Marcel Dekker:
   New York, **2002**.
10. Van Krevelen, D. W.; te Nijenhuis, K. *Properties of Polymers: Their
    Correlation with Chemical Structure; Their Numerical Estimation and
    Prediction from Additive Group Contributions*, 4th ed.; Elsevier:
    Amsterdam, **2009**. — the additive form *T*~g~ = *Y*~g~/*M* used for the
    group-contribution baseline in Section 2.6.
11. Askadskii, A. A. *Computational Materials Science of Polymers*; Cambridge
    International Science Publishing: Cambridge, **2003**. — an alternative
    additive scheme based on repeat-unit volume increments.
12. Otsuka, S.; Kuwajima, I.; Hosoya, J.; Xu, Y.; Yamazaki, M. PoLyInfo: Polymer
    Database for Polymeric Materials Design. In *2011 International Conference on
    Emerging Intelligent Data and Web Technologies*; IEEE, **2011**; pp 22–29.

## Uncertainty quantification and applicability domain

13. **[verified]** Akdoğan, M. Trust Beyond Accuracy: Conformal Uncertainty
    Quantification Reveals the Generalization Gap in Polymer Glass Transition
    Temperature Prediction. *ACS Omega* **2026**, *11* (30), 44926–44941.
    — the closest prior work; 410 simulation-derived samples, marginal split
    conformal only.
14. **[verified]** Agrawal, V. Polymer property prediction using ensemble of
    directed message-passing neural networks with uncertainty estimation and
    applicability domain analysis. *MRS Commun.* **2026**.
    DOI: 10.1557/s43579-026-00948-5
15. Papadopoulos, H.; Proedrou, K.; Vovk, V.; Gammerman, A. Inductive Confidence
    Machines for Regression. In *Machine Learning: ECML 2002*; LNCS 2430;
    Springer, **2002**; pp 345–356.
16. Vovk, V.; Gammerman, A.; Shafer, G. *Algorithmic Learning in a Random World*;
    Springer: New York, **2005**.
17. Vovk, V. Conditional validity of inductive conformal predictors.
    *Mach. Learn.* **2013**, *92*, 349–376. — Mondrian conformal prediction;
    cited as [Vovk2013].
18. Lei, J.; G'Sell, M.; Rinaldo, A.; Tibshirani, R. J.; Wasserman, L.
    Distribution-Free Predictive Inference for Regression.
    *J. Am. Stat. Assoc.* **2018**, *113* (523), 1094–1111.
19. Papadopoulos, H.; Haralambous, H. Reliable prediction intervals with
    regression neural networks. *Neural Netw.* **2011**, *24* (8), 842–851.
    — normalised nonconformity.
20. Boström, H.; Johansson, U. Mondrian conformal regressors.
    *Proc. Mach. Learn. Res. (COPA)* **2020**, *128*, 114–133.
21. Barber, R. F.; Candès, E. J.; Ramdas, A.; Tibshirani, R. J. Predictive
    inference with the jackknife+. *Ann. Stat.* **2021**, *49* (1), 486–507.
    — the data-reusing alternative to the split construction, noted in
    Section 2.7.

## Representations and evaluation protocol

22. Bemis, G. W.; Murcko, M. A. The Properties of Known Drugs. 1. Molecular
    Frameworks. *J. Med. Chem.* **1996**, *39* (15), 2887–2893.
23. Rogers, D.; Hahn, M. Extended-Connectivity Fingerprints.
    *J. Chem. Inf. Model.* **2010**, *50* (5), 742–754.
24. Sheridan, R. P. Time-Split Cross-Validation as a Method for Estimating the
    Effectiveness of QSAR Models. *J. Chem. Inf. Model.* **2013**, *53* (4),
    783–790.
25. Wu, Z.; et al. MoleculeNet: a benchmark for molecular machine learning.
    *Chem. Sci.* **2018**, *9*, 513–530. — scaffold splits.

## Software

26. Pedregosa, F.; et al. Scikit-learn: Machine Learning in Python.
    *J. Mach. Learn. Res.* **2011**, *12*, 2825–2830. — cited as [Pedregosa2011].
27. RDKit: Open-source cheminformatics. https://www.rdkit.org — cited as [RDKit].
28. Teh, S. X. Polymer Glass-Transition Temperature Prediction and Polymer-Family
    Generalisation. Software and report, GitHub repository.
    https://github.com/tehsongxuan/Polymer-Tg-Prediction-and-Polymer-Family-Generalisation
    (accessed 14 September 2026). — the closest prior analysis of this dataset,
    cited as [Teh]: it reports the 7,208 → 7,174 curation audit, the 31/34/28
    aggregation counts and the 115 K maximum within-structure range reproduced in
    Section 3.1, and states the family-holdout confound that Section 2.5 sets out
    to resolve. Bibliographic details were not re-verified from this environment;
    **check the repository, its authorship and its state before submission**.
