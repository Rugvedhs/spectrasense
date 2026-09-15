# References

## Polymer Tg prediction and polymer informatics

1. Tao, L.; Varshney, V.; Li, Y. Benchmarking Machine Learning
   Models for Polymer Informatics: An Example of Glass Transition Temperature.
   *J. Chem. Inf. Model.* **2021**, *61* (11), 5395–5413.
   DOI: 10.1021/acs.jcim.1c01031
2. Casanola-Martin, G. M.; Karuth, A.; Pham-The, H.; et al. Machine
   learning analysis of a large set of homopolymers to predict glass transition
   temperatures. *Commun. Chem.* **2024**, *7*, 226.
   DOI: 10.1038/s42004-024-01305-0
3. Kuenneth, C.; Ramprasad, R. polyBERT: a chemical language model
   to enable fully machine-driven ultrafast polymer informatics.
   *Nat. Commun.* **2023**, *14*, 4099. DOI: 10.1038/s41467-023-39868-6
4. Xu, C.; Wang, Y.; Barati Farimani, A. TransPolymer: a
   Transformer-based language model for polymer property predictions.
   *npj Comput. Mater.* **2023**, *9*, 64. DOI: 10.1038/s41524-023-01016-5
5. Kunchapu, S.; Jablonka, K. M. PolyMetriX: an ecosystem for
   digital polymer chemistry. *npj Comput. Mater.* **2025**, *11*, 312.
6. Xu, J.; et al. POINT²: A Polymer Informatics Training and
   Testing Database. *arXiv* **2025**, arXiv:2503.23491.
   — source of the *T*~g~ collection used here.
7. Babbar, A.; et al. Explainability and Transferability of Machine
   Learning Models for Predicting the Glass Transition Temperature of Polymers.
   *arXiv* **2023**, arXiv:2308.09898.
8. Qiu, H.; et al. Design of Polyimides with Targeted Glass
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

13. Akdoğan, M. Trust Beyond Accuracy: Conformal Uncertainty
    Quantification Reveals the Generalization Gap in Polymer Glass Transition
    Temperature Prediction. *ACS Omega* **2026**, *11* (30), 44926–44941.
    — the closest prior work; 410 simulation-derived samples, marginal split
    conformal only.
14. Agrawal, V. Polymer property prediction using ensemble of
    directed message-passing neural networks with uncertainty estimation and
    applicability domain analysis. *MRS Commun.* **2026**.
    DOI: 10.1557/s43579-026-00948-5
15. Norinder, U.; Carlsson, L.; Boyer, S.; Eklund, M. Introducing Conformal
    Prediction in Predictive Modeling. A Transparent and Flexible Alternative to
    Applicability Domain Determination. *J. Chem. Inf. Model.* **2014**, *54*
    (6), 1596–1603. DOI: 10.1021/ci5001168 — conformal prediction used as the
    applicability-domain statement itself; cited as [Norinder2014].
16. Papadopoulos, H.; Proedrou, K.; Vovk, V.; Gammerman, A. Inductive Confidence
    Machines for Regression. In *Machine Learning: ECML 2002*; LNCS 2430;
    Springer, **2002**; pp 345–356.
17. Vovk, V.; Gammerman, A.; Shafer, G. *Algorithmic Learning in a Random World*;
    Springer: New York, **2005**.
18. Vovk, V. Conditional validity of inductive conformal predictors.
    *Mach. Learn.* **2013**, *92*, 349–376. — Mondrian conformal prediction;
    cited as [Vovk2013].
19. Lei, J.; G'Sell, M.; Rinaldo, A.; Tibshirani, R. J.; Wasserman, L.
    Distribution-Free Predictive Inference for Regression.
    *J. Am. Stat. Assoc.* **2018**, *113* (523), 1094–1111.
20. Papadopoulos, H.; Haralambous, H. Reliable prediction intervals with
    regression neural networks. *Neural Netw.* **2011**, *24* (8), 842–851.
    — normalised nonconformity.
21. Boström, H.; Johansson, U. Mondrian conformal regressors.
    *Proc. Mach. Learn. Res. (COPA)* **2020**, *128*, 114–133.
22. Barber, R. F.; Candès, E. J.; Ramdas, A.; Tibshirani, R. J. Predictive
    inference with the jackknife+. *Ann. Stat.* **2021**, *49* (1), 486–507.
    — the data-reusing alternative to the split construction, noted in
    Section 2.7.
23. Tibshirani, R. J.; Barber, R. F.; Candès, E. J.; Ramdas, A. Conformal
    Prediction Under Covariate Shift. In *Advances in Neural Information
    Processing Systems 32 (NeurIPS 2019)*; **2019**. arXiv:1904.06019
    — weighted conformal prediction; cited as [Tibshirani2019] and discussed in
    Sections 2.7 and 4.3.

## Representations and evaluation protocol

24. Bemis, G. W.; Murcko, M. A. The Properties of Known Drugs. 1. Molecular
    Frameworks. *J. Med. Chem.* **1996**, *39* (15), 2887–2893.
25. Rogers, D.; Hahn, M. Extended-Connectivity Fingerprints.
    *J. Chem. Inf. Model.* **2010**, *50* (5), 742–754.
26. Sheridan, R. P. Time-Split Cross-Validation as a Method for Estimating the
    Effectiveness of QSAR Models. *J. Chem. Inf. Model.* **2013**, *53* (4),
    783–790.
27. Wu, Z.; et al. MoleculeNet: a benchmark for molecular machine learning.
    *Chem. Sci.* **2018**, *9*, 513–530. — scaffold splits.

## Software

28. Pedregosa, F.; et al. Scikit-learn: Machine Learning in Python.
    *J. Mach. Learn. Res.* **2011**, *12*, 2825–2830. — cited as [Pedregosa2011].
29. RDKit: Open-source cheminformatics. https://www.rdkit.org — cited as [RDKit].
30. Teh, S. X. Polymer Glass-Transition Temperature Prediction and Polymer-Family
    Generalisation. Software and report, GitHub repository.
    https://github.com/tehsongxuan/Polymer-Tg-Prediction-and-Polymer-Family-Generalisation
    (accessed 14 September 2026). — a prior public analysis of this dataset,
    cited as [Teh]: it reports the 7,208 → 7,174 curation audit, the 31/34/28
    aggregation counts and the 115 K maximum within-structure range reproduced in
    Section 3.1, and states the family-holdout confound that Section 2.5 sets out
    to resolve.
