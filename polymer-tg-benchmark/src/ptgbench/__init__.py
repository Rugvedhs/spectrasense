"""Benchmarking uncertainty-aware polymer glass-transition-temperature models.

The package is organised as a linear pipeline:

``data``      load / canonicalise / aggregate the raw repeat-unit records
``families``  backbone-aware assignment of each repeat unit to a polymer family
``features``  RDKit descriptors and Morgan fingerprints
``splits``    random, scaffold, cluster and family-holdout partitions
``models``    the regressor zoo used across every split regime
``conformal`` split / normalised / Mondrian conformal predictors
``metrics``   point-accuracy, interval and subgroup-coverage metrics
"""

__version__ = "0.1.0"
