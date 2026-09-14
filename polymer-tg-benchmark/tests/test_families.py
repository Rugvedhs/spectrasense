"""The family classifier decides what the holdout experiments mean, so its
chemistry is pinned down case by case."""

import pytest

from ptgbench.families import assign_family

CASES = [
    ("*CC(*)c1ccccc1", "Polystyrenes"),                       # polystyrene
    ("*CC(*)C(=O)OC", "Polyacrylics"),                        # PMA
    ("*CC(C)(*)C(=O)OC", "Polyacrylics"),                     # PMMA
    ("*CC(*)C#N", "Polyacrylics"),                            # PAN
    ("*CC(*)Cl", "Polyhalo-olefins"),                         # PVC
    ("*C(F)(F)C(F)(F)*", "Polyhalo-olefins"),                 # PTFE
    ("*CC*", "Polyolefins"),                                  # PE
    ("*CC(C)*", "Polyolefins"),                               # PP
    ("*OCC*", "Polyethers"),                                  # PEO
    ("*CC(*)O", "Polyvinyls"),                                # PVA
    ("*CC(*)OC(C)=O", "Polyvinyls"),                          # PVAc
    ("*OC(=O)c1ccc(C(=O)OCC*)cc1", "Polyesters"),             # PET-like
    ("*NC(=O)CCCCC(=O)N*", "Polyamides"),                     # nylon-like
    ("*NC(=O)OCCO*", "Polyurethanes"),
    ("*NC(=O)N*", "Polyureas"),
    ("*OC(=O)Oc1ccc(C(C)(C)c2ccc(*)cc2)cc1", "Polycarbonates"),
    ("*O[Si](C)(C)*", "Polysiloxanes"),                       # PDMS
    ("*N(*)C(=O)c1ccc2c(c1)C(=O)N(*)C2=O", "Polyimides"),
    ("*c1ccc(*)cc1", "Polyphenylenes"),
    ("*CC=CC*", "Polydienes"),                                # polybutadiene
    ("*SCC*", "Polysulfides"),
    ("*C(=O)OC(=O)CC*", "Polyanhydrides"),
]


@pytest.mark.parametrize("psmiles,expected", CASES)
def test_known_chemistry(psmiles, expected):
    assert assign_family(psmiles) == expected


def test_side_chain_ester_is_not_a_polyester():
    """The distinction the whole family analysis rests on."""
    assert assign_family("*CC(*)C(=O)OCCCC") == "Polyacrylics"
    assert assign_family("*OC(=O)CCCCC(=O)OCC*") == "Polyesters"


def test_unparseable_input_is_flagged_not_guessed():
    assert assign_family("not a smiles") == "Unassigned"


def test_priority_imide_beats_amide():
    """A cyclic imide contains two amide-like motifs; specificity must win."""
    assert assign_family("*N(*)C(=O)c1ccc2c(c1)C(=O)N(*)C2=O") == "Polyimides"


# Adversarial cases: repeat units whose characteristic linkage is split across
# the PSMILES boundary. Classifying the open fragment sees no linkage at all and
# reads these as hydrocarbon or aromatic chains, which is how nylon-6 and a
# polypeptide came to be labelled Polyvinyls.
BOUNDARY_CASES = [
    ("*NCCCCCC(=O)*", "Polyamides"),                  # nylon-6, cut at the amide
    ("*NC(C)CC(*)=O", "Polyamides"),                   # nylon-3
    ("*NC(CCC(=O)OC)C(*)=O", "Polyamides"),            # poly(gamma-methyl glutamate)
    ("*C(=O)c1ccc(N(*)CCC)cc1", "Polyamides"),         # N-alkyl aramid
    ("*OC(=O)CCCCCCCCC(*)=O", "Polyanhydrides"),       # sebacic polyanhydride
]


@pytest.mark.parametrize("psmiles,expected", BOUNDARY_CASES)
def test_linkage_split_across_the_repeat_unit_boundary(psmiles, expected):
    assert assign_family(psmiles) == expected


def test_backbone_heteroatoms_abstain_rather_than_becoming_polyolefins():
    """A silicon-backbone chain is not an olefin.

    _carbon_backbone_family used to absorb anything the linkage patterns missed,
    so poly(dimethylsilane) and polycarbosilanes were reported as Polyolefins and
    the taxonomy showed no unassigned structures at all.
    """
    assert assign_family("*[Si](C)(C)*") == "Other backbone"
    assert assign_family("*CC*") == "Polyolefins"


def test_dimerisation_does_not_fire_on_a_para_aromatic_self_closure():
    """Closing one unit onto itself would bridge a para-substituted ring.

    Joining the two ends of *c1ccc(*)cc1 within a single unit makes a bicyclic
    that aromatises into something the polymer does not contain; the dimer keeps
    the chemistry real.
    """
    assert assign_family("*c1ccc(*)cc1") == "Polyphenylenes"
