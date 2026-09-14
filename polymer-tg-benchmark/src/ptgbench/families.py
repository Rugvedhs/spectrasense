"""Backbone-aware assignment of polymer repeat units to chemical families.

Why backbone-aware?  A poly(alkyl acrylate) carries an ester group, but that ester
sits on a *side chain*: the chain that determines segmental mobility is an
all-carbon vinyl backbone.  Classifying on whole-molecule substructure matches
would label such a repeat unit a polyester, which is chemically wrong and would
contaminate any family-holdout experiment built on those labels.  We therefore
locate the backbone first and only then look for characteristic linkages on it.

The backbone is the union of shortest paths between the polymerisation
attachment points (the ``*`` dummy atoms of a PSMILES).  Linkage patterns are
tested in a fixed priority order because the groups are nested: every imide
contains two amide-like C(=O)N motifs, every urethane contains both an ester-like
and an amide-like fragment, and so on.  The first pattern whose *core* atoms lie
entirely on the backbone wins.
"""

from __future__ import annotations

from typing import Iterable, Sequence

from rdkit import Chem

# --------------------------------------------------------------------------
# Linkage patterns, most specific first.
#
# ``core`` lists the indices *within the SMARTS match* that must lie on the
# backbone.  Carbonyl oxygens are deliberately excluded from the core: they hang
# off the backbone rather than forming part of the chain path.
# --------------------------------------------------------------------------
_LINKAGE_PATTERNS: Sequence[tuple[str, str, tuple[int, ...]]] = (
    # Cyclic imide: the N and both carbonyl carbons share a ring.
    ("Polyimides", "[#6R](=O)[#7R]([#6R]=O)", (0, 2, 3)),
    ("Polyanhydrides", "[CX3](=O)[OX2][CX3](=O)", (0, 2, 3)),
    ("Polyurethanes", "[NX3][CX3](=O)[OX2]", (0, 1, 3)),
    ("Polyureas", "[NX3][CX3](=O)[NX3]", (0, 1, 3)),
    ("Polycarbonates", "[OX2][CX3](=O)[OX2]", (0, 1, 3)),
    ("Polyamides", "[NX3][CX3]=O", (0, 1)),
    ("Polyesters", "[OX2][CX3]=O", (0, 1)),
    ("Polysiloxanes", "[Si][OX2]", (0, 1)),
    ("Polyphosphazenes", "[P]=[NX2]", (0, 1)),
    ("Polysulfones", "[SX4](=O)(=O)", (0,)),
    ("Polyimines", "[CX3]=[NX2]", (0, 1)),
    ("Polysulfides", "[#6,#0][SX2][#6,#0]", (0, 1, 2)),
    ("Polyethers", "[#6,#0][OX2][#6,#0]", (0, 1, 2)),
)

_COMPILED: list[tuple[str, Chem.Mol, tuple[int, ...]]] = [
    (name, Chem.MolFromSmarts(sma), core) for name, sma, core in _LINKAGE_PATTERNS
]

# Carbon-backbone subfamilies, distinguished by what decorates the backbone.
# Pendant-group patterns for carbon-backbone polymers.  Index 0 of each match is
# the backbone atom carrying the group, which lets us insist that the
# substituent really hangs off the chain rather than off a side chain.
_PENDANT_PATTERNS: Sequence[tuple[str, str]] = (
    ("Polyacrylics", "[CX4][CX3](=O)[OX2,NX3]"),
    ("Polyacrylics", "[CX4][CX2]#[NX1]"),
    ("Polystyrenes", "[CX4][c]"),
    ("Polyhalo-olefins", "[CX4][F,Cl,Br,I]"),
    ("Polyvinyls", "[CX4][O,N,S]"),
)
_COMPILED_PENDANT: list[tuple[str, Chem.Mol]] = [
    (name, Chem.MolFromSmarts(sma)) for name, sma in _PENDANT_PATTERNS
]

UNASSIGNED = "Unassigned"
OTHER_BACKBONE = "Other backbone"


def _dummy_indices(mol: Chem.Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]


def dimerised(mol: Chem.Mol) -> tuple[Chem.Mol, list[int]] | None:
    """Join two copies of the repeat unit head-to-tail.

    A PSMILES repeat unit is a chain fragment cut at an arbitrary point, and that
    point is very often the characteristic linkage itself. Nylon-6 written as
    ``*NCCCCCC(=O)*`` has its amide split between the two ends, so the path
    between the attachment points contains no amide at all and the unit reads as
    a plain hydrocarbon chain.

    Joining *two* copies restores the linkage the cut destroyed, and does so
    without inventing strain: closing a single unit onto itself would bond the two
    ends of, say, a para-substituted aromatic into a bridged bicycle that
    aromatises into something the original polymer does not contain. The dimer's
    junction is a real bond in the real chain.

    Returns the dimer and its two surviving attachment points, or ``None`` when
    the unit does not present exactly two terminal attachment points.
    """
    stars = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(stars) != 2:
        return None
    for star in stars:
        if mol.GetAtomWithIdx(star).GetDegree() != 1:
            return None

    n = mol.GetNumAtoms()
    combo = Chem.RWMol(Chem.CombineMols(mol, mol))

    # Bond the tail of the first copy to the head of the second.
    tail_star, head_star = stars[1], stars[0] + n
    tail_anchor = combo.GetAtomWithIdx(tail_star).GetNeighbors()[0].GetIdx()
    head_anchor = combo.GetAtomWithIdx(head_star).GetNeighbors()[0].GetIdx()
    if combo.GetBondBetweenAtoms(tail_anchor, head_anchor) is not None:
        return None

    try:
        combo.AddBond(tail_anchor, head_anchor, Chem.BondType.SINGLE)
        for star in sorted((tail_star, head_star), reverse=True):
            combo.RemoveAtom(star)
        dimer = combo.GetMol()
        Chem.SanitizeMol(dimer)
    except Exception:
        return None

    remaining = [a.GetIdx() for a in dimer.GetAtoms() if a.GetAtomicNum() == 0]
    if len(remaining) != 2:
        return None
    return dimer, remaining


def backbone_atoms(mol: Chem.Mol) -> set[int]:
    """Atom indices lying on a shortest path between two attachment points.

    Repeat units with more than two ``*`` atoms (branch points) contribute every
    pairwise path.  Units with fewer than two attachment points have no
    well-defined backbone, so the whole molecule is returned and the caller
    degrades to whole-molecule matching.
    """
    stars = _dummy_indices(mol)
    if len(stars) < 2:
        return {a.GetIdx() for a in mol.GetAtoms()}

    atoms: set[int] = set()
    for i, start in enumerate(stars):
        for end in stars[i + 1 :]:
            path = Chem.GetShortestPath(mol, start, end)
            atoms.update(path)
    return atoms or {a.GetIdx() for a in mol.GetAtoms()}


def _backbone_ring_closure(mol: Chem.Mol, backbone: set[int]) -> set[int]:
    """Extend the backbone with whole rings it passes through.

    A backbone path that clips two atoms of an aromatic ring still means the ring
    itself is in the chain; imide and phenylene recognition both depend on the
    full ring being counted as backbone.
    """
    extended = set(backbone)
    for ring in mol.GetRingInfo().AtomRings():
        if len(set(ring) & backbone) >= 2:
            extended.update(ring)
    return extended


def _match_on_backbone(
    mol: Chem.Mol, patt: Chem.Mol, core: Iterable[int], backbone: set[int]
) -> bool:
    core = tuple(core)
    for match in mol.GetSubstructMatches(patt):
        if all(match[i] in backbone for i in core):
            return True
    return False


def _carbon_backbone_family(mol: Chem.Mol, backbone: set[int]) -> str:
    """Classify an all-carbon backbone by the groups pendant to it."""
    aromatic_backbone = sum(
        1 for i in backbone if mol.GetAtomWithIdx(i).GetIsAromatic()
    )
    heavy_backbone = sum(
        1 for i in backbone if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
    )
    # A chain that is mostly aromatic ring carbon is a phenylene-type polymer.
    if heavy_backbone and aromatic_backbone / heavy_backbone >= 0.6:
        return "Polyphenylenes"

    for name, patt in _COMPILED_PENDANT:
        if patt is None:
            continue
        # match[0] is the substituted atom: require it to sit on the backbone.
        if any(match[0] in backbone for match in mol.GetSubstructMatches(patt)):
            return name

    # Abstain rather than absorb. Anything the linkage patterns miss but that
    # still carries a chain heteroatom is not an olefin: poly(dimethylsilane),
    # polycarbosilanes and backbone phosphate esters all reached "Polyolefins"
    # by default, which flattered the taxonomy into reporting no unassigned
    # structures at all.
    chain_heteroatoms = {
        mol.GetAtomWithIdx(i).GetSymbol() for i in backbone
    } - {"C", "H", "*"}
    if chain_heteroatoms:
        return OTHER_BACKBONE

    has_backbone_unsaturation = any(
        bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetBeginAtomIdx() in backbone
        and bond.GetEndAtomIdx() in backbone
        and not bond.GetIsAromatic()
        for bond in mol.GetBonds()
    )
    if has_backbone_unsaturation:
        return "Polydienes"
    return "Polyolefins"


def assign_family(psmiles: str) -> str:
    """Return the polymer family of a PSMILES repeat unit.

    Linkages are sought on a head-to-tail *dimer* wherever one can be formed, so a
    unit cut at its own characteristic linkage is still recognised; the open form
    is used only as a fallback. Returns :data:`UNASSIGNED` when the string
    cannot be parsed.
    """
    mol = Chem.MolFromSmiles(psmiles)
    if mol is None:
        return UNASSIGNED

    dimer = dimerised(mol)
    if dimer is not None:
        dimer_mol, _ = dimer
        backbone = _backbone_ring_closure(dimer_mol, backbone_atoms(dimer_mol))
        for name, patt, core in _COMPILED:
            if patt is not None and _match_on_backbone(dimer_mol, patt, core, backbone):
                return name
        return _carbon_backbone_family(dimer_mol, backbone)

    backbone = _backbone_ring_closure(mol, backbone_atoms(mol))
    for name, patt, core in _COMPILED:
        if patt is not None and _match_on_backbone(mol, patt, core, backbone):
            return name
    return _carbon_backbone_family(mol, backbone)


def assign_families(psmiles: Sequence[str]) -> list[str]:
    """Vectorised convenience wrapper around :func:`assign_family`."""
    return [assign_family(s) for s in psmiles]
