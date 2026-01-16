"""
Shared fixtures for moleculariq_core tests.

This module provides common test fixtures including:
- Sample SMILES strings for various molecule types
- Pre-configured solver instances
- Common test data structures
"""

import pytest
from moleculariq_core import SymbolicSolver, NaturalLanguageFormatter


# =============================================================================
# Solver Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def solver():
    """Shared SymbolicSolver instance for all tests."""
    return SymbolicSolver()


@pytest.fixture(scope="module")
def formatter():
    """Shared NaturalLanguageFormatter instance."""
    return NaturalLanguageFormatter(seed=42)


# =============================================================================
# Sample SMILES - Simple Molecules
# =============================================================================

@pytest.fixture
def ethanol():
    """Ethanol: C2H5OH - simple alcohol."""
    return "CCO"


@pytest.fixture
def methane():
    """Methane: CH4 - simplest hydrocarbon."""
    return "C"


@pytest.fixture
def water():
    """Water: H2O."""
    return "O"


@pytest.fixture
def benzene():
    """Benzene: C6H6 - aromatic ring."""
    return "c1ccccc1"


@pytest.fixture
def cyclohexane():
    """Cyclohexane: C6H12 - saturated ring."""
    return "C1CCCCC1"


# =============================================================================
# Sample SMILES - Functional Groups
# =============================================================================

@pytest.fixture
def acetone():
    """Acetone: CH3COCH3 - ketone."""
    return "CC(=O)C"


@pytest.fixture
def acetic_acid():
    """Acetic acid: CH3COOH - carboxylic acid."""
    return "CC(=O)O"


@pytest.fixture
def acetaldehyde():
    """Acetaldehyde: CH3CHO - aldehyde."""
    return "CC=O"


@pytest.fixture
def ethylamine():
    """Ethylamine: C2H5NH2 - primary amine."""
    return "CCN"


@pytest.fixture
def chloroethane():
    """Chloroethane: C2H5Cl - halogenated."""
    return "CCCl"


# =============================================================================
# Sample SMILES - Complex Molecules
# =============================================================================

@pytest.fixture
def aspirin():
    """Aspirin (acetylsalicylic acid): C9H8O4."""
    return "CC(=O)Oc1ccccc1C(=O)O"


@pytest.fixture
def caffeine():
    """Caffeine: C8H10N4O2."""
    return "Cn1cnc2c1c(=O)n(c(=O)n2C)C"


@pytest.fixture
def ibuprofen():
    """Ibuprofen: C13H18O2."""
    return "CC(C)Cc1ccc(cc1)C(C)C(=O)O"


@pytest.fixture
def glucose():
    """Glucose: C6H12O6."""
    return "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"


# =============================================================================
# Sample SMILES - Stereochemistry
# =============================================================================

@pytest.fixture
def r_alanine():
    """R-Alanine: chiral amino acid."""
    return "C[C@H](N)C(=O)O"


@pytest.fixture
def s_alanine():
    """S-Alanine: chiral amino acid."""
    return "C[C@@H](N)C(=O)O"


@pytest.fixture
def e_butene():
    """E-2-butene: trans double bond."""
    return "C/C=C/C"


@pytest.fixture
def z_butene():
    """Z-2-butene: cis double bond."""
    return "C/C=C\\C"


# =============================================================================
# Sample SMILES - Ring Systems
# =============================================================================

@pytest.fixture
def naphthalene():
    """Naphthalene: fused aromatic rings."""
    return "c1ccc2ccccc2c1"


@pytest.fixture
def adamantane():
    """Adamantane: bridged ring system."""
    return "C1C2CC3CC1CC(C2)C3"


@pytest.fixture
def spiro_compound():
    """Spiro[4.5]decane: spiro ring system."""
    return "C1CCC2(CC1)CCCCC2"


# =============================================================================
# Sample SMILES - Branched and Chain Molecules
# =============================================================================

@pytest.fixture
def isobutane():
    """Isobutane: branched alkane with branch point."""
    return "CC(C)C"


@pytest.fixture
def neopentane():
    """Neopentane: highly branched alkane."""
    return "CC(C)(C)C"


@pytest.fixture
def hexane():
    """n-Hexane: linear alkane with chain termini."""
    return "CCCCCC"


@pytest.fixture
def isopropanol():
    """Isopropanol: secondary alcohol."""
    return "CC(O)C"


# =============================================================================
# Sample SMILES - Heterocycles
# =============================================================================

@pytest.fixture
def pyridine():
    """Pyridine: aromatic heterocycle with nitrogen."""
    return "c1ccncc1"


@pytest.fixture
def furan():
    """Furan: aromatic heterocycle with oxygen."""
    return "c1ccoc1"


@pytest.fixture
def thiophene():
    """Thiophene: aromatic heterocycle with sulfur."""
    return "c1ccsc1"


@pytest.fixture
def imidazole():
    """Imidazole: aromatic heterocycle with 2 nitrogens."""
    return "c1c[nH]cn1"


@pytest.fixture
def piperidine():
    """Piperidine: saturated heterocycle."""
    return "C1CCNCC1"


@pytest.fixture
def morpholine():
    """Morpholine: heterocycle with N and O."""
    return "C1COCCN1"


# =============================================================================
# Sample SMILES - E/Z Stereochemistry
# =============================================================================

@pytest.fixture
def stilbene_e():
    """E-Stilbene: trans double bond."""
    return "C(/C=C/c1ccccc1)c1ccccc1"


@pytest.fixture
def stilbene_z():
    """Z-Stilbene: cis double bond."""
    return r"C(\C=C/c1ccccc1)c1ccccc1"


# =============================================================================
# Sample SMILES - Polycyclic Systems
# =============================================================================

@pytest.fixture
def anthracene():
    """Anthracene: three fused aromatic rings."""
    return "c1ccc2cc3ccccc3cc2c1"


@pytest.fixture
def phenanthrene():
    """Phenanthrene: three fused aromatic rings (angular)."""
    return "c1ccc2c(c1)ccc1ccccc12"


@pytest.fixture
def decalin():
    """Decalin: fused saturated rings."""
    return "C1CCC2CCCCC2C1"


@pytest.fixture
def bicyclo_heptane():
    """Bicyclo[2.2.1]heptane (norbornane): bridged ring system."""
    return "C1CC2CCC1C2"


# =============================================================================
# Sample SMILES - Functional Group Rich
# =============================================================================

@pytest.fixture
def para_aminobenzoic_acid():
    """p-Aminobenzoic acid: amine and carboxylic acid."""
    return "Nc1ccc(C(=O)O)cc1"


@pytest.fixture
def nitrobenzene():
    """Nitrobenzene: nitro group."""
    return "c1ccc(cc1)[N+](=O)[O-]"


@pytest.fixture
def dimethyl_sulfoxide():
    """DMSO: sulfoxide."""
    return "CS(=O)C"


@pytest.fixture
def ethyl_acetate():
    """Ethyl acetate: ester."""
    return "CCOC(=O)C"


@pytest.fixture
def acetamide():
    """Acetamide: amide."""
    return "CC(=O)N"


@pytest.fixture
def acetonitrile():
    """Acetonitrile: nitrile."""
    return "CC#N"


@pytest.fixture
def diethyl_ether():
    """Diethyl ether: ether."""
    return "CCOCC"


# =============================================================================
# Sample SMILES - Drug-like Molecules
# =============================================================================

@pytest.fixture
def paracetamol():
    """Paracetamol (acetaminophen): common analgesic."""
    return "CC(=O)Nc1ccc(O)cc1"


@pytest.fixture
def naproxen():
    """Naproxen: NSAID with stereocenter."""
    return "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"


@pytest.fixture
def omeprazole():
    """Omeprazole: proton pump inhibitor."""
    return "COc1ccc2nc([nH]c2c1)S(=O)Cc1ncc(C)c(OC)c1C"


@pytest.fixture
def atorvastatin():
    """Atorvastatin: statin with multiple stereocenters."""
    return "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"


# =============================================================================
# Sample SMILES - Oxidation State Testing
# =============================================================================

@pytest.fixture
def formic_acid():
    """Formic acid: carbon at high oxidation state."""
    return "C(=O)O"


@pytest.fixture
def formaldehyde():
    """Formaldehyde: aldehyde carbon."""
    return "C=O"


@pytest.fixture
def methanol():
    """Methanol: carbon at low oxidation state."""
    return "CO"


@pytest.fixture
def carbon_dioxide():
    """Carbon dioxide: carbon at maximum oxidation."""
    return "O=C=O"


@pytest.fixture
def urea():
    """Urea: carbonyl with nitrogens."""
    return "NC(=O)N"


# =============================================================================
# Sample SMILES - Large/Complex Molecules
# =============================================================================

@pytest.fixture
def cholesterol():
    """Cholesterol: complex steroid."""
    return "CC(C)CCCC(C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"


@pytest.fixture
def taxol_core():
    """Taxol core scaffold (simplified)."""
    return "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)C)O)OC(=O)C)(C)[C@@H](C4(C)C)O)O)C)OC(=O)C"


# =============================================================================
# Invalid SMILES for Edge Case Testing
# =============================================================================

@pytest.fixture
def invalid_smiles_list():
    """Collection of invalid SMILES strings."""
    return [
        "",              # Empty
        "   ",           # Whitespace only
        "XYZ",           # Invalid atoms
        "C(C(C",         # Unbalanced parentheses
        "C1CC",          # Unclosed ring
        None,            # None value
        123,             # Wrong type
        "c1ccccc",       # Invalid aromatic (unclosed)
    ]


# =============================================================================
# Ground Truth Data for Regression Tests
# =============================================================================

@pytest.fixture
def known_molecule_properties():
    """
    Dictionary of molecules with known properties for regression testing.
    Format: {smiles: {property: expected_value}}
    """
    return {
        # Ethanol
        "CCO": {
            "ring_count": 0,
            "carbon_atom_count": 2,
            "heavy_atom_count": 3,
            "hbd_count": 1,  # OH
            "hba_count": 1,  # O
            "rotatable_bond_count": 0,
        },
        # Benzene
        "c1ccccc1": {
            "ring_count": 1,
            "aromatic_ring_count": 1,
            "aliphatic_ring_count": 0,
            "carbon_atom_count": 6,
            "heavy_atom_count": 6,
        },
        # Cyclohexane
        "C1CCCCC1": {
            "ring_count": 1,
            "aromatic_ring_count": 0,
            "aliphatic_ring_count": 1,
            "saturated_ring_count": 1,
            "carbon_atom_count": 6,
        },
        # Naphthalene
        "c1ccc2ccccc2c1": {
            "ring_count": 2,
            "aromatic_ring_count": 2,
            "fused_ring_count": 1,  # Counts fused ring systems, not individual rings
            "carbon_atom_count": 10,
        },
        # Aspirin
        "CC(=O)Oc1ccccc1C(=O)O": {
            "ring_count": 1,
            "aromatic_ring_count": 1,
            "carbon_atom_count": 9,
            "hbd_count": 1,  # COOH
            "hba_count": 3,  # RDKit HBA count for aspirin
        },
    }
