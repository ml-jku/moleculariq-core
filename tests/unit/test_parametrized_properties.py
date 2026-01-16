"""
Parametrized property tests for SymbolicSolver.

Uses pytest parametrize to systematically test properties across many molecules.
"""

import pytest
from moleculariq_core import SymbolicSolver


@pytest.fixture(scope="module")
def solver():
    """Create solver instance for all tests."""
    return SymbolicSolver()


# =============================================================================
# Test Data - Molecules with known properties
# =============================================================================

RING_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("CC", 0, "ethane"),
    ("CCO", 0, "ethanol"),
    ("CCCCCC", 0, "hexane"),
    ("c1ccccc1", 1, "benzene"),
    ("C1CCCCC1", 1, "cyclohexane"),
    ("c1ccncc1", 1, "pyridine"),
    ("C1CCOC1", 1, "tetrahydrofuran"),
    ("c1ccc2ccccc2c1", 2, "naphthalene"),
    ("c1ccc2cc3ccccc3cc2c1", 3, "anthracene"),
]

AROMATIC_RING_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("CCO", 0, "ethanol"),
    ("C1CCCCC1", 0, "cyclohexane"),
    ("c1ccccc1", 1, "benzene"),
    ("c1ccncc1", 1, "pyridine"),
    ("c1ccc2ccccc2c1", 2, "naphthalene"),
    ("c1ccc2cc3ccccc3cc2c1", 3, "anthracene"),
]

CARBON_COUNT_TEST_DATA = [
    ("C", 1, "methane"),
    ("CC", 2, "ethane"),
    ("CCO", 2, "ethanol"),
    ("C(=O)O", 1, "formic acid"),
    ("CC(=O)O", 2, "acetic acid"),
    ("c1ccccc1", 6, "benzene"),
    ("CCCCCC", 6, "hexane"),
    ("c1ccc2ccccc2c1", 10, "naphthalene"),
]

HETERO_ATOM_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("CC", 0, "ethane"),
    ("c1ccccc1", 0, "benzene"),
    ("CCO", 1, "ethanol"),
    ("CCN", 1, "ethylamine"),
    ("c1ccncc1", 1, "pyridine"),
    ("OC=O", 2, "formic acid"),
    ("CC(=O)O", 2, "acetic acid"),
    ("CC(=O)NC", 2, "N-methylacetamide"),
]

HBD_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("c1ccccc1", 0, "benzene"),
    ("CC(=O)C", 0, "acetone"),
    ("CCO", 1, "ethanol"),
    ("CCN", 1, "ethylamine"),
    ("OC=O", 1, "formic acid"),
]

HBA_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("c1ccccc1", 0, "benzene"),
    ("CCO", 1, "ethanol"),
    ("CC(=O)C", 1, "acetone"),
    ("CCN", 1, "ethylamine"),
    ("c1ccncc1", 1, "pyridine"),
]

HALOGEN_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("CCO", 0, "ethanol"),
    ("c1ccccc1", 0, "benzene"),
    ("CCCl", 1, "chloroethane"),
    ("CCBr", 1, "bromoethane"),
    ("CCI", 1, "iodoethane"),
    ("CCF", 1, "fluoroethane"),
    ("ClCCCl", 2, "1,2-dichloroethane"),
    ("FC(F)(F)C", 3, "1,1,1-trifluoroethane"),
]

STEREOCENTER_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("CCO", 0, "ethanol"),
    ("c1ccccc1", 0, "benzene"),
    ("C[C@H](N)C(=O)O", 1, "alanine"),
]

CSP3_CARBON_COUNT_TEST_DATA = [
    ("c1ccccc1", 0, "benzene - all sp2"),
    ("C=C", 0, "ethene - sp2 carbons"),
    ("C", 1, "methane"),
    ("CC", 2, "ethane"),
    ("CCO", 2, "ethanol"),
    ("CC(=O)C", 2, "acetone - 2 methyls are sp3"),
    ("CCCCCC", 6, "hexane"),
    ("C1CCCCC1", 6, "cyclohexane"),
]

HEAVY_ATOM_COUNT_TEST_DATA = [
    ("C", 1, "methane"),
    ("CC", 2, "ethane"),
    ("CCO", 3, "ethanol"),
    ("c1ccccc1", 6, "benzene"),
    ("CCCl", 3, "chloroethane"),
    ("O", 1, "water"),
]

HYDROGEN_COUNT_TEST_DATA = [
    ("C", 4, "methane"),
    ("CC", 6, "ethane"),
    ("CCO", 6, "ethanol"),
    ("c1ccccc1", 6, "benzene"),
    ("C1CCCCC1", 12, "cyclohexane"),
    ("O", 2, "water"),
]

ROTATABLE_BOND_COUNT_TEST_DATA = [
    ("C", 0, "methane"),
    ("CC", 0, "ethane"),
    ("CCO", 0, "ethanol"),
    ("c1ccccc1", 0, "benzene"),
    ("C1CCCCC1", 0, "cyclohexane"),
    ("CCCCC", 2, "pentane"),
    ("CCCCCC", 3, "hexane"),
]


# =============================================================================
# Parametrized Tests
# =============================================================================

class TestParametrizedRingCount:
    """Parametrized tests for ring count."""

    @pytest.mark.parametrize("smiles,expected,name", RING_COUNT_TEST_DATA)
    def test_ring_count(self, solver, smiles, expected, name):
        """Ring count for {name}."""
        assert solver.get_ring_count(smiles) == expected


class TestParametrizedAromaticRingCount:
    """Parametrized tests for aromatic ring count."""

    @pytest.mark.parametrize("smiles,expected,name", AROMATIC_RING_COUNT_TEST_DATA)
    def test_aromatic_ring_count(self, solver, smiles, expected, name):
        """Aromatic ring count for {name}."""
        assert solver.get_aromatic_ring_count(smiles) == expected


class TestParametrizedCarbonCount:
    """Parametrized tests for carbon count."""

    @pytest.mark.parametrize("smiles,expected,name", CARBON_COUNT_TEST_DATA)
    def test_carbon_count(self, solver, smiles, expected, name):
        """Carbon count for {name}."""
        assert solver.get_carbon_atom_count(smiles) == expected


class TestParametrizedHeteroAtomCount:
    """Parametrized tests for hetero atom count."""

    @pytest.mark.parametrize("smiles,expected,name", HETERO_ATOM_COUNT_TEST_DATA)
    def test_hetero_atom_count(self, solver, smiles, expected, name):
        """Hetero atom count for {name}."""
        assert solver.get_hetero_atom_count(smiles) == expected


class TestParametrizedHBDCount:
    """Parametrized tests for hydrogen bond donor count."""

    @pytest.mark.parametrize("smiles,expected,name", HBD_COUNT_TEST_DATA)
    def test_hbd_count(self, solver, smiles, expected, name):
        """HBD count for {name}."""
        assert solver.get_hbd_count(smiles) == expected


class TestParametrizedHBACount:
    """Parametrized tests for hydrogen bond acceptor count."""

    @pytest.mark.parametrize("smiles,expected,name", HBA_COUNT_TEST_DATA)
    def test_hba_count(self, solver, smiles, expected, name):
        """HBA count for {name}."""
        assert solver.get_hba_count(smiles) == expected


class TestParametrizedHalogenCount:
    """Parametrized tests for halogen count."""

    @pytest.mark.parametrize("smiles,expected,name", HALOGEN_COUNT_TEST_DATA)
    def test_halogen_count(self, solver, smiles, expected, name):
        """Halogen count for {name}."""
        assert solver.get_halogen_atom_count(smiles) == expected


class TestParametrizedStereocenterCount:
    """Parametrized tests for stereocenter count."""

    @pytest.mark.parametrize("smiles,expected,name", STEREOCENTER_COUNT_TEST_DATA)
    def test_stereocenter_count(self, solver, smiles, expected, name):
        """Stereocenter count for {name}."""
        assert solver.get_stereocenter_count(smiles) == expected


class TestParametrizedCsp3CarbonCount:
    """Parametrized tests for sp3 carbon count."""

    @pytest.mark.parametrize("smiles,expected,name", CSP3_CARBON_COUNT_TEST_DATA)
    def test_csp3_carbon_count(self, solver, smiles, expected, name):
        """sp3 carbon count for {name}."""
        assert solver.get_csp3_carbon_count(smiles) == expected


class TestParametrizedHeavyAtomCount:
    """Parametrized tests for heavy atom count."""

    @pytest.mark.parametrize("smiles,expected,name", HEAVY_ATOM_COUNT_TEST_DATA)
    def test_heavy_atom_count(self, solver, smiles, expected, name):
        """Heavy atom count for {name}."""
        assert solver.get_heavy_atom_count(smiles) == expected


class TestParametrizedHydrogenCount:
    """Parametrized tests for hydrogen count."""

    @pytest.mark.parametrize("smiles,expected,name", HYDROGEN_COUNT_TEST_DATA)
    def test_hydrogen_count(self, solver, smiles, expected, name):
        """Hydrogen count for {name}."""
        assert solver.get_hydrogen_count(smiles) == expected


class TestParametrizedRotatableBondCount:
    """Parametrized tests for rotatable bond count."""

    @pytest.mark.parametrize("smiles,expected,name", ROTATABLE_BOND_COUNT_TEST_DATA)
    def test_rotatable_bond_count(self, solver, smiles, expected, name):
        """Rotatable bond count for {name}."""
        assert solver.get_rotatable_bond_count(smiles) == expected


# =============================================================================
# Multi-property tests for specific molecules
# =============================================================================

MOLECULE_PROPERTY_MATRIX = [
    # (smiles, ring, aromatic_ring, carbon, hetero, hbd, hba, heavy, hydrogen)
    ("C", 0, 0, 1, 0, 0, 0, 1, 4),
    ("CCO", 0, 0, 2, 1, 1, 1, 3, 6),
    ("c1ccccc1", 1, 1, 6, 0, 0, 0, 6, 6),
    ("c1ccncc1", 1, 1, 5, 1, 0, 1, 6, 5),
    ("CC(=O)C", 0, 0, 3, 1, 0, 1, 4, 6),
]


class TestParametrizedMultiProperty:
    """Test multiple properties for each molecule."""

    @pytest.mark.parametrize(
        "smiles,ring,aromatic,carbon,hetero,hbd,hba,heavy,hydrogen",
        MOLECULE_PROPERTY_MATRIX
    )
    def test_all_properties(
        self, solver, smiles, ring, aromatic, carbon, hetero, hbd, hba, heavy, hydrogen
    ):
        """Test all basic properties for molecule."""
        assert solver.get_ring_count(smiles) == ring
        assert solver.get_aromatic_ring_count(smiles) == aromatic
        assert solver.get_carbon_atom_count(smiles) == carbon
        assert solver.get_hetero_atom_count(smiles) == hetero
        assert solver.get_hbd_count(smiles) == hbd
        assert solver.get_hba_count(smiles) == hba
        assert solver.get_heavy_atom_count(smiles) == heavy
        assert solver.get_hydrogen_count(smiles) == hydrogen


# =============================================================================
# Index consistency tests
# =============================================================================

INDEX_COUNT_CONSISTENCY = [
    ("c1ccccc1", "ring_indices", "ring_count"),
    ("c1ccccc1", "aromatic_ring_indices", "aromatic_ring_count"),
    ("C1CCCCC1", "aliphatic_ring_indices", "aliphatic_ring_count"),
    ("CCO", "carbon_atom_indices", "carbon_atom_count"),
    ("CCO", "hetero_atom_indices", "hetero_atom_count"),
]


class TestIndexCountConsistency:
    """Verify indices lists have correct length for counts."""

    @pytest.mark.parametrize("smiles,index_method,count_method", INDEX_COUNT_CONSISTENCY)
    def test_index_count_match(self, solver, smiles, index_method, count_method):
        """Index list length should match count for non-overlapping properties."""
        index_method_name = f"get_{index_method}"
        count_method_name = f"get_{count_method}"

        if hasattr(solver, index_method_name) and hasattr(solver, count_method_name):
            indices = getattr(solver, index_method_name)(smiles)
            count = getattr(solver, count_method_name)(smiles)

            # For ring indices, each ring's atoms are returned
            # For atom indices, count should match
            if "atom" in index_method:
                assert len(indices) == count
