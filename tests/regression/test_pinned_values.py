"""
Strict regression tests with pinned expected values.

These tests use exact expected values to catch any unintended changes.
If a test fails, either:
1. The solver behavior changed (investigate the change)
2. The expected value was wrong (update with correct value)

DO NOT make these tests flexible/lenient - that defeats the purpose.
"""

import pytest
from moleculariq_core import SymbolicSolver


@pytest.fixture(scope="module")
def solver():
    return SymbolicSolver()


class TestEthanolPinnedValues:
    """Pinned values for ethanol (CCO)."""

    SMILES = "CCO"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 0

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 0

    def test_aliphatic_ring_count(self, solver):
        assert solver.get_aliphatic_ring_count(self.SMILES) == 0

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 2

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 3

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 1

    def test_hydrogen_count(self, solver):
        assert solver.get_hydrogen_count(self.SMILES) == 6

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 1

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 1

    def test_rotatable_bond_count(self, solver):
        assert solver.get_rotatable_bond_count(self.SMILES) == 0

    def test_stereocenter_count(self, solver):
        assert solver.get_stereocenter_count(self.SMILES) == 0

    def test_halogen_atom_count(self, solver):
        assert solver.get_halogen_atom_count(self.SMILES) == 0

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 2


class TestBenzenePinnedValues:
    """Pinned values for benzene (c1ccccc1)."""

    SMILES = "c1ccccc1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 1

    def test_aliphatic_ring_count(self, solver):
        assert solver.get_aliphatic_ring_count(self.SMILES) == 0

    def test_saturated_ring_count(self, solver):
        assert solver.get_saturated_ring_count(self.SMILES) == 0

    def test_fused_ring_count(self, solver):
        assert solver.get_fused_ring_count(self.SMILES) == 0

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 6

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 6

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 0

    def test_hydrogen_count(self, solver):
        assert solver.get_hydrogen_count(self.SMILES) == 6

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 0

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 0

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 0

    def test_smallest_ring_size(self, solver):
        assert solver.get_smallest_or_largest_ring_count(self.SMILES, smallest=True) == 6

    def test_largest_ring_size(self, solver):
        assert solver.get_smallest_or_largest_ring_count(self.SMILES, smallest=False) == 6

    def test_ring_indices(self, solver):
        indices = solver.get_ring_indices(self.SMILES)
        assert set(indices) == {0, 1, 2, 3, 4, 5}

    def test_bridgehead_count(self, solver):
        assert solver.get_bridgehead_count(self.SMILES) == 0

    def test_spiro_count(self, solver):
        assert solver.get_spiro_count(self.SMILES) == 0

    def test_heterocycle_count(self, solver):
        assert solver.get_heterocycle_count(self.SMILES) == 0


class TestCyclohexanePinnedValues:
    """Pinned values for cyclohexane (C1CCCCC1)."""

    SMILES = "C1CCCCC1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 0

    def test_aliphatic_ring_count(self, solver):
        assert solver.get_aliphatic_ring_count(self.SMILES) == 1

    def test_saturated_ring_count(self, solver):
        assert solver.get_saturated_ring_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 6

    def test_hydrogen_count(self, solver):
        assert solver.get_hydrogen_count(self.SMILES) == 12

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 6


class TestNaphthalenePinnedValues:
    """Pinned values for naphthalene (c1ccc2ccccc2c1)."""

    SMILES = "c1ccc2ccccc2c1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 2

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 2

    def test_aliphatic_ring_count(self, solver):
        assert solver.get_aliphatic_ring_count(self.SMILES) == 0

    def test_fused_ring_count(self, solver):
        assert solver.get_fused_ring_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 10

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 10

    def test_hydrogen_count(self, solver):
        assert solver.get_hydrogen_count(self.SMILES) == 8


class TestAspirinPinnedValues:
    """Pinned values for aspirin (CC(=O)Oc1ccccc1C(=O)O)."""

    SMILES = "CC(=O)Oc1ccccc1C(=O)O"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 9

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 13

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 4

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 1

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 3


class TestPyridinePinnedValues:
    """Pinned values for pyridine (c1ccncc1)."""

    SMILES = "c1ccncc1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 1

    def test_heterocycle_count(self, solver):
        assert solver.get_heterocycle_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 5

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 1

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 1

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 0


class TestHexanePinnedValues:
    """Pinned values for n-hexane (CCCCCC)."""

    SMILES = "CCCCCC"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 0

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 6

    def test_hydrogen_count(self, solver):
        assert solver.get_hydrogen_count(self.SMILES) == 14

    def test_chain_termini_count(self, solver):
        assert solver.get_chain_termini_count(self.SMILES) == 2

    def test_branch_point_count(self, solver):
        assert solver.get_branch_point_count(self.SMILES) == 0

    def test_longest_carbon_chain(self, solver):
        assert solver.get_longest_carbon_chain_count(self.SMILES) == 6

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 6


class TestIsobutanePinnedValues:
    """Pinned values for isobutane (CC(C)C)."""

    SMILES = "CC(C)C"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 0

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 4

    def test_branch_point_count(self, solver):
        assert solver.get_branch_point_count(self.SMILES) == 1

    def test_longest_carbon_chain(self, solver):
        assert solver.get_longest_carbon_chain_count(self.SMILES) == 3

    def test_chain_termini_count(self, solver):
        assert solver.get_chain_termini_count(self.SMILES) == 3


class TestEButenePinnedValues:
    """Pinned values for E-2-butene (C/C=C/C)."""

    SMILES = "C/C=C/C"

    def test_e_double_bond_count(self, solver):
        assert solver.get_e_z_stereochemistry_double_bond_count(self.SMILES, e_count=True) == 1

    def test_z_double_bond_count(self, solver):
        assert solver.get_e_z_stereochemistry_double_bond_count(self.SMILES, e_count=False) == 0

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 4


class TestZButenePinnedValues:
    """Pinned values for Z-2-butene (C/C=C\\C)."""

    SMILES = "C/C=C\\C"

    def test_e_double_bond_count(self, solver):
        assert solver.get_e_z_stereochemistry_double_bond_count(self.SMILES, e_count=True) == 0

    def test_z_double_bond_count(self, solver):
        assert solver.get_e_z_stereochemistry_double_bond_count(self.SMILES, e_count=False) == 1


class TestRAlaninePinnedValues:
    """Pinned values for R-Alanine (C[C@H](N)C(=O)O)."""

    SMILES = "C[C@H](N)C(=O)O"

    def test_stereocenter_count(self, solver):
        assert solver.get_stereocenter_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 3

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 3

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 2

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 2


class TestAdamantanePinnedValues:
    """Pinned values for adamantane (C1C2CC3CC1CC(C2)C3)."""

    SMILES = "C1C2CC3CC1CC(C2)C3"

    def test_ring_count(self, solver):
        count = solver.get_ring_count(self.SMILES)
        assert count == 4  # Adamantane has 4 rings in its cage structure

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 10

    def test_bridgehead_count(self, solver):
        count = solver.get_bridgehead_count(self.SMILES)
        assert count == 4  # Adamantane has 4 bridgehead carbons


class TestSpiroCompoundPinnedValues:
    """Pinned values for spiro[5.5]undecane (C1CCC2(CC1)CCCCC2)."""

    SMILES = "C1CCC2(CC1)CCCCC2"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 2

    def test_spiro_count(self, solver):
        assert solver.get_spiro_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 11


class TestChloroethanePinnedValues:
    """Pinned values for chloroethane (CCCl)."""

    SMILES = "CCCl"

    def test_halogen_atom_count(self, solver):
        assert solver.get_halogen_atom_count(self.SMILES) == 1

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 2

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 3


class TestAcetonePinnedValues:
    """Pinned values for acetone (CC(=O)C)."""

    SMILES = "CC(=O)C"

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 3

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 4

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 0

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 1

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 2


class TestCaffeinePinnedValues:
    """Pinned values for caffeine (Cn1cnc2c1c(=O)n(c(=O)n2C)C)."""

    SMILES = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 2

    def test_carbon_atom_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 8

    def test_hetero_atom_count(self, solver):
        count = solver.get_hetero_atom_count(self.SMILES)
        assert count == 6  # 4N + 2O


class TestMolecularFormulaPinned:
    """Tests for molecular formula consistency."""

    @pytest.mark.parametrize("smiles,expected_parts", [
        ("C", ["C", "H4"]),
        ("CCO", ["C2", "H6", "O"]),
        ("c1ccccc1", ["C6", "H6"]),
        ("CC(=O)O", ["C2", "H4", "O2"]),
        ("CCN", ["C2", "H7", "N"]),
    ])
    def test_molecular_formula_parts(self, solver, smiles, expected_parts):
        """Molecular formula contains expected elements."""
        formula = solver.get_molecular_formula(smiles)
        for part in expected_parts:
            assert part in formula, f"{part} not in {formula} for {smiles}"
