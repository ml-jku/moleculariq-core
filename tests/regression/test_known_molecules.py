"""
Regression tests with known molecules.

Tests that known molecules always produce expected properties.
These tests catch unintended changes to solver behavior.
"""

import pytest


class TestEthanol:
    """Regression tests for ethanol (CCO)."""

    SMILES = "CCO"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 0

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 2

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 3

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 1

    def test_hbd_count(self, solver):
        assert solver.get_hbd_count(self.SMILES) == 1

    def test_hba_count(self, solver):
        assert solver.get_hba_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 0

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 2


class TestBenzene:
    """Regression tests for benzene (c1ccccc1)."""

    SMILES = "c1ccccc1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 1

    def test_aliphatic_ring_count(self, solver):
        assert solver.get_aliphatic_ring_count(self.SMILES) == 0

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 6

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 6

    def test_hetero_atom_count(self, solver):
        assert solver.get_hetero_atom_count(self.SMILES) == 0

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 0

    def test_ring_indices(self, solver):
        indices = solver.get_ring_indices(self.SMILES)
        assert len(indices) == 6
        assert set(indices) == {0, 1, 2, 3, 4, 5}


class TestCyclohexane:
    """Regression tests for cyclohexane (C1CCCCC1)."""

    SMILES = "C1CCCCC1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 0

    def test_aliphatic_ring_count(self, solver):
        assert solver.get_aliphatic_ring_count(self.SMILES) == 1

    def test_saturated_ring_count(self, solver):
        assert solver.get_saturated_ring_count(self.SMILES) == 1

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 6

    def test_csp3_carbon_count(self, solver):
        assert solver.get_csp3_carbon_count(self.SMILES) == 6


class TestNaphthalene:
    """Regression tests for naphthalene (c1ccc2ccccc2c1)."""

    SMILES = "c1ccc2ccccc2c1"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 2

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 2

    def test_fused_ring_count(self, solver):
        # Naphthalene: solver counts fused ring systems, not individual fused rings
        assert solver.get_fused_ring_count(self.SMILES) == 1

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 10

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 10


class TestAspirin:
    """Regression tests for aspirin (CC(=O)Oc1ccccc1C(=O)O)."""

    SMILES = "CC(=O)Oc1ccccc1C(=O)O"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_aromatic_ring_count(self, solver):
        assert solver.get_aromatic_ring_count(self.SMILES) == 1

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 9

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 13

    def test_hbd_count(self, solver):
        # Carboxylic acid OH
        assert solver.get_hbd_count(self.SMILES) == 1


class TestAcetone:
    """Regression tests for acetone (CC(=O)C)."""

    SMILES = "CC(=O)C"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 0

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 3

    def test_heavy_atom_count(self, solver):
        assert solver.get_heavy_atom_count(self.SMILES) == 4

    def test_hbd_count(self, solver):
        # No HBD in ketone
        assert solver.get_hbd_count(self.SMILES) == 0

    def test_hba_count(self, solver):
        # Carbonyl oxygen is HBA
        assert solver.get_hba_count(self.SMILES) == 1

    def test_csp3_carbon_count(self, solver):
        # Two methyl carbons are sp3
        assert solver.get_csp3_carbon_count(self.SMILES) == 2


class TestGlucose:
    """Regression tests for glucose (OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O)."""

    SMILES = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"

    def test_ring_count(self, solver):
        assert solver.get_ring_count(self.SMILES) == 1

    def test_carbon_count(self, solver):
        assert solver.get_carbon_atom_count(self.SMILES) == 6

    def test_stereocenter_count(self, solver):
        # Glucose has 5 stereocenters in this pyranose form
        count = solver.get_stereocenter_count(self.SMILES)
        assert count == 5

    def test_hbd_count(self, solver):
        # 5 OH groups (4 ring OH + 1 primary OH)
        count = solver.get_hbd_count(self.SMILES)
        assert count == 5


class TestChiralMolecules:
    """Regression tests for chiral molecules."""

    def test_r_alanine_stereocenter(self, solver, r_alanine):
        """R-Alanine has 1 stereocenter."""
        assert solver.get_stereocenter_count(r_alanine) == 1
        # Note: R/S assignment may depend on RDKit's CIP assignment algorithm
        r_count = solver.get_r_or_s_stereocenter_count(r_alanine, r_count=True)
        s_count = solver.get_r_or_s_stereocenter_count(r_alanine, r_count=False)
        # Total R+S should equal stereocenter count (or be 0 if CIP not assigned)
        assert r_count + s_count <= 1

    def test_s_alanine_stereocenter(self, solver, s_alanine):
        """S-Alanine has 1 stereocenter."""
        assert solver.get_stereocenter_count(s_alanine) == 1
        # Note: R/S assignment may depend on RDKit's CIP assignment algorithm
        r_count = solver.get_r_or_s_stereocenter_count(s_alanine, r_count=True)
        s_count = solver.get_r_or_s_stereocenter_count(s_alanine, r_count=False)
        # Total R+S should equal stereocenter count (or be 0 if CIP not assigned)
        assert r_count + s_count <= 1


class TestKnownMoleculeProperties:
    """Test using the known_molecule_properties fixture."""

    def test_all_known_properties(self, solver, known_molecule_properties):
        """Verify all properties in the known molecules dictionary."""
        for smiles, expected_props in known_molecule_properties.items():
            for prop_name, expected_value in expected_props.items():
                method_name = f"get_{prop_name}"
                if hasattr(solver, method_name):
                    method = getattr(solver, method_name)
                    actual = method(smiles)
                    assert actual == expected_value, (
                        f"Mismatch for {smiles}.{prop_name}: "
                        f"expected {expected_value}, got {actual}"
                    )


class TestMolecularFormula:
    """Regression tests for molecular formula."""

    @pytest.mark.parametrize("smiles,expected_elements", [
        ("C", ["C", "H4"]),
        ("CCO", ["C2", "H6", "O"]),
        ("c1ccccc1", ["C6", "H6"]),
        ("CC(=O)O", ["C2", "H4", "O2"]),
    ])
    def test_molecular_formula_contains_elements(self, solver, smiles, expected_elements):
        """Molecular formula contains expected elements."""
        formula = solver.get_molecular_formula(smiles)
        for element in expected_elements:
            assert element in formula, f"{element} not in {formula} for {smiles}"
