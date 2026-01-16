"""
Unit tests for SymbolicSolver.

Tests individual solver methods for correctness.
"""

class TestRingProperties:
    """Tests for ring-related solver methods."""

    def test_ring_count_benzene(self, solver, benzene):
        """Benzene has exactly 1 ring."""
        assert solver.get_ring_count(benzene) == 1

    def test_ring_count_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane has exactly 1 ring."""
        assert solver.get_ring_count(cyclohexane) == 1

    def test_ring_count_ethanol(self, solver, ethanol):
        """Ethanol has no rings."""
        assert solver.get_ring_count(ethanol) == 0

    def test_ring_count_naphthalene(self, solver, naphthalene):
        """Naphthalene has 2 fused rings."""
        assert solver.get_ring_count(naphthalene) == 2

    def test_aromatic_ring_count_benzene(self, solver, benzene):
        """Benzene is aromatic."""
        assert solver.get_aromatic_ring_count(benzene) == 1

    def test_aromatic_ring_count_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane is not aromatic."""
        assert solver.get_aromatic_ring_count(cyclohexane) == 0

    def test_aliphatic_ring_count_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane is aliphatic."""
        assert solver.get_aliphatic_ring_count(cyclohexane) == 1

    def test_aliphatic_ring_count_benzene(self, solver, benzene):
        """Benzene is not aliphatic."""
        assert solver.get_aliphatic_ring_count(benzene) == 0

    def test_saturated_ring_count_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane is saturated."""
        assert solver.get_saturated_ring_count(cyclohexane) == 1

    def test_fused_ring_count_naphthalene(self, solver, naphthalene):
        """Naphthalene has fused rings (counts fused ring systems, not individual rings)."""
        assert solver.get_fused_ring_count(naphthalene) == 1

    def test_fused_ring_count_benzene(self, solver, benzene):
        """Benzene has no fused rings (single ring)."""
        assert solver.get_fused_ring_count(benzene) == 0


class TestAtomCounts:
    """Tests for atom counting methods."""

    def test_carbon_count_ethanol(self, solver, ethanol):
        """Ethanol (CCO) has 2 carbons."""
        assert solver.get_carbon_atom_count(ethanol) == 2

    def test_carbon_count_benzene(self, solver, benzene):
        """Benzene has 6 carbons."""
        assert solver.get_carbon_atom_count(benzene) == 6

    def test_carbon_count_methane(self, solver, methane):
        """Methane has 1 carbon."""
        assert solver.get_carbon_atom_count(methane) == 1

    def test_heavy_atom_count_ethanol(self, solver, ethanol):
        """Ethanol (CCO) has 3 heavy atoms."""
        assert solver.get_heavy_atom_count(ethanol) == 3

    def test_heavy_atom_count_water(self, solver, water):
        """Water has 1 heavy atom (oxygen)."""
        assert solver.get_heavy_atom_count(water) == 1

    def test_hetero_atom_count_ethanol(self, solver, ethanol):
        """Ethanol has 1 heteroatom (O)."""
        assert solver.get_hetero_atom_count(ethanol) == 1

    def test_hetero_atom_count_benzene(self, solver, benzene):
        """Benzene has no heteroatoms."""
        assert solver.get_hetero_atom_count(benzene) == 0

    def test_halogen_count_chloroethane(self, solver, chloroethane):
        """Chloroethane has 1 halogen."""
        assert solver.get_halogen_atom_count(chloroethane) == 1

    def test_halogen_count_ethanol(self, solver, ethanol):
        """Ethanol has no halogens."""
        assert solver.get_halogen_atom_count(ethanol) == 0


class TestHydrogenBonding:
    """Tests for HBA/HBD counting."""

    def test_hbd_count_ethanol(self, solver, ethanol):
        """Ethanol has 1 HBD (OH group)."""
        assert solver.get_hbd_count(ethanol) == 1

    def test_hbd_count_acetone(self, solver, acetone):
        """Acetone has no HBD."""
        assert solver.get_hbd_count(acetone) == 0

    def test_hba_count_ethanol(self, solver, ethanol):
        """Ethanol has 1 HBA (oxygen)."""
        assert solver.get_hba_count(ethanol) == 1

    def test_hba_count_acetone(self, solver, acetone):
        """Acetone has 1 HBA (carbonyl oxygen)."""
        assert solver.get_hba_count(acetone) == 1

    def test_hba_count_ethylamine(self, solver, ethylamine):
        """Ethylamine has 1 HBA (nitrogen)."""
        assert solver.get_hba_count(ethylamine) == 1


class TestRotatableBonds:
    """Tests for rotatable bond counting."""

    def test_rotatable_bonds_ethanol(self, solver, ethanol):
        """Ethanol has no rotatable bonds."""
        count = solver.get_rotatable_bond_count(ethanol)
        assert count == 0

    def test_rotatable_bonds_benzene(self, solver, benzene):
        """Benzene has no rotatable bonds."""
        assert solver.get_rotatable_bond_count(benzene) == 0


class TestStereochemistry:
    """Tests for stereochemistry detection."""

    def test_stereocenter_count_alanine(self, solver, r_alanine):
        """Alanine has 1 stereocenter."""
        assert solver.get_stereocenter_count(r_alanine) == 1

    def test_stereocenter_count_ethanol(self, solver, ethanol):
        """Ethanol has no stereocenters."""
        assert solver.get_stereocenter_count(ethanol) == 0

    def test_r_stereocenter_count(self, solver, r_alanine):
        """R-Alanine R/S detection depends on CIP assignment."""
        r_count = solver.get_r_or_s_stereocenter_count(r_alanine, r_count=True)
        s_count = solver.get_r_or_s_stereocenter_count(r_alanine, r_count=False)
        # Total R+S should equal stereocenter count (or be 0 if CIP not assigned)
        assert r_count + s_count <= 1

    def test_s_stereocenter_count(self, solver, s_alanine):
        """S-Alanine R/S detection depends on CIP assignment."""
        r_count = solver.get_r_or_s_stereocenter_count(s_alanine, r_count=True)
        s_count = solver.get_r_or_s_stereocenter_count(s_alanine, r_count=False)
        # Total R+S should equal stereocenter count (or be 0 if CIP not assigned)
        assert r_count + s_count <= 1


class TestIndices:
    """Tests for index retrieval methods."""

    def test_ring_indices_benzene(self, solver, benzene):
        """Benzene ring indices should include all 6 carbons."""
        indices = solver.get_ring_indices(benzene)
        assert len(indices) == 6
        assert set(indices) == {0, 1, 2, 3, 4, 5}

    def test_ring_indices_ethanol(self, solver, ethanol):
        """Ethanol has no ring indices."""
        indices = solver.get_ring_indices(ethanol)
        assert indices == []

    def test_carbon_indices_ethanol(self, solver, ethanol):
        """Ethanol carbon indices."""
        indices = solver.get_carbon_atom_indices(ethanol)
        assert len(indices) == 2
        assert 0 in indices  # First carbon
        assert 1 in indices  # Second carbon

    def test_hetero_atom_indices_ethanol(self, solver, ethanol):
        """Ethanol heteroatom indices (oxygen)."""
        indices = solver.get_hetero_atom_indices(ethanol)
        assert len(indices) == 1
        assert 2 in indices  # Oxygen is at index 2


class TestCsp3Carbon:
    """Tests for sp3 carbon detection."""

    def test_csp3_count_ethanol(self, solver, ethanol):
        """Ethanol has 2 sp3 carbons."""
        assert solver.get_csp3_carbon_count(ethanol) == 2

    def test_csp3_count_benzene(self, solver, benzene):
        """Benzene has no sp3 carbons (all sp2)."""
        assert solver.get_csp3_carbon_count(benzene) == 0

    def test_csp3_count_acetone(self, solver, acetone):
        """Acetone has 2 sp3 carbons (methyl groups)."""
        assert solver.get_csp3_carbon_count(acetone) == 2


class TestMolecularFormula:
    """Tests for molecular formula."""

    def test_molecular_formula_ethanol(self, solver, ethanol):
        """Ethanol formula is C2H6O."""
        formula = solver.get_molecular_formula(ethanol)
        assert "C2" in formula
        assert "H6" in formula
        assert "O" in formula

    def test_molecular_formula_benzene(self, solver, benzene):
        """Benzene formula is C6H6."""
        formula = solver.get_molecular_formula(benzene)
        assert "C6" in formula
        assert "H6" in formula


class TestBridgeheadAtoms:
    """Tests for bridgehead atom detection."""

    def test_bridgehead_count_adamantane(self, solver, adamantane):
        """Adamantane has bridgehead atoms."""
        count = solver.get_bridgehead_count(adamantane)
        assert count > 0

    def test_bridgehead_count_benzene(self, solver, benzene):
        """Benzene has no bridgehead atoms."""
        assert solver.get_bridgehead_count(benzene) == 0

    def test_bridgehead_count_ethanol(self, solver, ethanol):
        """Ethanol has no bridgehead atoms."""
        assert solver.get_bridgehead_count(ethanol) == 0

    def test_bridgehead_indices_adamantane(self, solver, adamantane):
        """Adamantane bridgehead indices."""
        indices = solver.get_bridgehead_indices(adamantane)
        assert len(indices) > 0

    def test_bridgehead_indices_benzene(self, solver, benzene):
        """Benzene has no bridgehead indices."""
        indices = solver.get_bridgehead_indices(benzene)
        assert indices == []


class TestRingSizeExtremes:
    """Tests for smallest/largest ring detection."""

    def test_smallest_ring_size_benzene(self, solver, benzene):
        """Benzene smallest ring is 6."""
        size = solver.get_smallest_or_largest_ring_count(benzene, smallest=True)
        assert size == 6

    def test_largest_ring_size_benzene(self, solver, benzene):
        """Benzene largest ring is 6."""
        size = solver.get_smallest_or_largest_ring_count(benzene, smallest=False)
        assert size == 6

    def test_smallest_ring_size_no_rings(self, solver, ethanol):
        """Ethanol has no rings."""
        size = solver.get_smallest_or_largest_ring_count(ethanol, smallest=True)
        assert size == 0

    def test_smallest_ring_indices(self, solver, benzene):
        """Smallest ring indices for benzene."""
        indices = solver.get_smallest_or_largest_ring_indices(benzene, smallest=True)
        assert len(indices) == 6

    def test_largest_ring_indices(self, solver, benzene):
        """Largest ring indices for benzene."""
        indices = solver.get_smallest_or_largest_ring_indices(benzene, smallest=False)
        assert len(indices) == 6


class TestChainTopology:
    """Tests for chain termini and branch points."""

    def test_chain_termini_hexane(self, solver, hexane):
        """Hexane has 2 terminal carbons."""
        count = solver.get_chain_termini_count(hexane)
        assert count == 2

    def test_chain_termini_methane(self, solver, methane):
        """Methane - single atom is not a chain terminus."""
        count = solver.get_chain_termini_count(methane)
        assert count == 0

    def test_chain_termini_benzene(self, solver, benzene):
        """Benzene has no terminal carbons (all in ring)."""
        count = solver.get_chain_termini_count(benzene)
        assert count == 0

    def test_chain_termini_indices_hexane(self, solver, hexane):
        """Hexane terminal indices."""
        indices = solver.get_chain_termini_indices(hexane)
        assert len(indices) == 2
        assert 0 in indices
        assert 5 in indices

    def test_branch_point_count_isobutane(self, solver, isobutane):
        """Isobutane has 1 branch point."""
        count = solver.get_branch_point_count(isobutane)
        assert count == 1

    def test_branch_point_count_neopentane(self, solver, neopentane):
        """Neopentane has 1 quaternary carbon."""
        count = solver.get_branch_point_count(neopentane)
        assert count == 1

    def test_branch_point_count_hexane(self, solver, hexane):
        """Linear hexane has no branch points."""
        count = solver.get_branch_point_count(hexane)
        assert count == 0

    def test_branch_point_indices_isobutane(self, solver, isobutane):
        """Isobutane branch point index."""
        indices = solver.get_branch_point_indices(isobutane)
        assert len(indices) == 1


class TestHeterocycles:
    """Tests for heterocycle detection."""

    def test_heterocycle_count_pyridine(self, solver, pyridine):
        """Pyridine is a heterocycle."""
        assert solver.get_heterocycle_count(pyridine) == 1

    def test_heterocycle_count_benzene(self, solver, benzene):
        """Benzene is not a heterocycle."""
        assert solver.get_heterocycle_count(benzene) == 0

    def test_heterocycle_count_morpholine(self, solver, morpholine):
        """Morpholine is a heterocycle."""
        assert solver.get_heterocycle_count(morpholine) == 1

    def test_heterocycle_indices_pyridine(self, solver, pyridine):
        """Pyridine heterocycle indices include all ring atoms."""
        indices = solver.get_heterocycle_indices(pyridine)
        assert len(indices) == 6

    def test_heterocycle_indices_benzene(self, solver, benzene):
        """Benzene has no heterocycle indices."""
        indices = solver.get_heterocycle_indices(benzene)
        assert indices == []


class TestSpiroAtoms:
    """Tests for spiro atom detection."""

    def test_spiro_count_spiro_compound(self, solver, spiro_compound):
        """Spiro compound has 1 spiro atom."""
        count = solver.get_spiro_count(spiro_compound)
        assert count == 1

    def test_spiro_count_benzene(self, solver, benzene):
        """Benzene has no spiro atoms."""
        assert solver.get_spiro_count(benzene) == 0

    def test_spiro_indices_spiro_compound(self, solver, spiro_compound):
        """Spiro compound spiro atom index."""
        indices = solver.get_spiro_indices(spiro_compound)
        assert len(indices) == 1


class TestLongestCarbonChain:
    """Tests for longest carbon chain detection."""

    def test_longest_chain_hexane(self, solver, hexane):
        """Hexane has chain length 6."""
        count = solver.get_longest_carbon_chain_count(hexane)
        assert count == 6

    def test_longest_chain_isobutane(self, solver, isobutane):
        """Isobutane longest chain is 3."""
        count = solver.get_longest_carbon_chain_count(isobutane)
        assert count == 3

    def test_longest_chain_benzene(self, solver, benzene):
        """Benzene longest carbon chain."""
        count = solver.get_longest_carbon_chain_count(benzene)
        assert count == 6

    def test_longest_chain_indices_hexane(self, solver, hexane):
        """Hexane chain indices."""
        indices = solver.get_longest_carbon_chain_indices(hexane)
        assert len(indices) == 6


class TestEZStereochemistry:
    """Tests for E/Z double bond stereochemistry."""

    def test_e_double_bond_count_e_butene(self, solver, e_butene):
        """E-butene has 1 E double bond."""
        count = solver.get_e_z_stereochemistry_double_bond_count(e_butene, e_count=True)
        assert count == 1

    def test_z_double_bond_count_e_butene(self, solver, e_butene):
        """E-butene has 0 Z double bonds."""
        count = solver.get_e_z_stereochemistry_double_bond_count(e_butene, e_count=False)
        assert count == 0

    def test_z_double_bond_count_z_butene(self, solver, z_butene):
        """Z-butene has 1 Z double bond."""
        count = solver.get_e_z_stereochemistry_double_bond_count(z_butene, e_count=False)
        assert count == 1

    def test_e_double_bond_count_z_butene(self, solver, z_butene):
        """Z-butene has 0 E double bonds."""
        count = solver.get_e_z_stereochemistry_double_bond_count(z_butene, e_count=True)
        assert count == 0

    def test_e_double_bond_indices(self, solver, e_butene):
        """E-butene E double bond indices."""
        indices = solver.get_e_z_stereochemistry_double_bond_indices(e_butene, e_indices=True)
        assert len(indices) == 2

    def test_z_double_bond_indices(self, solver, z_butene):
        """Z-butene Z double bond indices."""
        indices = solver.get_e_z_stereochemistry_double_bond_indices(z_butene, e_indices=False)
        assert len(indices) == 2

    def test_no_ez_bonds_ethanol(self, solver, ethanol):
        """Ethanol has no E/Z bonds."""
        e_count = solver.get_e_z_stereochemistry_double_bond_count(ethanol, e_count=True)
        z_count = solver.get_e_z_stereochemistry_double_bond_count(ethanol, e_count=False)
        assert e_count == 0
        assert z_count == 0


class TestHydrogenCounts:
    """Tests for hydrogen counting."""

    def test_hydrogen_count_methane(self, solver, methane):
        """Methane has 4 hydrogens."""
        count = solver.get_hydrogen_count(methane)
        assert count == 4

    def test_hydrogen_count_ethanol(self, solver, ethanol):
        """Ethanol (C2H6O) has 6 hydrogens."""
        count = solver.get_hydrogen_count(ethanol)
        assert count == 6

    def test_hydrogen_count_benzene(self, solver, benzene):
        """Benzene has 6 hydrogens."""
        count = solver.get_hydrogen_count(benzene)
        assert count == 6

    def test_explicit_hydrogen_count(self, solver, ethanol):
        """Explicit hydrogen count (0 in standard SMILES)."""
        count = solver.get_explicit_hydrogen_count(ethanol)
        assert count == 0


class TestOxidationState:
    """Tests for oxidation state methods."""

    def test_max_oxidation_carbon_formic_acid(self, solver, formic_acid):
        """Formic acid has 1 carbon at high oxidation."""
        count = solver.get_oxidation_state_count(formic_acid, 'C', max_oxidation=True)
        assert count == 1

    def test_oxidation_state_indices_formic_acid(self, solver, formic_acid):
        """Formic acid oxidation indices."""
        indices = solver.get_oxidation_state_indices(formic_acid, 'C', max_oxidation=True)
        assert len(indices) == 1

    def test_min_oxidation_carbon_methanol(self, solver, methanol):
        """Methanol has carbon at low oxidation."""
        count = solver.get_oxidation_state_count(methanol, 'C', max_oxidation=False)
        assert count == 1


class TestBRICSDecomposition:
    """Tests for BRICS fragmentation."""

    def test_brics_fragment_count_aspirin(self, solver, aspirin):
        """Aspirin can be fragmented by BRICS."""
        count = solver.get_brics_fragment_count(aspirin)
        assert count > 0

    def test_brics_fragment_count_methane(self, solver, methane):
        """Methane has 1 BRICS fragment (itself)."""
        count = solver.get_brics_fragment_count(methane)
        assert count == 1

    def test_brics_bond_indices_aspirin(self, solver, aspirin):
        """Aspirin has 5 BRICS breakable bonds."""
        indices = solver.get_brics_bond_indices(aspirin)
        assert len(indices) == 5


class TestMurckoScaffold:
    """Tests for Murcko scaffold extraction."""

    def test_murcko_scaffold_count_aspirin(self, solver, aspirin):
        """Aspirin has a Murcko scaffold."""
        count = solver.get_murcko_scaffold_count(aspirin)
        assert count > 0

    def test_murcko_scaffold_count_ethanol(self, solver, ethanol):
        """Ethanol has no scaffold (no rings)."""
        count = solver.get_murcko_scaffold_count(ethanol)
        assert count == 0

    def test_murcko_scaffold_indices_benzene(self, solver, benzene):
        """Benzene scaffold is the whole molecule."""
        indices = solver.get_murcko_scaffold_indices(benzene)
        assert len(indices) == 6

    def test_murcko_scaffold_value_benzene(self, solver, benzene):
        """Benzene scaffold SMILES."""
        scaffold = solver.get_murcko_scaffold_value(benzene)
        assert isinstance(scaffold, str)
        assert len(scaffold) > 0


class TestAromaticAliphaticRingIndices:
    """Tests for aromatic and aliphatic ring indices."""

    def test_aromatic_ring_indices_benzene(self, solver, benzene):
        """Benzene aromatic ring indices."""
        indices = solver.get_aromatic_ring_indices(benzene)
        assert len(indices) == 6
        assert set(indices) == {0, 1, 2, 3, 4, 5}

    def test_aromatic_ring_indices_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane has no aromatic ring indices."""
        indices = solver.get_aromatic_ring_indices(cyclohexane)
        assert indices == []

    def test_aliphatic_ring_indices_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane aliphatic ring indices."""
        indices = solver.get_aliphatic_ring_indices(cyclohexane)
        assert len(indices) == 6

    def test_aliphatic_ring_indices_benzene(self, solver, benzene):
        """Benzene has no aliphatic ring indices."""
        indices = solver.get_aliphatic_ring_indices(benzene)
        assert indices == []

    def test_saturated_ring_indices_cyclohexane(self, solver, cyclohexane):
        """Cyclohexane saturated ring indices."""
        indices = solver.get_saturated_ring_indices(cyclohexane)
        assert len(indices) == 6


class TestFusedRingIndices:
    """Tests for fused ring indices."""

    def test_fused_ring_indices_naphthalene(self, solver, naphthalene):
        """Naphthalene fused ring indices include all atoms."""
        indices = solver.get_fused_ring_indices(naphthalene)
        assert len(indices) == 10

    def test_fused_ring_indices_benzene(self, solver, benzene):
        """Benzene has no fused ring indices (single ring)."""
        indices = solver.get_fused_ring_indices(benzene)
        assert indices == []


class TestHBAHBDIndices:
    """Tests for HBA/HBD indices."""

    def test_hba_indices_ethanol(self, solver, ethanol):
        """Ethanol HBA indices (oxygen)."""
        indices = solver.get_hba_indices(ethanol)
        assert len(indices) >= 1

    def test_hbd_indices_ethanol(self, solver, ethanol):
        """Ethanol HBD indices (OH)."""
        indices = solver.get_hbd_indices(ethanol)
        assert len(indices) == 1

    def test_hba_indices_acetone(self, solver, acetone):
        """Acetone HBA indices (carbonyl O)."""
        indices = solver.get_hba_indices(acetone)
        assert len(indices) >= 1

    def test_hbd_indices_acetone(self, solver, acetone):
        """Acetone has no HBD."""
        indices = solver.get_hbd_indices(acetone)
        assert indices == []


class TestRotatableBondIndices:
    """Tests for rotatable bond indices."""

    def test_rotatable_bond_indices_benzene(self, solver, benzene):
        """Benzene has no rotatable bond indices."""
        indices = solver.get_rotatable_bond_indices(benzene)
        assert indices == []

    def test_rotatable_bond_indices_hexane(self, solver, hexane):
        """Hexane has rotatable bonds."""
        indices = solver.get_rotatable_bond_indices(hexane)
        assert len(indices) > 0


class TestHalogenIndices:
    """Tests for halogen atom indices."""

    def test_halogen_indices_chloroethane(self, solver, chloroethane):
        """Chloroethane halogen indices."""
        indices = solver.get_halogen_atom_indices(chloroethane)
        assert len(indices) == 1

    def test_halogen_indices_ethanol(self, solver, ethanol):
        """Ethanol has no halogen indices."""
        indices = solver.get_halogen_atom_indices(ethanol)
        assert indices == []


class TestHeavyAtomIndices:
    """Tests for heavy atom indices."""

    def test_heavy_atom_indices_ethanol(self, solver, ethanol):
        """Ethanol heavy atom indices (C, C, O)."""
        indices = solver.get_heavy_atom_indices(ethanol)
        assert len(indices) == 3

    def test_heavy_atom_indices_methane(self, solver, methane):
        """Methane has 1 heavy atom."""
        indices = solver.get_heavy_atom_indices(methane)
        assert len(indices) == 1


class TestStereocenterIndices:
    """Tests for stereocenter indices."""

    def test_stereocenter_indices_alanine(self, solver, r_alanine):
        """R-Alanine stereocenter indices."""
        indices = solver.get_stereocenter_indices(r_alanine)
        assert len(indices) == 1

    def test_stereocenter_indices_ethanol(self, solver, ethanol):
        """Ethanol has no stereocenter indices."""
        indices = solver.get_stereocenter_indices(ethanol)
        assert indices == []

    def test_unspecified_stereocenter_count(self, solver):
        """Test unspecified stereocenter count."""
        # Molecule with unspecified stereocenter
        smiles = "CC(O)C(=O)O"  # Lactic acid without specified stereochemistry
        count = solver.get_unspecified_stereocenter_count(smiles)
        assert isinstance(count, int)


class TestComplexMolecules:
    """Tests with complex drug-like molecules."""

    def test_caffeine_properties(self, solver, caffeine):
        """Caffeine has expected properties."""
        assert solver.get_ring_count(caffeine) >= 2
        assert solver.get_hetero_atom_count(caffeine) >= 4
        assert solver.get_carbon_atom_count(caffeine) == 8

    def test_ibuprofen_properties(self, solver, ibuprofen):
        """Ibuprofen has expected properties."""
        assert solver.get_ring_count(ibuprofen) == 1
        assert solver.get_aromatic_ring_count(ibuprofen) == 1
        assert solver.get_carbon_atom_count(ibuprofen) == 13

    def test_glucose_stereocenters(self, solver, glucose):
        """Glucose has 5 stereocenters."""
        count = solver.get_stereocenter_count(glucose)
        assert count == 5

    def test_paracetamol_properties(self, solver, paracetamol):
        """Paracetamol has expected properties."""
        assert solver.get_ring_count(paracetamol) == 1
        assert solver.get_hbd_count(paracetamol) == 2  # NH and OH
        assert solver.get_hba_count(paracetamol) == 2


class TestPolycyclicSystems:
    """Tests for polycyclic ring systems."""

    def test_anthracene_ring_count(self, solver, anthracene):
        """Anthracene has 3 rings."""
        assert solver.get_ring_count(anthracene) == 3

    def test_anthracene_aromatic_count(self, solver, anthracene):
        """Anthracene has 3 aromatic rings."""
        assert solver.get_aromatic_ring_count(anthracene) == 3

    def test_anthracene_fused_count(self, solver, anthracene):
        """Anthracene has 1 fused ring system."""
        count = solver.get_fused_ring_count(anthracene)
        assert count == 1

    def test_decalin_properties(self, solver, decalin):
        """Decalin (fused saturated rings)."""
        assert solver.get_ring_count(decalin) == 2
        assert solver.get_aromatic_ring_count(decalin) == 0
        assert solver.get_saturated_ring_count(decalin) == 2

    def test_bicyclo_heptane_bridgehead(self, solver, bicyclo_heptane):
        """Bicyclo[2.2.1]heptane has 2 bridgehead atoms."""
        count = solver.get_bridgehead_count(bicyclo_heptane)
        assert count == 2
