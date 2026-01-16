"""
Unit tests for TemplateBasedReactionSolver.

Tests SMIRKS-based reaction template evaluation.
"""

import pytest
from moleculariq_core import TemplateBasedReactionSolver


@pytest.fixture(scope="module")
def reaction_solver():
    """Create a TemplateBasedReactionSolver instance."""
    return TemplateBasedReactionSolver()


class TestReactionSolverInitialization:
    """Tests for reaction solver initialization."""

    def test_initialization(self, reaction_solver):
        """Solver initializes correctly."""
        assert reaction_solver is not None

    def test_templates_loaded(self, reaction_solver):
        """Reaction templates are loaded."""
        assert len(reaction_solver.templates) > 0

    def test_templates_compiled(self, reaction_solver):
        """Reaction templates are compiled."""
        assert len(reaction_solver.compiled_templates) > 0

    def test_center_patterns_compiled(self, reaction_solver):
        """Reaction center patterns are compiled."""
        assert len(reaction_solver.compiled_center_patterns) > 0


class TestGetReactionData:
    """Tests for get_reaction_data method."""

    def test_returns_dict(self, reaction_solver):
        """get_reaction_data returns a dictionary."""
        result = reaction_solver.get_reaction_data("CCO")
        assert isinstance(result, dict)

    def test_result_has_keys(self, reaction_solver):
        """Result contains expected keys."""
        result = reaction_solver.get_reaction_data("CCO")
        # Should have success, products, count, and index keys
        has_success = any("_success" in k for k in result.keys())
        has_products = any("_products" in k for k in result.keys())
        has_count = any("_count" in k for k in result.keys())
        has_index = any("_index" in k for k in result.keys())
        assert has_success
        assert has_products
        assert has_count
        assert has_index


class TestOxidationReactions:
    """Tests for oxidation reactions."""

    def test_primary_alcohol_to_aldehyde(self, reaction_solver):
        """Primary alcohol oxidation to aldehyde."""
        # Ethanol (CCO) has a primary alcohol
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_primary_alcohol_to_aldehyde_success"
        assert key in result
        assert result[key] == 1  # Should succeed

    def test_primary_alcohol_to_aldehyde_products(self, reaction_solver):
        """Primary alcohol oxidation produces aldehyde."""
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_primary_alcohol_to_aldehyde_products"
        assert key in result
        assert result[key] is not None
        # Product should be acetaldehyde-like

    def test_secondary_alcohol_to_ketone(self, reaction_solver):
        """Secondary alcohol oxidation to ketone."""
        # Isopropanol (CC(O)C) has a secondary alcohol
        result = reaction_solver.get_reaction_data("CC(O)C")
        key = "template_based_reaction_prediction_secondary_alcohol_to_ketone_success"
        assert key in result
        assert result[key] == 1  # Should succeed

    def test_aldehyde_to_carboxylic_acid(self, reaction_solver):
        """Aldehyde oxidation to carboxylic acid."""
        # Acetaldehyde (CC=O)
        result = reaction_solver.get_reaction_data("CC=O")
        key = "template_based_reaction_prediction_aldehyde_to_carboxylic_acid_success"
        assert key in result
        assert result[key] == 1


class TestReductionReactions:
    """Tests for reduction reactions."""

    def test_ketone_to_secondary_alcohol(self, reaction_solver):
        """Ketone reduction to secondary alcohol."""
        # Acetone (CC(=O)C)
        result = reaction_solver.get_reaction_data("CC(=O)C")
        key = "template_based_reaction_prediction_ketone_to_secondary_alcohol_success"
        assert key in result
        assert result[key] == 1

    def test_aldehyde_to_primary_alcohol(self, reaction_solver):
        """Aldehyde reduction to primary alcohol."""
        # Acetaldehyde (CC=O)
        result = reaction_solver.get_reaction_data("CC=O")
        key = "template_based_reaction_prediction_aldehyde_to_primary_alcohol_success"
        assert key in result
        assert result[key] == 1

    def test_nitro_to_amine(self, reaction_solver):
        """Nitro reduction to amine."""
        # Nitrobenzene
        result = reaction_solver.get_reaction_data("c1ccc(cc1)[N+](=O)[O-]")
        key = "template_based_reaction_prediction_nitro_to_amine_success"
        assert key in result
        assert result[key] == 1

    def test_alkene_to_alkane(self, reaction_solver):
        """Alkene hydrogenation to alkane."""
        # Ethene (C=C)
        result = reaction_solver.get_reaction_data("C=C")
        key = "template_based_reaction_prediction_alkene_to_alkane_success"
        assert key in result
        assert result[key] == 1


class TestSubstitutionReactions:
    """Tests for substitution reactions."""

    def test_alcohol_to_alkyl_chloride(self, reaction_solver):
        """Alcohol to alkyl chloride."""
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_alcohol_to_alkyl_halide_Cl_success"
        assert key in result
        assert result[key] == 1

    def test_alcohol_to_alkyl_bromide(self, reaction_solver):
        """Alcohol to alkyl bromide."""
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_alcohol_to_alkyl_halide_Br_success"
        assert key in result
        assert result[key] == 1

    def test_alkyl_halide_to_alcohol(self, reaction_solver):
        """Alkyl halide hydrolysis to alcohol."""
        # Chloroethane
        result = reaction_solver.get_reaction_data("CCCl")
        key = "template_based_reaction_prediction_alkyl_halide_to_alcohol_success"
        assert key in result
        assert result[key] == 1


class TestHydrolysisReactions:
    """Tests for hydrolysis reactions."""

    def test_ester_hydrolysis(self, reaction_solver):
        """Ester hydrolysis to carboxylic acid."""
        # Methyl acetate (CC(=O)OC)
        result = reaction_solver.get_reaction_data("CC(=O)OC")
        key = "template_based_reaction_prediction_ester_hydrolysis_to_acid_success"
        assert key in result
        assert result[key] == 1

    def test_amide_hydrolysis(self, reaction_solver):
        """Amide hydrolysis to carboxylic acid."""
        # Acetamide (CC(=O)N)
        result = reaction_solver.get_reaction_data("CC(=O)NC")
        key = "template_based_reaction_prediction_amide_hydrolysis_to_acid_success"
        assert key in result
        assert result[key] == 1

    def test_nitrile_hydrolysis(self, reaction_solver):
        """Nitrile hydrolysis."""
        # Acetonitrile (CC#N)
        result = reaction_solver.get_reaction_data("CC#N")
        key = "template_based_reaction_prediction_nitrile_hydrolysis_success"
        assert key in result
        assert result[key] == 1


class TestAdditionReactions:
    """Tests for addition reactions."""

    def test_hydration_of_alkene(self, reaction_solver):
        """Alkene hydration."""
        result = reaction_solver.get_reaction_data("C=C")
        key = "template_based_reaction_prediction_hydration_of_alkene_success"
        assert key in result
        assert result[key] == 1

    def test_hydrohalogenation(self, reaction_solver):
        """Hydrohalogenation of alkene."""
        result = reaction_solver.get_reaction_data("C=C")
        key_cl = "template_based_reaction_prediction_hydrohalogenation_HCl_success"
        key_br = "template_based_reaction_prediction_hydrohalogenation_HBr_success"
        assert key_cl in result
        assert key_br in result
        assert result[key_cl] == 1
        assert result[key_br] == 1

    def test_epoxidation(self, reaction_solver):
        """Alkene epoxidation."""
        result = reaction_solver.get_reaction_data("C=C")
        key = "template_based_reaction_prediction_epoxidation_success"
        assert key in result
        assert result[key] == 1


class TestAromaticSubstitution:
    """Tests for aromatic substitution reactions."""

    def test_nitration(self, reaction_solver):
        """Benzene nitration."""
        result = reaction_solver.get_reaction_data("c1ccccc1")
        key = "template_based_reaction_prediction_nitration_success"
        assert key in result
        assert result[key] == 1

    def test_bromination(self, reaction_solver):
        """Benzene bromination."""
        result = reaction_solver.get_reaction_data("c1ccccc1")
        key = "template_based_reaction_prediction_bromination_success"
        assert key in result
        assert result[key] == 1


class TestProtectingGroups:
    """Tests for protecting group reactions."""

    def test_boc_protection(self, reaction_solver):
        """BOC protection of amine."""
        # Aniline (c1ccc(N)cc1)
        result = reaction_solver.get_reaction_data("c1ccc(N)cc1")
        key = "template_based_reaction_prediction_boc_protection_success"
        assert key in result
        assert result[key] == 1

    def test_acetylation(self, reaction_solver):
        """Acetylation of alcohol."""
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_acetylation_success"
        assert key in result
        assert result[key] == 1


class TestRingOpening:
    """Tests for ring opening reactions."""

    def test_epoxide_to_diol(self, reaction_solver):
        """Epoxide opening to diol."""
        # Ethylene oxide (C1OC1)
        result = reaction_solver.get_reaction_data("C1OC1")
        key = "template_based_reaction_prediction_epoxide_to_diol_success"
        assert key in result
        assert result[key] == 1


class TestReactionCenterCounts:
    """Tests for reaction center count/index methods."""

    def test_primary_alcohol_count(self, reaction_solver):
        """Count of primary alcohol reaction centers."""
        # 1,4-butanediol has two primary alcohols
        result = reaction_solver.get_reaction_data("OCCCCO")
        key = "template_based_reaction_prediction_primary_alcohol_to_aldehyde_count"
        assert key in result
        assert result[key] >= 1

    def test_primary_alcohol_indices(self, reaction_solver):
        """Indices of primary alcohol reaction centers."""
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_primary_alcohol_to_aldehyde_index"
        assert key in result
        assert isinstance(result[key], list)

    def test_aromatic_h_count(self, reaction_solver):
        """Count of aromatic H sites."""
        # Benzene has 6 aromatic H
        result = reaction_solver.get_reaction_data("c1ccccc1")
        key = "template_based_reaction_prediction_nitration_count"
        assert key in result
        assert result[key] == 6


class TestNoMatchCases:
    """Tests for molecules that don't match reaction patterns."""

    def test_no_alcohol_no_oxidation(self, reaction_solver):
        """Molecule without alcohol doesn't undergo alcohol oxidation."""
        # Hexane (CCCCCC) - no alcohol
        result = reaction_solver.get_reaction_data("CCCCCC")
        key = "template_based_reaction_prediction_primary_alcohol_to_aldehyde_success"
        assert result[key] == 0

    def test_no_alkene_no_hydrogenation(self, reaction_solver):
        """Molecule without alkene doesn't undergo hydrogenation."""
        # Ethanol (CCO) - no alkene
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_alkene_to_alkane_success"
        assert result[key] == 0

    def test_no_nitro_no_reduction(self, reaction_solver):
        """Molecule without nitro doesn't undergo nitro reduction."""
        result = reaction_solver.get_reaction_data("CCO")
        key = "template_based_reaction_prediction_nitro_to_amine_success"
        assert result[key] == 0


class TestInvalidInputs:
    """Tests for invalid inputs."""

    def test_invalid_smiles(self, reaction_solver):
        """Invalid SMILES returns zeros."""
        result = reaction_solver.get_reaction_data("invalid_xyz")
        # All counts should be 0
        count_keys = [k for k in result if "_count" in k]
        for key in count_keys:
            assert result[key] == 0

    def test_invalid_smiles_no_success(self, reaction_solver):
        """Invalid SMILES has no successful reactions."""
        result = reaction_solver.get_reaction_data("invalid_xyz")
        success_keys = [k for k in result if "_success" in k]
        for key in success_keys:
            assert result[key] == 0

    def test_invalid_smiles_empty_indices(self, reaction_solver):
        """Invalid SMILES has empty indices."""
        result = reaction_solver.get_reaction_data("invalid_xyz")
        index_keys = [k for k in result if "_index" in k]
        for key in index_keys:
            assert result[key] == []

    def test_empty_smiles(self, reaction_solver):
        """Empty SMILES."""
        result = reaction_solver.get_reaction_data("")
        assert isinstance(result, dict)


class TestProductStructure:
    """Tests for reaction product structure."""

    def test_products_are_list(self, reaction_solver):
        """Products are returned as list."""
        result = reaction_solver.get_reaction_data("CCO")
        product_keys = [k for k in result if "_products" in k and result[k] is not None]
        for key in product_keys:
            assert isinstance(result[key], list)

    def test_products_are_smiles(self, reaction_solver):
        """Products are valid SMILES strings."""
        from rdkit import Chem
        result = reaction_solver.get_reaction_data("CCO")
        product_keys = [k for k in result if "_products" in k and result[k] is not None]
        for key in product_keys:
            for product in result[key]:
                mol = Chem.MolFromSmiles(product)
                assert mol is not None, f"Invalid product SMILES: {product}"


class TestComplexMolecules:
    """Tests with complex molecules."""

    def test_aspirin(self, reaction_solver):
        """Aspirin (CC(=O)Oc1ccccc1C(=O)O)."""
        result = reaction_solver.get_reaction_data("CC(=O)Oc1ccccc1C(=O)O")
        assert isinstance(result, dict)
        # Aspirin has ester - should match ester hydrolysis
        key = "template_based_reaction_prediction_ester_hydrolysis_to_acid_success"
        assert result[key] == 1

    def test_paracetamol(self, reaction_solver):
        """Paracetamol (CC(=O)Nc1ccc(O)cc1)."""
        result = reaction_solver.get_reaction_data("CC(=O)Nc1ccc(O)cc1")
        assert isinstance(result, dict)
        # Has phenol - should match alcohol reactions
        key = "template_based_reaction_prediction_acetylation_success"
        assert result[key] == 1

    def test_caffeine(self, reaction_solver):
        """Caffeine."""
        caffeine = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
        result = reaction_solver.get_reaction_data(caffeine)
        assert isinstance(result, dict)


class TestMultipleReactionSites:
    """Tests for molecules with multiple reaction sites."""

    def test_diol_multiple_alcohols(self, reaction_solver):
        """Diol has multiple oxidizable alcohols."""
        # 1,4-butanediol
        result = reaction_solver.get_reaction_data("OCCCCO")
        count_key = "template_based_reaction_prediction_primary_alcohol_to_aldehyde_count"
        assert result[count_key] == 2

    def test_polyene_multiple_alkenes(self, reaction_solver):
        """Polyene has multiple hydrogenatable sites."""
        # 1,3-butadiene (C=CC=C)
        result = reaction_solver.get_reaction_data("C=CC=C")
        count_key = "template_based_reaction_prediction_alkene_to_alkane_count"
        assert result[count_key] == 2


class TestReactionCategories:
    """Tests that verify reaction categories are represented."""

    def test_oxidation_reactions_present(self, reaction_solver):
        """Oxidation reactions are available."""
        templates = reaction_solver.templates
        oxidation = [k for k in templates if "to_aldehyde" in k or "to_ketone" in k]
        assert len(oxidation) > 0

    def test_reduction_reactions_present(self, reaction_solver):
        """Reduction reactions are available."""
        templates = reaction_solver.templates
        reduction = [k for k in templates if "to_alcohol" in k or "to_amine" in k or "to_alkane" in k]
        assert len(reduction) > 0

    def test_substitution_reactions_present(self, reaction_solver):
        """Substitution reactions are available."""
        templates = reaction_solver.templates
        substitution = [k for k in templates if "halide" in k or "to_alkyl" in k]
        assert len(substitution) > 0
