"""
Tests for the property_solver_mapping module.

Validates property-to-solver mappings and helper functions.
"""

import pytest

from moleculariq_core.rewards.property_solver_mapping import (
    PROPERTY_TO_SOLVER_MAP,
    get_functional_group_mapping,
    get_reaction_template_mapping,
    is_string_valued_property,
    get_solver_mapping,
)
from moleculariq_core import SymbolicSolver
from moleculariq_core.properties import FG_KEYS


class TestPropertyToSolverMap:
    """Tests for PROPERTY_TO_SOLVER_MAP dictionary."""

    def test_map_not_empty(self):
        """PROPERTY_TO_SOLVER_MAP should contain entries."""
        assert len(PROPERTY_TO_SOLVER_MAP) > 0

    def test_map_has_core_properties(self):
        """Should have core property mappings."""
        core_properties = [
            "ring_count",
            "ring_index",
            "aromatic_ring_count",
            "carbon_atom_count",
            "hetero_atom_count",
            "hba_count",
            "hbd_count",
            "stereocenter_count",
            "molecular_formula",
        ]
        for prop in core_properties:
            assert prop in PROPERTY_TO_SOLVER_MAP, f"Missing core property: {prop}"

    def test_map_values_are_tuples(self):
        """All values should be (method_name, params) tuples."""
        for key, value in PROPERTY_TO_SOLVER_MAP.items():
            assert isinstance(value, tuple), f"Value for {key} is not a tuple"
            assert len(value) == 2, f"Value for {key} should have 2 elements"
            method_name, params = value
            assert isinstance(method_name, str), f"Method name for {key} is not a string"
            assert isinstance(params, dict), f"Params for {key} is not a dict"

    def test_method_names_start_with_get(self):
        """All method names should start with 'get_'."""
        for key, (method_name, _) in PROPERTY_TO_SOLVER_MAP.items():
            assert method_name.startswith("get_"), f"Method {method_name} for {key} doesn't start with 'get_'"

    def test_count_properties_map_to_count_methods(self):
        """Properties ending in _count should map to count methods."""
        for key, (method_name, _) in PROPERTY_TO_SOLVER_MAP.items():
            if key.endswith("_count") and not key.startswith("molecular_formula"):
                assert "count" in method_name.lower() or "formula" in method_name.lower(), \
                    f"Count property {key} maps to non-count method {method_name}"

    def test_index_properties_map_to_index_methods(self):
        """Properties ending in _index should map to index methods."""
        for key, (method_name, _) in PROPERTY_TO_SOLVER_MAP.items():
            if key.endswith("_index"):
                assert "indic" in method_name.lower() or "bond" in method_name.lower(), \
                    f"Index property {key} maps to non-index method {method_name}"


class TestMappingMethodsExist:
    """Tests that mapped methods exist on SymbolicSolver."""

    @pytest.fixture(scope="class")
    def solver(self):
        """Create solver instance for method existence checks."""
        return SymbolicSolver()

    def test_all_mapped_methods_exist(self, solver):
        """All methods referenced in mapping should exist on solver."""
        missing_methods = []
        for prop_name, (method_name, _) in PROPERTY_TO_SOLVER_MAP.items():
            if not hasattr(solver, method_name):
                missing_methods.append(f"{prop_name} -> {method_name}")

        assert len(missing_methods) == 0, f"Missing methods: {missing_methods[:10]}"

    def test_mapped_methods_are_callable(self, solver):
        """All mapped methods should be callable."""
        for prop_name, (method_name, _) in PROPERTY_TO_SOLVER_MAP.items():
            method = getattr(solver, method_name, None)
            if method is not None:
                assert callable(method), f"Method {method_name} for {prop_name} is not callable"


class TestGetFunctionalGroupMapping:
    """Tests for get_functional_group_mapping function."""

    def test_returns_tuple(self):
        """Should return a (method_name, params) tuple."""
        result = get_functional_group_mapping("alcohol", "count")
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_method_name_is_correct(self):
        """Should return correct method name."""
        method_name, _ = get_functional_group_mapping("ketone", "count")
        assert method_name == "get_functional_group_count_and_indices"

    def test_params_contain_group_name(self):
        """Params should contain the group name."""
        _, params = get_functional_group_mapping("aldehyde", "index")
        assert "group_name" in params
        assert params["group_name"] == "aldehyde"

    def test_different_groups_same_method(self):
        """All functional groups should use the same method."""
        groups = ["alcohol", "ketone", "ester", "amine"]
        methods = set()
        for group in groups:
            method_name, _ = get_functional_group_mapping(group, "count")
            methods.add(method_name)
        assert len(methods) == 1, "All functional groups should use same method"

    def test_all_fg_keys_work(self):
        """Should work for all functional group keys."""
        for fg_key in FG_KEYS[:10]:  # Test first 10
            result = get_functional_group_mapping(fg_key, "count")
            assert result is not None


class TestGetReactionTemplateMapping:
    """Tests for get_reaction_template_mapping function."""

    def test_returns_tuple(self):
        """Should return a (method_name, params) tuple."""
        result = get_reaction_template_mapping("bromination")
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_method_name_is_correct(self):
        """Should return correct method name."""
        method_name, _ = get_reaction_template_mapping("nitration")
        assert method_name == "get_reaction_counts_and_indices"

    def test_params_contain_template_name(self):
        """Params should contain the template name."""
        _, params = get_reaction_template_mapping("epoxidation")
        assert "template_name" in params
        assert params["template_name"] == "epoxidation"

    def test_different_templates_same_method(self):
        """All templates should use the same method."""
        templates = ["bromination", "nitration", "epoxidation", "chlorination"]
        methods = set()
        for template in templates:
            method_name, _ = get_reaction_template_mapping(template)
            methods.add(method_name)
        assert len(methods) == 1, "All templates should use same method"


class TestIsStringValuedProperty:
    """Tests for is_string_valued_property function."""

    def test_molecular_formula_is_string(self):
        """molecular_formula should be string-valued."""
        assert is_string_valued_property("molecular_formula") is True
        assert is_string_valued_property("molecular_formula_count") is True
        assert is_string_valued_property("molecular_formula_value") is True

    def test_murcko_scaffold_value_is_string(self):
        """murcko_scaffold_value should be string-valued."""
        assert is_string_valued_property("murcko_scaffold_value") is True

    def test_count_properties_are_not_string(self):
        """Count properties should not be string-valued."""
        assert is_string_valued_property("ring_count") is False
        assert is_string_valued_property("carbon_atom_count") is False
        assert is_string_valued_property("hba_count") is False

    def test_index_properties_are_not_string(self):
        """Index properties should not be string-valued."""
        assert is_string_valued_property("ring_index") is False
        assert is_string_valued_property("carbon_atom_index") is False

    def test_murcko_scaffold_count_is_not_string(self):
        """murcko_scaffold_count should not be string-valued."""
        assert is_string_valued_property("murcko_scaffold_count") is False


class TestGetSolverMapping:
    """Tests for get_solver_mapping function."""

    def test_explicit_mapping(self):
        """Should return mapping for explicitly mapped properties."""
        result = get_solver_mapping("ring_count")
        assert result is not None
        method_name, params = result
        assert method_name == "get_ring_count"
        assert params == {}

    def test_parametrized_mapping(self):
        """Should return mapping with parameters for parametrized properties."""
        result = get_solver_mapping("r_s_stereocenter_r_count")
        assert result is not None
        method_name, params = result
        assert method_name == "get_r_or_s_stereocenter_count"
        assert "r_count" in params
        assert params["r_count"] is True

    def test_functional_group_mapping(self):
        """Should return mapping for functional group properties."""
        result = get_solver_mapping("functional_group_alcohol_count")
        assert result is not None
        method_name, params = result
        assert method_name == "get_functional_group_count_and_indices"
        assert params["group_name"] == "alcohol"

    def test_functional_group_index_mapping(self):
        """Should return mapping for functional group index properties."""
        result = get_solver_mapping("functional_group_ketone_index")
        assert result is not None
        method_name, params = result
        assert method_name == "get_functional_group_count_and_indices"
        assert params["group_name"] == "ketone"

    def test_reaction_template_mapping(self):
        """Should return mapping for reaction template properties."""
        result = get_solver_mapping("template_based_reaction_prediction_bromination_success")
        assert result is not None
        method_name, params = result
        assert method_name == "get_reaction_counts_and_indices"

    def test_legacy_reaction_format(self):
        """Should handle legacy reaction format."""
        result = get_solver_mapping("reaction_bromination_count")
        assert result is not None
        method_name, _ = result
        assert method_name == "get_reaction_counts_and_indices"

    def test_unknown_property_returns_none(self):
        """Should return None for unknown properties."""
        result = get_solver_mapping("completely_unknown_property")
        assert result is None

    def test_oxidation_state_case_handling(self):
        """Should handle oxidation state element case variations."""
        # Uppercase
        result = get_solver_mapping("oxidation_state_C_max_count")
        assert result is not None
        assert result[1]["element"] == "C"

        # Lowercase should still work (normalized)
        result_lower = get_solver_mapping("oxidation_state_c_max_count")
        assert result_lower is not None


class TestMappingCompleteness:
    """Tests for mapping completeness."""

    def test_ring_properties_complete(self):
        """All ring properties should have mappings."""
        ring_properties = [
            "ring_count", "ring_index",
            "aromatic_ring_count", "aromatic_ring_index",
            "aliphatic_ring_count", "aliphatic_ring_index",
            "fused_ring_count", "fused_ring_index",
            "saturated_ring_count", "saturated_ring_index",
        ]
        for prop in ring_properties:
            assert get_solver_mapping(prop) is not None, f"Missing mapping for {prop}"

    def test_stereochemistry_properties_complete(self):
        """All stereochemistry properties should have mappings."""
        stereo_properties = [
            "stereocenter_count", "stereocenter_index",
            "r_s_stereocenter_r_count", "r_s_stereocenter_s_count",
            "r_s_stereocenter_r_index", "r_s_stereocenter_s_index",
            "e_z_stereochemistry_double_bond_e_count",
            "e_z_stereochemistry_double_bond_z_count",
        ]
        for prop in stereo_properties:
            assert get_solver_mapping(prop) is not None, f"Missing mapping for {prop}"

    def test_atom_properties_complete(self):
        """All atom count/index properties should have mappings."""
        atom_properties = [
            "carbon_atom_count", "carbon_atom_index",
            "hetero_atom_count", "hetero_atom_index",
            "halogen_atom_count", "halogen_atom_index",
            "heavy_atom_count", "heavy_atom_index",
            "hydrogen_atom_count",
        ]
        for prop in atom_properties:
            assert get_solver_mapping(prop) is not None, f"Missing mapping for {prop}"

    def test_oxidation_state_properties_complete(self):
        """All oxidation state properties should have mappings."""
        elements = ["C", "N", "O", "P", "S"]
        for element in elements:
            for suffix in ["max_count", "min_count", "max_index", "min_index"]:
                prop = f"oxidation_state_{element}_{suffix}"
                assert get_solver_mapping(prop) is not None, f"Missing mapping for {prop}"


class TestMappingIntegration:
    """Integration tests that verify mappings work with actual solver calls."""

    @pytest.fixture(scope="class")
    def solver(self):
        """Create solver instance for integration tests."""
        return SymbolicSolver()

    def test_ring_count_mapping_works(self, solver):
        """Ring count mapping should work end-to-end."""
        mapping = get_solver_mapping("ring_count")
        method_name, params = mapping
        method = getattr(solver, method_name)
        result = method("c1ccccc1", **params)
        assert result == 1

    def test_carbon_count_mapping_works(self, solver):
        """Carbon count mapping should work end-to-end."""
        mapping = get_solver_mapping("carbon_atom_count")
        method_name, params = mapping
        method = getattr(solver, method_name)
        result = method("CCO", **params)
        assert result == 2

    def test_stereocenter_mapping_works(self, solver):
        """Stereocenter mapping should work end-to-end."""
        mapping = get_solver_mapping("stereocenter_count")
        method_name, params = mapping
        method = getattr(solver, method_name)
        result = method("C[C@H](N)C(=O)O", **params)
        assert result == 1

    def test_parametrized_mapping_works(self, solver):
        """Parametrized mapping should pass parameters correctly."""
        mapping = get_solver_mapping("smallest_largest_ring_size_smallest_count")
        method_name, params = mapping
        method = getattr(solver, method_name)
        result = method("c1ccccc1", **params)
        assert result == 6  # Benzene has 6 atoms in smallest ring

    def test_molecular_formula_mapping_works(self, solver):
        """Molecular formula mapping should work end-to-end."""
        mapping = get_solver_mapping("molecular_formula")
        method_name, params = mapping
        method = getattr(solver, method_name)
        result = method("CCO", **params)
        assert "C" in result
        assert "H" in result
        assert "O" in result
