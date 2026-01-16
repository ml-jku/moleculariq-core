"""
Tests for the properties module.

Validates property maps, alias maps, and helper functions.
"""

from moleculariq_core.properties import (
    INDEX_MAP,
    COUNT_MAP,
    CONSTRAINT_MAP,
    KEY_ALIAS_MAP,
    FG_KEYS,
    COUNT_TO_INDEX_MAP,
    SUBGROUP_DEFINITIONS,
    SUBGROUP_COUNT_MAP,
    SUBGROUP_INDEX_MAP,
    SUBGROUP_CONSTRAINT_MAP,
    ALIAS_TO_KEY_MAP,
    get_alias,
    canonicalize_property_name,
)


class TestIndexMap:
    """Tests for INDEX_MAP dictionary."""

    def test_index_map_not_empty(self):
        """INDEX_MAP should contain entries."""
        assert len(INDEX_MAP) > 0

    def test_index_map_values_are_lists(self):
        """All values should be lists."""
        for key, value in INDEX_MAP.items():
            assert isinstance(value, list), f"Value for {key} is not a list"

    def test_index_map_has_core_categories(self):
        """Should have core property categories."""
        core_categories = [
            "ring",
            "aromatic_ring",
            "carbon_atom",
            "hetero_atom",
            "hba",
            "hbd",
            "stereocenter",
        ]
        for category in core_categories:
            assert category in INDEX_MAP, f"Missing category: {category}"

    def test_index_map_values_end_with_index(self):
        """All values in lists should end with '_index'."""
        for category, props in INDEX_MAP.items():
            for prop in props:
                assert prop.endswith("_index"), f"Property {prop} in {category} doesn't end with '_index'"

    def test_index_map_functional_group_populated(self):
        """functional_group category should be populated from FG_KEYS."""
        fg_indices = INDEX_MAP.get("functional_group", [])
        assert len(fg_indices) > 0, "functional_group indices not populated"
        assert len(fg_indices) == len(FG_KEYS), "functional_group indices count mismatch"


class TestCountMap:
    """Tests for COUNT_MAP dictionary."""

    def test_count_map_not_empty(self):
        """COUNT_MAP should contain entries."""
        assert len(COUNT_MAP) > 0

    def test_count_map_values_are_lists(self):
        """All values should be lists."""
        for key, value in COUNT_MAP.items():
            assert isinstance(value, list), f"Value for {key} is not a list"

    def test_count_map_has_core_categories(self):
        """Should have core property categories."""
        core_categories = [
            "ring",
            "aromatic_ring",
            "carbon_atom",
            "hetero_atom",
            "hba",
            "hbd",
            "stereocenter",
        ]
        for category in core_categories:
            assert category in COUNT_MAP, f"Missing category: {category}"

    def test_count_map_values_end_with_count(self):
        """All values in lists should end with '_count'."""
        for category, props in COUNT_MAP.items():
            for prop in props:
                assert prop.endswith("_count"), f"Property {prop} in {category} doesn't end with '_count'"

    def test_count_map_functional_group_populated(self):
        """functional_group category should be populated from FG_KEYS."""
        fg_counts = COUNT_MAP.get("functional_group", [])
        assert len(fg_counts) > 0, "functional_group counts not populated"
        assert len(fg_counts) == len(FG_KEYS), "functional_group counts count mismatch"


class TestConstraintMap:
    """Tests for CONSTRAINT_MAP dictionary."""

    def test_constraint_map_not_empty(self):
        """CONSTRAINT_MAP should contain entries."""
        assert len(CONSTRAINT_MAP) > 0

    def test_constraint_map_has_reaction_success(self):
        """Should have reaction_success category."""
        assert "reaction_success" in CONSTRAINT_MAP

    def test_constraint_map_reaction_success_populated(self):
        """reaction_success should have multiple entries."""
        reaction_success = CONSTRAINT_MAP.get("reaction_success", [])
        assert len(reaction_success) > 10, "reaction_success should have many entries"

    def test_constraint_map_functional_group_uses_nbrInstances(self):
        """functional_group should use _nbrInstances suffix."""
        fg_constraints = CONSTRAINT_MAP.get("functional_group", [])
        for prop in fg_constraints:
            assert "_nbrInstances" in prop, f"Constraint {prop} should use _nbrInstances"


class TestKeyAliasMap:
    """Tests for KEY_ALIAS_MAP dictionary."""

    def test_key_alias_map_not_empty(self):
        """KEY_ALIAS_MAP should contain entries."""
        assert len(KEY_ALIAS_MAP) > 0

    def test_key_alias_map_values_are_strings(self):
        """All values should be strings."""
        for key, value in KEY_ALIAS_MAP.items():
            assert isinstance(value, str), f"Value for {key} is not a string"

    def test_key_alias_map_has_ring_size_aliases(self):
        """Should have ring size composite aliases."""
        ring_size_keys = [
            "smallest_largest_ring_size_smallest_count",
            "smallest_largest_ring_size_largest_count",
        ]
        for key in ring_size_keys:
            assert key in KEY_ALIAS_MAP, f"Missing alias for {key}"

    def test_key_alias_map_has_stereocenter_aliases(self):
        """Should have stereocenter aliases."""
        stereo_keys = [
            "r_s_stereocenter_r_count",
            "r_s_stereocenter_s_count",
            "stereocenter_count",
        ]
        for key in stereo_keys:
            assert key in KEY_ALIAS_MAP, f"Missing alias for {key}"

    def test_key_alias_map_has_oxidation_state_aliases(self):
        """Should have oxidation state aliases."""
        oxidation_keys = [
            "oxidation_state_C_max_count",
            "oxidation_state_C_min_count",
        ]
        for key in oxidation_keys:
            assert key in KEY_ALIAS_MAP, f"Missing alias for {key}"

    def test_functional_group_aliases_generated(self):
        """Functional group aliases should be generated."""
        for fg in FG_KEYS[:5]:  # Check first 5
            count_key = f"functional_group_{fg}_count"
            assert count_key in KEY_ALIAS_MAP, f"Missing alias for {count_key}"


class TestFGKeys:
    """Tests for FG_KEYS list."""

    def test_fg_keys_not_empty(self):
        """FG_KEYS should contain entries."""
        assert len(FG_KEYS) > 0

    def test_fg_keys_are_strings(self):
        """All entries should be strings."""
        for key in FG_KEYS:
            assert isinstance(key, str)

    def test_fg_keys_has_common_groups(self):
        """Should have common functional groups."""
        # Note: "amine" is split into primary_amine, secondary_amine, tertiary_amine
        common_groups = [
            "alcohol",
            "ketone",
            "aldehyde",
            "primary_amine",
            "carboxylic_acid",
            "ester",
            "amide",
            "ether",
            "thiol",
        ]
        for group in common_groups:
            assert group in FG_KEYS, f"Missing common group: {group}"

    def test_fg_keys_no_duplicates(self):
        """Should have no duplicate entries."""
        assert len(FG_KEYS) == len(set(FG_KEYS))

    def test_fg_keys_uses_underscores(self):
        """Multi-word groups should use underscores."""
        multi_word = [k for k in FG_KEYS if len(k.split('_')) > 1]
        assert len(multi_word) > 0, "Should have multi-word functional groups"
        for key in FG_KEYS:
            assert ' ' not in key, f"Group {key} uses spaces instead of underscores"


class TestCountToIndexMap:
    """Tests for COUNT_TO_INDEX_MAP dictionary."""

    def test_count_to_index_map_not_empty(self):
        """COUNT_TO_INDEX_MAP should contain entries."""
        assert len(COUNT_TO_INDEX_MAP) > 0

    def test_count_to_index_map_keys_end_with_count(self):
        """All keys should end with '_count'."""
        for key in COUNT_TO_INDEX_MAP.keys():
            assert key.endswith("_count"), f"Key {key} doesn't end with '_count'"

    def test_count_to_index_map_values_end_with_index(self):
        """All values should end with '_index'."""
        for key, value in COUNT_TO_INDEX_MAP.items():
            assert value.endswith("_index"), f"Value {value} for {key} doesn't end with '_index'"

    def test_count_to_index_map_has_core_mappings(self):
        """Should have core count-to-index mappings."""
        core_mappings = [
            "ring_count",
            "aromatic_ring_count",
            "carbon_atom_count",
            "hba_count",
            "hbd_count",
        ]
        for key in core_mappings:
            assert key in COUNT_TO_INDEX_MAP, f"Missing mapping for {key}"

    def test_count_to_index_map_has_functional_group_mappings(self):
        """Should have functional group mappings."""
        for fg in FG_KEYS[:5]:  # Check first 5
            count_key = f"functional_group_{fg}_count"
            assert count_key in COUNT_TO_INDEX_MAP, f"Missing mapping for {count_key}"

    def test_count_to_index_map_consistency(self):
        """Count key prefix should match index value prefix."""
        for count_key, index_value in COUNT_TO_INDEX_MAP.items():
            count_prefix = count_key[:-6]  # Remove "_count"
            index_prefix = index_value[:-6]  # Remove "_index"
            assert count_prefix == index_prefix, f"Prefix mismatch: {count_key} -> {index_value}"


class TestSubgroupDefinitions:
    """Tests for SUBGROUP_DEFINITIONS dictionary."""

    def test_subgroup_definitions_not_empty(self):
        """SUBGROUP_DEFINITIONS should contain entries."""
        assert len(SUBGROUP_DEFINITIONS) > 0

    def test_subgroup_definitions_has_expected_groups(self):
        """Should have expected subgroups."""
        expected_groups = [
            "graph_topology",
            "chemistry_typed_topology",
            "composition",
            "chemical_perception",
            "functional_groups",
            "synthesis",
        ]
        for group in expected_groups:
            assert group in SUBGROUP_DEFINITIONS, f"Missing subgroup: {group}"

    def test_subgroup_definitions_values_are_lists(self):
        """All values should be lists."""
        for key, value in SUBGROUP_DEFINITIONS.items():
            assert isinstance(value, list), f"Value for {key} is not a list"

    def test_subgroup_definitions_values_not_empty(self):
        """All subgroups should have at least one category."""
        for key, value in SUBGROUP_DEFINITIONS.items():
            assert len(value) > 0, f"Empty subgroup: {key}"


class TestSubgroupMaps:
    """Tests for dynamically built SUBGROUP_*_MAP dictionaries."""

    def test_subgroup_count_map_not_empty(self):
        """SUBGROUP_COUNT_MAP should contain entries."""
        assert len(SUBGROUP_COUNT_MAP) > 0

    def test_subgroup_index_map_not_empty(self):
        """SUBGROUP_INDEX_MAP should contain entries."""
        assert len(SUBGROUP_INDEX_MAP) > 0

    def test_subgroup_constraint_map_not_empty(self):
        """SUBGROUP_CONSTRAINT_MAP should contain entries."""
        assert len(SUBGROUP_CONSTRAINT_MAP) > 0

    def test_subgroup_maps_values_are_lists(self):
        """All values in subgroup maps should be lists."""
        for map_name, map_dict in [
            ("count", SUBGROUP_COUNT_MAP),
            ("index", SUBGROUP_INDEX_MAP),
            ("constraint", SUBGROUP_CONSTRAINT_MAP),
        ]:
            for key, value in map_dict.items():
                assert isinstance(value, list), f"Value for {key} in {map_name} is not a list"


class TestAliasToKeyMap:
    """Tests for ALIAS_TO_KEY_MAP dictionary."""

    def test_alias_to_key_map_not_empty(self):
        """ALIAS_TO_KEY_MAP should contain entries."""
        assert len(ALIAS_TO_KEY_MAP) > 0

    def test_alias_to_key_map_values_are_strings(self):
        """All values should be strings."""
        for key, value in ALIAS_TO_KEY_MAP.items():
            assert isinstance(value, str), f"Value for {key} is not a string"

    def test_alias_to_key_map_contains_lowercase_variants(self):
        """Should contain lowercase variants."""
        # Check if both original and lowercase versions exist for some keys
        sample_keys = ["ring_count", "aromatic_ring_count", "carbon_atom_count"]
        for key in sample_keys:
            assert key in ALIAS_TO_KEY_MAP or key.lower() in ALIAS_TO_KEY_MAP, f"Missing {key}"


class TestGetAlias:
    """Tests for get_alias function."""

    def test_get_alias_returns_alias_for_known_key(self):
        """Should return alias for known technical key."""
        result = get_alias("smallest_largest_ring_size_smallest_count")
        assert result == "smallest_ring_atom_count"

    def test_get_alias_returns_original_for_unknown_key(self):
        """Should return original for unknown key."""
        result = get_alias("unknown_property_key")
        assert result == "unknown_property_key"

    def test_get_alias_handles_empty_string(self):
        """Should handle empty string."""
        result = get_alias("")
        assert result == ""

    def test_get_alias_returns_string(self):
        """Should always return a string."""
        test_cases = [
            "ring_count",
            "unknown",
            "",
            "stereocenter_count",
        ]
        for key in test_cases:
            result = get_alias(key)
            assert isinstance(result, str)


class TestCanonicalizePropertyName:
    """Tests for canonicalize_property_name function."""

    def test_canonicalize_returns_canonical_for_alias(self):
        """Should return canonical key for alias."""
        result = canonicalize_property_name("smallest_ring_atom_count")
        assert result == "smallest_largest_ring_size_smallest_count"

    def test_canonicalize_returns_original_for_canonical(self):
        """Should return same value for canonical key."""
        result = canonicalize_property_name("ring_count")
        assert result == "ring_count"

    def test_canonicalize_case_insensitive(self):
        """Should be case insensitive."""
        result1 = canonicalize_property_name("ring_count")
        result2 = canonicalize_property_name("RING_COUNT")
        result3 = canonicalize_property_name("Ring_Count")
        assert result1 == result2 == result3

    def test_canonicalize_strips_whitespace(self):
        """Should strip whitespace."""
        result = canonicalize_property_name("  ring_count  ")
        assert result == "ring_count"

    def test_canonicalize_returns_original_for_unknown(self):
        """Should return original for unknown key."""
        result = canonicalize_property_name("completely_unknown_key")
        assert result == "completely_unknown_key"

    def test_canonicalize_handles_non_string(self):
        """Should handle non-string input gracefully."""
        result = canonicalize_property_name(123)
        assert result == 123


class TestMapConsistency:
    """Tests for consistency between different maps."""

    def test_index_count_map_categories_match(self):
        """INDEX_MAP and COUNT_MAP should have same categories."""
        index_categories = set(INDEX_MAP.keys())
        count_categories = set(COUNT_MAP.keys())
        # Most categories should be in both (some exceptions allowed)
        common = index_categories & count_categories
        assert len(common) >= 15, "Too few common categories between INDEX_MAP and COUNT_MAP"

    def test_count_constraint_map_overlap(self):
        """COUNT_MAP and CONSTRAINT_MAP should have significant overlap."""
        count_categories = set(COUNT_MAP.keys())
        constraint_categories = set(CONSTRAINT_MAP.keys())
        common = count_categories & constraint_categories
        assert len(common) >= 15, "Too few common categories between COUNT_MAP and CONSTRAINT_MAP"

    def test_key_alias_map_bidirectional_consistency(self):
        """KEY_ALIAS_MAP should have consistent bidirectional mapping in ALIAS_TO_KEY_MAP."""
        for technical_key, alias in KEY_ALIAS_MAP.items():
            # The alias should map back to a valid technical key
            # Note: some aliases may map to uppercase or lowercase variants
            # (e.g., oxidation_state_c_max_count vs oxidation_state_C_max_count)
            if alias in ALIAS_TO_KEY_MAP:
                mapped_key = ALIAS_TO_KEY_MAP[alias]
                # Accept either exact match or case-insensitive match
                assert mapped_key.lower() == technical_key.lower() or mapped_key == technical_key, \
                    f"Bidirectional mismatch for {alias}: {mapped_key} vs {technical_key}"


class TestFunctionalGroupIntegration:
    """Tests for functional group integration across maps."""

    def test_all_fg_keys_in_index_map(self):
        """All FG_KEYS should have entries in INDEX_MAP functional_group."""
        fg_indices = INDEX_MAP.get("functional_group", [])
        for fg in FG_KEYS:
            expected = f"functional_group_{fg}_index"
            assert expected in fg_indices, f"Missing {expected} in INDEX_MAP"

    def test_all_fg_keys_in_count_map(self):
        """All FG_KEYS should have entries in COUNT_MAP functional_group."""
        fg_counts = COUNT_MAP.get("functional_group", [])
        for fg in FG_KEYS:
            expected = f"functional_group_{fg}_count"
            assert expected in fg_counts, f"Missing {expected} in COUNT_MAP"

    def test_all_fg_keys_in_constraint_map(self):
        """All FG_KEYS should have entries in CONSTRAINT_MAP functional_group."""
        fg_constraints = CONSTRAINT_MAP.get("functional_group", [])
        for fg in FG_KEYS:
            expected = f"functional_group_{fg}_nbrInstances"
            assert expected in fg_constraints, f"Missing {expected} in CONSTRAINT_MAP"

    def test_all_fg_keys_in_count_to_index_map(self):
        """All FG_KEYS should have entries in COUNT_TO_INDEX_MAP."""
        for fg in FG_KEYS:
            count_key = f"functional_group_{fg}_count"
            assert count_key in COUNT_TO_INDEX_MAP, f"Missing {count_key} in COUNT_TO_INDEX_MAP"

    def test_all_fg_keys_have_aliases(self):
        """All FG_KEYS should have aliases in KEY_ALIAS_MAP."""
        for fg in FG_KEYS:
            count_key = f"functional_group_{fg}_count"
            index_key = f"functional_group_{fg}_index"
            assert count_key in KEY_ALIAS_MAP, f"Missing alias for {count_key}"
            assert index_key in KEY_ALIAS_MAP, f"Missing alias for {index_key}"


class TestReactionSuccessIntegration:
    """Tests for reaction success integration in CONSTRAINT_MAP."""

    def test_reaction_success_entries_have_correct_format(self):
        """reaction_success entries should have correct naming format."""
        reaction_success = CONSTRAINT_MAP.get("reaction_success", [])
        for entry in reaction_success:
            assert entry.startswith("template_based_reaction_prediction_"), f"Wrong prefix: {entry}"
            assert entry.endswith("_success"), f"Wrong suffix: {entry}"

    def test_reaction_success_entries_have_aliases(self):
        """reaction_success entries should have aliases."""
        reaction_success = CONSTRAINT_MAP.get("reaction_success", [])
        for entry in reaction_success[:5]:  # Check first 5
            assert entry in KEY_ALIAS_MAP, f"Missing alias for {entry}"
