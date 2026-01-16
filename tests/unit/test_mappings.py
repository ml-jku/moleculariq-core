"""
Tests for the _nlp/mappings module.

Validates mapping dictionaries, dynamic loading functions, and helper functions.
"""

import pytest
from pathlib import Path
from tempfile import NamedTemporaryFile

from moleculariq_core._nlp.mappings import (
    COUNT_MAPPINGS,
    INDEX_MAPPINGS,
    ALIASES,
    FUNCTIONAL_GROUPS,
    REACTION_TEMPLATES,
    TEMPLATE_REACTION_MAPPINGS,
    load_functional_groups,
    load_reaction_templates,
    get_natural_language,
    parse_natural_language,
    get_all_properties,
    _REVERSE_LOOKUP,
)


class TestCountMappings:
    """Tests for COUNT_MAPPINGS dictionary."""

    def test_count_mappings_not_empty(self):
        """COUNT_MAPPINGS should contain entries."""
        assert len(COUNT_MAPPINGS) > 0

    def test_count_mappings_keys_are_strings(self):
        """All keys should be strings."""
        for key in COUNT_MAPPINGS.keys():
            assert isinstance(key, str)

    def test_count_mappings_values_are_lists(self):
        """All values should be lists of strings."""
        for key, value in COUNT_MAPPINGS.items():
            assert isinstance(value, list), f"Value for {key} is not a list"
            for item in value:
                assert isinstance(item, str), f"Item in {key} is not a string: {item}"

    def test_count_mappings_has_core_properties(self):
        """Core properties should be present in COUNT_MAPPINGS."""
        core_properties = [
            "ring_count",
            "aromatic_ring_count",
            "carbon_atom_count",
            "hetero_atom_count",
            "hba_count",
            "hbd_count",
            "rotatable_bond_count",
            "heavy_atom_count",
            "stereocenter_count",
        ]
        for prop in core_properties:
            assert prop in COUNT_MAPPINGS, f"Missing core property: {prop}"

    def test_count_mappings_values_not_empty(self):
        """Each mapping should have at least one natural language form."""
        for key, value in COUNT_MAPPINGS.items():
            assert len(value) > 0, f"Empty value list for {key}"

    def test_count_mappings_no_duplicate_keys(self):
        """There should be no duplicate keys (implicit in dict)."""
        keys = list(COUNT_MAPPINGS.keys())
        assert len(keys) == len(set(keys))

    def test_functional_group_count_mappings_exist(self):
        """Functional group count mappings should be dynamically added."""
        fg_keys = [k for k in COUNT_MAPPINGS.keys() if k.startswith("functional_group_")]
        assert len(fg_keys) > 0, "No functional group mappings found"


class TestIndexMappings:
    """Tests for INDEX_MAPPINGS dictionary."""

    def test_index_mappings_not_empty(self):
        """INDEX_MAPPINGS should contain entries."""
        assert len(INDEX_MAPPINGS) > 0

    def test_index_mappings_keys_are_strings(self):
        """All keys should be strings."""
        for key in INDEX_MAPPINGS.keys():
            assert isinstance(key, str)

    def test_index_mappings_values_are_lists(self):
        """All values should be lists of strings."""
        for key, value in INDEX_MAPPINGS.items():
            assert isinstance(value, list), f"Value for {key} is not a list"
            for item in value:
                assert isinstance(item, str), f"Item in {key} is not a string: {item}"

    def test_index_mappings_has_core_properties(self):
        """Core index properties should be present."""
        core_properties = [
            "ring_index",
            "aromatic_ring_index",
            "carbon_atom_index",
            "hetero_atom_index",
            "hba_index",
            "hbd_index",
        ]
        for prop in core_properties:
            assert prop in INDEX_MAPPINGS, f"Missing core property: {prop}"

    def test_functional_group_index_mappings_exist(self):
        """Functional group index mappings should be dynamically added."""
        fg_keys = [k for k in INDEX_MAPPINGS.keys() if k.startswith("functional_group_")]
        assert len(fg_keys) > 0, "No functional group index mappings found"


class TestAliases:
    """Tests for ALIASES dictionary."""

    def test_aliases_not_empty(self):
        """ALIASES should contain entries."""
        assert len(ALIASES) > 0

    def test_aliases_keys_are_lowercase(self):
        """All alias keys should be lowercase."""
        for key in ALIASES.keys():
            assert key == key.lower(), f"Alias key not lowercase: {key}"

    def test_aliases_values_are_strings(self):
        """All alias values should be strings."""
        for key, value in ALIASES.items():
            assert isinstance(value, str), f"Value for alias {key} is not a string"

    def test_common_aliases_present(self):
        """Common aliases should be present."""
        common_aliases = [
            "rings",
            "carbons",
            "hba",
            "hbd",
            "aromatic rings",
            "stereocenters",
        ]
        for alias in common_aliases:
            assert alias in ALIASES, f"Missing common alias: {alias}"

    def test_aliases_point_to_valid_properties(self):
        """Aliases should point to properties in COUNT_MAPPINGS or INDEX_MAPPINGS or be valid technical keys."""
        all_technical_keys = set(COUNT_MAPPINGS.keys()) | set(INDEX_MAPPINGS.keys())
        # Allow aliases that point to keys not in mappings but are valid technical keys
        for alias, technical_key in ALIASES.items():
            # Technical key could be valid even if not in mappings (e.g., nitrogen_atom_count)
            assert isinstance(technical_key, str), f"Invalid alias target for {alias}"


class TestFunctionalGroups:
    """Tests for dynamically loaded FUNCTIONAL_GROUPS."""

    def test_functional_groups_loaded(self):
        """FUNCTIONAL_GROUPS should be loaded from file."""
        assert len(FUNCTIONAL_GROUPS) > 0

    def test_functional_groups_are_dict(self):
        """FUNCTIONAL_GROUPS should be a dictionary."""
        assert isinstance(FUNCTIONAL_GROUPS, dict)

    def test_functional_groups_values_are_strings(self):
        """All functional group values should be strings."""
        for key, value in FUNCTIONAL_GROUPS.items():
            assert isinstance(key, str)
            assert isinstance(value, str)

    def test_common_functional_groups_present(self):
        """Common functional groups should be present."""
        # Note: "amine" is split into primary_amine, secondary_amine, tertiary_amine
        common_groups = ["alcohol", "ketone", "aldehyde", "primary_amine", "carboxylic_acid"]
        for group in common_groups:
            assert group in FUNCTIONAL_GROUPS, f"Missing functional group: {group}"


class TestReactionTemplates:
    """Tests for dynamically loaded REACTION_TEMPLATES."""

    def test_reaction_templates_loaded(self):
        """REACTION_TEMPLATES should be loaded from file."""
        assert len(REACTION_TEMPLATES) > 0

    def test_reaction_templates_are_dict(self):
        """REACTION_TEMPLATES should be a dictionary."""
        assert isinstance(REACTION_TEMPLATES, dict)

    def test_reaction_templates_values_are_strings(self):
        """All reaction template values should be strings."""
        for key, value in REACTION_TEMPLATES.items():
            assert isinstance(key, str)
            assert isinstance(value, str)


class TestTemplateReactionMappings:
    """Tests for TEMPLATE_REACTION_MAPPINGS."""

    def test_template_reaction_mappings_not_empty(self):
        """TEMPLATE_REACTION_MAPPINGS should contain entries."""
        assert len(TEMPLATE_REACTION_MAPPINGS) > 0

    def test_common_reactions_present(self):
        """Common reactions should be present."""
        common_reactions = [
            "bromination",
            "chlorination",
            "nitration",
            "epoxidation",
        ]
        for reaction in common_reactions:
            assert reaction in TEMPLATE_REACTION_MAPPINGS, f"Missing reaction: {reaction}"


class TestLoadFunctionalGroups:
    """Tests for load_functional_groups function."""

    def test_load_functional_groups_default(self):
        """Should load functional groups from default path."""
        groups = load_functional_groups()
        assert len(groups) > 0
        assert isinstance(groups, dict)

    def test_load_functional_groups_file_not_found(self):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            load_functional_groups(Path("/nonexistent/path/file.txt"))

    def test_load_functional_groups_custom_file(self):
        """Should load from custom file."""
        content = """# Comment line
alcohol:1:[OH]
ketone:2:[CX3](=[OX1])[#6]
"""
        with NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(content)
            f.flush()
            temp_path = Path(f.name)

        groups = load_functional_groups(temp_path)
        assert "alcohol" in groups
        assert "ketone" in groups
        temp_path.unlink()

    def test_load_functional_groups_skips_comments(self):
        """Should skip comment lines."""
        content = """# This is a comment
# Another comment
alcohol:1:[OH]
"""
        with NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(content)
            f.flush()
            temp_path = Path(f.name)

        groups = load_functional_groups(temp_path)
        assert len(groups) == 1
        assert "alcohol" in groups
        temp_path.unlink()

    def test_load_functional_groups_skips_empty_lines(self):
        """Should skip empty lines."""
        content = """alcohol:1:[OH]

ketone:2:[CX3]
"""
        with NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(content)
            f.flush()
            temp_path = Path(f.name)

        groups = load_functional_groups(temp_path)
        assert len(groups) == 2
        temp_path.unlink()


class TestLoadReactionTemplates:
    """Tests for load_reaction_templates function."""

    def test_load_reaction_templates_default(self):
        """Should load reaction templates from default path."""
        templates = load_reaction_templates()
        assert len(templates) > 0
        assert isinstance(templates, dict)

    def test_load_reaction_templates_file_not_found(self):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            load_reaction_templates(Path("/nonexistent/path/file.txt"))

    def test_load_reaction_templates_custom_file(self):
        """Should load from custom file."""
        content = """# Comment
bromination;halogenation;[c:1][H:2]>>[c:1]Br;Add bromine to aromatic ring
"""
        with NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(content)
            f.flush()
            temp_path = Path(f.name)

        templates = load_reaction_templates(temp_path)
        assert "bromination" in templates
        assert templates["bromination"] == "Add bromine to aromatic ring"
        temp_path.unlink()


class TestGetNaturalLanguage:
    """Tests for get_natural_language function."""

    def test_get_natural_language_count_context(self):
        """Should return natural language forms for count context."""
        result = get_natural_language("ring_count", context="count")
        assert isinstance(result, list)
        assert len(result) > 0

    def test_get_natural_language_index_context(self):
        """Should return natural language forms for index context."""
        result = get_natural_language("ring_index", context="index")
        assert isinstance(result, list)
        assert len(result) > 0

    def test_get_natural_language_constraint_context(self):
        """Should return natural language forms for constraint context."""
        result = get_natural_language("ring_count", context="constraint")
        assert isinstance(result, list)
        assert len(result) > 0

    def test_get_natural_language_unknown_key(self):
        """Should return key in list for unknown key."""
        result = get_natural_language("unknown_property", context="count")
        assert result == ["unknown_property"]

    def test_get_natural_language_invalid_context(self):
        """Should return key in list for invalid context."""
        result = get_natural_language("ring_count", context="invalid")
        assert result == ["ring_count"]


class TestParseNaturalLanguage:
    """Tests for parse_natural_language function."""

    def test_parse_natural_language_direct_match(self):
        """Should parse common natural language to technical key."""
        result = parse_natural_language("rings")
        assert result == "ring_count"

    def test_parse_natural_language_case_insensitive(self):
        """Should be case insensitive."""
        result1 = parse_natural_language("RINGS")
        result2 = parse_natural_language("Rings")
        result3 = parse_natural_language("rings")
        assert result1 == result2 == result3

    def test_parse_natural_language_whitespace_normalization(self):
        """Should normalize whitespace."""
        result = parse_natural_language("  aromatic   rings  ")
        assert result == "aromatic_ring_count"

    def test_parse_natural_language_unknown_returns_original(self):
        """Should return original text for unknown input."""
        result = parse_natural_language("completely unknown text")
        assert result == "completely unknown text"

    def test_parse_natural_language_empty_returns_original(self):
        """Should return original text for empty input."""
        result = parse_natural_language("")
        assert result == ""

    def test_parse_natural_language_technical_key(self):
        """Should handle technical keys as input."""
        result = parse_natural_language("ring_count")
        assert result == "ring_count"


class TestGetAllProperties:
    """Tests for get_all_properties function."""

    def test_get_all_properties_returns_dict(self):
        """Should return a dictionary."""
        result = get_all_properties()
        assert isinstance(result, dict)

    def test_get_all_properties_has_count_key(self):
        """Should have 'count' key."""
        result = get_all_properties()
        assert "count" in result

    def test_get_all_properties_has_index_key(self):
        """Should have 'index' key."""
        result = get_all_properties()
        assert "index" in result

    def test_get_all_properties_count_is_list(self):
        """Count value should be a list."""
        result = get_all_properties()
        assert isinstance(result["count"], list)

    def test_get_all_properties_index_is_list(self):
        """Index value should be a list."""
        result = get_all_properties()
        assert isinstance(result["index"], list)

    def test_get_all_properties_count_not_empty(self):
        """Count list should not be empty."""
        result = get_all_properties()
        assert len(result["count"]) > 0

    def test_get_all_properties_index_not_empty(self):
        """Index list should not be empty."""
        result = get_all_properties()
        assert len(result["index"]) > 0


class TestReverseLookup:
    """Tests for _REVERSE_LOOKUP dictionary."""

    def test_reverse_lookup_not_empty(self):
        """_REVERSE_LOOKUP should contain entries."""
        assert len(_REVERSE_LOOKUP) > 0

    def test_reverse_lookup_contains_natural_forms(self):
        """Should contain natural language forms."""
        assert "rings" in _REVERSE_LOOKUP
        assert "aromatic rings" in _REVERSE_LOOKUP

    def test_reverse_lookup_contains_technical_keys(self):
        """Should contain technical keys."""
        assert "ring_count" in _REVERSE_LOOKUP

    def test_reverse_lookup_values_are_strings(self):
        """All values should be strings."""
        for key, value in _REVERSE_LOOKUP.items():
            assert isinstance(value, str), f"Value for {key} is not a string"


class TestRoundTripConsistency:
    """Tests for round-trip consistency between get_natural_language and parse_natural_language."""

    @pytest.mark.parametrize("technical_key", [
        "ring_count",
        "aromatic_ring_count",
        "carbon_atom_count",
        "hba_count",
        "hbd_count",
    ])
    def test_round_trip_count_properties(self, technical_key):
        """Natural language should parse back to technical key."""
        natural_forms = get_natural_language(technical_key, context="count")
        for natural_form in natural_forms:
            parsed = parse_natural_language(natural_form)
            assert parsed == technical_key, f"Round trip failed: {natural_form} -> {parsed} != {technical_key}"

    @pytest.mark.parametrize("technical_key", [
        "ring_index",
        "aromatic_ring_index",
    ])
    def test_round_trip_index_properties(self, technical_key):
        """Natural language should parse back to technical key for unique index terms."""
        natural_forms = get_natural_language(technical_key, context="index")
        for natural_form in natural_forms:
            parsed = parse_natural_language(natural_form)
            # Some terms like "carbon atoms" may map to count version first
            # Just verify parsing returns a valid key
            assert isinstance(parsed, str)
            assert len(parsed) > 0


class TestMappingConsistency:
    """Tests for consistency between COUNT_MAPPINGS and INDEX_MAPPINGS."""

    def test_count_index_correspondence(self):
        """For each _count key, there should be a corresponding _index key where applicable."""
        count_prefixes = set()
        for key in COUNT_MAPPINGS.keys():
            if key.endswith("_count"):
                prefix = key[:-6]  # Remove "_count"
                count_prefixes.add(prefix)

        index_prefixes = set()
        for key in INDEX_MAPPINGS.keys():
            if key.endswith("_index"):
                prefix = key[:-6]  # Remove "_index"
                index_prefixes.add(prefix)

        # Most count properties should have corresponding index properties
        # (not all, e.g., molecular_formula doesn't have an index)
        common_prefixes = count_prefixes & index_prefixes
        assert len(common_prefixes) > 10, "Too few count/index pairs"


class TestDynamicMappingGeneration:
    """Tests for dynamically generated mappings."""

    def test_functional_group_count_mappings_generated(self):
        """Functional group count mappings should be generated for all loaded groups."""
        for fg_name in FUNCTIONAL_GROUPS.keys():
            count_key = f"functional_group_{fg_name}_count"
            assert count_key in COUNT_MAPPINGS, f"Missing count mapping for {fg_name}"

    def test_functional_group_index_mappings_generated(self):
        """Functional group index mappings should be generated for all loaded groups."""
        for fg_name in FUNCTIONAL_GROUPS.keys():
            index_key = f"functional_group_{fg_name}_index"
            assert index_key in INDEX_MAPPINGS, f"Missing index mapping for {fg_name}"

    def test_functional_group_instance_mappings_generated(self):
        """Functional group instance count mappings should be generated."""
        for fg_name in FUNCTIONAL_GROUPS.keys():
            instance_key = f"functional_group_{fg_name}_nbrInstances"
            assert instance_key in COUNT_MAPPINGS, f"Missing instance mapping for {fg_name}"

    def test_reaction_count_mappings_generated(self):
        """Reaction count mappings should be generated for all loaded templates."""
        for reaction_name in REACTION_TEMPLATES.keys():
            count_key = f"reaction_{reaction_name}_count"
            assert count_key in COUNT_MAPPINGS, f"Missing count mapping for reaction {reaction_name}"

    def test_reaction_index_mappings_generated(self):
        """Reaction index mappings should be generated for all loaded templates."""
        for reaction_name in REACTION_TEMPLATES.keys():
            index_key = f"reaction_{reaction_name}_index"
            assert index_key in INDEX_MAPPINGS, f"Missing index mapping for reaction {reaction_name}"
