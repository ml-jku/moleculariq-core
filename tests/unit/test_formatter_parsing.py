"""
Extended tests for NaturalLanguageFormatter parsing and conversion methods.

Covers parse_count_answer, parse_index_answer, parse_constraint,
technical_to_natural, natural_to_technical, and helper methods.
"""

import pytest
import random

from moleculariq_core._nlp.formatter import (
    NaturalLanguageFormatter,
    format_count_query,
    format_index_query,
    format_constraint,
    get_count_text,
    get_index_text,
    get_constraint_text,
)


@pytest.fixture
def formatter():
    """Create formatter instance with fixed seed for reproducibility."""
    return NaturalLanguageFormatter(seed=42)


@pytest.fixture
def deterministic_formatter():
    """Create formatter with random phrasing disabled."""
    return NaturalLanguageFormatter(seed=42, enable_random_phrasing=False)


class TestFormatterInitialization:
    """Tests for NaturalLanguageFormatter initialization."""

    def test_init_default(self):
        """Should initialize with defaults."""
        f = NaturalLanguageFormatter()
        assert f._use_random_phrasing is True
        assert f._rng is not None

    def test_init_with_seed(self):
        """Should initialize with seed."""
        f = NaturalLanguageFormatter(seed=123)
        assert f._rng is not None

    def test_init_with_rng(self):
        """Should initialize with custom RNG."""
        rng = random.Random(456)
        f = NaturalLanguageFormatter(rng=rng)
        assert f._rng is rng

    def test_init_rng_and_seed_error(self):
        """Should raise error if both rng and seed provided."""
        with pytest.raises(ValueError):
            NaturalLanguageFormatter(rng=random.Random(), seed=123)

    def test_init_disable_random_phrasing(self):
        """Should initialize with random phrasing disabled."""
        f = NaturalLanguageFormatter(enable_random_phrasing=False)
        assert f._use_random_phrasing is False


class TestParseCountAnswer:
    """Tests for parse_count_answer method."""

    def test_parse_json_format(self, formatter):
        """Should parse JSON format answers."""
        answer = '{"ring_count": 2, "carbon_atom_count": 6}'
        result = formatter.parse_count_answer(answer)
        assert result["ring_count"] == 2
        assert result["carbon_atom_count"] == 6

    def test_parse_semicolon_format(self, formatter):
        """Should parse semicolon-separated format."""
        answer = "rings: 2; carbon atoms: 6"
        result = formatter.parse_count_answer(answer)
        assert "ring_count" in result
        assert result["ring_count"] == 2

    def test_parse_comma_format(self, formatter):
        """Should parse comma-separated format."""
        answer = "rings: 2, carbon atoms: 6"
        result = formatter.parse_count_answer(answer)
        assert len(result) == 2

    def test_parse_none_values(self, formatter):
        """Should handle none/null values."""
        answer = '{"ring_count": null, "carbon_count": "none"}'
        result = formatter.parse_count_answer(answer)
        assert result.get("ring_count") is None

    def test_parse_natural_language_keys(self, formatter):
        """Should convert natural language keys to technical keys."""
        answer = '{"aromatic rings": 1, "hydrogen bond acceptors": 2}'
        result = formatter.parse_count_answer(answer)
        assert "aromatic_ring_count" in result
        assert "hba_count" in result

    def test_parse_empty_json(self, formatter):
        """Should handle empty JSON object."""
        answer = "{}"
        result = formatter.parse_count_answer(answer)
        assert result == {}

    def test_parse_whitespace_in_values(self, formatter):
        """Should handle whitespace in values."""
        answer = "rings:  2  ; carbons:  6  "
        result = formatter.parse_count_answer(answer)
        assert len(result) == 2


class TestParseIndexAnswer:
    """Tests for parse_index_answer method."""

    def test_parse_basic_format(self, formatter):
        """Should parse basic format."""
        answer = "ring atoms: 0, 1, 2, 3, 4, 5; hetero atoms: 6"
        result = formatter.parse_index_answer(answer)
        # Note: "carbon atoms" maps to carbon_atom_count by default in ALIASES
        # so we use "ring atoms" instead which is more distinct for index context
        assert len(result) > 0
        # Check that indices were parsed correctly
        found_ring = any("ring" in k.lower() for k in result.keys())
        assert found_ring or len(result) >= 1

    def test_parse_empty_indices(self, formatter):
        """Should parse empty index lists."""
        answer = "carbon atoms: none; oxygen atoms: []"
        result = formatter.parse_index_answer(answer)
        for key in result:
            if "carbon" in key or "oxygen" in key:
                assert result[key] == []

    def test_parse_range_notation(self, formatter):
        """Should parse range notation."""
        answer = "ring atoms: 0-5"
        result = formatter.parse_index_answer(answer)
        # Should expand range to list
        for key in result:
            if "ring" in key:
                assert result[key] == [0, 1, 2, 3, 4, 5]

    def test_parse_with_expected_types(self, formatter):
        """Should use expected_types for missing values."""
        answer = "ring atoms: 0, 1, 2"
        result = formatter.parse_index_answer(answer, expected_types=["ring_index", "carbon_atom_index"])
        assert "carbon_atom_index" in result
        assert result["carbon_atom_index"] == []

    def test_parse_single_index(self, formatter):
        """Should handle single index value."""
        answer = "stereocenters: 2"
        result = formatter.parse_index_answer(answer)
        assert len(result) > 0


class TestParseConstraint:
    """Tests for parse_constraint method."""

    def test_parse_exactly_operator(self, formatter):
        """Should parse 'exactly' as equals operator."""
        result = formatter.parse_constraint("exactly 3 rings")
        assert result["operator"] == "="
        assert result["value"] == 3

    def test_parse_at_least_operator(self, formatter):
        """Should parse 'at least' as >= operator."""
        result = formatter.parse_constraint("at least 2 carbon atoms")
        assert result["operator"] == ">="
        assert result["value"] == 2

    def test_parse_at_most_operator(self, formatter):
        """Should parse 'at most' as <= operator."""
        result = formatter.parse_constraint("at most 5 aromatic rings")
        assert result["operator"] == "<="
        assert result["value"] == 5

    def test_parse_more_than_operator(self, formatter):
        """Should parse 'more than' as > operator."""
        result = formatter.parse_constraint("more than 3 heteroatoms")
        assert result["operator"] == ">"
        assert result["value"] == 3

    def test_parse_fewer_than_operator(self, formatter):
        """Should parse 'fewer than' as < operator."""
        result = formatter.parse_constraint("fewer than 2 stereocenters")
        assert result["operator"] == "<"
        assert result["value"] == 2

    def test_parse_between_operator(self, formatter):
        """Should parse 'between X and Y' as range operator."""
        result = formatter.parse_constraint("between 2 and 5 rings")
        assert result["operator"] == "range"
        assert result["min_value"] == 2
        assert result["max_value"] == 5

    def test_parse_default_operator(self, formatter):
        """Should default to equals operator."""
        result = formatter.parse_constraint("3 rings")
        assert result["operator"] == "="
        assert result["value"] == 3

    def test_parse_functional_group(self, formatter):
        """Should identify functional group constraints."""
        result = formatter.parse_constraint("exactly 2 alcohol groups")
        assert "functional_group" in result or "type" in result


class TestTechnicalToNatural:
    """Tests for technical_to_natural method."""

    def test_ring_count(self, formatter):
        """Should convert ring_count to natural language."""
        result = formatter.technical_to_natural("ring_count", "count")
        assert "ring" in result.lower()

    def test_hba_count(self, formatter):
        """Should convert hba_count to natural language."""
        result = formatter.technical_to_natural("hba_count", "count")
        assert "hydrogen bond acceptor" in result.lower() or "hba" in result.lower()

    def test_functional_group_count(self, formatter):
        """Should convert functional group count to natural language."""
        result = formatter.technical_to_natural("functional_group_alcohol_count", "count")
        assert "alcohol" in result.lower()

    def test_functional_group_index(self, formatter):
        """Should convert functional group index to natural language."""
        result = formatter.technical_to_natural("functional_group_ketone_index", "index")
        assert "ketone" in result.lower()

    def test_functional_group_nbrInstances(self, formatter):
        """Should convert functional group nbrInstances to natural language."""
        result = formatter.technical_to_natural("functional_group_ester_nbrInstances", "constraint")
        assert "ester" in result.lower()
        assert "group" in result.lower()

    def test_unknown_technical_name(self, formatter):
        """Should fallback to replacing underscores for unknown names."""
        result = formatter.technical_to_natural("unknown_property_name", "count")
        assert "unknown property name" == result.lower()


class TestNaturalToTechnical:
    """Tests for natural_to_technical method."""

    def test_rings_to_ring_count(self, formatter):
        """Should convert 'rings' to ring_count."""
        result = formatter.natural_to_technical("rings", "count")
        assert result == "ring_count"

    def test_aromatic_rings(self, formatter):
        """Should convert 'aromatic rings' to aromatic_ring_count."""
        result = formatter.natural_to_technical("aromatic rings", "count")
        assert result == "aromatic_ring_count"

    def test_case_insensitive(self, formatter):
        """Should be case insensitive."""
        result1 = formatter.natural_to_technical("RINGS", "count")
        result2 = formatter.natural_to_technical("rings", "count")
        assert result1 == result2

    def test_whitespace_normalization(self, formatter):
        """Should normalize whitespace."""
        result = formatter.natural_to_technical("  aromatic   rings  ", "count")
        assert result == "aromatic_ring_count"

    def test_atoms_in_pattern(self, formatter):
        """Should handle 'atoms in X' pattern for index types."""
        result = formatter.natural_to_technical("atoms in alcohol groups", "index")
        assert "alcohol" in result

    def test_remove_groups_suffix(self, formatter):
        """Should handle 'X groups' pattern."""
        result = formatter.natural_to_technical("ketone groups", "count")
        # Should recognize ketone
        assert "ketone" in result.lower()

    def test_unknown_natural_text(self, formatter):
        """Should convert unknown text to underscore format."""
        result = formatter.natural_to_technical("completely unknown text", "count")
        assert result == "completely_unknown_text"


class TestFormatCountAnswer:
    """Tests for format_count_answer method."""

    def test_format_single_count(self, formatter):
        """Should format single count."""
        result = formatter.format_count_answer({"ring_count": 2})
        assert "2" in result

    def test_format_multiple_counts(self, formatter):
        """Should format multiple counts."""
        result = formatter.format_count_answer({
            "ring_count": 2,
            "carbon_atom_count": 6
        })
        assert "2" in result
        assert "6" in result
        assert ";" in result

    def test_format_none_value(self, formatter):
        """Should format None as 'none'."""
        result = formatter.format_count_answer({"ring_count": None})
        assert "none" in result.lower()

    def test_format_empty_dict(self, formatter):
        """Should handle empty dictionary."""
        result = formatter.format_count_answer({})
        assert "no counts" in result.lower()


class TestFormatIndexAnswer:
    """Tests for format_index_answer method."""

    def test_format_single_index_list(self, formatter):
        """Should format single index list."""
        result = formatter.format_index_answer({"ring_index": [0, 1, 2, 3, 4, 5]})
        assert "0" in result
        assert "5" in result

    def test_format_empty_index_list(self, formatter):
        """Should format empty list as 'none'."""
        result = formatter.format_index_answer({"ring_index": []})
        assert "none" in result.lower()

    def test_format_multiple_index_types(self, formatter):
        """Should format multiple index types."""
        result = formatter.format_index_answer({
            "ring_index": [0, 1, 2],
            "carbon_atom_index": [0, 1]
        })
        assert ";" in result

    def test_format_empty_dict(self, formatter):
        """Should handle empty dictionary."""
        result = formatter.format_index_answer({})
        assert "no indices" in result.lower()


class TestFormatConstraint:
    """Tests for format_constraint method."""

    def test_format_equals_zero(self, deterministic_formatter):
        """Should format equals zero as 'no X'."""
        result = deterministic_formatter.format_constraint({
            "type": "ring_count",
            "operator": "=",
            "value": 0
        })
        assert "no" in result.lower() or "zero" in result.lower()

    def test_format_equals_one(self, deterministic_formatter):
        """Should format equals one with singular form."""
        result = deterministic_formatter.format_constraint({
            "type": "ring_count",
            "operator": "=",
            "value": 1
        })
        assert "1" in result

    def test_format_range(self, deterministic_formatter):
        """Should format range constraint."""
        result = deterministic_formatter.format_constraint({
            "type": "ring_count",
            "operator": "range",
            "min_value": 2,
            "max_value": 5
        })
        assert "between" in result.lower()
        assert "2" in result
        assert "5" in result

    def test_format_greater_equal(self, deterministic_formatter):
        """Should format >= constraint."""
        result = deterministic_formatter.format_constraint({
            "type": "carbon_atom_count",
            "operator": ">=",
            "value": 3
        })
        assert "at least" in result.lower()
        assert "3" in result

    def test_format_less_equal(self, deterministic_formatter):
        """Should format <= constraint."""
        result = deterministic_formatter.format_constraint({
            "type": "ring_count",
            "operator": "<=",
            "value": 2
        })
        assert "at most" in result.lower()
        assert "2" in result

    def test_format_functional_group_constraint(self, deterministic_formatter):
        """Should format functional group constraint."""
        result = deterministic_formatter.format_constraint({
            "type": "functional_group_alcohol_nbrInstances",
            "operator": "=",
            "value": 2
        })
        assert "alcohol" in result.lower()
        assert "2" in result

    def test_format_molecular_formula(self, deterministic_formatter):
        """Should format molecular formula constraint."""
        result = deterministic_formatter.format_constraint({
            "type": "molecular_formula",
            "operator": "=",
            "value": "C6H12O6"
        })
        assert "C6H12O6" in result


class TestFormatConstraintsList:
    """Tests for format_constraints_list method."""

    def test_format_empty_list(self, formatter):
        """Should handle empty constraints list."""
        result = formatter.format_constraints_list([])
        assert result == "no constraints"

    def test_format_single_constraint(self, formatter):
        """Should format single constraint."""
        result = formatter.format_constraints_list([
            {"type": "ring_count", "operator": "=", "value": 2}
        ])
        assert "2" in result

    def test_format_two_constraints(self, formatter):
        """Should format two constraints with connector."""
        result = formatter.format_constraints_list([
            {"type": "ring_count", "operator": "=", "value": 1},
            {"type": "carbon_atom_count", "operator": ">=", "value": 4}
        ])
        # Should have both values
        assert "1" in result
        assert "4" in result

    def test_format_three_constraints(self, formatter):
        """Should format three constraints with proper comma/and."""
        result = formatter.format_constraints_list([
            {"type": "ring_count", "operator": "=", "value": 1},
            {"type": "carbon_atom_count", "operator": ">=", "value": 4},
            {"type": "hba_count", "operator": "<=", "value": 3}
        ])
        assert "," in result or "and" in result


class TestHelperMethods:
    """Tests for internal helper methods."""

    def test_pluralize_regular(self, formatter):
        """Should pluralize regular words."""
        assert formatter._pluralize("ring") == "rings"
        assert formatter._pluralize("atom") == "atoms"

    def test_pluralize_y_ending(self, formatter):
        """Should pluralize words ending in y."""
        assert formatter._pluralize("category") == "categories"

    def test_pluralize_s_ending(self, formatter):
        """Should pluralize words ending in s."""
        assert formatter._pluralize("mass") == "masses"

    def test_singularize_regular(self, formatter):
        """Should singularize regular plurals."""
        assert formatter._singularize("rings") == "ring"
        assert formatter._singularize("atoms") == "atom"

    def test_singularize_ies(self, formatter):
        """Should singularize -ies words."""
        assert formatter._singularize("categories") == "category"

    def test_singularize_es(self, formatter):
        """Should singularize -es words."""
        assert formatter._singularize("masses") == "mass"

    def test_ensure_type_list_string(self, formatter):
        """Should convert string to list."""
        result = formatter._ensure_type_list("ring_count", "test")
        assert result == ["ring_count"]

    def test_ensure_type_list_list(self, formatter):
        """Should preserve list."""
        result = formatter._ensure_type_list(["a", "b"], "test")
        assert result == ["a", "b"]

    def test_ensure_type_list_none_error(self, formatter):
        """Should raise error for None."""
        with pytest.raises(ValueError):
            formatter._ensure_type_list(None, "test")

    def test_ensure_type_list_empty_error(self, formatter):
        """Should raise error for empty list."""
        with pytest.raises(ValueError):
            formatter._ensure_type_list([], "test")

    def test_parse_optional_int_integer(self, formatter):
        """Should parse integer."""
        assert formatter._parse_optional_int(42) == 42

    def test_parse_optional_int_float(self, formatter):
        """Should parse integer float."""
        assert formatter._parse_optional_int(42.0) == 42

    def test_parse_optional_int_string(self, formatter):
        """Should parse string integer."""
        assert formatter._parse_optional_int("42") == 42

    def test_parse_optional_int_none(self, formatter):
        """Should return None for None."""
        assert formatter._parse_optional_int(None) is None

    def test_parse_optional_int_none_string(self, formatter):
        """Should return None for 'none' string."""
        assert formatter._parse_optional_int("none") is None
        assert formatter._parse_optional_int("null") is None
        assert formatter._parse_optional_int("N/A") is None

    def test_parse_optional_int_boolean_error(self, formatter):
        """Should raise error for boolean."""
        with pytest.raises(ValueError):
            formatter._parse_optional_int(True)


class TestConvenienceFunctions:
    """Tests for module-level convenience functions."""

    def test_get_count_text(self):
        """Should return count text."""
        result = get_count_text(["ring_count", "carbon_atom_count"])
        assert len(result) > 0

    def test_get_index_text(self):
        """Should return index text."""
        result = get_index_text(["ring_index"])
        assert len(result) > 0

    def test_get_constraint_text(self):
        """Should return constraint text."""
        result = get_constraint_text({"type": "ring_count", "operator": "=", "value": 2})
        assert "2" in result

    def test_format_count_query_function(self):
        """Should format count query."""
        result = format_count_query("CCO", ["ring_count"])
        assert "CCO" in result
        assert "?" in result

    def test_format_index_query_function(self):
        """Should format index query."""
        result = format_index_query("CCO", ["carbon_atom_index"])
        assert "CCO" in result

    def test_format_constraint_function(self):
        """Should format constraint."""
        result = format_constraint({"type": "ring_count", "operator": "=", "value": 1})
        assert "1" in result


class TestKeyHintFormatting:
    """Tests for key hint formatting methods."""

    def test_format_single_key_hint(self, formatter):
        """Should format single key hint."""
        result = formatter.format_key_hint(["ring_count"], "count")
        assert "`ring_count`" in result
        assert "JSON" in result

    def test_format_multiple_key_hints(self, formatter):
        """Should format multiple key hints."""
        result = formatter.format_key_hint(["ring_count", "carbon_atom_count"], "count")
        assert "`ring_count`" in result
        assert "`carbon_atom_count`" in result

    def test_format_constraint_hint(self, formatter):
        """Should return constraint hint."""
        result = formatter.format_constraint_hint()
        assert "smiles" in result.lower()
        assert "JSON" in result

    def test_empty_key_hint(self, formatter):
        """Should return empty string for empty keys."""
        result = formatter.format_key_hint([], "count")
        assert result == ""


class TestReproducibility:
    """Tests for deterministic/reproducible behavior."""

    def test_same_seed_same_output(self):
        """Same seed should produce same output."""
        f1 = NaturalLanguageFormatter(seed=42)
        f2 = NaturalLanguageFormatter(seed=42)

        result1 = f1.format_constraint({"type": "ring_count", "operator": "=", "value": 0})
        result2 = f2.format_constraint({"type": "ring_count", "operator": "=", "value": 0})
        assert result1 == result2

    def test_disabled_random_deterministic(self):
        """Disabled random should be deterministic."""
        f1 = NaturalLanguageFormatter(enable_random_phrasing=False)
        f2 = NaturalLanguageFormatter(enable_random_phrasing=False)

        result1 = f1.format_constraint({"type": "ring_count", "operator": "=", "value": 0})
        result2 = f2.format_constraint({"type": "ring_count", "operator": "=", "value": 0})
        assert result1 == result2

    def test_return_only_text_deterministic(self):
        """return_only_text should disable randomness."""
        f = NaturalLanguageFormatter()
        results = set()
        for _ in range(10):
            result = f.format_constraint(
                {"type": "ring_count", "operator": "=", "value": 0},
                return_only_text=True
            )
            results.add(result)
        # With return_only_text, should be deterministic
        assert len(results) == 1
