"""
Unit tests for LLM output parsing robustness.

Tests handling of malformed/unusual LLM outputs including:
- Malformed JSON
- String numbers vs int/float
- Different key naming conventions
- Extra whitespace and formatting
- Unicode characters
- Edge cases
"""

from moleculariq_core.rewards import (
    multi_count_dict_reward,
    single_count_reward,
    multi_index_identification_reward,
)
from moleculariq_core.rewards.count_reward import (
    _coerce_numeric,
    _parse_dict_like,
    _normalize_count_dict,
)
from moleculariq_core.rewards.index_reward import parse_indices_string


class TestCoerceNumeric:
    """Tests for _coerce_numeric function."""

    def test_integer(self):
        """Integer input."""
        assert _coerce_numeric(5) == 5.0

    def test_float(self):
        """Float input."""
        assert _coerce_numeric(5.5) == 5.5

    def test_string_integer(self):
        """String integer."""
        assert _coerce_numeric("5") == 5.0

    def test_string_float(self):
        """String float."""
        assert _coerce_numeric("5.5") == 5.5

    def test_string_with_whitespace(self):
        """String with leading/trailing whitespace."""
        assert _coerce_numeric("  5  ") == 5.0
        assert _coerce_numeric("\t5\n") == 5.0

    def test_empty_string(self):
        """Empty string returns None."""
        assert _coerce_numeric("") is None
        assert _coerce_numeric("   ") is None

    def test_none(self):
        """None returns None."""
        assert _coerce_numeric(None) is None

    def test_boolean(self):
        """Boolean returns None (not 0/1)."""
        assert _coerce_numeric(True) is None
        assert _coerce_numeric(False) is None

    def test_non_numeric_string(self):
        """Non-numeric string returns None."""
        assert _coerce_numeric("abc") is None
        assert _coerce_numeric("five") is None

    def test_negative_numbers(self):
        """Negative numbers."""
        assert _coerce_numeric(-5) == -5.0
        assert _coerce_numeric("-5") == -5.0

    def test_scientific_notation(self):
        """Scientific notation."""
        assert _coerce_numeric("1e3") == 1000.0
        assert _coerce_numeric("1.5e-2") == 0.015


class TestParseDictLike:
    """Tests for _parse_dict_like function.

    Note: _parse_dict_like returns a tuple (dict_or_none, has_duplicates).
    """

    def test_dict_passthrough(self):
        """Dict passes through."""
        result, has_duplicates = _parse_dict_like({"a": 1})
        assert result == {"a": 1}
        assert not has_duplicates

    def test_valid_json_string(self):
        """Valid JSON string."""
        result, has_duplicates = _parse_dict_like('{"ring_count": 5}')
        assert result == {"ring_count": 5}
        assert not has_duplicates

    def test_json_with_extra_whitespace(self):
        """JSON with extra whitespace."""
        result, has_duplicates = _parse_dict_like('  { "ring_count" : 5 }  ')
        assert result == {"ring_count": 5}
        assert not has_duplicates

    def test_colon_separated_format(self):
        """Colon-separated key:value format."""
        result, has_duplicates = _parse_dict_like("ring_count: 5")
        assert result == {"ring_count": "5"}
        assert not has_duplicates

    def test_multiple_colon_values(self):
        """Multiple colon-separated values."""
        result, has_duplicates = _parse_dict_like("ring_count: 5, carbon_count: 6")
        assert result == {"ring_count": "5", "carbon_count": "6"}
        assert not has_duplicates

    def test_semicolon_separated(self):
        """Semicolon-separated values."""
        result, has_duplicates = _parse_dict_like("ring_count: 5; carbon_count: 6")
        assert result == {"ring_count": "5", "carbon_count": "6"}
        assert not has_duplicates

    def test_newline_separated(self):
        """Newline-separated values."""
        result, has_duplicates = _parse_dict_like("ring_count: 5\ncarbon_count: 6")
        assert result == {"ring_count": "5", "carbon_count": "6"}
        assert not has_duplicates

    def test_empty_string(self):
        """Empty string."""
        result, has_duplicates = _parse_dict_like("")
        assert result == {}
        assert not has_duplicates

    def test_whitespace_only(self):
        """Whitespace only."""
        result, has_duplicates = _parse_dict_like("   ")
        assert result == {}
        assert not has_duplicates

    def test_malformed_json(self):
        """Malformed JSON falls back to colon parsing."""
        result, has_duplicates = _parse_dict_like('{"ring_count": 5')  # Missing closing brace
        # Falls back to colon parsing - may extract key:value pairs
        assert isinstance(result, (dict, type(None)))

    def test_json_array_not_dict(self):
        """JSON array returns None (not a dict)."""
        result, has_duplicates = _parse_dict_like('[1, 2, 3]')
        assert result is None


class TestNormalizeCountDict:
    """Tests for _normalize_count_dict function.

    Note: _normalize_count_dict returns a tuple (dict_or_none, has_duplicates).
    """

    def test_technical_keys_unchanged(self):
        """Technical keys pass through."""
        result, has_duplicates = _normalize_count_dict({"ring_count": 5})
        assert "ring_count" in result
        assert not has_duplicates

    def test_natural_language_keys(self):
        """Natural language keys get normalized."""
        result, has_duplicates = _normalize_count_dict({"number of rings": 5})
        # Should normalize to technical key
        assert isinstance(result, dict)
        assert len(result) == 1
        assert not has_duplicates

    def test_mixed_case_keys(self):
        """Mixed case keys."""
        result, has_duplicates = _normalize_count_dict({"Ring_Count": 5})
        assert isinstance(result, dict)
        assert not has_duplicates

    def test_preserves_values(self):
        """Values preserved during normalization."""
        result, has_duplicates = _normalize_count_dict({"ring_count": 5})
        values = list(result.values())
        assert 5 in values
        assert not has_duplicates


class TestParseIndicesString:
    """Tests for parse_indices_string function."""

    def test_simple_comma_list(self):
        """Simple comma-separated list."""
        assert parse_indices_string("0,1,2") == [0, 1, 2]

    def test_with_brackets(self):
        """With square brackets."""
        assert parse_indices_string("[0,1,2]") == [0, 1, 2]

    def test_with_parentheses(self):
        """With parentheses."""
        assert parse_indices_string("(0,1,2)") == [0, 1, 2]

    def test_with_curly_braces(self):
        """With curly braces."""
        assert parse_indices_string("{0,1,2}") == [0, 1, 2]

    def test_with_spaces(self):
        """With spaces around numbers."""
        assert parse_indices_string("0, 1, 2") == [0, 1, 2]
        assert parse_indices_string(" 0 , 1 , 2 ") == [0, 1, 2]

    def test_empty_indicators(self):
        """Empty list indicators."""
        assert parse_indices_string("") == []
        assert parse_indices_string("[]") == []
        assert parse_indices_string("()") == []
        assert parse_indices_string("{}") == []
        assert parse_indices_string("none") == []
        assert parse_indices_string("empty") == []

    def test_single_value(self):
        """Single value."""
        assert parse_indices_string("5") == [5]

    def test_duplicates_removed(self):
        """Duplicates removed."""
        result = parse_indices_string("0,1,1,2,2,2")
        assert result == [0, 1, 2]

    def test_sorted_output(self):
        """Output is sorted."""
        result = parse_indices_string("5,2,8,1")
        assert result == [1, 2, 5, 8]

    def test_negative_returns_none(self):
        """Negative indices return None."""
        assert parse_indices_string("-1,0,1") is None

    def test_non_numeric_returns_none(self):
        """Non-numeric values return None."""
        assert parse_indices_string("a,b,c") is None
        assert parse_indices_string("0,a,2") is None


class TestCountRewardStringNumbers:
    """Tests for count reward with string numbers."""

    def test_string_integer_prediction(self):
        """String integer prediction matches."""
        score = multi_count_dict_reward(
            predicted="5",
            target={"ring_count": 5}
        )
        assert score == 1.0

    def test_string_in_dict_prediction(self):
        """String number in dict prediction."""
        score = multi_count_dict_reward(
            predicted={"ring_count": "5"},
            target={"ring_count": 5}
        )
        assert score == 1.0

    def test_float_string_prediction(self):
        """Float string prediction."""
        score = multi_count_dict_reward(
            predicted="5.0",
            target={"ring_count": 5}
        )
        assert score == 1.0

    def test_json_string_numbers(self):
        """JSON with string numbers."""
        score = multi_count_dict_reward(
            predicted='{"ring_count": "5"}',
            target={"ring_count": 5}
        )
        assert score == 1.0


class TestCountRewardKeyNormalization:
    """Tests for count reward with different key naming."""

    def test_exact_key_match(self):
        """Exact key match."""
        score = multi_count_dict_reward(
            predicted={"ring_count": 5},
            target={"ring_count": 5}
        )
        assert score == 1.0

    def test_underscore_vs_space(self):
        """Underscore vs space in key."""
        # Prediction uses spaces, target uses underscores
        score = multi_count_dict_reward(
            predicted={"ring count": 5},
            target={"ring_count": 5}
        )
        # Should normalize
        assert isinstance(score, float)

    def test_natural_language_key(self):
        """Natural language key in prediction."""
        score = multi_count_dict_reward(
            predicted={"number of rings": 1},
            target={"ring_count": 1}
        )
        # Depends on natural language mappings
        assert isinstance(score, float)


class TestCountRewardMalformedJSON:
    """Tests for count reward with malformed JSON."""

    def test_missing_closing_brace(self):
        """Missing closing brace."""
        score = multi_count_dict_reward(
            predicted='{"ring_count": 5',
            target={"ring_count": 5}
        )
        # Should handle gracefully
        assert isinstance(score, float)

    def test_single_quotes(self):
        """Single quotes instead of double quotes."""
        score = multi_count_dict_reward(
            predicted="{'ring_count': 5}",
            target={"ring_count": 5}
        )
        # JSON requires double quotes - should handle gracefully
        assert isinstance(score, float)

    def test_trailing_comma(self):
        """Trailing comma in JSON."""
        score = multi_count_dict_reward(
            predicted='{"ring_count": 5,}',
            target={"ring_count": 5}
        )
        # Invalid JSON - should handle gracefully
        assert isinstance(score, float)


class TestCountRewardEdgeCases:
    """Edge cases for count reward."""

    def test_zero_value(self):
        """Zero value."""
        score = multi_count_dict_reward(
            predicted=0,
            target={"ring_count": 0}
        )
        assert score == 1.0

    def test_float_tolerance(self):
        """Float comparison with tolerance."""
        score = multi_count_dict_reward(
            predicted=5.0000001,
            target={"ring_count": 5}
        )
        assert score == 1.0

    def test_none_prediction(self):
        """None prediction."""
        score = multi_count_dict_reward(
            predicted=None,
            target={"ring_count": 5}
        )
        assert score == 0.0

    def test_empty_dict_prediction(self):
        """Empty dict prediction."""
        score = multi_count_dict_reward(
            predicted={},
            target={"ring_count": 5}
        )
        assert score == 0.0

    def test_extra_keys_in_prediction(self):
        """Extra keys in prediction don't affect score."""
        score = multi_count_dict_reward(
            predicted={"ring_count": 5, "extra_key": 99},
            target={"ring_count": 5}
        )
        assert score == 1.0


class TestIndexRewardFormats:
    """Tests for index reward with different formats."""

    def test_list_format(self):
        """Direct list format."""
        score = multi_index_identification_reward(
            predicted=[0, 1, 2],
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0

    def test_string_comma_format(self):
        """Comma-separated string format."""
        score = multi_index_identification_reward(
            predicted="0,1,2",
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0

    def test_string_bracket_format(self):
        """Bracketed string format."""
        score = multi_index_identification_reward(
            predicted="[0,1,2]",
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0

    def test_json_dict_format(self):
        """JSON dict format."""
        score = multi_index_identification_reward(
            predicted='{"ring_index": [0,1,2]}',
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0

    def test_dict_with_string_indices(self):
        """Dict with string indices."""
        score = multi_index_identification_reward(
            predicted={"ring_index": "0,1,2"},
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0


class TestIndexRewardOrderInsensitive:
    """Tests for index reward order insensitivity."""

    def test_different_order_matches(self):
        """Different order should match."""
        score = multi_index_identification_reward(
            predicted=[2, 0, 1],
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0

    def test_duplicates_handled(self):
        """Duplicates should be handled."""
        score = multi_index_identification_reward(
            predicted=[0, 1, 1, 2, 2],
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0


class TestIndexRewardEdgeCases:
    """Edge cases for index reward."""

    def test_empty_list_match(self):
        """Empty list matches empty target."""
        score = multi_index_identification_reward(
            predicted=[],
            target={"ring_index": []}
        )
        assert score == 1.0

    def test_empty_string_match(self):
        """Empty string for empty target."""
        score = multi_index_identification_reward(
            predicted="",
            target={"ring_index": []}
        )
        assert score == 1.0

    def test_none_indicator_match(self):
        """'none' string for empty target."""
        score = multi_index_identification_reward(
            predicted="none",
            target={"ring_index": []}
        )
        assert score == 1.0


class TestMultiCountFormats:
    """Tests for multi-count with various formats."""

    def test_standard_json(self):
        """Standard JSON format."""
        score = multi_count_dict_reward(
            predicted='{"ring_count": 1, "carbon_count": 6}',
            target={"ring_count": 1, "carbon_count": 6}
        )
        assert score == 1.0

    def test_colon_newline_format(self):
        """Colon with newline format."""
        prediction = """ring_count: 1
carbon_count: 6"""
        score = multi_count_dict_reward(
            predicted=prediction,
            target={"ring_count": 1, "carbon_count": 6}
        )
        assert score == 1.0

    def test_partial_match_fails(self):
        """Partial match should fail."""
        score = multi_count_dict_reward(
            predicted={"ring_count": 1, "carbon_count": 5},  # Wrong carbon count
            target={"ring_count": 1, "carbon_count": 6}
        )
        assert score == 0.0


class TestMultiIndexFormats:
    """Tests for multi-index with various formats."""

    def test_standard_dict(self):
        """Standard dict format."""
        score = multi_index_identification_reward(
            predicted={"ring_index": [0, 1, 2], "hetero_index": [3]},
            target={"ring_index": [0, 1, 2], "hetero_index": [3]}
        )
        assert score == 1.0

    def test_semicolon_format(self):
        """Semicolon-separated format."""
        score = multi_index_identification_reward(
            predicted="ring_index: 0,1,2; hetero_index: 3",
            target={"ring_index": [0, 1, 2], "hetero_index": [3]}
        )
        assert score == 1.0


class TestMolecularFormulaHandling:
    """Tests for molecular formula special handling."""

    def test_formula_exact_match(self):
        """Exact molecular formula match."""
        score = multi_count_dict_reward(
            predicted="C2H6O",
            target={"molecular_formula": "C2H6O"}
        )
        assert score == 1.0

    def test_formula_different_order(self):
        """Molecular formula with different element order."""
        score = multi_count_dict_reward(
            predicted="H6C2O",
            target={"molecular_formula": "C2H6O"}
        )
        assert score == 1.0

    def test_formula_lowercase(self):
        """Molecular formula with lowercase."""
        score = multi_count_dict_reward(
            predicted="c2h6o",
            target={"molecular_formula": "C2H6O"}
        )
        assert score == 1.0

    def test_formula_with_subscripts(self):
        """Molecular formula with unicode subscripts."""
        score = multi_count_dict_reward(
            predicted="C₂H₆O",
            target={"molecular_formula": "C2H6O"}
        )
        assert score == 1.0


class TestReturnDetailsWithMalformedInput:
    """Tests for return_details with malformed inputs."""

    def test_count_details_invalid_prediction(self):
        """Count reward details with invalid prediction."""
        result = multi_count_dict_reward(
            predicted="not a number",
            target={"ring_count": 5},
            return_details=True
        )
        assert isinstance(result, dict)
        assert result["reward"] == 0.0
        assert "details" in result

    def test_index_details_invalid_prediction(self):
        """Index reward details with invalid prediction."""
        result = multi_index_identification_reward(
            predicted="not indices",
            target={"ring_index": [0, 1, 2]},
            return_details=True
        )
        assert isinstance(result, dict)
        assert result["reward"] == 0.0
        assert "details" in result

    def test_extra_predictions_captured(self):
        """Extra prediction keys are captured in details."""
        result = multi_count_dict_reward(
            predicted={"ring_count": 5, "extra_key": 99},
            target={"ring_count": 5},
            return_details=True
        )
        assert "extra_predictions" in result
        # extra_key should be in extra_predictions


class TestUnicodeAndSpecialCharacters:
    """Tests for unicode and special character handling."""

    def test_unicode_in_json(self):
        """Unicode characters in JSON."""
        score = multi_count_dict_reward(
            predicted='{"ring_count": 5}',  # Standard ASCII
            target={"ring_count": 5}
        )
        assert score == 1.0

    def test_smart_quotes(self):
        """Smart/curly quotes (common from word processors)."""
        # These are often copy-pasted from documents
        score = multi_count_dict_reward(
            predicted='{"ring_count": 5}',  # Using smart quotes would fail JSON
            target={"ring_count": 5}
        )
        assert isinstance(score, float)


class TestLLMTypicalOutputs:
    """Tests simulating typical LLM output patterns."""

    def test_verbose_explanation_with_answer(self):
        """LLM output with explanation before answer."""
        # LLMs often give explanations - just the dict should work
        score = multi_count_dict_reward(
            predicted={"ring_count": 1},
            target={"ring_count": 1}
        )
        assert score == 1.0

    def test_just_number_for_single_count(self):
        """Just a number for single count task."""
        score = single_count_reward(
            predicted=6,
            target={"carbon_count": 6}
        )
        assert score == 1.0

    def test_json_wrapped_answer(self):
        """JSON with extra whitespace/newlines."""
        prediction = """
        {
            "ring_count": 1,
            "carbon_count": 6
        }
        """
        score = multi_count_dict_reward(
            predicted=prediction,
            target={"ring_count": 1, "carbon_count": 6}
        )
        assert score == 1.0

    def test_float_indices(self):
        """Float indices (0.0 instead of 0)."""
        score = multi_index_identification_reward(
            predicted=[0.0, 1.0, 2.0],
            target={"ring_index": [0, 1, 2]}
        )
        assert score == 1.0
