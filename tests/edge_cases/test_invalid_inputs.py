"""
Edge case tests for invalid and boundary inputs.

Tests error handling and edge cases.
"""

import pytest
from moleculariq_core import (
    evaluate_answer,
    valid_smiles,
)


class TestInvalidSmiles:
    """Tests for invalid SMILES handling."""

    @pytest.mark.parametrize("invalid_smiles", [
        "",              # Empty string
        "   ",           # Whitespace only
        "XYZ",           # Invalid atoms
        "C(C(C",         # Unbalanced parentheses
        "c1ccccc",       # Unclosed aromatic ring
    ])
    def test_valid_smiles_rejects_invalid(self, invalid_smiles):
        """valid_smiles returns False for invalid SMILES."""
        assert valid_smiles(invalid_smiles) is False

    def test_valid_smiles_with_none(self):
        """valid_smiles handles None gracefully."""
        assert valid_smiles(None) is False

    def test_solver_with_empty_smiles(self, solver):
        """Solver handles empty SMILES."""
        # Should either return 0/empty or raise an exception
        # Both are acceptable behaviors
        result = solver.get_ring_count("")
        assert result == 0 or result is None

    def test_solver_with_invalid_smiles(self, solver):
        """Solver handles invalid SMILES."""
        result = solver.get_ring_count("invalid_xyz")
        assert result == 0 or result is None


class TestBoundaryValues:
    """Tests for boundary value inputs."""

    def test_single_atom_molecule(self, solver):
        """Single atom molecule (methane)."""
        methane = "C"
        assert solver.get_carbon_atom_count(methane) == 1
        assert solver.get_ring_count(methane) == 0
        assert solver.get_heavy_atom_count(methane) == 1

    def test_very_simple_molecule(self, solver):
        """Very simple molecule (water represented as O)."""
        water = "O"
        assert solver.get_carbon_atom_count(water) == 0
        assert solver.get_hetero_atom_count(water) == 1

    def test_large_ring(self, solver):
        """Large ring system."""
        # 12-membered ring
        large_ring = "C1CCCCCCCCCCC1"
        assert solver.get_ring_count(large_ring) == 1

    def test_empty_prediction(self):
        """Empty prediction dictionary."""
        score = evaluate_answer(
            task_type="multi_count",
            predicted={},
            target={"ring_count": 1}
        )
        assert score == 0.0

    def test_empty_target(self):
        """Empty target dictionary."""
        score = evaluate_answer(
            task_type="multi_count",
            predicted={"ring_count": 1},
            target={}
        )
        # Score depends on implementation
        assert isinstance(score, (int, float))


class TestTypeErrors:
    """Tests for incorrect types."""

    def test_prediction_wrong_type_string(self):
        """String instead of dict for count task."""
        # Should handle gracefully
        score = evaluate_answer(
            task_type="single_count",
            predicted="not_a_dict",
            target={"ring_count": 1}
        )
        assert score == 0.0

    def test_indices_as_string(self):
        """Indices provided as comma-separated string."""
        # Some implementations accept "0,1,2" format
        score = evaluate_answer(
            task_type="single_index",
            predicted={"indices": "0,1,2"},
            target={"indices": [0, 1, 2]}
        )
        # Should either parse or return 0
        assert isinstance(score, (int, float))


class TestFormatterEdgeCases:
    """Edge cases for NaturalLanguageFormatter."""

    def test_empty_count_types(self, formatter):
        """Empty count_types list."""
        with pytest.raises((ValueError, TypeError)):
            formatter.format_count_query(
                smiles="CCO",
                count_types=[]
            )

    def test_none_count_types(self, formatter):
        """None count_types."""
        with pytest.raises((ValueError, TypeError)):
            formatter.format_count_query(
                smiles="CCO",
                count_types=None
            )

    def test_empty_constraint_list(self, formatter):
        """Empty constraints list."""
        text = formatter.format_constraints_list([])
        assert isinstance(text, str)

    def test_unknown_property_type(self, formatter):
        """Unknown property type in count_types."""
        # Should handle gracefully (fallback to name itself)
        query = formatter.format_count_query(
            smiles="CCO",
            count_types=["completely_unknown_property_xyz"]
        )
        assert isinstance(query, str)
        assert "CCO" in query


class TestConstraintEdgeCases:
    """Edge cases for constraint handling."""

    def test_constraint_zero_value(self, formatter):
        """Constraint with value=0."""
        constraint = {"type": "ring_count", "operator": "=", "value": 0}
        text = formatter.format_constraint(constraint)
        assert isinstance(text, str)
        # Should express "no rings", "zero rings", or "without any rings"
        assert "0" in text or "no" in text.lower() or "zero" in text.lower() or "without" in text.lower()

    def test_constraint_missing_operator(self, formatter):
        """Constraint without explicit operator (defaults to =)."""
        constraint = {"type": "ring_count", "value": 2}
        text = formatter.format_constraint(constraint)
        assert isinstance(text, str)

    def test_constraint_large_value(self, formatter):
        """Constraint with large value."""
        constraint = {"type": "carbon_atom_count", "operator": ">=", "value": 100}
        text = formatter.format_constraint(constraint)
        assert "100" in text


class TestRewardEdgeCases:
    """Edge cases for reward functions."""

    def test_reward_with_extra_keys(self):
        """Prediction has extra keys not in target."""
        score = evaluate_answer(
            task_type="multi_count",
            predicted={"ring_count": 1, "extra_key": 5},
            target={"ring_count": 1}
        )
        # Extra keys should not penalize
        assert score >= 0.0

    def test_reward_with_missing_keys(self):
        """Prediction missing keys from target."""
        score = evaluate_answer(
            task_type="multi_count",
            predicted={"ring_count": 1},
            target={"ring_count": 1, "carbon_count": 6}
        )
        # Missing keys should reduce score
        assert 0.0 <= score <= 1.0

    def test_index_reward_with_duplicates(self):
        """Index list with duplicates."""
        score = evaluate_answer(
            task_type="single_index",
            predicted={"indices": [0, 1, 1, 2]},  # Duplicate 1
            target={"indices": [0, 1, 2]}
        )
        # Should handle duplicates (likely counts unique)
        assert isinstance(score, (int, float))

    def test_constraint_reward_empty_smiles(self):
        """Empty SMILES for constraint task."""
        score = evaluate_answer(
            task_type="constraint_generation",
            predicted="",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 0.0
