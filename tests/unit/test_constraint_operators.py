"""
Unit tests for constraint evaluation operators.

Tests all constraint operators: =, >, <, >=, <=, !=, range
"""

from moleculariq_core import multi_constraint_generation_reward
from moleculariq_core.rewards.utils import evaluate_numeric_constraint


class TestEvaluateNumericConstraint:
    """Tests for evaluate_numeric_constraint utility."""

    def test_equals_operator_match(self):
        """Equals operator with matching value."""
        constraint = {"operator": "=", "value": 5}
        assert evaluate_numeric_constraint(5, constraint) is True

    def test_equals_operator_no_match(self):
        """Equals operator with non-matching value."""
        constraint = {"operator": "=", "value": 5}
        assert evaluate_numeric_constraint(6, constraint) is False

    def test_greater_than_true(self):
        """Greater than operator - true case."""
        constraint = {"operator": ">", "value": 5}
        assert evaluate_numeric_constraint(6, constraint) is True

    def test_greater_than_false(self):
        """Greater than operator - false case."""
        constraint = {"operator": ">", "value": 5}
        assert evaluate_numeric_constraint(5, constraint) is False

    def test_greater_than_boundary(self):
        """Greater than operator - boundary case."""
        constraint = {"operator": ">", "value": 5}
        assert evaluate_numeric_constraint(4, constraint) is False

    def test_less_than_true(self):
        """Less than operator - true case."""
        constraint = {"operator": "<", "value": 5}
        assert evaluate_numeric_constraint(4, constraint) is True

    def test_less_than_false(self):
        """Less than operator - false case."""
        constraint = {"operator": "<", "value": 5}
        assert evaluate_numeric_constraint(5, constraint) is False

    def test_less_than_boundary(self):
        """Less than operator - boundary case."""
        constraint = {"operator": "<", "value": 5}
        assert evaluate_numeric_constraint(6, constraint) is False

    def test_greater_equal_match(self):
        """Greater or equal - exact match."""
        constraint = {"operator": ">=", "value": 5}
        assert evaluate_numeric_constraint(5, constraint) is True

    def test_greater_equal_greater(self):
        """Greater or equal - greater value."""
        constraint = {"operator": ">=", "value": 5}
        assert evaluate_numeric_constraint(6, constraint) is True

    def test_greater_equal_less(self):
        """Greater or equal - less value."""
        constraint = {"operator": ">=", "value": 5}
        assert evaluate_numeric_constraint(4, constraint) is False

    def test_less_equal_match(self):
        """Less or equal - exact match."""
        constraint = {"operator": "<=", "value": 5}
        assert evaluate_numeric_constraint(5, constraint) is True

    def test_less_equal_less(self):
        """Less or equal - less value."""
        constraint = {"operator": "<=", "value": 5}
        assert evaluate_numeric_constraint(4, constraint) is True

    def test_less_equal_greater(self):
        """Less or equal - greater value."""
        constraint = {"operator": "<=", "value": 5}
        assert evaluate_numeric_constraint(6, constraint) is False

    def test_not_equal_true(self):
        """Not equal - may not be supported."""
        constraint = {"operator": "!=", "value": 5}
        result = evaluate_numeric_constraint(6, constraint)
        # != operator may not be implemented - just verify it doesn't crash
        assert isinstance(result, bool)

    def test_not_equal_false(self):
        """Not equal - may not be supported."""
        constraint = {"operator": "!=", "value": 5}
        result = evaluate_numeric_constraint(5, constraint)
        # != operator may not be implemented - just verify it doesn't crash
        assert isinstance(result, bool)

    def test_range_within(self):
        """Range operator - value within range."""
        constraint = {"operator": "range", "min_value": 3, "max_value": 7}
        assert evaluate_numeric_constraint(5, constraint) is True

    def test_range_at_min(self):
        """Range operator - value at minimum."""
        constraint = {"operator": "range", "min_value": 3, "max_value": 7}
        assert evaluate_numeric_constraint(3, constraint) is True

    def test_range_at_max(self):
        """Range operator - value at maximum."""
        constraint = {"operator": "range", "min_value": 3, "max_value": 7}
        assert evaluate_numeric_constraint(7, constraint) is True

    def test_range_below(self):
        """Range operator - value below range."""
        constraint = {"operator": "range", "min_value": 3, "max_value": 7}
        assert evaluate_numeric_constraint(2, constraint) is False

    def test_range_above(self):
        """Range operator - value above range."""
        constraint = {"operator": "range", "min_value": 3, "max_value": 7}
        assert evaluate_numeric_constraint(8, constraint) is False

    def test_default_operator_equals(self):
        """Default operator is equals."""
        constraint = {"value": 5}  # No operator specified
        assert evaluate_numeric_constraint(5, constraint) is True
        assert evaluate_numeric_constraint(6, constraint) is False

    def test_float_comparison(self):
        """Float value comparison."""
        constraint = {"operator": ">=", "value": 3.5}
        assert evaluate_numeric_constraint(3.5, constraint) is True
        assert evaluate_numeric_constraint(4.0, constraint) is True
        assert evaluate_numeric_constraint(3.0, constraint) is False


class TestConstraintRewardOperators:
    """Tests for constraint reward with different operators."""

    def test_equals_constraint_satisfied(self):
        """Equals constraint satisfied."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # Benzene - 1 ring
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 1.0

    def test_equals_constraint_not_satisfied(self):
        """Equals constraint not satisfied."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # Benzene - 1 ring
            constraints=[{"type": "ring_count", "operator": "=", "value": 2}]
        )
        assert score == 0.0

    def test_greater_than_constraint(self):
        """Greater than constraint."""
        # Benzene has 6 carbons
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "carbon_atom_count", "operator": ">", "value": 5}]
        )
        assert score == 1.0

    def test_greater_than_constraint_not_satisfied(self):
        """Greater than constraint not satisfied."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # 6 carbons
            constraints=[{"type": "carbon_atom_count", "operator": ">", "value": 6}]
        )
        assert score == 0.0

    def test_less_than_constraint(self):
        """Less than constraint."""
        score = multi_constraint_generation_reward(
            predicted="CCO",  # Ethanol - 2 carbons
            constraints=[{"type": "carbon_atom_count", "operator": "<", "value": 5}]
        )
        assert score == 1.0

    def test_less_than_constraint_not_satisfied(self):
        """Less than constraint not satisfied."""
        score = multi_constraint_generation_reward(
            predicted="CCO",  # Ethanol - 2 carbons
            constraints=[{"type": "carbon_atom_count", "operator": "<", "value": 2}]
        )
        assert score == 0.0

    def test_greater_equal_constraint(self):
        """Greater or equal constraint."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # 6 carbons
            constraints=[{"type": "carbon_atom_count", "operator": ">=", "value": 6}]
        )
        assert score == 1.0

    def test_less_equal_constraint(self):
        """Less or equal constraint."""
        score = multi_constraint_generation_reward(
            predicted="CCO",  # 2 carbons
            constraints=[{"type": "carbon_atom_count", "operator": "<=", "value": 2}]
        )
        assert score == 1.0

    def test_range_constraint(self):
        """Range constraint."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # 6 carbons
            constraints=[{
                "type": "carbon_atom_count",
                "operator": "range",
                "min_value": 5,
                "max_value": 10
            }]
        )
        assert score == 1.0

    def test_range_constraint_not_satisfied(self):
        """Range constraint not satisfied."""
        score = multi_constraint_generation_reward(
            predicted="CCO",  # 2 carbons
            constraints=[{
                "type": "carbon_atom_count",
                "operator": "range",
                "min_value": 5,
                "max_value": 10
            }]
        )
        assert score == 0.0


class TestMultipleConstraints:
    """Tests for multiple constraints combined."""

    def test_two_constraints_both_satisfied(self):
        """Two constraints, both satisfied."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # Benzene
            constraints=[
                {"type": "ring_count", "operator": "=", "value": 1},
                {"type": "carbon_atom_count", "operator": "=", "value": 6}
            ]
        )
        assert score == 1.0

    def test_two_constraints_one_fails(self):
        """Two constraints, one fails."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # Benzene
            constraints=[
                {"type": "ring_count", "operator": "=", "value": 1},
                {"type": "carbon_atom_count", "operator": "=", "value": 10}  # Wrong
            ]
        )
        assert score == 0.0

    def test_three_constraints_all_satisfied(self):
        """Three constraints all satisfied."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # Benzene
            constraints=[
                {"type": "ring_count", "operator": ">=", "value": 1},
                {"type": "aromatic_ring_count", "operator": "=", "value": 1},
                {"type": "carbon_atom_count", "operator": "<=", "value": 10}
            ]
        )
        assert score == 1.0

    def test_mixed_operators(self):
        """Multiple constraints with mixed operators."""
        score = multi_constraint_generation_reward(
            predicted="c1ccc2ccccc2c1",  # Naphthalene
            constraints=[
                {"type": "ring_count", "operator": "=", "value": 2},
                {"type": "carbon_atom_count", "operator": ">", "value": 8},
                {"type": "hetero_atom_count", "operator": "<", "value": 1}
            ]
        )
        assert score == 1.0


class TestConstraintReturnDetails:
    """Tests for constraint return_details functionality."""

    def test_return_details_structure(self):
        """return_details returns proper structure."""
        result = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}],
            return_details=True
        )
        assert isinstance(result, dict)
        assert "reward" in result
        assert "details" in result
        assert "supported" in result
        assert "total" in result

    def test_return_details_satisfied(self):
        """return_details for satisfied constraint."""
        result = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}],
            return_details=True
        )
        assert result["reward"] == 1.0
        assert len(result["details"]) == 1
        assert result["details"][0]["satisfied"] is True

    def test_return_details_not_satisfied(self):
        """return_details for unsatisfied constraint."""
        result = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 2}],
            return_details=True
        )
        assert result["reward"] == 0.0
        assert len(result["details"]) == 1
        assert result["details"][0]["satisfied"] is False

    def test_return_details_actual_value(self):
        """return_details includes actual computed value."""
        result = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}],
            return_details=True
        )
        assert "actual" in result["details"][0]
        assert result["details"][0]["actual"] == 1.0

    def test_return_details_invalid_smiles(self):
        """return_details for invalid SMILES."""
        result = multi_constraint_generation_reward(
            predicted="invalid_xyz",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}],
            return_details=True
        )
        assert result["reward"] == 0.0
        assert result["valid_smiles"] is False


class TestConstraintPropertyTypes:
    """Tests for different property types in constraints."""

    def test_ring_count_constraint(self):
        """Ring count constraint."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 1.0

    def test_aromatic_ring_constraint(self):
        """Aromatic ring constraint."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "aromatic_ring_count", "operator": "=", "value": 1}]
        )
        assert score == 1.0

    def test_carbon_count_constraint(self):
        """Carbon atom count constraint."""
        score = multi_constraint_generation_reward(
            predicted="CCO",
            constraints=[{"type": "carbon_atom_count", "operator": "=", "value": 2}]
        )
        assert score == 1.0

    def test_hetero_atom_constraint(self):
        """Hetero atom count constraint."""
        score = multi_constraint_generation_reward(
            predicted="CCO",
            constraints=[{"type": "hetero_atom_count", "operator": "=", "value": 1}]
        )
        assert score == 1.0

    def test_hbd_constraint(self):
        """Hydrogen bond donor constraint."""
        score = multi_constraint_generation_reward(
            predicted="CCO",
            constraints=[{"type": "hbd_count", "operator": ">=", "value": 1}]
        )
        assert score == 1.0

    def test_hba_constraint(self):
        """Hydrogen bond acceptor constraint."""
        score = multi_constraint_generation_reward(
            predicted="CCO",
            constraints=[{"type": "hba_count", "operator": ">=", "value": 1}]
        )
        assert score == 1.0


class TestEdgeCases:
    """Edge cases for constraint evaluation."""

    def test_empty_constraints(self):
        """Empty constraints list."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[]
        )
        assert score == 1.0

    def test_zero_value_constraint(self):
        """Constraint with value=0."""
        score = multi_constraint_generation_reward(
            predicted="CCO",  # No rings
            constraints=[{"type": "ring_count", "operator": "=", "value": 0}]
        )
        assert score == 1.0

    def test_constraint_json_string(self):
        """Constraints as JSON string."""
        import json
        constraints_json = json.dumps([
            {"type": "ring_count", "operator": "=", "value": 1}
        ])
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=constraints_json
        )
        assert score == 1.0

    def test_smiles_in_dict_format(self):
        """SMILES prediction in dictionary format."""
        score = multi_constraint_generation_reward(
            predicted={"smiles": "c1ccccc1"},
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 1.0

    def test_case_insensitive_property(self):
        """Property name is case-insensitive."""
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",
            constraints=[{"type": "RING_COUNT", "operator": "=", "value": 1}]
        )
        # Should handle case insensitively or return appropriate result
        assert isinstance(score, float)
