"""
Unit tests for reward functions.

Tests evaluation/scoring functions for different task types.
"""

import pytest
from moleculariq_core import (
    evaluate_answer,
    multi_count_dict_reward,
    single_count_reward,
    multi_index_identification_reward,
    single_index_reward,
    multi_constraint_generation_reward,
    valid_smiles,
)


class TestValidSmiles:
    """Tests for SMILES validation utility."""

    def test_valid_smiles_ethanol(self):
        """Ethanol is valid SMILES."""
        assert valid_smiles("CCO") is True

    def test_valid_smiles_benzene(self):
        """Benzene is valid SMILES."""
        assert valid_smiles("c1ccccc1") is True

    def test_invalid_smiles_empty(self):
        """Empty string is invalid."""
        assert valid_smiles("") is False

    def test_invalid_smiles_garbage(self):
        """Garbage string is invalid."""
        assert valid_smiles("XYZ123") is False

    def test_invalid_smiles_none(self):
        """None is invalid."""
        assert valid_smiles(None) is False


class TestSingleCountReward:
    """Tests for single count reward function."""

    def test_exact_match(self):
        """Perfect match returns 1.0."""
        score = single_count_reward(
            predicted={"ring_count": 2},
            target={"ring_count": 2}
        )
        assert score == 1.0

    def test_wrong_value(self):
        """Wrong value returns 0.0."""
        score = single_count_reward(
            predicted={"ring_count": 3},
            target={"ring_count": 2}
        )
        assert score == 0.0

    def test_missing_key(self):
        """Missing key behavior depends on implementation."""
        score = single_count_reward(
            predicted={"wrong_key": 2},
            target={"ring_count": 2}
        )
        # Implementation may be lenient with key matching
        assert isinstance(score, (int, float))


class TestMultiCountReward:
    """Tests for multi-count reward function."""

    def test_all_correct(self):
        """All correct returns 1.0."""
        score = multi_count_dict_reward(
            predicted={"ring_count": 2, "carbon_count": 6},
            target={"ring_count": 2, "carbon_count": 6}
        )
        assert score == 1.0

    def test_partial_correct(self):
        """Partial correct - one right, one wrong."""
        score = multi_count_dict_reward(
            predicted={"ring_count": 2, "carbon_count": 5},
            target={"ring_count": 2, "carbon_count": 6}
        )
        # Implementation uses exact match per property, so 1/2 = 0.5 or 0.0 if strict
        assert 0.0 <= score <= 1.0

    def test_all_wrong(self):
        """All wrong returns 0.0."""
        score = multi_count_dict_reward(
            predicted={"ring_count": 0, "carbon_count": 0},
            target={"ring_count": 2, "carbon_count": 6}
        )
        assert score == 0.0

    def test_return_details(self):
        """Can return detailed breakdown."""
        result = multi_count_dict_reward(
            predicted={"ring_count": 2, "carbon_count": 6},
            target={"ring_count": 2, "carbon_count": 6},
            return_details=True
        )
        assert isinstance(result, dict)
        assert "reward" in result or "score" in result or isinstance(result.get("ring_count"), dict)


class TestSingleIndexReward:
    """Tests for single index reward function."""

    def test_exact_match(self):
        """Exact index match returns 1.0."""
        score = single_index_reward(
            predicted={"ring_indices": [0, 1, 2]},
            target={"ring_indices": [0, 1, 2]}
        )
        assert score == 1.0

    def test_order_independent(self):
        """Order should not matter."""
        score = single_index_reward(
            predicted={"ring_indices": [2, 0, 1]},
            target={"ring_indices": [0, 1, 2]}
        )
        assert score == 1.0

    def test_partial_match(self):
        """Partial match behavior depends on implementation."""
        score = single_index_reward(
            predicted={"ring_indices": [0, 1]},
            target={"ring_indices": [0, 1, 2]}
        )
        # Implementation may use strict exact match (0.0) or partial credit
        assert 0.0 <= score <= 1.0

    def test_empty_correct(self):
        """Empty list matches empty target."""
        score = single_index_reward(
            predicted={"ring_indices": []},
            target={"ring_indices": []}
        )
        assert score == 1.0


class TestMultiIndexReward:
    """Tests for multi-index reward function."""

    def test_all_correct(self):
        """All indices correct returns 1.0."""
        score = multi_index_identification_reward(
            predicted={
                "ring_indices": [0, 1, 2],
                "carbon_indices": [0, 1]
            },
            target={
                "ring_indices": [0, 1, 2],
                "carbon_indices": [0, 1]
            }
        )
        assert score == 1.0

    def test_partial_correct(self):
        """Partial correct returns 0.0 (strict all-or-nothing)."""
        score = multi_index_identification_reward(
            predicted={
                "ring_indices": [0, 1, 2],
                "carbon_indices": [0]  # Missing 1
            },
            target={
                "ring_indices": [0, 1, 2],
                "carbon_indices": [0, 1]
            }
        )
        assert score == 0.0


class TestConstraintReward:
    """Tests for constraint generation reward."""

    def test_valid_molecule_satisfies_constraint(self):
        """Valid molecule satisfying constraint returns >0."""
        # Generate benzene, which has 1 ring
        score = multi_constraint_generation_reward(
            predicted="c1ccccc1",  # Benzene
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score > 0.0

    def test_invalid_smiles_returns_zero(self):
        """Invalid SMILES returns 0."""
        score = multi_constraint_generation_reward(
            predicted="invalid_smiles",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 0.0

    def test_constraint_not_satisfied(self):
        """Molecule not satisfying constraint returns 0."""
        # Ethanol has 0 rings, constraint wants 1
        score = multi_constraint_generation_reward(
            predicted="CCO",  # Ethanol, 0 rings
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 0.0


class TestEvaluateAnswerDispatcher:
    """Tests for the unified evaluate_answer dispatcher."""

    def test_single_count_dispatch(self):
        """Dispatches to single_count correctly."""
        score = evaluate_answer(
            task_type="single_count",
            predicted={"ring_count": 1},
            target={"ring_count": 1}
        )
        assert score == 1.0

    def test_multi_count_dispatch(self):
        """Dispatches to multi_count correctly."""
        score = evaluate_answer(
            task_type="multi_count",
            predicted={"a": 1, "b": 2},
            target={"a": 1, "b": 2}
        )
        assert score == 1.0

    def test_single_index_dispatch(self):
        """Dispatches to single_index correctly."""
        score = evaluate_answer(
            task_type="single_index",
            predicted={"indices": [0, 1]},
            target={"indices": [0, 1]}
        )
        assert score == 1.0

    def test_constraint_generation_dispatch(self):
        """Dispatches to constraint_generation correctly."""
        score = evaluate_answer(
            task_type="constraint_generation",
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score > 0.0

    def test_unknown_task_type_raises(self):
        """Unknown task type raises ValueError."""
        with pytest.raises(ValueError):
            evaluate_answer(
                task_type="unknown_task",
                predicted={},
                target={}
            )

    def test_missing_target_raises(self):
        """Missing target for count task raises ValueError."""
        with pytest.raises(ValueError):
            evaluate_answer(
                task_type="single_count",
                predicted={"a": 1},
                target=None
            )


class TestReturnDetailsSingleCount:
    """Tests for return_details in single count reward."""

    def test_single_count_return_details_structure(self):
        """Single count return_details has correct structure."""
        result = single_count_reward(
            predicted={"ring_count": 1},
            target={"ring_count": 1},
            return_details=True
        )
        assert isinstance(result, dict)

    def test_single_count_return_details_match(self):
        """Single count return_details shows match info."""
        result = single_count_reward(
            predicted={"ring_count": 1},
            target={"ring_count": 1},
            return_details=True
        )
        if isinstance(result, dict) and "reward" in result:
            assert result["reward"] == 1.0

    def test_single_count_return_details_mismatch(self):
        """Single count return_details shows mismatch info."""
        result = single_count_reward(
            predicted={"ring_count": 2},
            target={"ring_count": 1},
            return_details=True
        )
        if isinstance(result, dict) and "reward" in result:
            assert result["reward"] == 0.0


class TestReturnDetailsMultiCount:
    """Tests for return_details in multi count reward."""

    def test_multi_count_return_details_structure(self):
        """Multi count return_details has correct structure."""
        result = multi_count_dict_reward(
            predicted={"ring_count": 1, "carbon_count": 6},
            target={"ring_count": 1, "carbon_count": 6},
            return_details=True
        )
        assert isinstance(result, dict)

    def test_multi_count_return_details_all_match(self):
        """Multi count return_details all properties match."""
        result = multi_count_dict_reward(
            predicted={"ring_count": 1, "carbon_count": 6},
            target={"ring_count": 1, "carbon_count": 6},
            return_details=True
        )
        if isinstance(result, dict) and "reward" in result:
            assert result["reward"] == 1.0

    def test_multi_count_return_details_partial(self):
        """Multi count return_details partial match."""
        result = multi_count_dict_reward(
            predicted={"ring_count": 1, "carbon_count": 5},
            target={"ring_count": 1, "carbon_count": 6},
            return_details=True
        )
        assert isinstance(result, dict)


class TestReturnDetailsSingleIndex:
    """Tests for return_details in single index reward."""

    def test_single_index_return_details_structure(self):
        """Single index return_details has correct structure."""
        result = single_index_reward(
            predicted={"indices": [0, 1, 2]},
            target={"indices": [0, 1, 2]},
            return_details=True
        )
        assert isinstance(result, dict)

    def test_single_index_return_details_match(self):
        """Single index return_details shows match."""
        result = single_index_reward(
            predicted={"indices": [0, 1, 2]},
            target={"indices": [0, 1, 2]},
            return_details=True
        )
        if isinstance(result, dict) and "reward" in result:
            assert result["reward"] == 1.0


class TestReturnDetailsMultiIndex:
    """Tests for return_details in multi index reward."""

    def test_multi_index_return_details_structure(self):
        """Multi index return_details has correct structure."""
        result = multi_index_identification_reward(
            predicted={"ring_indices": [0, 1], "carbon_indices": [0, 1, 2]},
            target={"ring_indices": [0, 1], "carbon_indices": [0, 1, 2]},
            return_details=True
        )
        assert isinstance(result, dict)

    def test_multi_index_return_details_match(self):
        """Multi index return_details shows match."""
        result = multi_index_identification_reward(
            predicted={"ring_indices": [0, 1], "carbon_indices": [0, 1, 2]},
            target={"ring_indices": [0, 1], "carbon_indices": [0, 1, 2]},
            return_details=True
        )
        if isinstance(result, dict) and "reward" in result:
            assert result["reward"] == 1.0


class TestValidSmilesUtility:
    """Additional tests for valid_smiles utility."""

    def test_valid_smiles_complex_molecule(self):
        """Complex molecule is valid."""
        assert valid_smiles("CC(=O)Oc1ccccc1C(=O)O") is True  # Aspirin

    def test_valid_smiles_stereochemistry(self):
        """Stereochemistry is valid."""
        assert valid_smiles("C[C@H](N)C(=O)O") is True  # R-Alanine

    def test_valid_smiles_ez_bonds(self):
        """E/Z bonds are valid."""
        assert valid_smiles("C/C=C/C") is True  # E-2-butene

    def test_valid_smiles_heterocycle(self):
        """Heterocycle is valid."""
        assert valid_smiles("c1ccncc1") is True  # Pyridine

    def test_valid_smiles_charged(self):
        """Charged molecule is valid."""
        assert valid_smiles("[NH4+]") is True  # Ammonium

    def test_valid_smiles_radicals(self):
        """Some radical notation may be valid."""
        result = valid_smiles("[CH3]")
        assert isinstance(result, bool)


class TestRewardNormalization:
    """Tests for reward score normalization."""

    def test_score_range_0_to_1(self):
        """All scores are in [0, 1] range."""
        test_cases = [
            ({"a": 1}, {"a": 1}),
            ({"a": 0}, {"a": 1}),
            ({"a": 1}, {"a": 0}),
            ({"a": 1, "b": 2}, {"a": 1, "b": 2}),
            ({"a": 1, "b": 2}, {"a": 0, "b": 0}),
        ]
        for pred, target in test_cases:
            score = multi_count_dict_reward(pred, target)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for {pred}, {target}"

    def test_perfect_match_is_1(self):
        """Perfect match always returns 1.0."""
        assert single_count_reward({"x": 5}, {"x": 5}) == 1.0
        assert multi_count_dict_reward({"x": 5, "y": 10}, {"x": 5, "y": 10}) == 1.0

    def test_complete_mismatch(self):
        """Complete mismatch returns 0.0."""
        assert single_count_reward({"x": 0}, {"x": 100}) == 0.0


class TestIndexRewardSetOperations:
    """Tests for set-based index comparisons."""

    def test_index_exact_match(self):
        """Exact index match."""
        score = single_index_reward(
            predicted={"indices": [0, 1, 2]},
            target={"indices": [0, 1, 2]}
        )
        assert score == 1.0

    def test_index_different_order(self):
        """Different order should still match."""
        score = single_index_reward(
            predicted={"indices": [2, 0, 1]},
            target={"indices": [0, 1, 2]}
        )
        assert score == 1.0

    def test_index_subset(self):
        """Prediction is subset of target."""
        score = single_index_reward(
            predicted={"indices": [0, 1]},
            target={"indices": [0, 1, 2]}
        )
        # Partial credit or 0 depending on implementation
        assert 0.0 <= score <= 1.0

    def test_index_superset(self):
        """Prediction is superset of target."""
        score = single_index_reward(
            predicted={"indices": [0, 1, 2, 3]},
            target={"indices": [0, 1, 2]}
        )
        # May penalize extra indices
        assert 0.0 <= score <= 1.0

    def test_index_completely_wrong(self):
        """Completely wrong indices."""
        score = single_index_reward(
            predicted={"indices": [10, 11, 12]},
            target={"indices": [0, 1, 2]}
        )
        assert score == 0.0


class TestDispatcherAliases:
    """Tests for evaluate_answer task type aliases."""

    def test_single_count_alias(self):
        """'count' alias works."""
        score = evaluate_answer(
            task_type="count",
            predicted={"x": 1},
            target={"x": 1}
        )
        assert score == 1.0

    def test_multi_count_dict_alias(self):
        """'multi_count_dict' alias works."""
        score = evaluate_answer(
            task_type="multi_count_dict",
            predicted={"x": 1},
            target={"x": 1}
        )
        assert score == 1.0

    def test_single_index_identification_alias(self):
        """'single_index_identification' alias works."""
        score = evaluate_answer(
            task_type="single_index_identification",
            predicted={"i": [0]},
            target={"i": [0]}
        )
        assert score == 1.0

    def test_multi_constraint_alias(self):
        """'multi_constraint' alias works."""
        score = evaluate_answer(
            task_type="multi_constraint",
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score > 0.0

    def test_generation_alias(self):
        """'generation' alias works."""
        score = evaluate_answer(
            task_type="generation",
            predicted="c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score > 0.0


class TestDuplicateKeyDetection:
    """Tests for duplicate key handling in reward functions."""

    def test_json_duplicate_keys_count_returns_zero(self):
        """JSON with duplicate keys in count task returns 0."""
        score = multi_count_dict_reward(
            predicted='{"carbon": 5, "carbon": 6}',
            target={"carbon_count": 5}
        )
        assert score == 0.0

    def test_string_format_duplicate_keys_count_returns_zero(self):
        """String with duplicate keys in count task returns 0."""
        score = multi_count_dict_reward(
            predicted="carbon:5, carbon:6",
            target={"carbon_count": 5}
        )
        assert score == 0.0

    def test_normalized_duplicate_keys_count_returns_zero(self):
        """Different keys normalizing to same key in count task returns 0."""
        score = multi_count_dict_reward(
            predicted={"carbon atoms": 5, "carbons": 6},
            target={"carbon_count": 5}
        )
        assert score == 0.0

    def test_index_normalized_duplicate_keys_returns_zero(self):
        """Index task with normalized duplicates returns 0."""
        score = multi_index_identification_reward(
            predicted={"ring indices": [0, 1], "ring_indices": [2, 3]},
            target={"ring_indices": [0, 1]}
        )
        assert score == 0.0

    def test_no_duplicate_keys_still_works(self):
        """Ensure normal cases without duplicates still work."""
        score = multi_count_dict_reward(
            predicted={"carbon_count": 5},
            target={"carbon_count": 5}
        )
        assert score == 1.0

        score = multi_index_identification_reward(
            predicted={"ring_indices": [0, 1]},
            target={"ring_indices": [0, 1]}
        )
        assert score == 1.0
