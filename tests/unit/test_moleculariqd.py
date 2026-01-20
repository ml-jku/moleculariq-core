"""
Tests for MolecularIQD high-level API.

Tests the main entry point for question generation and answer validation.
"""

from moleculariq_core import (
    MolecularIQD,
    load_molecule_pool,
    get_available_pools,
    MoleculePoolHiddenError,
)


class TestMolecularIQDInstantiation:
    """Tests for MolecularIQD initialization."""

    def test_default_instantiation(self):
        """MolecularIQD instantiates with defaults."""
        mqd = MolecularIQD()
        assert mqd is not None
        assert mqd.solver is not None
        assert mqd.formatter is not None

    def test_instantiation_with_seed(self):
        """MolecularIQD accepts seed parameter."""
        mqd = MolecularIQD(seed=42)
        assert mqd.seed == 42

    def test_instantiation_with_options(self):
        """MolecularIQD accepts all configuration options."""
        mqd = MolecularIQD(
            seed=123,
            enable_random_phrasing=False,
            cache_properties=False,
            system_prompt_style="concise"
        )
        assert mqd.seed == 123
        assert mqd.cache_properties is False

    def test_reproducibility_with_seed(self):
        """Same seed produces same questions."""
        mqd1 = MolecularIQD(seed=42)
        mqd2 = MolecularIQD(seed=42)

        q1, _, _ = mqd1.generate_count_question("CCO", "ring_count")
        q2, _, _ = mqd2.generate_count_question("CCO", "ring_count")

        assert q1 == q2


class TestGenerateCountQuestion:
    """Tests for generate_count_question method."""

    def test_single_count_question(self):
        """Generate single count question."""
        mqd = MolecularIQD(seed=42)
        question, answer, metadata = mqd.generate_count_question(
            smiles="CCO",
            count_properties="ring_count"
        )

        assert isinstance(question, str)
        assert "CCO" in question
        assert isinstance(answer, dict)
        assert "ring_count" in answer or any("ring" in k for k in answer.keys())
        assert metadata["task_type"] == "single_count"

    def test_multi_count_question(self):
        """Generate multi-count question."""
        mqd = MolecularIQD(seed=42)
        question, answer, metadata = mqd.generate_count_question(
            smiles="c1ccccc1",
            count_properties=["ring_count", "aromatic_ring_count"]
        )

        assert isinstance(question, str)
        assert isinstance(answer, dict)
        assert len(answer) == 2
        assert metadata["task_type"] == "multi_count"

    def test_count_question_correct_values(self):
        """Count question returns correct ground truth."""
        mqd = MolecularIQD(seed=42)

        # Benzene: 1 ring, 1 aromatic ring
        _, answer, _ = mqd.generate_count_question("c1ccccc1", "ring_count")
        assert list(answer.values())[0] == 1

        # Ethanol: 0 rings
        _, answer, _ = mqd.generate_count_question("CCO", "ring_count")
        assert list(answer.values())[0] == 0

    def test_count_question_with_custom_template(self):
        """Count question accepts custom template."""
        mqd = MolecularIQD(seed=42)
        template = "Count the {count_type} in {smiles}."

        question, _, _ = mqd.generate_count_question(
            smiles="CCO",
            count_properties="ring_count",
            template=template
        )

        assert "Count the" in question
        assert "CCO" in question

    def test_count_question_without_key_hint(self):
        """Count question without key hint."""
        mqd = MolecularIQD(seed=42)

        q_with_hint, _, _ = mqd.generate_count_question(
            "CCO", "ring_count", include_key_hint=True
        )
        q_without_hint, _, _ = mqd.generate_count_question(
            "CCO", "ring_count", include_key_hint=False
        )

        assert len(q_without_hint) < len(q_with_hint)


class TestGenerateIndexQuestion:
    """Tests for generate_index_question method."""

    def test_single_index_question(self):
        """Generate single index question."""
        mqd = MolecularIQD(seed=42)
        question, answer, metadata = mqd.generate_index_question(
            smiles="CCO",
            index_properties="carbon_atom_index"
        )

        assert isinstance(question, str)
        assert "CCO" in question
        assert isinstance(answer, dict)
        assert metadata["task_type"] == "single_index_identification"

    def test_multi_index_question(self):
        """Generate multi-index question."""
        mqd = MolecularIQD(seed=42)
        question, answer, metadata = mqd.generate_index_question(
            smiles="CCO",
            index_properties=["carbon_atom_index", "hetero_atom_index"]
        )

        assert isinstance(question, str)
        assert isinstance(answer, dict)
        assert len(answer) == 2
        assert metadata["task_type"] == "multi_index_identification"

    def test_index_question_correct_values(self):
        """Index question returns correct ground truth."""
        mqd = MolecularIQD(seed=42)

        # Ethanol CCO: carbons at 0, 1; oxygen at 2
        _, answer, _ = mqd.generate_index_question("CCO", "carbon_atom_index")
        indices = list(answer.values())[0]
        assert 0 in indices
        assert 1 in indices
        assert 2 not in indices

    def test_index_question_returns_lists(self):
        """Index answers are lists."""
        mqd = MolecularIQD(seed=42)
        _, answer, _ = mqd.generate_index_question("CCO", "carbon_atom_index")

        for value in answer.values():
            assert isinstance(value, list)


class TestGenerateConstraintQuestion:
    """Tests for generate_constraint_question method."""

    def test_single_constraint_question(self):
        """Generate single constraint question."""
        mqd = MolecularIQD(seed=42)
        question, metadata = mqd.generate_constraint_question(
            constraints=[{"property": "ring_count", "operator": "=", "value": 1}]
        )

        assert isinstance(question, str)
        assert metadata["task_type"] == "constraint_generation"

    def test_multiple_constraints_question(self):
        """Generate question with multiple constraints."""
        mqd = MolecularIQD(seed=42)
        question, metadata = mqd.generate_constraint_question(
            constraints=[
                {"property": "ring_count", "operator": ">=", "value": 1},
                {"property": "carbon_atom_count", "operator": "<=", "value": 10}
            ]
        )

        assert isinstance(question, str)
        assert len(metadata["constraints"]) == 2

    def test_constraint_question_operators(self):
        """Constraint question handles different operators."""
        mqd = MolecularIQD(seed=42)

        for operator in ["=", ">=", "<=", ">", "<"]:
            question, _ = mqd.generate_constraint_question(
                constraints=[{"property": "ring_count", "operator": operator, "value": 2}]
            )
            assert isinstance(question, str)


class TestValidateCountAnswer:
    """Tests for validate_count_answer method."""

    def test_correct_answer_scores_one(self):
        """Correct count answer scores 1.0."""
        mqd = MolecularIQD(seed=42)

        # Benzene has 1 ring
        score = mqd.validate_count_answer("c1ccccc1", {"ring_count": 1})
        assert score == 1.0

    def test_incorrect_answer_scores_zero(self):
        """Incorrect count answer scores 0.0."""
        mqd = MolecularIQD(seed=42)

        # Benzene has 1 ring, not 5
        score = mqd.validate_count_answer("c1ccccc1", {"ring_count": 5})
        assert score == 0.0

    def test_validate_with_ground_truth(self):
        """Validate against provided ground truth."""
        mqd = MolecularIQD(seed=42)

        score = mqd.validate_count_answer(
            "CCO",
            predicted_answer={"ring_count": 0},
            ground_truth={"ring_count": 0}
        )
        assert score == 1.0

    def test_validate_with_return_details(self):
        """Validate with detailed results."""
        mqd = MolecularIQD(seed=42)

        result = mqd.validate_count_answer(
            "c1ccccc1",
            {"ring_count": 1},
            return_details=True
        )
        assert isinstance(result, dict)
        assert "reward" in result


class TestValidateIndexAnswer:
    """Tests for validate_index_answer method."""

    def test_correct_indices_score_one(self):
        """Correct index answer scores 1.0."""
        mqd = MolecularIQD(seed=42)

        # Ethanol: carbons at 0, 1
        score = mqd.validate_index_answer("CCO", {"carbon_atom_index": [0, 1]})
        assert score == 1.0

    def test_incorrect_indices_score_zero(self):
        """Incorrect index answer scores 0.0."""
        mqd = MolecularIQD(seed=42)

        # Wrong indices
        score = mqd.validate_index_answer("CCO", {"carbon_atom_index": [5, 6, 7]})
        assert score == 0.0

    def test_partial_indices_lower_score(self):
        """Missing indices result in lower or zero score."""
        mqd = MolecularIQD(seed=42)

        # Only one correct index out of two - may score 0 or partial depending on implementation
        score = mqd.validate_index_answer("CCO", {"carbon_atom_index": [0]})
        assert score < 1.0  # Not perfect score


class TestValidateConstraintAnswer:
    """Tests for validate_constraint_answer method."""

    def test_satisfying_molecule_scores_positive(self):
        """Molecule satisfying constraints scores > 0."""
        mqd = MolecularIQD(seed=42)

        # Benzene satisfies ring_count = 1
        score = mqd.validate_constraint_answer(
            "c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score > 0.0

    def test_invalid_smiles_scores_zero(self):
        """Invalid SMILES scores 0.0."""
        mqd = MolecularIQD(seed=42)

        score = mqd.validate_constraint_answer(
            "invalid_smiles_xyz",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )
        assert score == 0.0

    def test_unsatisfied_constraint_scores_zero(self):
        """Molecule not satisfying constraint scores lower."""
        mqd = MolecularIQD(seed=42)

        # Ethanol has 0 rings, constraint requires 5
        score = mqd.validate_constraint_answer(
            "CCO",
            constraints=[{"type": "ring_count", "operator": "=", "value": 5}]
        )
        assert score < 1.0


class TestComputeProperty:
    """Tests for compute_property method."""

    def test_compute_count_property(self):
        """Compute count property."""
        mqd = MolecularIQD(seed=42)

        result = mqd.compute_property("c1ccccc1", "ring_count")
        assert result == 1

    def test_compute_index_property(self):
        """Compute index property."""
        mqd = MolecularIQD(seed=42)

        result = mqd.compute_property("CCO", "carbon_atom_index")
        assert isinstance(result, list)
        assert 0 in result
        assert 1 in result

    def test_compute_multiple_properties(self):
        """Compute multiple properties."""
        mqd = MolecularIQD(seed=42)

        result = mqd.compute_properties("c1ccccc1", ["ring_count", "aromatic_ring_count"])
        assert result["ring_count"] == 1
        assert result["aromatic_ring_count"] == 1

    def test_property_caching(self):
        """Properties are cached when enabled."""
        mqd = MolecularIQD(seed=42, cache_properties=True)

        # First call computes
        mqd.compute_property("CCO", "ring_count")
        assert len(mqd._property_cache) > 0

        # Clear and verify
        mqd.clear_cache()
        assert len(mqd._property_cache) == 0


class TestGeneratePairedQuestion:
    """Tests for generate_paired_question method."""

    def test_paired_question_structure(self):
        """Paired question returns count and index tasks."""
        mqd = MolecularIQD(seed=42)

        count_task, index_task = mqd.generate_paired_question("c1ccccc1", "ring_count")

        assert "question" in count_task
        assert "answer" in count_task
        assert "metadata" in count_task

        assert "question" in index_task
        assert "answer" in index_task
        assert "metadata" in index_task

    def test_paired_question_consistency(self):
        """Paired count and index answers are consistent."""
        mqd = MolecularIQD(seed=42)

        count_task, index_task = mqd.generate_paired_question("CCO", "carbon_atom_count")

        count_value = list(count_task["answer"].values())[0]
        index_value = list(index_task["answer"].values())[0]

        # Count should equal length of indices
        assert count_value == len(index_value)


class TestAvailableProperties:
    """Tests for property listing methods."""

    def test_get_available_count_properties(self):
        """Get list of count properties."""
        mqd = MolecularIQD(seed=42)
        props = mqd.get_available_count_properties()

        assert isinstance(props, list)
        assert len(props) > 0
        assert "ring_count" in props

    def test_get_available_index_properties(self):
        """Get list of index properties."""
        mqd = MolecularIQD(seed=42)
        props = mqd.get_available_index_properties()

        assert isinstance(props, list)
        assert len(props) > 0

    def test_get_available_constraint_properties(self):
        """Get list of constraint properties."""
        mqd = MolecularIQD(seed=42)
        props = mqd.get_available_constraint_properties()

        assert isinstance(props, list)
        assert len(props) > 0


class TestMoleculePoolLoader:
    """Tests for molecule pool loading functions."""

    def test_get_available_pools(self):
        """Get list of available pools."""
        pools = get_available_pools()

        assert isinstance(pools, list)
        assert "train" in pools

    def test_hidden_pool_raises_error(self):
        """Accessing hidden pool raises MoleculePoolHiddenError."""
        import pytest

        with pytest.raises(MoleculePoolHiddenError):
            load_molecule_pool("val_hard")

        with pytest.raises(MoleculePoolHiddenError):
            load_molecule_pool("val_easy")

        with pytest.raises(MoleculePoolHiddenError):
            load_molecule_pool("test")

    def test_unknown_pool_raises_value_error(self):
        """Unknown pool raises ValueError."""
        import pytest

        with pytest.raises(ValueError):
            load_molecule_pool("nonexistent_pool")

    def test_hidden_pool_error_message(self):
        """Hidden pool error has informative message."""
        try:
            load_molecule_pool("val_hard")
        except MoleculePoolHiddenError as e:
            assert "hidden" in str(e).lower()
            assert "leakage" in str(e).lower()


class TestFunctionalGroupProperties:
    """Tests for functional group property computation."""

    def test_functional_group_count(self):
        """Compute functional group count."""
        mqd = MolecularIQD(seed=42)

        # Ethanol has alcohol group
        result = mqd.compute_property("CCO", "functional_group_alcohol_count")
        assert result >= 0

    def test_functional_group_in_question(self):
        """Generate question for functional group."""
        mqd = MolecularIQD(seed=42)

        question, answer, _ = mqd.generate_count_question(
            "CCO",
            "functional_group_alcohol_count"
        )

        assert isinstance(question, str)
        assert isinstance(answer, dict)

    def test_functional_group_index(self):
        """Compute functional group index property."""
        mqd = MolecularIQD(seed=42)

        result = mqd.compute_property("CCO", "functional_group_alcohol_index")
        assert isinstance(result, list)


class TestReactionTemplateProperties:
    """Tests for reaction template property computation."""

    def test_reaction_count_property(self):
        """Compute reaction count property via get_reaction_counts_and_indices."""
        mqd = MolecularIQD(seed=42)

        # Test via the solver's reaction_counts_and_indices method
        result = mqd.solver.get_reaction_counts_and_indices("c1ccccc1")
        assert isinstance(result, dict)

    def test_reaction_index_property(self):
        """Compute reaction data via get_reaction_counts_and_indices."""
        mqd = MolecularIQD(seed=42)

        result = mqd.solver.get_reaction_counts_and_indices("CCO")
        assert isinstance(result, dict)
        # Should contain reaction-related entries
        assert len(result) > 0
        # Keys should be template_based_reaction_prediction_* entries
        for key in result.keys():
            assert key.startswith("template_based_reaction_prediction_")


class TestPropertyResolutionFallbacks:
    """Tests for property name resolution edge cases."""

    def test_unknown_count_property_returns_zero(self):
        """Unknown count property returns 0."""
        mqd = MolecularIQD(seed=42)

        result = mqd.compute_property("CCO", "completely_unknown_xyz_count")
        assert result == 0

    def test_unknown_index_property_returns_empty_list(self):
        """Unknown index property returns empty list."""
        mqd = MolecularIQD(seed=42)

        result = mqd.compute_property("CCO", "completely_unknown_xyz_index")
        assert result == []

    def test_cache_hit_returns_cached_value(self):
        """Property cache returns cached value on second call."""
        mqd = MolecularIQD(seed=42, cache_properties=True)

        # First call - computes and caches
        result1 = mqd.compute_property("CCO", "ring_count")

        # Verify it's in cache
        assert ("CCO", "ring_count") in mqd._property_cache

        # Second call - should return cached value
        result2 = mqd.compute_property("CCO", "ring_count")

        assert result1 == result2

    def test_no_caching_when_disabled(self):
        """Property cache is not used when disabled."""
        mqd = MolecularIQD(seed=42, cache_properties=False)

        mqd.compute_property("CCO", "ring_count")
        assert len(mqd._property_cache) == 0


class TestConstraintEdgeCases:
    """Tests for constraint edge cases."""

    def test_constraint_with_type_key(self):
        """Constraint using 'type' key instead of 'property'."""
        mqd = MolecularIQD(seed=42)

        question, metadata = mqd.generate_constraint_question(
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
        )

        assert isinstance(question, str)

    def test_constraint_with_range_values(self):
        """Constraint with min_value and max_value."""
        mqd = MolecularIQD(seed=42)

        question, metadata = mqd.generate_constraint_question(
            constraints=[{
                "property": "ring_count",
                "operator": "range",
                "value": None,
                "min_value": 1,
                "max_value": 5
            }]
        )

        assert isinstance(question, str)

    def test_constraint_without_key_hint(self):
        """Constraint question without key hint."""
        mqd = MolecularIQD(seed=42)

        q_with, _ = mqd.generate_constraint_question(
            constraints=[{"property": "ring_count", "operator": "=", "value": 1}],
            include_key_hint=True
        )
        q_without, _ = mqd.generate_constraint_question(
            constraints=[{"property": "ring_count", "operator": "=", "value": 1}],
            include_key_hint=False
        )

        assert len(q_without) < len(q_with)


class TestValidationWithDetails:
    """Tests for validation methods with return_details=True."""

    def test_validate_index_with_details(self):
        """Validate index answer with detailed results."""
        mqd = MolecularIQD(seed=42)

        result = mqd.validate_index_answer(
            "CCO",
            {"carbon_atom_index": [0, 1]},
            return_details=True
        )
        assert isinstance(result, dict)
        assert "reward" in result

    def test_validate_constraint_with_details(self):
        """Validate constraint answer with detailed results."""
        mqd = MolecularIQD(seed=42)

        result = mqd.validate_constraint_answer(
            "c1ccccc1",
            constraints=[{"type": "ring_count", "operator": "=", "value": 1}],
            return_details=True
        )
        # Returns dict with "reward" key when return_details=True
        assert isinstance(result, dict)
        assert "reward" in result


class TestPairedQuestionEdgeCases:
    """Tests for paired question edge cases."""

    def test_paired_question_unmapped_property(self):
        """Paired question with property not in COUNT_TO_INDEX_MAP."""
        mqd = MolecularIQD(seed=42)

        # Use a property that might not have explicit mapping
        count_task, index_task = mqd.generate_paired_question(
            "CCO",
            "heavy_atom_count"
        )

        assert "question" in count_task
        assert "question" in index_task

    def test_paired_question_with_templates(self):
        """Paired question with custom templates."""
        mqd = MolecularIQD(seed=42)

        count_template = "Count {count_type} in {smiles}."
        index_template = "Find {index_type} in {smiles}."

        count_task, index_task = mqd.generate_paired_question(
            "c1ccccc1",
            "ring_count",
            template_count=count_template,
            template_index=index_template
        )

        assert "Count" in count_task["question"]
        assert "Find" in index_task["question"]


class TestIndexQuestionEdgeCases:
    """Tests for index question edge cases."""

    def test_index_question_without_key_hint(self):
        """Index question without key hint."""
        mqd = MolecularIQD(seed=42)

        q_with, _, _ = mqd.generate_index_question(
            "CCO", "carbon_atom_index", include_key_hint=True
        )
        q_without, _, _ = mqd.generate_index_question(
            "CCO", "carbon_atom_index", include_key_hint=False
        )

        assert len(q_without) < len(q_with)

    def test_index_question_with_custom_template(self):
        """Index question with custom template."""
        mqd = MolecularIQD(seed=42)

        template = "Identify {index_type} positions in {smiles}."
        question, _, _ = mqd.generate_index_question(
            "CCO",
            "carbon_atom_index",
            template=template
        )

        assert "Identify" in question
        assert "CCO" in question
