"""
Integration tests for end-to-end workflows.

Tests complete pipelines: compute → format → evaluate.
"""

from moleculariq_core import (
    TASKS,
    SYSTEM_PROMPTS,
    evaluate_answer,
)


class TestCountWorkflow:
    """End-to-end tests for count task workflow."""

    def test_single_count_workflow(self, solver, formatter, benzene):
        """Complete single count workflow."""
        # 1. Compute ground truth
        ring_count = solver.get_ring_count(benzene)

        # 2. Format question (use default template - TASKS templates have different placeholders)
        question = formatter.format_count_query(
            smiles=benzene,
            count_types=["ring_count"]
        )

        # 3. Simulate perfect prediction
        predicted = {"ring_count": ring_count}
        target = {"ring_count": ring_count}

        # 4. Evaluate
        score = evaluate_answer(
            task_type="single_count",
            predicted=predicted,
            target=target
        )

        assert score == 1.0
        assert benzene in question

    def test_multi_count_workflow(self, solver, formatter, aspirin):
        """Complete multi-count workflow."""
        # 1. Compute ground truths
        ground_truth = {
            "ring_count": solver.get_ring_count(aspirin),
            "carbon_atom_count": solver.get_carbon_atom_count(aspirin),
            "hba_count": solver.get_hba_count(aspirin),
        }

        # 2. Format question
        template = TASKS["multi_count"]["question_templates"][0]
        question = formatter.format_count_query(
            smiles=aspirin,
            count_types=list(ground_truth.keys()),
            template=template
        )

        # 3. Evaluate perfect prediction
        score = evaluate_answer(
            task_type="multi_count",
            predicted=ground_truth,
            target=ground_truth
        )

        assert score == 1.0
        assert aspirin in question


class TestIndexWorkflow:
    """End-to-end tests for index identification workflow."""

    def test_single_index_workflow(self, solver, formatter, benzene):
        """Complete single index workflow."""
        # 1. Compute ground truth indices
        ring_indices = solver.get_ring_indices(benzene)

        # 2. Format question (use default template - TASKS templates have different placeholders)
        question = formatter.format_index_query(
            smiles=benzene,
            index_types=["ring_index"]
        )

        # 3. Evaluate perfect prediction
        predicted = {"ring_index": ring_indices}
        target = {"ring_index": ring_indices}

        score = evaluate_answer(
            task_type="single_index",
            predicted=predicted,
            target=target
        )

        assert score == 1.0
        assert benzene in question

    def test_multi_index_workflow(self, solver, formatter, ethanol):
        """Complete multi-index workflow."""
        # 1. Compute ground truth indices
        ground_truth = {
            "carbon_atom_index": solver.get_carbon_atom_indices(ethanol),
            "hetero_atom_index": solver.get_hetero_atom_indices(ethanol),
        }

        # 2. Format question
        question = formatter.format_index_query(
            smiles=ethanol,
            index_types=list(ground_truth.keys())
        )

        # 3. Evaluate perfect prediction
        score = evaluate_answer(
            task_type="multi_index",
            predicted=ground_truth,
            target=ground_truth
        )

        assert score == 1.0


class TestConstraintWorkflow:
    """End-to-end tests for constraint generation workflow."""

    def test_constraint_workflow_ring(self, solver, formatter):
        """Constraint workflow: generate molecule with 1 ring."""
        # 1. Define constraint
        constraint = {"type": "ring_count", "operator": "=", "value": 1}

        # 2. Format question
        constraint_text = formatter.format_constraint(constraint)
        template = TASKS["constraint_generation"]["question_templates"][0]
        question = template.format(constraint=constraint_text)

        # 3. Use benzene as "generated" molecule
        generated_smiles = "c1ccccc1"  # Benzene has 1 ring

        # 4. Verify it satisfies constraint
        actual_rings = solver.get_ring_count(generated_smiles)
        assert actual_rings == 1

        # 5. Evaluate
        score = evaluate_answer(
            task_type="constraint_generation",
            predicted=generated_smiles,
            constraints=[constraint]
        )

        assert score > 0.0

    def test_constraint_workflow_multiple(self, solver, formatter):
        """Constraint workflow: multiple constraints."""
        # 1. Define constraints
        constraints = [
            {"type": "ring_count", "operator": "=", "value": 1},
            {"type": "aromatic_ring_count", "operator": "=", "value": 1},
        ]

        # 2. Format constraints
        text = formatter.format_constraints_list(constraints)
        assert len(text) > 0

        # 3. Use benzene (satisfies both: 1 ring, 1 aromatic)
        generated_smiles = "c1ccccc1"

        # 4. Evaluate
        score = evaluate_answer(
            task_type="constraint_generation",
            predicted=generated_smiles,
            constraints=constraints
        )

        assert score > 0.0


class TestSystemPromptIntegration:
    """Tests for system prompt usage."""

    def test_system_prompt_with_question(self, formatter, benzene):
        """System prompt + question forms complete input."""
        system_prompt = SYSTEM_PROMPTS["with_key_hints"]
        question = formatter.format_count_query(
            smiles=benzene,
            count_types=["ring_count"],
            include_key_hint=True,
            key_names=["ring_count"]
        )

        # Combined should form valid LLM input
        full_input = f"{system_prompt}\n\n{question}"
        assert len(full_input) > 100
        assert "ring_count" in full_input
        assert benzene in full_input


class TestRandomTemplateSelection:
    """Tests for random template selection."""

    def test_different_templates_produce_different_questions(self, formatter, ethanol):
        """Different templates produce varied questions."""
        # Use formatter-compatible templates (with {count_types} placeholder)
        templates = [
            "How many {count_types} are in {smiles}?",
            "Count the {count_types} in molecule {smiles}.",
            "What is the count of {count_types} in {smiles}?",
            "For {smiles}, determine the {count_types}.",
            "Calculate {count_types} for {smiles}.",
        ]
        questions = set()

        for template in templates:
            q = formatter.format_count_query(
                smiles=ethanol,
                count_types=["ring_count"],
                template=template
            )
            questions.add(q)

        # Should have multiple unique questions
        assert len(questions) > 1


class TestDataConsistency:
    """Tests for data consistency across modules."""

    def test_solver_and_reward_agree(self, solver):
        """Solver output matches expected reward evaluation."""
        smiles = "CCO"

        # Compute with solver
        ring_count = solver.get_ring_count(smiles)

        # Build prediction matching ground truth
        predicted = {"ring_count": ring_count}
        target = {"ring_count": ring_count}

        # Reward should be perfect
        score = evaluate_answer(
            task_type="single_count",
            predicted=predicted,
            target=target
        )

        assert score == 1.0

    def test_multiple_properties_consistent(self, solver, aspirin):
        """Multiple properties computed consistently."""
        # Compute various properties
        props = {
            "ring_count": solver.get_ring_count(aspirin),
            "aromatic_ring_count": solver.get_aromatic_ring_count(aspirin),
            "carbon_atom_count": solver.get_carbon_atom_count(aspirin),
            "heavy_atom_count": solver.get_heavy_atom_count(aspirin),
        }

        # Sanity checks
        assert props["ring_count"] >= props["aromatic_ring_count"]
        assert props["heavy_atom_count"] >= props["carbon_atom_count"]
