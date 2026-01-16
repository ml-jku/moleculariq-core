"""
Tests for the questions module.

Validates SYSTEM_PROMPTS and TASKS dictionaries for structure and content.
"""

import re

from moleculariq_core.questions import SYSTEM_PROMPTS, TASKS


class TestSystemPrompts:
    """Tests for SYSTEM_PROMPTS dictionary."""

    def test_system_prompts_not_empty(self):
        """SYSTEM_PROMPTS should contain entries."""
        assert len(SYSTEM_PROMPTS) > 0

    def test_system_prompts_has_required_keys(self):
        """Should have expected prompt types."""
        expected_keys = ["with_key_hints", "concise"]
        for key in expected_keys:
            assert key in SYSTEM_PROMPTS, f"Missing system prompt: {key}"

    def test_system_prompts_values_are_strings(self):
        """All prompts should be non-empty strings."""
        for key, value in SYSTEM_PROMPTS.items():
            assert isinstance(value, str), f"Prompt {key} is not a string"
            assert len(value) > 0, f"Prompt {key} is empty"

    def test_system_prompts_contain_answer_tags(self):
        """Prompts should mention answer tags."""
        for key, value in SYSTEM_PROMPTS.items():
            assert "<answer>" in value, f"Prompt {key} missing <answer> tag reference"
            assert "</answer>" in value, f"Prompt {key} missing </answer> tag reference"

    def test_system_prompts_contain_json_reference(self):
        """Prompts should mention JSON format."""
        for key, value in SYSTEM_PROMPTS.items():
            assert "JSON" in value or "json" in value, f"Prompt {key} missing JSON reference"

    def test_system_prompts_contain_indexing_rules(self):
        """Prompts should contain indexing rules."""
        for key, value in SYSTEM_PROMPTS.items():
            # Check for 0-based indexing mention
            assert "0" in value, f"Prompt {key} should mention 0-based indexing"

    def test_with_key_hints_prompt_is_detailed(self):
        """with_key_hints prompt should be more detailed."""
        assert len(SYSTEM_PROMPTS["with_key_hints"]) > len(SYSTEM_PROMPTS["concise"])

    def test_concise_prompt_contains_examples(self):
        """Concise prompt should contain examples."""
        prompt = SYSTEM_PROMPTS["concise"]
        # Should have example answers
        assert "alcohol_count" in prompt or "ketone" in prompt or "ring_count" in prompt


class TestTasksStructure:
    """Tests for TASKS dictionary structure."""

    def test_tasks_not_empty(self):
        """TASKS should contain entries."""
        assert len(TASKS) > 0

    def test_tasks_has_expected_task_types(self):
        """Should have expected task types."""
        expected_tasks = [
            "constraint_generation",
            "single_count",
            "multi_count",
            "single_index_identification",
            "multi_index_identification",
        ]
        for task in expected_tasks:
            assert task in TASKS, f"Missing task type: {task}"

    def test_task_has_required_fields(self):
        """Each task should have required fields."""
        required_fields = [
            "task_name",
            "task_type",
            "description",
            "difficulty",
            "output_type",
            "input_fields",
            "output_fields",
            "question_templates",
        ]
        for task_name, task_config in TASKS.items():
            for field in required_fields:
                assert field in task_config, f"Task {task_name} missing field: {field}"


class TestTaskFields:
    """Tests for individual task field values."""

    def test_task_names_match_keys(self):
        """task_name field should match dictionary key."""
        for task_key, task_config in TASKS.items():
            assert task_config["task_name"] == task_key, f"task_name mismatch for {task_key}"

    def test_task_types_are_valid(self):
        """task_type should be one of expected values."""
        valid_types = ["count", "index", "generation"]
        for task_name, task_config in TASKS.items():
            assert task_config["task_type"] in valid_types, f"Invalid task_type for {task_name}"

    def test_difficulty_levels_are_valid(self):
        """difficulty should be one of expected values."""
        valid_difficulties = ["easy", "medium", "hard"]
        for task_name, task_config in TASKS.items():
            assert task_config["difficulty"] in valid_difficulties, f"Invalid difficulty for {task_name}"

    def test_output_types_are_valid(self):
        """output_type should be one of expected values."""
        valid_output_types = ["INTEGER", "DICT", "LIST", "SMILES"]
        for task_name, task_config in TASKS.items():
            assert task_config["output_type"] in valid_output_types, f"Invalid output_type for {task_name}"

    def test_input_fields_are_lists(self):
        """input_fields should be non-empty lists."""
        for task_name, task_config in TASKS.items():
            assert isinstance(task_config["input_fields"], list), f"input_fields not list for {task_name}"
            assert len(task_config["input_fields"]) > 0, f"Empty input_fields for {task_name}"

    def test_output_fields_are_lists(self):
        """output_fields should be non-empty lists."""
        for task_name, task_config in TASKS.items():
            assert isinstance(task_config["output_fields"], list), f"output_fields not list for {task_name}"
            assert len(task_config["output_fields"]) > 0, f"Empty output_fields for {task_name}"

    def test_descriptions_are_meaningful(self):
        """descriptions should be non-empty strings."""
        for task_name, task_config in TASKS.items():
            desc = task_config["description"]
            assert isinstance(desc, str), f"Description not string for {task_name}"
            assert len(desc) > 10, f"Description too short for {task_name}"


class TestQuestionTemplates:
    """Tests for question_templates in each task."""

    def test_question_templates_are_lists(self):
        """question_templates should be lists."""
        for task_name, task_config in TASKS.items():
            assert isinstance(task_config["question_templates"], list), f"Templates not list for {task_name}"

    def test_question_templates_not_empty(self):
        """question_templates should not be empty."""
        for task_name, task_config in TASKS.items():
            assert len(task_config["question_templates"]) > 0, f"Empty templates for {task_name}"

    def test_question_templates_have_minimum_count(self):
        """Each task should have at least 10 templates for variety."""
        for task_name, task_config in TASKS.items():
            count = len(task_config["question_templates"])
            assert count >= 10, f"Task {task_name} has only {count} templates, need at least 10"

    def test_question_templates_are_strings(self):
        """All templates should be non-empty strings."""
        for task_name, task_config in TASKS.items():
            for i, template in enumerate(task_config["question_templates"]):
                assert isinstance(template, str), f"Template {i} in {task_name} is not a string"
                assert len(template) > 0, f"Template {i} in {task_name} is empty"

    def test_question_templates_no_duplicates(self):
        """Templates should be unique within each task."""
        for task_name, task_config in TASKS.items():
            templates = task_config["question_templates"]
            unique_templates = set(templates)
            assert len(templates) == len(unique_templates), f"Duplicate templates in {task_name}"


class TestTemplatePlaceholders:
    """Tests for placeholder validation in templates."""

    def test_single_count_templates_have_required_placeholders(self):
        """single_count templates should have {smiles} and {count_type}."""
        templates = TASKS["single_count"]["question_templates"]
        for template in templates:
            assert "{smiles}" in template, f"Missing {{smiles}} in: {template[:50]}..."
            assert "{count_type}" in template, f"Missing {{count_type}} in: {template[:50]}..."

    def test_multi_count_templates_have_required_placeholders(self):
        """multi_count templates should have {smiles} and {count_types}."""
        templates = TASKS["multi_count"]["question_templates"]
        missing_placeholders = []
        for template in templates:
            assert "{smiles}" in template, f"Missing {{smiles}} in: {template[:50]}..."
            # Some templates may refer to "listed features" or similar without explicit {count_types}
            if "{count_types}" not in template and "listed" not in template.lower():
                missing_placeholders.append(template[:50])
        # Allow up to 2 templates with implicit feature references
        assert len(missing_placeholders) <= 2, f"Too many templates missing {{count_types}}: {missing_placeholders}"

    def test_single_index_templates_have_required_placeholders(self):
        """single_index_identification templates should have {smiles} and {index_type}."""
        templates = TASKS["single_index_identification"]["question_templates"]
        for template in templates:
            assert "{smiles}" in template, f"Missing {{smiles}} in: {template[:50]}..."
            assert "{index_type}" in template, f"Missing {{index_type}} in: {template[:50]}..."

    def test_multi_index_templates_have_required_placeholders(self):
        """multi_index_identification templates should have {smiles} and {index_types}."""
        templates = TASKS["multi_index_identification"]["question_templates"]
        for template in templates:
            assert "{smiles}" in template, f"Missing {{smiles}} in: {template[:50]}..."
            assert "{index_types}" in template, f"Missing {{index_types}} in: {template[:50]}..."

    def test_constraint_generation_templates_have_required_placeholders(self):
        """constraint_generation templates should have {constraint}."""
        templates = TASKS["constraint_generation"]["question_templates"]
        for template in templates:
            assert "{constraint}" in template, f"Missing {{constraint}} in: {template[:50]}..."

    def test_no_invalid_placeholders(self):
        """Templates should not have malformed placeholders."""
        placeholder_pattern = re.compile(r'\{[^}]*\}')
        valid_placeholders = {"{smiles}", "{count_type}", "{count_types}", "{index_type}", "{index_types}", "{constraint}"}

        for task_name, task_config in TASKS.items():
            for template in task_config["question_templates"]:
                found = placeholder_pattern.findall(template)
                for placeholder in found:
                    assert placeholder in valid_placeholders, f"Unknown placeholder {placeholder} in {task_name}"


class TestTaskTypeConsistency:
    """Tests for consistency between task_type and other fields."""

    def test_count_tasks_have_count_output_type(self):
        """Tasks with task_type 'count' should have appropriate output_type."""
        for task_name, task_config in TASKS.items():
            if task_config["task_type"] == "count":
                assert task_config["output_type"] in ["INTEGER", "DICT"], f"Count task {task_name} has wrong output_type"

    def test_index_tasks_have_appropriate_output_type(self):
        """Tasks with task_type 'index' should have LIST or DICT output_type."""
        for task_name, task_config in TASKS.items():
            if task_config["task_type"] == "index":
                assert task_config["output_type"] in ["LIST", "DICT"], f"Index task {task_name} has wrong output_type"

    def test_generation_tasks_have_smiles_output_type(self):
        """Tasks with task_type 'generation' should have SMILES output_type."""
        for task_name, task_config in TASKS.items():
            if task_config["task_type"] == "generation":
                assert task_config["output_type"] == "SMILES", f"Generation task {task_name} has wrong output_type"

    def test_single_tasks_have_integer_or_list_output(self):
        """Single tasks (not multi) should have INTEGER or LIST output."""
        for task_name, task_config in TASKS.items():
            if "single" in task_name:
                assert task_config["output_type"] in ["INTEGER", "LIST"], f"Single task {task_name} has DICT output"

    def test_multi_tasks_have_dict_output(self):
        """Multi tasks should have DICT output."""
        for task_name, task_config in TASKS.items():
            if "multi" in task_name:
                assert task_config["output_type"] == "DICT", f"Multi task {task_name} doesn't have DICT output"


class TestInputOutputFieldsConsistency:
    """Tests for input/output field consistency."""

    def test_smiles_in_input_fields_for_count_index_tasks(self):
        """Count and index tasks should have 'smiles' in input_fields."""
        for task_name, task_config in TASKS.items():
            if task_config["task_type"] in ["count", "index"]:
                assert "smiles" in task_config["input_fields"], f"Task {task_name} missing 'smiles' input"

    def test_constraint_in_input_fields_for_generation_tasks(self):
        """Generation tasks should have 'constraint' in input_fields."""
        for task_name, task_config in TASKS.items():
            if task_config["task_type"] == "generation":
                assert "constraint" in task_config["input_fields"], f"Task {task_name} missing 'constraint' input"

    def test_output_fields_match_task_type(self):
        """Output fields should be appropriate for task type."""
        expected_outputs = {
            "single_count": ["count"],
            "multi_count": ["counts"],
            "single_index_identification": ["indices"],
            "multi_index_identification": ["index_mappings"],
            "constraint_generation": ["smiles"],
        }
        for task_name, expected in expected_outputs.items():
            assert TASKS[task_name]["output_fields"] == expected, f"Wrong output_fields for {task_name}"


class TestTemplateFormatting:
    """Tests for template string formatting."""

    def test_templates_can_be_formatted(self):
        """Templates should be formattable with appropriate kwargs."""
        test_values = {
            "single_count": {"smiles": "CCO", "count_type": "carbon atoms"},
            "multi_count": {"smiles": "CCO", "count_types": "carbon atoms, oxygen atoms"},
            "single_index_identification": {"smiles": "CCO", "index_type": "carbon atoms"},
            "multi_index_identification": {"smiles": "CCO", "index_types": "carbon atoms, oxygen atoms"},
            "constraint_generation": {"constraint": "2 carbon atoms"},
        }

        for task_name, kwargs in test_values.items():
            templates = TASKS[task_name]["question_templates"]
            for template in templates:
                formatted = template.format(**kwargs)
                assert isinstance(formatted, str)
                assert len(formatted) > 0
                # Placeholders should be replaced
                assert "{" not in formatted, f"Unformatted placeholder in: {formatted}"

    def test_templates_end_with_punctuation(self):
        """Most templates should end with proper punctuation."""
        punctuation = {'.', '?', '!'}
        for task_name, task_config in TASKS.items():
            for template in task_config["question_templates"]:
                last_char = template.strip()[-1]
                assert last_char in punctuation, f"Template doesn't end with punctuation: {template[-50:]}"


class TestTemplateVariety:
    """Tests for template variety and quality."""

    def test_templates_have_varied_beginnings(self):
        """Templates should not all start the same way."""
        for task_name, task_config in TASKS.items():
            templates = task_config["question_templates"]
            first_words = [t.split()[0] for t in templates]
            unique_first_words = set(first_words)
            # Should have at least 5 different starting words
            assert len(unique_first_words) >= 5, f"Task {task_name} has too little variety in template starts"

    def test_templates_have_reasonable_length(self):
        """Templates should have reasonable length (not too short or long)."""
        for task_name, task_config in TASKS.items():
            for template in task_config["question_templates"]:
                length = len(template)
                assert length >= 20, f"Template too short in {task_name}: {template}"
                assert length <= 200, f"Template too long in {task_name}: {template[:50]}..."
