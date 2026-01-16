"""
Unit tests for NaturalLanguageFormatter.

Tests question formatting and natural language conversion.
"""

from moleculariq_core import (
    TASKS,
    SYSTEM_PROMPTS,
    get_natural_language,
    parse_natural_language,
)


class TestNaturalLanguageFormatter:
    """Tests for NaturalLanguageFormatter class."""

    def test_formatter_initialization(self, formatter):
        """Formatter initializes correctly."""
        assert formatter is not None
        assert hasattr(formatter, 'format_count_query')
        assert hasattr(formatter, 'format_index_query')

    def test_format_count_query_single(self, formatter):
        """Format single count query."""
        query = formatter.format_count_query(
            smiles="CCO",
            count_types=["ring_count"]
        )
        assert "CCO" in query
        assert isinstance(query, str)

    def test_format_count_query_multiple(self, formatter):
        """Format multiple count query."""
        query = formatter.format_count_query(
            smiles="CCO",
            count_types=["ring_count", "carbon_atom_count"]
        )
        assert "CCO" in query
        assert isinstance(query, str)

    def test_format_index_query_single(self, formatter):
        """Format single index query."""
        query = formatter.format_index_query(
            smiles="c1ccccc1",
            index_types=["ring_index"]
        )
        assert "c1ccccc1" in query
        assert isinstance(query, str)

    def test_format_constraint(self, formatter):
        """Format constraint into natural language."""
        constraint = {"type": "ring_count", "operator": "=", "value": 2}
        text = formatter.format_constraint(constraint)
        assert isinstance(text, str)
        assert "2" in text or "two" in text.lower()

    def test_format_constraint_range(self, formatter):
        """Format range constraint."""
        constraint = {
            "type": "carbon_atom_count",
            "operator": "range",
            "min_value": 5,
            "max_value": 10
        }
        text = formatter.format_constraint(constraint)
        assert "5" in text
        assert "10" in text

    def test_format_constraints_list(self, formatter):
        """Format list of constraints."""
        constraints = [
            {"type": "ring_count", "operator": "=", "value": 1},
            {"type": "carbon_atom_count", "operator": ">=", "value": 6}
        ]
        text = formatter.format_constraints_list(constraints)
        assert isinstance(text, str)
        assert len(text) > 0


class TestTechnicalNaturalConversion:
    """Tests for technical ↔ natural language conversion."""

    def test_technical_to_natural_ring(self, formatter):
        """Convert 'ring_count' to natural language."""
        natural = formatter.technical_to_natural("ring_count", "count")
        assert isinstance(natural, str)
        assert len(natural) > 0

    def test_technical_to_natural_hba(self, formatter):
        """Convert 'hba_count' to natural language."""
        natural = formatter.technical_to_natural("hba_count", "count")
        assert isinstance(natural, str)
        # Should expand HBA to something readable

    def test_natural_to_technical_round_trip(self, formatter):
        """Round trip conversion preserves meaning."""
        original = "ring_count"
        natural = formatter.technical_to_natural(original, "count")
        # The reverse should give something related
        technical = formatter.natural_to_technical(natural, "count")
        assert isinstance(technical, str)


class TestGetNaturalLanguage:
    """Tests for get_natural_language function."""

    def test_get_natural_language_known_key(self):
        """Known key returns natural language forms."""
        forms = get_natural_language("ring_count", "count")
        assert isinstance(forms, list)
        assert len(forms) > 0

    def test_get_natural_language_unknown_key(self):
        """Unknown key returns fallback."""
        forms = get_natural_language("unknown_xyz_property", "count")
        assert isinstance(forms, list)


class TestParseNaturalLanguage:
    """Tests for parse_natural_language function."""

    def test_parse_natural_language_basic(self):
        """Parse basic natural language to technical."""
        result = parse_natural_language("rings")
        assert isinstance(result, str)

    def test_parse_natural_language_passthrough(self):
        """Unknown text passes through."""
        result = parse_natural_language("completely_unknown_text")
        assert isinstance(result, str)


class TestTaskDefinitions:
    """Tests for TASKS dictionary structure."""

    def test_tasks_has_required_types(self):
        """TASKS contains required task types."""
        required_tasks = [
            "single_count",
            "multi_count",
            "single_index_identification",
            "multi_index_identification",
            "constraint_generation"
        ]
        for task in required_tasks:
            assert task in TASKS, f"Missing task: {task}"

    def test_task_has_templates(self):
        """Each task has question_templates."""
        for task_name, task_def in TASKS.items():
            assert "question_templates" in task_def, f"{task_name} missing templates"
            assert len(task_def["question_templates"]) > 0

    def test_task_has_metadata(self):
        """Each task has required metadata."""
        for task_name, task_def in TASKS.items():
            assert "task_name" in task_def
            assert "task_type" in task_def
            assert "output_type" in task_def


class TestSystemPrompts:
    """Tests for SYSTEM_PROMPTS dictionary."""

    def test_system_prompts_exist(self):
        """SYSTEM_PROMPTS dictionary exists and has content."""
        assert SYSTEM_PROMPTS is not None
        assert len(SYSTEM_PROMPTS) > 0

    def test_system_prompts_are_strings(self):
        """All prompts are non-empty strings."""
        for key, prompt in SYSTEM_PROMPTS.items():
            assert isinstance(prompt, str), f"{key} is not a string"
            assert len(prompt) > 100, f"{key} prompt too short"

    def test_with_key_hints_prompt(self):
        """with_key_hints prompt exists."""
        assert "with_key_hints" in SYSTEM_PROMPTS

    def test_concise_prompt(self):
        """concise prompt exists."""
        assert "concise" in SYSTEM_PROMPTS


class TestFormatterWithKeyHints:
    """Tests for key hint functionality."""

    def test_count_query_with_key_hint(self, formatter):
        """Count query can include key hints."""
        query = formatter.format_count_query(
            smiles="CCO",
            count_types=["ring_count"],
            include_key_hint=True,
            key_names=["ring_count"]
        )
        assert "ring_count" in query

    def test_index_query_with_key_hint(self, formatter):
        """Index query can include key hints."""
        query = formatter.format_index_query(
            smiles="CCO",
            index_types=["carbon_atom_index"],
            include_key_hint=True,
            key_names=["carbon_atom_index"]
        )
        assert "carbon_atom_index" in query
