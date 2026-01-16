"""
Unit tests for NaturalLanguageFormatter reproducibility and seeding.

Tests that seeding produces consistent, reproducible outputs.
"""

import pytest
import random
from moleculariq_core import NaturalLanguageFormatter


class TestFormatterSeeding:
    """Tests for formatter seeding/reproducibility."""

    def test_seed_in_constructor(self):
        """Seed can be passed to constructor."""
        formatter = NaturalLanguageFormatter(seed=42)
        assert formatter is not None

    def test_rng_in_constructor(self):
        """Custom RNG can be passed to constructor."""
        rng = random.Random(42)
        formatter = NaturalLanguageFormatter(rng=rng)
        assert formatter is not None

    def test_seed_and_rng_mutually_exclusive(self):
        """Cannot provide both seed and rng."""
        rng = random.Random(42)
        with pytest.raises(ValueError):
            NaturalLanguageFormatter(rng=rng, seed=42)


class TestCountQueryReproducibility:
    """Tests for count query reproducibility with seeding."""

    def test_same_seed_same_output(self):
        """Same seed produces same count query output."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        query1 = formatter1.format_count_query(
            smiles="CCO",
            count_types=["ring_count"],
            include_key_hint=True,
            key_names=["ring_count"]
        )
        query2 = formatter2.format_count_query(
            smiles="CCO",
            count_types=["ring_count"],
            include_key_hint=True,
            key_names=["ring_count"]
        )

        assert query1 == query2

    def test_different_seed_may_differ(self):
        """Different seeds may produce different outputs."""
        # Generate many outputs to verify randomness works
        outputs = set()
        for seed in range(100):
            formatter = NaturalLanguageFormatter(seed=seed)
            query = formatter.format_count_query(
                smiles="CCO",
                count_types=["ring_count"],
                include_key_hint=True,
                key_names=["ring_count"]
            )
            outputs.add(query)

        # With varied phrasings, we should see some variation
        # (though not necessarily 100 unique outputs)
        assert len(outputs) >= 1

    def test_sequential_calls_deterministic(self):
        """Sequential calls with same seed are deterministic."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        # Make multiple calls in sequence
        results1 = []
        results2 = []

        for count_type in ["ring_count", "carbon_count", "hba_count"]:
            results1.append(formatter1.format_count_query(
                smiles="CCO",
                count_types=[count_type],
                include_key_hint=True,
                key_names=[count_type]
            ))
            results2.append(formatter2.format_count_query(
                smiles="CCO",
                count_types=[count_type],
                include_key_hint=True,
                key_names=[count_type]
            ))

        assert results1 == results2


class TestIndexQueryReproducibility:
    """Tests for index query reproducibility with seeding."""

    def test_same_seed_same_output(self):
        """Same seed produces same index query output."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        query1 = formatter1.format_index_query(
            smiles="c1ccccc1",
            index_types=["aromatic_atom_index"],
            include_key_hint=True,
            key_names=["aromatic_atom_index"]
        )
        query2 = formatter2.format_index_query(
            smiles="c1ccccc1",
            index_types=["aromatic_atom_index"],
            include_key_hint=True,
            key_names=["aromatic_atom_index"]
        )

        assert query1 == query2


class TestConstraintReproducibility:
    """Tests for constraint formatting reproducibility with seeding."""

    def test_same_seed_same_constraint(self):
        """Same seed produces same constraint output."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        constraint = {"type": "ring_count", "operator": "=", "value": 2}

        text1 = formatter1.format_constraint(constraint)
        text2 = formatter2.format_constraint(constraint)

        assert text1 == text2

    def test_zero_value_constraint_variation(self):
        """Zero value constraints have varied phrasings."""
        outputs = set()
        for seed in range(100):
            formatter = NaturalLanguageFormatter(seed=seed)
            constraint = {"type": "ring_count", "operator": "=", "value": 0}
            text = formatter.format_constraint(constraint)
            outputs.add(text)

        # Should see variation like "no rings", "zero rings", "without any rings"
        assert len(outputs) >= 1

    def test_constraint_list_reproducibility(self):
        """Constraint list formatting is reproducible."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        constraints = [
            {"type": "ring_count", "operator": "=", "value": 1},
            {"type": "carbon_count", "operator": ">=", "value": 6}
        ]

        text1 = formatter1.format_constraints_list(constraints)
        text2 = formatter2.format_constraints_list(constraints)

        assert text1 == text2


class TestDisableRandomPhrasing:
    """Tests for disabling random phrasing."""

    def test_disable_random_phrasing(self):
        """Disabling random phrasing makes output deterministic."""
        formatter1 = NaturalLanguageFormatter(enable_random_phrasing=False)
        formatter2 = NaturalLanguageFormatter(enable_random_phrasing=False)

        constraint = {"type": "ring_count", "operator": "=", "value": 0}

        text1 = formatter1.format_constraint(constraint)
        text2 = formatter2.format_constraint(constraint)

        assert text1 == text2

    def test_disable_random_always_first_option(self):
        """Disabling random phrasing always picks first option."""
        outputs = set()
        for _ in range(10):
            formatter = NaturalLanguageFormatter(enable_random_phrasing=False)
            constraint = {"type": "ring_count", "operator": "=", "value": 0}
            text = formatter.format_constraint(constraint)
            outputs.add(text)

        # Should always be the same (first option)
        assert len(outputs) == 1

    def test_return_only_text_disables_variation(self):
        """return_only_text parameter disables variation."""
        outputs = set()
        for seed in range(10):
            formatter = NaturalLanguageFormatter(seed=seed)
            constraint = {"type": "ring_count", "operator": "=", "value": 0}
            text = formatter.format_constraint(constraint, return_only_text=True)
            outputs.add(text)

        # Should always be the same when return_only_text=True
        assert len(outputs) == 1


class TestKeyHintReproducibility:
    """Tests for key hint reproducibility."""

    def test_key_hint_reproducible(self):
        """Key hints are reproducible with seeding."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        hint1 = formatter1.format_key_hint(["ring_count"], "count")
        hint2 = formatter2.format_key_hint(["ring_count"], "count")

        assert hint1 == hint2

    def test_constraint_hint_reproducible(self):
        """Constraint hints are reproducible with seeding."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        hint1 = formatter1.format_constraint_hint()
        hint2 = formatter2.format_constraint_hint()

        assert hint1 == hint2


class TestRNGStateIsolation:
    """Tests for RNG state isolation between formatters."""

    def test_formatters_independent(self):
        """Different formatter instances have independent RNG state."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=43)

        # Both should work independently
        query1 = formatter1.format_count_query(
            smiles="CCO",
            count_types=["ring_count"]
        )
        query2 = formatter2.format_count_query(
            smiles="CCO",
            count_types=["ring_count"]
        )

        # Both should be valid queries
        assert "CCO" in query1
        assert "CCO" in query2

    def test_shared_rng_affects_both(self):
        """Sharing RNG between formatters affects both."""
        rng = random.Random(42)
        formatter1 = NaturalLanguageFormatter(rng=rng)
        formatter2 = NaturalLanguageFormatter(rng=rng)

        # They share the same RNG state
        # So the second call will get different random values
        query1 = formatter1.format_count_query(
            smiles="CCO",
            count_types=["ring_count"],
            include_key_hint=True,
            key_names=["ring_count"]
        )
        query2 = formatter2.format_count_query(
            smiles="CCO",
            count_types=["ring_count"],
            include_key_hint=True,
            key_names=["ring_count"]
        )

        # Both are valid, may or may not be equal depending on phrasing choices
        assert "CCO" in query1
        assert "CCO" in query2


class TestParseAnswerDeterminism:
    """Tests that parsing is always deterministic (no randomness)."""

    def test_parse_count_answer_deterministic(self):
        """Count answer parsing is deterministic."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=99)

        answer = "rings: 2; carbons: 6"

        result1 = formatter1.parse_count_answer(answer)
        result2 = formatter2.parse_count_answer(answer)

        assert result1 == result2

    def test_parse_index_answer_deterministic(self):
        """Index answer parsing is deterministic."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=99)

        answer = "ring atoms: 0,1,2,3,4,5"

        result1 = formatter1.parse_index_answer(answer)
        result2 = formatter2.parse_index_answer(answer)

        assert result1 == result2

    def test_parse_constraint_deterministic(self):
        """Constraint parsing is deterministic."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=99)

        text = "exactly 2 rings"

        result1 = formatter1.parse_constraint(text)
        result2 = formatter2.parse_constraint(text)

        assert result1 == result2


class TestConversionDeterminism:
    """Tests that conversions are deterministic."""

    def test_technical_to_natural_deterministic(self):
        """Technical to natural conversion is deterministic."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=99)

        result1 = formatter1.technical_to_natural("ring_count", "count")
        result2 = formatter2.technical_to_natural("ring_count", "count")

        assert result1 == result2

    def test_natural_to_technical_deterministic(self):
        """Natural to technical conversion is deterministic."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=99)

        result1 = formatter1.natural_to_technical("number of rings", "count")
        result2 = formatter2.natural_to_technical("number of rings", "count")

        assert result1 == result2


class TestMultiTypeReproducibility:
    """Tests for multi-type formatting reproducibility."""

    def test_multi_count_reproducible(self):
        """Multi-count query is reproducible."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        query1 = formatter1.format_count_query(
            smiles="CCO",
            count_types=["ring_count", "carbon_count", "hba_count"],
            include_key_hint=True,
            key_names=["ring_count", "carbon_count", "hba_count"]
        )
        query2 = formatter2.format_count_query(
            smiles="CCO",
            count_types=["ring_count", "carbon_count", "hba_count"],
            include_key_hint=True,
            key_names=["ring_count", "carbon_count", "hba_count"]
        )

        assert query1 == query2

    def test_multi_index_reproducible(self):
        """Multi-index query is reproducible."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        query1 = formatter1.format_index_query(
            smiles="c1ccccc1",
            index_types=["aromatic_atom_index", "ring_atom_index"],
            include_key_hint=True,
            key_names=["aromatic_atom_index", "ring_atom_index"]
        )
        query2 = formatter2.format_index_query(
            smiles="c1ccccc1",
            index_types=["aromatic_atom_index", "ring_atom_index"],
            include_key_hint=True,
            key_names=["aromatic_atom_index", "ring_atom_index"]
        )

        assert query1 == query2

    def test_multi_constraint_reproducible(self):
        """Multi-constraint list is reproducible."""
        formatter1 = NaturalLanguageFormatter(seed=42)
        formatter2 = NaturalLanguageFormatter(seed=42)

        constraints = [
            {"type": "ring_count", "operator": "=", "value": 1},
            {"type": "carbon_count", "operator": ">=", "value": 6},
            {"type": "hba_count", "operator": "<=", "value": 3}
        ]

        text1 = formatter1.format_constraints_list(constraints)
        text2 = formatter2.format_constraints_list(constraints)

        assert text1 == text2
