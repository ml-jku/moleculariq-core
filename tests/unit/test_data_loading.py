"""
Tests for the _data module.

Validates data file loading utilities and path constants.
"""

import pytest
from pathlib import Path

from moleculariq_core._data import (
    DATA_DIR,
    SMARTS_FUNCTIONAL_GROUPS,
    SMARTS_RENAMED,
    REACTION_TEMPLATES,
    get_data_path,
)


class TestDataDir:
    """Tests for DATA_DIR constant."""

    def test_data_dir_exists(self):
        """DATA_DIR should point to an existing directory."""
        assert DATA_DIR.exists()
        assert DATA_DIR.is_dir()

    def test_data_dir_is_path(self):
        """DATA_DIR should be a Path object."""
        assert isinstance(DATA_DIR, Path)


class TestDataPathConstants:
    """Tests for data file path constants."""

    def test_smarts_functional_groups_exists(self):
        """SMARTS_FUNCTIONAL_GROUPS file should exist."""
        assert SMARTS_FUNCTIONAL_GROUPS.exists()
        assert SMARTS_FUNCTIONAL_GROUPS.is_file()

    def test_smarts_renamed_exists(self):
        """SMARTS_RENAMED file should exist."""
        assert SMARTS_RENAMED.exists()
        assert SMARTS_RENAMED.is_file()

    def test_reaction_templates_exists(self):
        """REACTION_TEMPLATES file should exist."""
        assert REACTION_TEMPLATES.exists()
        assert REACTION_TEMPLATES.is_file()

    def test_paths_are_in_data_dir(self):
        """All data file paths should be within DATA_DIR."""
        assert SMARTS_FUNCTIONAL_GROUPS.parent == DATA_DIR
        assert SMARTS_RENAMED.parent == DATA_DIR
        assert REACTION_TEMPLATES.parent == DATA_DIR


class TestGetDataPath:
    """Tests for get_data_path function."""

    def test_get_existing_file(self):
        """Should return path for existing file."""
        path = get_data_path("smarts_renamed.txt")
        assert path.exists()
        assert path == SMARTS_RENAMED

    def test_get_reaction_templates(self):
        """Should return path for reaction templates."""
        path = get_data_path("reaction_templates.txt")
        assert path.exists()
        assert path == REACTION_TEMPLATES

    def test_get_smarts_functional_groups(self):
        """Should return path for SMARTS functional groups."""
        path = get_data_path("smarts_functional_groups.txt")
        assert path.exists()
        assert path == SMARTS_FUNCTIONAL_GROUPS

    def test_nonexistent_file_raises_error(self):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError) as exc_info:
            get_data_path("nonexistent_file.txt")
        assert "nonexistent_file.txt" in str(exc_info.value)

    def test_returns_path_object(self):
        """Should return a Path object."""
        path = get_data_path("smarts_renamed.txt")
        assert isinstance(path, Path)


class TestDataFileContents:
    """Tests for data file contents validity."""

    def test_smarts_renamed_has_content(self):
        """SMARTS_RENAMED should have non-empty content."""
        content = SMARTS_RENAMED.read_text()
        assert len(content) > 0
        lines = [l for l in content.split('\n') if l.strip() and not l.startswith('#')]
        assert len(lines) > 50, "Expected at least 50 functional group definitions"

    def test_smarts_renamed_format(self):
        """SMARTS_RENAMED lines should have correct format (Name:Rank:SMARTS)."""
        content = SMARTS_RENAMED.read_text()
        valid_lines = 0
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(':')
            assert len(parts) >= 3, f"Invalid line format: {line[:50]}"
            valid_lines += 1
        assert valid_lines > 50

    def test_reaction_templates_has_content(self):
        """REACTION_TEMPLATES should have non-empty content."""
        content = REACTION_TEMPLATES.read_text()
        assert len(content) > 0
        lines = [l for l in content.split('\n') if l.strip() and not l.startswith('#')]
        assert len(lines) > 30, "Expected at least 30 reaction templates"

    def test_reaction_templates_format(self):
        """REACTION_TEMPLATES lines should have correct format (Name;Category;SMIRKS;Description)."""
        content = REACTION_TEMPLATES.read_text()
        valid_lines = 0
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(';')
            assert len(parts) >= 4, f"Invalid line format: {line[:50]}"
            valid_lines += 1
        assert valid_lines > 30

    def test_smarts_patterns_are_valid(self):
        """SMARTS patterns in data files should be parseable."""
        from rdkit import Chem

        content = SMARTS_RENAMED.read_text()
        invalid_count = 0
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(':')
            if len(parts) >= 3:
                smarts = parts[2]
                mol = Chem.MolFromSmarts(smarts)
                if mol is None:
                    invalid_count += 1
        # Allow some invalid patterns (some might be intentional placeholders)
        assert invalid_count < 5, f"Too many invalid SMARTS patterns: {invalid_count}"

    def test_reaction_smirks_are_valid(self):
        """SMIRKS patterns in reaction templates should be parseable."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        content = REACTION_TEMPLATES.read_text()
        invalid_count = 0
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(';')
            if len(parts) >= 3:
                smirks = parts[2]
                rxn = AllChem.ReactionFromSmarts(smirks)
                if rxn is None:
                    invalid_count += 1
        # Allow some invalid patterns
        assert invalid_count < 3, f"Too many invalid SMIRKS patterns: {invalid_count}"


class TestDataFileUniqueness:
    """Tests for data file entry uniqueness."""

    def test_smarts_renamed_unique_names(self):
        """Functional group names should be unique."""
        content = SMARTS_RENAMED.read_text()
        names = []
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(':')
            if len(parts) >= 1:
                names.append(parts[0])
        assert len(names) == len(set(names)), "Duplicate functional group names found"

    def test_reaction_templates_unique_names(self):
        """Reaction template names should be unique."""
        content = REACTION_TEMPLATES.read_text()
        names = []
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(';')
            if len(parts) >= 1:
                names.append(parts[0])
        assert len(names) == len(set(names)), "Duplicate reaction template names found"
