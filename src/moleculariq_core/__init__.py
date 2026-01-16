"""
MolecularIQ Core
================

Core shared modules for molecular property calculations for different chemistry tasks
and .

Quick Start
-----------

Compute properties::

    from moleculariq_core import SymbolicSolver

    solver = SymbolicSolver()
    smiles = "CCO"  # ethanol

    rings = solver.get_ring_count(smiles)
    carbons = solver.get_carbon_atom_count(smiles)

Generate questions::

    from moleculariq_core import NaturalLanguageFormatter, TASKS
    import random

    formatter = NaturalLanguageFormatter()
    template = random.choice(TASKS["single_count"]["question_templates"])
    question = formatter.format_count_query("CCO", ["ring_count"], template)

Evaluate answers::

    from moleculariq_core import evaluate_answer

    score = evaluate_answer(
        task_type="single_count",
        predicted={"ring_count": 0},
        target={"ring_count": 0}
    )

Modules
-------
- solver: Molecular property computation (SymbolicSolver)
- _nlp: Natural language formatting (internal)
- rewards: Reward/evaluation functions
- questions: Task definitions and prompts
- properties: Property category mappings
"""

__version__ = "0.1.0"

# =============================================================================
# Primary API - What most users need
# =============================================================================

from .solver import SymbolicSolver
from ._nlp import NaturalLanguageFormatter
from .questions import TASKS, SYSTEM_PROMPTS
from .rewards import chemical_reward, valid_smiles

# Friendly alias for chemical_reward
evaluate_answer = chemical_reward

# =============================================================================
# Secondary API - For advanced usage
# =============================================================================

# Additional solver classes (rarely needed directly)
from .solver import FunctionalGroupSolver, TemplateBasedReactionSolver

# NLP utilities
from ._nlp import (
    COUNT_MAPPINGS,
    INDEX_MAPPINGS,
    get_natural_language,
    parse_natural_language,
)

# Additional reward functions (use evaluate_answer for most cases)
from .rewards import (
    multi_count_dict_reward,
    single_count_reward,
    multi_index_identification_reward,
    single_index_reward,
    multi_constraint_generation_reward,
    constraint_reward,
    is_reasonable_molecule,
)

# Property definitions
from .properties import (
    COUNT_MAP,
    INDEX_MAP,
    CONSTRAINT_MAP,
    COUNT_TO_INDEX_MAP,
    KEY_ALIAS_MAP,
    get_alias,
    canonicalize_property_name,
)

# =============================================================================
# Public API
# =============================================================================

__all__ = [
    # Primary API
    'SymbolicSolver',
    'NaturalLanguageFormatter',
    'TASKS',
    'SYSTEM_PROMPTS',
    'evaluate_answer',
    'valid_smiles',

    # Secondary API - Solver
    'FunctionalGroupSolver',
    'TemplateBasedReactionSolver',

    # Secondary API - NLP
    'COUNT_MAPPINGS',
    'INDEX_MAPPINGS',
    'get_natural_language',
    'parse_natural_language',

    # Secondary API - Rewards
    'chemical_reward',
    'multi_count_dict_reward',
    'single_count_reward',
    'multi_index_identification_reward',
    'single_index_reward',
    'multi_constraint_generation_reward',
    'constraint_reward',
    'is_reasonable_molecule',

    # Secondary API - Properties
    'COUNT_MAP',
    'INDEX_MAP',
    'CONSTRAINT_MAP',
    'COUNT_TO_INDEX_MAP',
    'KEY_ALIAS_MAP',
    'get_alias',
    'canonicalize_property_name',
]
