<div align="center">

# MolecularIQ Core

**Shared library for molecular property calculations and chemistry reasoning tasks**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![RDKit](https://img.shields.io/badge/Built%20with-RDKit-3838ff.svg)](https://www.rdkit.org/)
[![Tests](https://img.shields.io/badge/tests-1035%20passed-brightgreen.svg)]()

<p align="center">
  <em>Symbolic solvers, reward functions, and NLP formatting for the MolecularIQ ecosystem</em>
</p>

[Installation](#installation) • [Quick Start](#quick-start) • [API Overview](#api-overview) • [MolecularIQ Family](#moleculariq-family)

</div>

---

## 🎯 Overview

**MolecularIQ Core** provides the foundational components shared across the MolecularIQ ecosystem:

- **SymbolicSolver**: Compute molecular properties (rings, atoms, functional groups, stereochemistry, etc.)
- **Reward Functions**: Evaluate model predictions against ground truth
- **NaturalLanguageFormatter**: Convert between technical property keys and natural language
- **Task Definitions**: Standardized task configurations and question templates

## Installation

```bash
# From source (recommended for development)
pip install -e ".[dev]"

# From GitHub
pip install git+https://github.com/ml-jku/moleculariq-core.git
```

**Requirements**: Python 3.9+ and RDKit

## 🚀 Quick Start

### Compute Molecular Properties

```python
from moleculariq_core import SymbolicSolver

solver = SymbolicSolver()
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

# Counts
solver.get_ring_count(smiles)           # 1
solver.get_aromatic_ring_count(smiles)  # 1
solver.get_carbon_atom_count(smiles)    # 9
solver.get_hba_count(smiles)            # 3
solver.get_hbd_count(smiles)            # 1

# Indices (0-based atom positions)
solver.get_ring_indices(smiles)         # [[2, 3, 4, 5, 6, 7]]
solver.get_carbon_atom_indices(smiles)  # [0, 1, 2, 3, 4, 5, 6, 7, 8]

# Functional groups
solver.get_functional_group_count_and_indices(smiles, group_name="ester")
```

### Evaluate Predictions

```python
from moleculariq_core import evaluate_answer

# Count task
score = evaluate_answer(
    task_type="single_count",
    predicted={"ring_count": 1},
    target={"ring_count": 1}
)  # 1.0

# Index task
score = evaluate_answer(
    task_type="single_index",
    predicted={"ring_index": [0, 1, 2, 3, 4, 5]},
    target={"ring_index": [0, 1, 2, 3, 4, 5]}
)  # 1.0

# Constraint generation task
score = evaluate_answer(
    task_type="constraint_generation",
    predicted="c1ccccc1",  # Benzene
    constraints=[{"type": "ring_count", "operator": "=", "value": 1}]
)  # 1.0
```

### Format Questions

```python
from moleculariq_core import NaturalLanguageFormatter

formatter = NaturalLanguageFormatter(seed=42)

# Format count question
question = formatter.format_count_query(
    smiles="CCO",
    count_types=["ring_count", "carbon_atom_count"]
)
# "How many rings and carbon atoms are in CCO?"

# Format constraint
constraint_text = formatter.format_constraint({
    "type": "ring_count",
    "operator": ">=",
    "value": 2
})
# "at least 2 rings"
```

## API Overview

### Core Components

| Module | Description |
|--------|-------------|
| `SymbolicSolver` | Compute 100+ molecular properties from SMILES |
| `FunctionalGroupSolver` | SMARTS-based functional group detection |
| `TemplateBasedReactionSolver` | SMIRKS-based reaction site prediction |
| `NaturalLanguageFormatter` | Bidirectional NL ↔ technical key conversion |
| `evaluate_answer` | Unified reward function dispatcher |

### Property Categories

| Category | Examples |
|----------|----------|
| **Topology** | ring_count, fused_ring_count, bridgehead_atom_count |
| **Aromaticity** | aromatic_ring_count, aliphatic_ring_count |
| **Atoms** | carbon_atom_count, hetero_atom_count, halogen_atom_count |
| **Bonds** | rotatable_bond_count, double_bond_count |
| **H-bonding** | hba_count, hbd_count |
| **Stereochemistry** | stereocenter_count, r_s_stereocenter_r_count |
| **Functional Groups** | 60+ groups (alcohol, ketone, amine, etc.) |
| **Reactions** | 40+ templates (bromination, oxidation, etc.) |

### Task Types

| Task | Output | Description |
|------|--------|-------------|
| `single_count` | INTEGER | Count one property |
| `multi_count` | DICT | Count multiple properties |
| `single_index` | LIST | Identify atom indices for one property |
| `multi_index` | DICT | Identify indices for multiple properties |
| `constraint_generation` | SMILES | Generate molecule satisfying constraints |

## 👨‍👧‍👦 MolecularIQ Family

This package is part of the MolecularIQ ecosystem:

| Repository | Purpose | 
|------------|---------|
| **[moleculariq](https://github.com/ml-jku/moleculariq)** | Central hub for the MolecularIQ benchmark ecosystem|
| **[moleculariq-leaderboard](https://github.com/ml-jku/moleculariq-leaderboard)** | Leaderboard: HuggingFace space, displays results, handles submissions |
| 📍 **[moleculariq-core](#moleculariq-core)** | Shared library providing core functionality, e.g. symbolic verifiers and question formatting | 
| **[moleculariq-benchmark](https://github.com/ml-jku/moleculariq-benchmark)** | Dataset creation: task definitions, symbolic verifiers implementations, question generator| 
| **[moleculariq-eval](https://github.com/ml-jku/moleculariq-eval)** | Evaluation code: integration with lm-eval-harness, model configs, reward functions, extraction functions, and system prompts|
|**[moleculariqd](https://github.com/ml-jku/moleculariqd)**| Dynamic MolecularIQ version: includes molecule pools from which questions/samples can be created and evaluated on the fly| 

## License

MIT License - see [LICENSE](LICENSE) for details.

---

<div align="center">

**[Back to Top](#moleculariq-core)**

</div>
