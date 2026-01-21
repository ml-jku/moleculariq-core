<div align="center">

# MolecularIQ Core

**The complete library for molecular reasoning: question generation, property computation, and answer evaluation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![RDKit](https://img.shields.io/badge/Built%20with-RDKit-3838ff.svg)](https://www.rdkit.org/)
[![Tests](https://img.shields.io/badge/tests-1092%20passed-brightgreen.svg)]()

<p align="center">
  <em>Everything you need for molecular reasoning benchmarks in one package</em>
</p>

[Installation](#installation) • [Quick Start](#quick-start) • [API Overview](#api-overview) • [MolecularIQ Family](#moleculariq-family)

</div>

---

## Overview

**MolecularIQ Core** is the central library for the MolecularIQ benchmark ecosystem:

- **MolecularIQD**: High-level API for dynamic question generation and evaluation
- **Molecule Pools**: Access training molecules (validation pools hidden to prevent leakage)
- **SymbolicSolver**: Compute 100+ molecular properties from SMILES
- **Reward Functions**: Evaluate model predictions against ground truth
- **NaturalLanguageFormatter**: Convert between technical keys and natural language

## Installation

```bash
# From GitHub
pip install git+https://github.com/ml-jku/moleculariq-core.git

# From source (for development)
pip install -e ".[dev]"
```

**Requirements**: Python 3.9+ and RDKit

## Quick Start

### The Easy Way: MolecularIQD

For most users, `MolecularIQD` is the recommended entry point:

```python
from moleculariq_core import MolecularIQD

mqd = MolecularIQD(seed=42)

# Generate a count question
question, answer, metadata = mqd.generate_count_question(
    smiles="c1ccccc1",
    count_properties="aromatic_ring_count"
)
print(question)
# "How many aromatic rings are in c1ccccc1? Return the result as JSON with key `aromatic_ring_count`."
print(answer)
# {'aromatic_ring_count': 1}

# Generate an index question
question, answer, metadata = mqd.generate_index_question(
    smiles="CCO",
    index_properties="carbon_atom_index"
)
print(answer)
# {'carbon_atom_index': [0, 1]}

# Validate predictions
score = mqd.validate_count_answer("c1ccccc1", {"aromatic_ring_count": 1})
# 1.0

# Generate a constraint question (molecule generation task)
question, metadata = mqd.generate_constraint_question(
    constraints=[{"property": "ring_count", "operator": ">=", "value": 2}]
)
# "Generate a molecule with at least 2 rings..."

# Validate a generated molecule against constraints
score = mqd.validate_constraint_answer(
    "c1ccc2ccccc2c1",  # Naphthalene (2 fused rings)
    metadata["constraints"]
)
# 1.0
```

### Load Training Molecules

```python
from moleculariq_core import load_molecule_pool

# Load training pool for development
train_smiles = load_molecule_pool("train")
print(f"Loaded {len(train_smiles)} training molecules")

# Validation pools are hidden to prevent data leakage
load_molecule_pool("val_hard")  # Raises MoleculePoolHiddenError
```

### Training Loop Example

```python
from moleculariq_core import MolecularIQD, load_molecule_pool
import random

mqd = MolecularIQD(seed=42)
train_smiles = load_molecule_pool("train")

for epoch in range(num_epochs):
    smiles = random.choice(train_smiles)
    question, answer, _ = mqd.generate_count_question(smiles, "ring_count")

    prediction = model(question)
    reward = mqd.validate_count_answer(smiles, prediction, answer)
    # ... update model
```

### Low-Level Access (Advanced)

For custom pipelines, access the primitives directly:

```python
from moleculariq_core import SymbolicSolver, evaluate_answer

solver = SymbolicSolver()
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

# Compute properties
solver.get_ring_count(smiles)           # 1
solver.get_aromatic_ring_count(smiles)  # 1
solver.get_carbon_atom_count(smiles)    # 9

# Evaluate predictions
score = evaluate_answer(
    task_type="single_count",
    predicted={"ring_count": 1},
    target={"ring_count": 1}
)  # 1.0
```

## API Overview

### High-Level API

| Component | Description |
|-----------|-------------|
| `MolecularIQD` | All-in-one class for question generation and evaluation |
| `load_molecule_pool` | Load training molecules from HuggingFace |

### MolecularIQD Methods

| Method | Description |
|--------|-------------|
| `generate_count_question()` | Generate counting questions |
| `generate_index_question()` | Generate atom index identification questions |
| `generate_constraint_question()` | Generate molecule generation questions |
| `validate_count_answer()` | Validate count predictions |
| `validate_index_answer()` | Validate index predictions |
| `validate_constraint_answer()` | Validate generated molecules |
| `generate_paired_question()` | Generate matched count/index pairs |
| `compute_property()` | Compute any molecular property |

### Low-Level Primitives

| Module | Description |
|--------|-------------|
| `SymbolicSolver` | Compute 100+ molecular properties from SMILES |
| `FunctionalGroupSolver` | SMARTS-based functional group detection |
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

## MolecularIQ Family

This package is part of the MolecularIQ ecosystem:

| Repository | Purpose |
|------------|---------|
| **[moleculariq](https://anonymous.4open.science/r/moleculariq-EE40/README.md)** | Central hub for the MolecularIQ benchmark ecosystem |
| **moleculariq-leaderboard]** | Leaderboard: HuggingFace space for results and submissions; please see supplementary material |
| **[moleculariq-core](https://anonymous.4open.science/r/moleculariq-core-2F02/README.md)** | Core library: question generation, property computation, evaluation, molecule pools |
| **[moleculariq-benchmark](https://anonymous.4open.science/r/moleculariq-benchmark-5A40/README.md)** | Dataset creation pipeline |
| **moleculariq-eval** | Evaluation code: lm-eval-harness integration: please see supplementary material |

## License

MIT License - see [LICENSE](LICENSE) for details.

---

<div align="center">

**[Back to Top]([#moleculariq-core](https://anonymous.4open.science/r/moleculariq-core-2F02/README.md))**

</div>
