# Exact finite-assurance capsule

This directory reproduces the finite-atlas proposition (stable manuscript label `prop:finite-atlas`) without any private dependency.

## Contents

- `inputs/critical_ratio_skeleton.json`: exact `k=4,6` cells and four named
  algebraic weights;
- `schemas/`: strict Draft 2020-12 atlas and row schemas;
- `scripts/produce.py`: exact word-product implementation;
- `scripts/validate.py`: independent base-4-code reconstruction of both the
  critical-ratio skeleton and every atlas row;
- `results/certified_terminal_atlas.json`: canonical 30-row certificate;
- `results/certified_validation.json`: independent validation and mutation
  receipt.

## Stable payloads

```text
atlas      aca5d5fc81c839c218b6eefaaf00acdc733b6e86c342e3d13546c82b0ce07f1f
validation b79051b4190d6f13b9bb183fe1299fb67b34a1eef21bbed5bf09c516f243ce46
```

All mathematical decisions use exact integer, rational, or isolating-interval
arithmetic. Floating-point tolerances are forbidden. Runtime measurements do
not enter the scientific payload.
