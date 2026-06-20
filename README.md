# CDI Proteomics Analysis System

**Results-First Edition**

This project provides a reproducible system for moving from differential protein abundance result tables to biological interpretation.

## Core Workflow

```text
Proteomics result tables
        ↓
Quality control
        ↓
Differential protein filtering
        ↓
Protein ranking
        ↓
Identifier cleaning and annotation
        ↓
GO and pathway enrichment
        ↓
STRING network analysis
        ↓
Biological interpretation
        ↓
Reproducible report
```

## Render

```bash
quarto render
```

## First Test

```bash
python scripts/python/01-filter-differential-proteins.py
```
