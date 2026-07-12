#!/usr/bin/env python3
"""
Custom graphify extraction for the Animals Vax Atlas project.
Scans all directories and builds a knowledge graph of scripts, tables,
documents, figures, and their relationships.

Output: graphify-out/graph.json (graphify-compatible format)
"""

import json
import os
import re
from pathlib import Path

PROJECT_ROOT = Path("/home/wasim/Área de trabalho/Github/animals_vax_atlas").resolve()
OUT_DIR = PROJECT_ROOT / "AI" / "graphify-out"

INCLUDE_DIRS = [
    "scripts_notebooks", "tables", "DataCuration", "Quality control",
    "Structures", "Genomic", "VaxGO", "ArrayQM", "Figures", "example",
]

SKIP_DIRS = {".git", ".Rproj.user", "__pycache__"}
SKIP_EXTS = {".pdf", ".png", ".jpg", ".jpeg", ".svg", ".gif", ".zip", ".gz", ".xlsx", ".pse", ".cif"}

RE_LIBRARY = re.compile(r'(?:library|require)\s*\(\s*["\']?(\w+)["\']?\s*\)')
RE_SOURCE = re.compile(r'source\s*\(\s*["\']([^"\']+)["\']\s*\)')
RE_READRDS = re.compile(r'readRDS\s*\(\s*["\']([^"\']+)["\']')
RE_WRITERDS = re.compile(r'saveRDS\s*\(\s*(?:[^,]+,\s*)?["\']([^"\']+)["\']')
RE_READ_CSV = re.compile(r'(?:read_csv|read_tsv|read_delim|read_table|read_excel|read_xlsx)\s*\(\s*["\']([^"\']+)["\']')
RE_WRITE_CSV = re.compile(r'(?:write_csv|write_tsv|write_xlsx|write_excel)\s*\(\s*(?:[^,]+,\s*)?["\']([^"\']+)["\']')
RE_GGSAVE = re.compile(r'ggsave\s*\(\s*(?:[^,]+,\s*)?["\']([^"\']+)["\']')
RE_LOAD = re.compile(r'load\s*\(\s*["\']([^"\']+)["\']')
RE_RENDER_FULL = re.compile(r'rmarkdown\s*::\s*render\s*\(\s*["\']([^"\']+)["\']')

# Also check how many here() calls use string literals vs variables
here_literal = re.compile(r'here\s*\(\s*["\']')
here_multi = re.compile(r'here\s*\([^)]*\)')
here_single_literal = re.compile(r'here\s*\(\s*["\'][^"\']+["\']\s*\)')
here_multi_literal = re.compile(r'here\s*\(\s*["\'][^"\']+["\']\s*,\s*["\'][^"\']+["\']')

total_here = 0
literal_here = 0
multi_literal = 0
for f in files:
    try:
        text = f.read_text(errors='replace')
    except:
        continue
    total_here += len(re.findall(r'here\s*\(', text))
    literal_here += len(here_literal.findall(text))
    multi_literal += len(here_multi.findall(text))

print(f"Total here(): {total_here}")
print(f"Single literal: {literal_here}")
print(f"Multi literal: {multi_literal}")
