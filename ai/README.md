# AI — opencode configuration for animals_vax_atlas

## Directory structure

```
ai/
├── README.md          ← this file
├── SKILL.md           ← coding style guide (r-tidyverse-vaxgo skill)
└── opencode.json      ← opencode project configuration
```

A global copy of the skill is also installed at:
`~/.config/opencode/skills/r-tidyverse-vaxgo/SKILL.md`

## Purpose

This skill captures the **tidyverse-oriented coding style** used throughout this project. It enables opencode to:

- Annotate and refactor existing code following project conventions
- Write new code that is consistent with the existing style
- Avoid re-explaining conventions at every session

## Usage

```bash
# Load the skill manually (works in any project)
/skill r-tidyverse-vaxgo

# Or open the project — opencode.json loads it automatically
opencode
```

## Style summary (30 seconds)

| Topic | Convention |
|-------|-----------|
| Pipe | `%>%` (not `|>`) |
| Paths | `here::here("folder", "file")` |
| Data I/O | `readRDS()` / `saveRDS()` |
| Joins | `inner_join()`, `full_join()`, `anti_join()`, `bind_rows()` |
| Factors | `fct_relevel()`, `fct_reorder()` (forcats) |
| Conditionals | `case_when()`, `if_else()` |
| ggplot2 | `aes()` outside geom + `theme_vaxgo()` |
| Naming | `snake_case`, descriptive |
| Project | Numbered Rmds (`0_`, `1_`...), `required.R` as central config |
| Namespace | Explicit: `dplyr::first()`, `rstatix::t_test()` |

## How the skill was created

The skill was derived by analyzing all R scripts in `scripts_notebooks/` and `example/`. Common patterns across pipes, joins, factors, ggplot2, naming, bioinformatics workflow, and helper functions were extracted and organized into the `SKILL.md` document.

## Maintenance

- **To update** the skill with new patterns, edit `ai/SKILL.md` (and the global copy at `~/.config/opencode/skills/r-tidyverse-vaxgo/SKILL.md`)
- **To verify** the skill is applied, ask opencode to annotate a block of code
- The `opencode.json` file binds the skill to this project automatically
