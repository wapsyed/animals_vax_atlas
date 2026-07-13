#!/usr/bin/env python3
import json, os, re
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
RE_LOAD = re.compile(r'load\s*\(\s*["\']([^"\']+)["\']')
RE_RENDER_FULL = re.compile(r'rmarkdown\s*::\s*render\s*\(\s*["\']([^"\']+)["\']')

def make_node_id(fpath):
    rel = fpath.relative_to(PROJECT_ROOT)
    return f"file://{rel}"

def make_pkg_id(pkg):
    return "pkg://" + pkg

def classify_file(fpath):
    ext = fpath.suffix.lower()
    if ext in (".r", ".rmd"):
        return "notebook" if ext == ".rmd" else "r-script"
    if ext == ".rds":
        return "r-object"
    if ext == ".csv":
        return "table"
    if ext in (".tsv", ".txt"):
        return "data"
    if ext in (".md", ".html", ".qmd"):
        return "document"
    if ext in (".py",):
        return "python-script"
    if ext in (".yml", ".yaml", ".json", ".toml"):
        return "config"
    if ext in (".png", ".jpg", ".jpeg", ".svg", ".gif"):
        return "figure"
    if ext in (".pdb", ".cif", ".pse"):
        return "structure"
    if ext in (".bed",):
        return "genomic-data"
    if ext in (".xls",):
        return "spreadsheet"
    if ext in (".css", ".js"):
        return "web-asset"
    return "other"

def resolve_ref(ref, base_dir):
    ref = ref.strip().strip('"').strip("'")
    if not ref:
        return None
    if ref.startswith("/"):
        candidate = Path(ref)
    else:
        candidate = (base_dir / ref).resolve()
    try:
        candidate = candidate.resolve()
    except (OSError, RuntimeError):
        return None
    if candidate.exists():
        return candidate
    return None

def extract_markdown_links(filepath):
    edges = []
    try:
        text = filepath.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return edges
    for m in re.finditer(r'\[([^\]]*)\]\(([^)]+)\)', text):
        target_str = m.group(2).split("#")[0].split("?")[0].strip()
        if target_str and not target_str.startswith(("http://", "https://", "mailto:")):
            target = resolve_ref(target_str, filepath.parent)
            if target:
                edges.append({
                    "source": make_node_id(filepath),
                    "target": make_node_id(target),
                    "relation": "references",
                    "context": "markdown-link: " + m.group(2)[:60],
                })
    return edges


HERE_ARGS_RE = re.compile(r'here\s*::?\s*here\s*\(([^)]*)\)')
HERE_SIMPLE_RE = re.compile(r'here\s*\(([^)]*)\)')
STRING_ARG_RE = re.compile(r'"([^"]*)"')

def parse_here_args(text):
    paths = []
    for m in HERE_SIMPLE_RE.finditer(text):
        args_str = m.group(1)
        args = STRING_ARG_RE.findall(args_str)
        if args:
            joined = os.path.join(*args)
            paths.append(joined)
    for m in HERE_ARGS_RE.finditer(text):
        args_str = m.group(1)
        args = STRING_ARG_RE.findall(args_str)
        if args:
            joined = os.path.join(*args)
            paths.append(joined)
    return paths


FUNC_CALLS_WITH_HERE = [
    (r'readRDS\s*\(\s*', "reads", "readRDS"),
    (r'saveRDS\s*\(\s*(?:[^,]+,\s*)?', "writes", "saveRDS"),
    (r'read_csv\s*\(\s*', "reads", "read_csv"),
    (r'read_tsv\s*\(\s*', "reads", "read_tsv"),
    (r'read_delim\s*\(\s*', "reads", "read_delim"),
    (r'read_table\s*\(\s*', "reads", "read_table"),
    (r'read_excel\s*\(\s*', "reads", "read_excel"),
    (r'write_csv\s*\(\s*(?:[^,]+,\s*)?', "writes", "write_csv"),
    (r'write_tsv\s*\(\s*(?:[^,]+,\s*)?', "writes", "write_tsv"),
    (r'ggsave\s*\(\s*(?:[^,]+,\s*)?', "generates", "ggsave"),
]

def extract_r_rmd(filepath):
    edges = []
    packages = set()
    try:
        text = filepath.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return edges, []
    script_dir = filepath.parent

    for m in RE_LIBRARY.finditer(text):
        packages.add(m.group(1))
    for m in RE_SOURCE.finditer(text):
        target = resolve_ref(m.group(1), script_dir)
        if target:
            edges.append({"source": make_node_id(filepath), "target": make_node_id(target), "relation": "sources", "context": "source(" + m.group(1) + ")"})
    for m in RE_LOAD.finditer(text):
        target = resolve_ref(m.group(1), script_dir)
        if target:
            edges.append({"source": make_node_id(filepath), "target": make_node_id(target), "relation": "loads", "context": "load(" + m.group(1) + ")"})
    for m in RE_RENDER_FULL.finditer(text):
        target = resolve_ref(m.group(1), script_dir)
        if target:
            edges.append({"source": make_node_id(filepath), "target": make_node_id(target), "relation": "renders", "context": "rmarkdown::render(" + m.group(1) + ")"})

    for m in re.finditer(r'(?:readRDS|saveRDS|read_csv|read_tsv|read_delim|read_table|read_excel|read_xlsx|write_csv|write_tsv|ggsave|load|source)\s*\(\s*["\']([^"\']+)["\']', text):
        target = resolve_ref(m.group(1), script_dir)
        if target:
            func = m.group(0).split("(")[0]
            if func in ("saveRDS", "write_csv", "write_tsv", "write_xlsx", "write_excel"):
                rel = "writes"
            elif func == "ggsave":
                rel = "generates"
            else:
                rel = "reads"
            edges.append({"source": make_node_id(filepath), "target": make_node_id(target), "relation": rel, "context": m.group(0)[:60]})

    here_paths = parse_here_args(text)
    for hp in here_paths:
        candidate = (PROJECT_ROOT / hp).resolve()
        if candidate.exists():
            edges.append({"source": make_node_id(filepath), "target": make_node_id(candidate), "relation": "references", "context": "here(" + hp + ")"})

    READ_FUNCS = r'(?:readRDS|read_csv|read_tsv|read_delim|read_table|read_excel|read_xlsx)'
    SAVE_FUNCS = r'(?:saveRDS|write_csv|write_tsv|write_xlsx|write_excel)'
    GSAVE = r'ggsave'

    for m in re.finditer(READ_FUNCS + r'\s*\(\s*here\s*\(', text):
        nid = make_node_id(filepath)
        start = m.end()
        depth = 1
        i = start
        while i < len(text) and depth > 0:
            if text[i] == '(':
                depth += 1
            elif text[i] == ')':
                depth -= 1
            i += 1
        inner = text[start:i-1]
        args = STRING_ARG_RE.findall(inner)
        if args:
            joined = os.path.join(*args)
            candidate = (PROJECT_ROOT / joined).resolve()
            if candidate.exists():
                edges.append({"source": nid, "target": make_node_id(candidate), "relation": "reads", "context": "read(here(" + joined + "))"})

    for m in re.finditer(SAVE_FUNCS + r'\s*\(\s*(?:[^,]+,\s*)?here\s*\(', text):
        nid = make_node_id(filepath)
        start = m.end()
        depth = 1
        i = start
        while i < len(text) and depth > 0:
            if text[i] == '(':
                depth += 1
            elif text[i] == ')':
                depth -= 1
            i += 1
        inner = text[start:i-1]
        args = STRING_ARG_RE.findall(inner)
        if args:
            joined = os.path.join(*args)
            candidate = (PROJECT_ROOT / joined).resolve()
            if candidate.exists():
                edges.append({"source": nid, "target": make_node_id(candidate), "relation": "writes", "context": "save(here(" + joined + "))"})

    for m in re.finditer(GSAVE + r'\s*\(\s*(?:[^,]+,\s*)?here\s*\(', text):
        nid = make_node_id(filepath)
        start = m.end()
        depth = 1
        i = start
        while i < len(text) and depth > 0:
            if text[i] == '(':
                depth += 1
            elif text[i] == ')':
                depth -= 1
            i += 1
        inner = text[start:i-1]
        args = STRING_ARG_RE.findall(inner)
        if args:
            joined = os.path.join(*args)
            candidate = (PROJECT_ROOT / joined).resolve()
            if candidate.exists():
                edges.append({"source": nid, "target": make_node_id(candidate), "relation": "generates", "context": "ggsave(here(" + joined + "))"})

    return edges, sorted(set(packages))


def scan_project():
    nodes = {}
    edges = []
    for dir_name in INCLUDE_DIRS:
        dir_path = PROJECT_ROOT / dir_name
        if not dir_path.exists():
            print("  [SKIP] " + dir_name + " -- not found")
            continue
        print("  Scanning " + dir_name + "/ ...")
        for root, dirs, files in os.walk(dir_path):
            root_path = Path(root)
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.startswith(".")]
            for fname in sorted(files):
                fpath = root_path / fname
                ext = fpath.suffix.lower()
                if ext in SKIP_EXTS:
                    continue
                rel = fpath.relative_to(PROJECT_ROOT)
                nid = make_node_id(fpath)
                ftype = classify_file(fpath)
                try:
                    size = fpath.stat().st_size
                except OSError:
                    size = 0
                nodes[nid] = {
                    "id": nid,
                    "label": str(rel),
                    "type": ftype,
                    "file_path": str(fpath),
                    "size_bytes": size,
                    "properties": {"extension": ext, "directory": str(rel.parent)},
                }

    print("  Extracting dependencies from R/Rmd files ...")
    for nid, node in list(nodes.items()):
        fpath = Path(node["file_path"])
        ext = fpath.suffix.lower()
        if ext in (".r", ".rmd"):
            r_edges, pkgs = extract_r_rmd(fpath)
            edges.extend(r_edges)
            for pkg in pkgs:
                pkg_id = make_pkg_id(pkg)
                if pkg_id not in nodes:
                    nodes[pkg_id] = {
                        "id": pkg_id,
                        "label": pkg,
                        "type": "r-package",
                        "file_path": "",
                        "size_bytes": 0,
                        "properties": {"package": pkg},
                    }
                edges.append({
                    "source": make_node_id(fpath),
                    "target": pkg_id,
                    "relation": "imports",
                    "context": "library(" + pkg + ")",
                })

    print("  Extracting markdown links ...")
    for nid, node in list(nodes.items()):
        fpath = Path(node["file_path"])
        if fpath.suffix.lower() in (".md", ".mdx", ".qmd"):
            md_edges = extract_markdown_links(fpath)
            edges.extend(md_edges)

    notebook_order = [
        "0_Data_Curation.Rmd",
        "1_QualityControl.Rmd",
        "2_Preprocessing_and_DGE.Rmd",
        "3.1_Comparing_Human_Mouse_DGE_analyses.Rmd",
        "3.2_Comparing_Human_Mouse_ByCondition_GSEA.Rmd",
        "3.3_Comparing_Human_Mouse_ByCondition_Functional_Analyses.Rmd",
        "4_Performance.Rmd",
        "5_MachineLearning.Rmd",
    ]
    for i in range(len(notebook_order) - 1):
        src = "scripts_notebooks/" + notebook_order[i]
        tgt = "scripts_notebooks/" + notebook_order[i + 1]
        src_id = make_node_id(PROJECT_ROOT / src)
        tgt_id = make_node_id(PROJECT_ROOT / tgt)
        if src_id in nodes and tgt_id in nodes:
            edges.append({
                "source": src_id,
                "target": tgt_id,
                "relation": "precedes",
                "context": "analysis pipeline order",
            })

    return nodes, edges


def build_graph(nodes, edges):
    node_list = []
    for nid, node in nodes.items():
        node_list.append({
            "id": nid,
            "label": node["label"],
            "type": node["type"],
            "file_path": node.get("file_path", ""),
            "size_bytes": node.get("size_bytes", 0),
            "properties": node.get("properties", {}),
        })
    edge_list = []
    for e in edges:
        edge_list.append({
            "source": e["source"],
            "target": e["target"],
            "relation": e["relation"],
            "context": e.get("context", ""),
        })
    return {
        "directed": True,
        "nodes": node_list,
        "edges": edge_list,
        "metadata": {
            "project": "Animals Vax Atlas",
            "description": "Cross-species transcriptomic analysis of vaccine and infection responses",
            "node_count": len(node_list),
            "edge_count": len(edge_list),
            "generated_by": "graphify-custom-extractor",
        },
    }


def main():
    print("=" * 60)
    print("Animals Vax Atlas - Graphify Knowledge Graph Extraction")
    print("=" * 60)
    print("Project root: " + str(PROJECT_ROOT))
    print("Output: " + str(OUT_DIR))
    print()

    nodes, edges = scan_project()
    graph = build_graph(nodes, edges)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "graph.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2, ensure_ascii=False)

    print()
    print("=" * 60)
    print("Graph built successfully!")
    print("  Nodes: " + str(graph["metadata"]["node_count"]))
    print("  Edges: " + str(graph["metadata"]["edge_count"]))
    print("  Output: " + str(out_path))

    type_counts = {}
    for n in graph["nodes"]:
        t = n["type"]
        type_counts[t] = type_counts.get(t, 0) + 1
    print()
    print("  Node types:")
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print("    " + t + ": " + str(c))

    rel_counts = {}
    for e in graph["edges"]:
        r = e["relation"]
        rel_counts[r] = rel_counts.get(r, 0) + 1
    print()
    print("  Edge relations:")
    for r, c in sorted(rel_counts.items(), key=lambda x: -x[1]):
        print("    " + r + ": " + str(c))


if __name__ == "__main__":
    main()
