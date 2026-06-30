#!/usr/bin/env python3
"""Mechanical lint checks for the markdown knowledge base."""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote, urlparse


ROOT = Path(__file__).resolve().parents[1]
WIKI = ROOT / "wiki"
RAW = ROOT / "raw"
AGENTS = ROOT / "AGENTS.md"
CATALOG = WIKI / "catalog.md"
GLOSSARY = WIKI / "glossary.md"
LOG = WIKI / "log.md"

REQUIRED_LEAF_SECTIONS = ["Summary", "Detail", "Sources", "See Also"]
ROUTER_NAMES = {"index.md"}
SUPPORT_PAGES = {
    WIKI / "SCHEMA.md",
    WIKI / "workflows.md",
    WIKI / "glossary.md",
    WIKI / "catalog.md",
    WIKI / "log.md",
    WIKI / "source-map.md",
    WIKI / "contradictions.md",
    WIKI / "lint" / "CHECKS.md",
}
EXEMPT_INDEX_DIRS = {WIKI / "lint"}
SKIP_LINK_SCHEMES = {"http", "https", "mailto", "tel"}
LOG_HEAD_RE = re.compile(
    r"^## \[\d{4}-\d{2}-\d{2}\] "
    r"(ingest|query-filed|lint|source-drift|reconcile|page-add|page-move|page-delete) \| .+"
)
LOG_FORBIDDEN_RE = re.compile(
    r"\b(chat transcript|token report|routine lint dump|agent-status|scratch output|planning dump)\b",
    re.IGNORECASE,
)
ROOT_KB_ARTIFACT_RE = re.compile(r"^kb_audit_.*[.]md$")


@dataclass(frozen=True)
class Link:
    source: Path
    line: int
    label: str
    target: str
    target_path: Path | None
    anchor: str
    is_source_section: bool


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def fail(failures: list[str], path: Path, line: int | None, message: str) -> None:
    prefix = rel(path)
    if line is not None:
        prefix += f":{line}"
    failures.append(f"{prefix}: {message}")


def iter_md_files() -> list[Path]:
    files = [AGENTS]
    if WIKI.exists():
        files.extend(sorted(WIKI.rglob("*.md")))
    sources = RAW / "SOURCES.md"
    if sources.exists():
        files.append(sources)
    return sorted({p.resolve() for p in files})


def wiki_pages() -> list[Path]:
    if not WIKI.exists():
        return []
    return sorted(p.resolve() for p in WIKI.rglob("*.md"))


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def line_for_offset(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def mask_code_fences(text: str) -> str:
    def repl(match: re.Match[str]) -> str:
        return "".join("\n" if ch == "\n" else " " for ch in match.group(0))

    return re.sub(r"```.*?```", repl, text, flags=re.DOTALL)


def slugify_heading(title: str) -> str:
    slug = re.sub(r"[^\w\s-]", "", title.strip().lower())
    return re.sub(r"\s+", "-", slug)


def heading_anchors(path: Path) -> set[str]:
    anchors: set[str] = set()
    for line in read(path).splitlines():
        match = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
        if not match:
            continue
        title = re.sub(r"\s+#+\s*$", "", match.group(2)).strip()
        anchors.add(slugify_heading(title))
    return anchors


def headings(text: str) -> list[tuple[int, str]]:
    result: list[tuple[int, str]] = []
    for match in re.finditer(r"^(#{2,6})\s+(.+?)\s*$", text, re.MULTILINE):
        title = re.sub(r"\s+#+\s*$", "", match.group(2)).strip()
        result.append((len(match.group(1)), title))
    return result


def section_body(text: str, section: str) -> str:
    match = re.search(rf"^## {re.escape(section)}\s*$", text, re.MULTILINE)
    if not match:
        return ""
    start = match.end()
    next_match = re.search(r"^##\s+", text[start:], re.MULTILINE)
    end = start + next_match.start() if next_match else len(text)
    return text[start:end].strip()


def path_from_link(source: Path, target: str) -> tuple[Path | None, str]:
    target = target.strip()
    if not target:
        return None, ""
    parsed = urlparse(target)
    if parsed.scheme in SKIP_LINK_SCHEMES or parsed.netloc:
        return None, ""
    raw_path = parsed.path
    anchor = unquote(parsed.fragment or "")
    if not raw_path:
        return source, anchor
    if raw_path.startswith("/"):
        candidate = ROOT / raw_path.lstrip("/")
    else:
        candidate = source.parent / unquote(raw_path)
    return candidate.resolve(), anchor


def in_sources_section(text: str, offset: int) -> bool:
    sources_start = re.search(r"^## Sources\s*$", text, re.MULTILINE)
    if not sources_start:
        return False
    next_match = re.search(r"^##\s+", text[sources_start.end() :], re.MULTILINE)
    sources_end = sources_start.end() + next_match.start() if next_match else len(text)
    return sources_start.start() <= offset < sources_end


def extract_links(path: Path) -> list[Link]:
    text = read(path)
    search_text = mask_code_fences(text)
    pattern = re.compile(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)")
    links: list[Link] = []
    for match in pattern.finditer(search_text):
        label, raw_target = match.groups()
        target = raw_target.strip()
        path_part, anchor = path_from_link(path, target)
        links.append(
            Link(
                source=path,
                line=line_for_offset(search_text, match.start()),
                label=label.strip(),
                target=target,
                target_path=path_part,
                anchor=anchor,
                is_source_section=in_sources_section(search_text, match.start()),
            )
        )
    return links


def is_router(path: Path) -> bool:
    return path.name in ROUTER_NAMES or path == AGENTS


def is_support_page(path: Path) -> bool:
    return path in SUPPORT_PAGES


def is_normal_leaf(path: Path) -> bool:
    return path.is_relative_to(WIKI) and not is_router(path) and not is_support_page(path)


def check_wikilinks(files: list[Path], failures: list[str]) -> None:
    pattern = re.compile(r"\[\[[^\]]+\]\]")
    for path in files:
        text = read(path)
        for match in pattern.finditer(text):
            fail(failures, path, line_for_offset(text, match.start()), "wikilink found")


def check_links(files: list[Path], failures: list[str]) -> dict[Path, list[Path]]:
    graph: dict[Path, list[Path]] = defaultdict(list)
    anchor_cache: dict[Path, set[str]] = {}
    for path in files:
        for link in extract_links(path):
            target = link.target_path
            if target is None:
                continue
            if not target.exists():
                message = f"link target missing: {link.target}"
                if link.is_source_section:
                    message = f"source link target missing: {link.target}"
                fail(failures, path, link.line, message)
                continue
            if link.anchor:
                anchor_cache.setdefault(target, heading_anchors(target))
                if slugify_heading(link.anchor) not in anchor_cache[target]:
                    fail(failures, path, link.line, f"link anchor missing: {link.target}")
                    continue
            graph[path].append(target)
    return graph


def check_router_detail(files: list[Path], failures: list[str]) -> None:
    for path in files:
        if path == AGENTS or path.is_relative_to(WIKI):
            if is_router(path) and re.search(r"^## Detail\s*$", read(path), re.MULTILINE):
                fail(failures, path, None, "router has Detail section")


def check_leaf_contract(pages: list[Path], failures: list[str]) -> None:
    for path in pages:
        if not is_normal_leaf(path):
            continue
        text = read(path)
        lines = text.splitlines()
        if not lines or not lines[0].startswith("# "):
            fail(failures, path, 1, "leaf missing H1 title")
        header_lines = [line for line in lines[1:4] if line.startswith("> ")]
        if len(header_lines) < 2:
            fail(failures, path, 2, "leaf missing two-line blockquote header")

        found = [title for level, title in headings(text) if level == 2]
        positions = []
        for section in REQUIRED_LEAF_SECTIONS:
            try:
                positions.append(found.index(section))
            except ValueError:
                fail(failures, path, None, f"leaf missing section: {section}")
                positions.append(-1)
        if all(pos >= 0 for pos in positions) and positions != sorted(positions):
            fail(
                failures,
                path,
                None,
                "leaf sections out of order; want Summary, Detail, Sources, See Also",
            )
        for section in ("Sources", "See Also"):
            if not section_body(text, section):
                fail(failures, path, None, f"leaf has empty {section}")


def check_reachability(graph: dict[Path, list[Path]], failures: list[str]) -> None:
    if not AGENTS.exists():
        return
    reachable: set[Path] = set()
    queue: deque[Path] = deque([AGENTS.resolve()])
    while queue:
        path = queue.popleft()
        if path in reachable:
            continue
        reachable.add(path)
        for target in graph.get(path, []):
            if target.is_file() and (
                target == AGENTS.resolve() or target.is_relative_to(WIKI)
            ):
                queue.append(target)

    for page in wiki_pages():
        if page not in reachable:
            fail(failures, page, None, "wiki page not reachable from AGENTS.md")


def check_catalog(failures: list[str]) -> None:
    if not CATALOG.exists():
        return
    header, rows = parse_table(CATALOG)
    try:
        page_i = header.index("Page")
    except ValueError:
        fail(failures, CATALOG, None, "catalog lacks Page column")
        return
    links: list[Path] = []
    for row in rows:
        if len(row) <= page_i:
            continue
        for link in re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", row[page_i]):
            target, _ = path_from_link(CATALOG, link.group(2))
            if target is not None and target.is_relative_to(WIKI):
                links.append(target)
    counts = Counter(links)
    for page in wiki_pages():
        count = counts.get(page, 0)
        if count == 0:
            fail(failures, CATALOG, None, f"missing wiki page: {rel(page)}")
        elif count > 1:
            fail(failures, CATALOG, None, f"duplicate wiki page: {rel(page)}")


def check_content_dir_indexes(failures: list[str]) -> None:
    for directory, _, files in os.walk(WIKI):
        path = Path(directory).resolve()
        if path in EXEMPT_INDEX_DIRS:
            continue
        if path == WIKI:
            continue
        if any(name.endswith(".md") for name in files) and "index.md" not in files:
            fail(failures, path, None, "wiki content directory lacks index.md")


def parse_table(path: Path) -> tuple[list[str], list[list[str]]]:
    if not path.exists():
        return [], []
    lines = read(path).splitlines()
    header: list[str] = []
    rows: list[list[str]] = []
    for line in lines:
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells):
            continue
        if not header:
            header = cells
        else:
            rows.append(cells)
    return header, rows


def check_glossary_catalog_signal(failures: list[str]) -> None:
    if not (CATALOG.exists() and GLOSSARY.exists()):
        return
    catalog_header, catalog_rows = parse_table(CATALOG)
    try:
        page_i = catalog_header.index("Page")
        candidate_i = catalog_header.index("Concept hub candidate")
    except ValueError:
        fail(failures, CATALOG, None, "catalog lacks Concept hub candidate column")
        return

    candidates: set[str] = set()
    for row in catalog_rows:
        if len(row) <= max(page_i, candidate_i):
            continue
        marker = row[candidate_i].strip().lower()
        if marker and marker not in {"-", "no", "n/a", "none"}:
            candidates.add(row[page_i].lower())

    glossary_header, glossary_rows = parse_table(GLOSSARY)
    try:
        term_i = glossary_header.index("Term")
        meaning_i = glossary_header.index("Meaning")
    except ValueError:
        fail(failures, GLOSSARY, None, "glossary lacks Term/Meaning table")
        return

    for row in glossary_rows:
        if len(row) <= max(term_i, meaning_i):
            continue
        term = row[term_i].strip()
        meaning = row[meaning_i].strip()
        has_owner = bool(re.search(r"\[[^\]]+\]\([^)]+\)", meaning))
        external = "external" in meaning.lower() or "background" in meaning.lower()
        candidate = any(term.lower() in entry for entry in candidates)
        if not (has_owner or external or candidate):
            fail(
                failures,
                GLOSSARY,
                None,
                f"glossary term lacks owner/external/catalog candidate signal: {term}",
            )


def check_log_format(failures: list[str]) -> None:
    if not LOG.exists():
        fail(failures, LOG, None, "missing KB log")
        return
    text = read(LOG)
    for idx, line in enumerate(text.splitlines(), start=1):
        if line.startswith("## [") and not LOG_HEAD_RE.fullmatch(line):
            fail(failures, LOG, idx, "unparseable log heading")

    required = [
        "- `Sources:`",
        "- `Pages touched:`",
        "- `Decision:`",
        "- `Checks:`",
        "- `Follow-up:`",
    ]
    chunks = re.split(r"(?=^## \[\d{4}-\d{2}-\d{2}\] )", text, flags=re.MULTILINE)
    for chunk in chunks:
        if not chunk.startswith("## "):
            continue
        heading = chunk.splitlines()[0]
        if LOG_FORBIDDEN_RE.search(chunk):
            fail(failures, LOG, line_for_offset(text, text.find(chunk)), "forbidden transient log content")
        missing = [field for field in required if field not in chunk]
        if missing:
            fail(
                failures,
                LOG,
                line_for_offset(text, text.find(heading)),
                f"log entry missing fields: {', '.join(missing)}",
            )


def check_root_and_generated_artifacts(failures: list[str]) -> None:
    allowed_root_docs = {"AGENTS.md", "README.md", "CITATION.md", "coding_style.md"}
    for path in ROOT.glob("*.md"):
        if path.name in allowed_root_docs:
            continue
        if ROOT_KB_ARTIFACT_RE.match(path.name):
            fail(failures, path, None, "root-level KB maintenance artifact")
    for path in ROOT.rglob("__pycache__"):
        fail(failures, path, None, "generated Python artifact")
    for path in ROOT.rglob("*.pyc"):
        fail(failures, path, None, "generated Python artifact")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all",
        action="store_true",
        help="accepted for future expansion; current default already checks all",
    )
    parser.parse_args()

    failures: list[str] = []
    files = iter_md_files()
    pages = wiki_pages()

    check_wikilinks(files, failures)
    graph = check_links(files, failures)
    check_router_detail(files, failures)
    check_leaf_contract(pages, failures)
    check_reachability(graph, failures)
    check_catalog(failures)
    check_content_dir_indexes(failures)
    check_glossary_catalog_signal(failures)
    check_log_format(failures)
    check_root_and_generated_artifacts(failures)

    if failures:
        print("KB lint failed:")
        for item in sorted(failures):
            print(f"- {item}")
        return 1

    print("KB lint passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
