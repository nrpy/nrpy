#!/usr/bin/env python3
"""
Mechanical lint checks for the markdown knowledge base.

Author: NRPy contributors
"""

import argparse
import os

# Regexes are needed for Markdown link, heading, table, and log-pattern checks.
import re
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, Dict, List, Match, Optional, Set, Tuple
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
class _Link:
    """Markdown link with resolved local target metadata."""

    source: Path
    line: int
    label: str
    target: str
    target_path: Optional[Path]
    anchor: str
    is_source_section: bool


def _rel(path: Path) -> str:
    """
    Return a repository-relative POSIX path.

    :param path: Absolute or relative path inside the repository.
    :return: POSIX path relative to the repository root.
    """
    return path.resolve().relative_to(ROOT).as_posix()


def _fail(failures: List[str], path: Path, line: Optional[int], message: str) -> None:
    """
    Append a formatted lint failure.

    :param failures: Mutable failure list.
    :param path: File or directory associated with the failure.
    :param line: Optional one-based line number.
    :param message: Failure detail.
    """
    prefix = _rel(path)
    if line is not None:
        prefix += f":{line}"
    failures.append(f"{prefix}: {message}")


def _iter_md_files() -> List[Path]:
    """
    Return markdown files governed by compiled-KB lint rules.

    :return: Markdown paths included in normal lint checks.
    """
    files = [AGENTS]
    if WIKI.exists():
        files.extend(sorted(WIKI.rglob("*.md")))
    sources = RAW / "SOURCES.md"
    if sources.exists():
        files.append(sources)
    return sorted({p.resolve() for p in files})


def _wiki_pages() -> List[Path]:
    """
    Return all wiki markdown pages.

    :return: Markdown files under ``wiki/``.
    """
    if not WIKI.exists():
        return []
    return sorted(p.resolve() for p in WIKI.rglob("*.md"))


def _read(path: Path) -> str:
    """
    Read a UTF-8 text file.

    :param path: File to read.
    :return: File contents.
    """
    return path.read_text(encoding="utf-8")


def _line_for_offset(text: str, offset: int) -> int:
    """
    Return one-based line number for character offset.

    :param text: Text containing the offset.
    :param offset: Zero-based character offset.
    :return: One-based line number.
    """
    return text.count("\n", 0, offset) + 1


def _mask_code_fences(text: str) -> str:
    """
    Mask fenced code while preserving offsets and line numbers.

    :param text: Markdown text.
    :return: Text with fenced-code characters replaced by spaces.
    """

    def repl(match: Match[str]) -> str:
        return "".join("\n" if ch == "\n" else " " for ch in match.group(0))

    return re.sub(r"```.*?```", repl, text, flags=re.DOTALL)


def _slugify_heading(title: str) -> str:
    """
    Convert a Markdown heading to GitHub-style anchor slug.

    :param title: Heading title.
    :return: Anchor slug.
    """
    slug = re.sub(r"[^\w\s-]", "", title.strip().lower())
    return re.sub(r"\s+", "-", slug)


def _heading_anchors(path: Path) -> Set[str]:
    """
    Return Markdown heading anchors defined in a file.

    :param path: Markdown file to scan.
    :return: Set of heading anchor slugs.
    """
    anchors: Set[str] = set()
    for line in _read(path).splitlines():
        match = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
        if not match:
            continue
        title = re.sub(r"\s+#+\s*$", "", match.group(2)).strip()
        anchors.add(_slugify_heading(title))
    return anchors


def _headings(text: str) -> List[Tuple[int, str]]:
    """
    Return Markdown headings as ``(level, title)`` tuples.

    :param text: Markdown text.
    :return: Heading level and title pairs.
    """
    result: List[Tuple[int, str]] = []
    for match in re.finditer(r"^(#{2,6})\s+(.+?)\s*$", text, re.MULTILINE):
        title = re.sub(r"\s+#+\s*$", "", match.group(2)).strip()
        result.append((len(match.group(1)), title))
    return result


def _section_body(text: str, section: str) -> str:
    """
    Return content under a second-level Markdown section.

    :param text: Markdown text.
    :param section: Section title without hashes.
    :return: Section body, or an empty string when absent.
    """
    match = re.search(rf"^## {re.escape(section)}\s*$", text, re.MULTILINE)
    if not match:
        return ""
    start = match.end()
    next_match = re.search(r"^##\s+", text[start:], re.MULTILINE)
    end = start + next_match.start() if next_match else len(text)
    return text[start:end].strip()


def _path_from_link(source: Path, target: str) -> Tuple[Optional[Path], str]:
    """
    Resolve a Markdown link target relative to a source file.

    :param source: Markdown file containing the link.
    :param target: Raw link target.
    :return: Resolved local path, if any, and decoded anchor.
    """
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


def _in_sources_section(text: str, offset: int) -> bool:
    """
    Return whether an offset lies inside a ``Sources`` section.

    :param text: Markdown text.
    :param offset: Character offset to test.
    :return: Whether the offset is inside ``## Sources``.
    """
    sources_start = re.search(r"^## Sources\s*$", text, re.MULTILINE)
    if not sources_start:
        return False
    next_match = re.search(r"^##\s+", text[sources_start.end() :], re.MULTILINE)
    sources_end = sources_start.end() + next_match.start() if next_match else len(text)
    return sources_start.start() <= offset < sources_end


def _extract_links(path: Path) -> List[_Link]:
    """
    Extract Markdown links outside fenced code blocks.

    :param path: Markdown file to scan.
    :return: Links with local resolution metadata.
    """
    text = _read(path)
    search_text = _mask_code_fences(text)
    pattern = re.compile(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)")
    links: List[_Link] = []
    for match in pattern.finditer(search_text):
        label, raw_target = match.groups()
        target = raw_target.strip()
        path_part, anchor = _path_from_link(path, target)
        links.append(
            _Link(
                source=path,
                line=_line_for_offset(search_text, match.start()),
                label=label.strip(),
                target=target,
                target_path=path_part,
                anchor=anchor,
                is_source_section=_in_sources_section(search_text, match.start()),
            )
        )
    return links


def _is_router(path: Path) -> bool:
    """
    Return whether a path is a KB router page.

    :param path: Candidate page path.
    :return: Whether the path is a router.
    """
    return path.name in ROUTER_NAMES or path == AGENTS


def _is_support_page(path: Path) -> bool:
    """
    Return whether a path is a governance/support page.

    :param path: Candidate page path.
    :return: Whether the path is a support page.
    """
    return path in SUPPORT_PAGES


def _is_normal_leaf(path: Path) -> bool:
    """
    Return whether a path is a normal sourced wiki leaf.

    :param path: Candidate page path.
    :return: Whether the path is a normal leaf.
    """
    return (
        path.is_relative_to(WIKI)
        and not _is_router(path)
        and not _is_support_page(path)
    )


def _check_wikilinks(files: List[Path], failures: List[str]) -> None:
    """
    Check governed files for Obsidian-style wikilinks.

    :param files: Markdown files to scan.
    :param failures: Mutable failure list.
    """
    pattern = re.compile(r"\[\[[^\]]+\]\]")
    for path in files:
        text = _read(path)
        for match in pattern.finditer(text):
            _fail(
                failures, path, _line_for_offset(text, match.start()), "wikilink found"
            )


def _check_links(files: List[Path], failures: List[str]) -> Dict[Path, List[Path]]:
    """
    Validate Markdown links and return local link graph.

    :param files: Markdown files to scan.
    :param failures: Mutable failure list.
    :return: Directed graph of local links.
    """
    graph: Dict[Path, List[Path]] = defaultdict(list)
    anchor_cache: Dict[Path, Set[str]] = {}
    for path in files:
        for link in _extract_links(path):
            target = link.target_path
            if target is None:
                continue
            if not target.exists():
                message = f"link target missing: {link.target}"
                if link.is_source_section:
                    message = f"source link target missing: {link.target}"
                _fail(failures, path, link.line, message)
                continue
            if link.anchor:
                anchor_cache.setdefault(target, _heading_anchors(target))
                if _slugify_heading(link.anchor) not in anchor_cache[target]:
                    _fail(
                        failures, path, link.line, f"link anchor missing: {link.target}"
                    )
                    continue
            graph[path].append(target)
    return graph


def _check_router_detail(files: List[Path], failures: List[str]) -> None:
    """
    Check router pages do not contain detail sections.

    :param files: Markdown files to scan.
    :param failures: Mutable failure list.
    """
    for path in files:
        if path == AGENTS or path.is_relative_to(WIKI):
            if _is_router(path) and re.search(
                r"^## Detail\s*$", _read(path), re.MULTILINE
            ):
                _fail(failures, path, None, "router has Detail section")


def _check_leaf_contract(pages: List[Path], failures: List[str]) -> None:
    """
    Check normal leaves satisfy the page contract.

    :param pages: Wiki pages to scan.
    :param failures: Mutable failure list.
    """
    for path in pages:
        if not _is_normal_leaf(path):
            continue
        text = _read(path)
        lines = text.splitlines()
        if not lines or not lines[0].startswith("# "):
            _fail(failures, path, 1, "leaf missing H1 title")
        header_lines = [line for line in lines[1:4] if line.startswith("> ")]
        if len(header_lines) < 2:
            _fail(failures, path, 2, "leaf missing two-line blockquote header")

        found = [title for level, title in _headings(text) if level == 2]
        positions = []
        for section in REQUIRED_LEAF_SECTIONS:
            try:
                positions.append(found.index(section))
            except ValueError:
                _fail(failures, path, None, f"leaf missing section: {section}")
                positions.append(-1)
        if all(pos >= 0 for pos in positions) and positions != sorted(positions):
            _fail(
                failures,
                path,
                None,
                "leaf sections out of order; want Summary, Detail, Sources, See Also",
            )
        for section in ("Sources", "See Also"):
            if not _section_body(text, section):
                _fail(failures, path, None, f"leaf has empty {section}")


def _check_reachability(graph: Dict[Path, List[Path]], failures: List[str]) -> None:
    """
    Check all wiki pages are reachable from ``AGENTS.md``.

    :param graph: Directed graph of local links.
    :param failures: Mutable failure list.
    """
    if not AGENTS.exists():
        return
    reachable: Set[Path] = set()
    queue: Deque[Path] = deque([AGENTS.resolve()])
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

    for page in _wiki_pages():
        if page not in reachable:
            _fail(failures, page, None, "wiki page not reachable from AGENTS.md")


def _check_catalog(failures: List[str]) -> None:
    """
    Check the global catalog lists each wiki page once.

    :param failures: Mutable failure list.
    """
    if not CATALOG.exists():
        return
    header, rows = _parse_table(CATALOG)
    try:
        page_i = header.index("Page")
    except ValueError:
        _fail(failures, CATALOG, None, "catalog lacks Page column")
        return
    links: List[Path] = []
    for row in rows:
        if len(row) <= page_i:
            continue
        for link in re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", row[page_i]):
            target, _ = _path_from_link(CATALOG, link.group(2))
            if target is not None and target.is_relative_to(WIKI):
                links.append(target)
    counts = Counter(links)
    for page in _wiki_pages():
        count = counts.get(page, 0)
        if count == 0:
            _fail(failures, CATALOG, None, f"missing wiki page: {_rel(page)}")
        elif count > 1:
            _fail(failures, CATALOG, None, f"duplicate wiki page: {_rel(page)}")


def _check_content_dir_indexes(failures: List[str]) -> None:
    """
    Check content-bearing wiki directories have routers.

    :param failures: Mutable failure list.
    """
    for directory, _, files in os.walk(WIKI):
        path = Path(directory).resolve()
        if path in EXEMPT_INDEX_DIRS:
            continue
        if path == WIKI:
            continue
        if any(name.endswith(".md") for name in files) and "index.md" not in files:
            _fail(failures, path, None, "wiki content directory lacks index.md")


def _parse_table(path: Path) -> Tuple[List[str], List[List[str]]]:
    """
    Parse the first simple Markdown table in a file.

    :param path: Markdown file containing a table.
    :return: Header cells and row cells.
    """
    if not path.exists():
        return [], []
    lines = _read(path).splitlines()
    header: List[str] = []
    rows: List[List[str]] = []
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


def _check_glossary_catalog_signal(failures: List[str]) -> None:
    """
    Check glossary terms have owner or concept-candidate signals.

    :param failures: Mutable failure list.
    """
    if not (CATALOG.exists() and GLOSSARY.exists()):
        return
    catalog_header, catalog_rows = _parse_table(CATALOG)
    try:
        page_i = catalog_header.index("Page")
        candidate_i = catalog_header.index("Concept hub candidate")
    except ValueError:
        _fail(failures, CATALOG, None, "catalog lacks Concept hub candidate column")
        return

    candidates: Set[str] = set()
    for row in catalog_rows:
        if len(row) <= max(page_i, candidate_i):
            continue
        marker = row[candidate_i].strip().lower()
        if marker and marker not in {"-", "no", "n/a", "none"}:
            candidates.add(row[page_i].lower())

    glossary_header, glossary_rows = _parse_table(GLOSSARY)
    try:
        term_i = glossary_header.index("Term")
        meaning_i = glossary_header.index("Meaning")
    except ValueError:
        _fail(failures, GLOSSARY, None, "glossary lacks Term/Meaning table")
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
            _fail(
                failures,
                GLOSSARY,
                None,
                f"glossary term lacks owner/external/catalog candidate signal: {term}",
            )


def _check_log_format(failures: List[str]) -> None:
    """
    Check KB log headings, required fields, and hygiene.

    :param failures: Mutable failure list.
    """
    if not LOG.exists():
        _fail(failures, LOG, None, "missing KB log")
        return
    text = _read(LOG)
    for idx, line in enumerate(text.splitlines(), start=1):
        if line.startswith("## [") and not LOG_HEAD_RE.fullmatch(line):
            _fail(failures, LOG, idx, "unparseable log heading")

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
            _fail(
                failures,
                LOG,
                _line_for_offset(text, text.find(chunk)),
                "forbidden transient log content",
            )
        missing = [field for field in required if field not in chunk]
        if missing:
            _fail(
                failures,
                LOG,
                _line_for_offset(text, text.find(heading)),
                f"log entry missing fields: {', '.join(missing)}",
            )


def _check_root_and_generated_artifacts(failures: List[str]) -> None:
    """
    Check root maintenance artifacts and generated Python files.

    :param failures: Mutable failure list.
    """
    allowed_root_docs = {"AGENTS.md", "README.md", "CITATION.md", "coding_style.md"}
    for path in ROOT.glob("*.md"):
        if path.name in allowed_root_docs:
            continue
        if ROOT_KB_ARTIFACT_RE.match(path.name):
            _fail(failures, path, None, "root-level KB maintenance artifact")
    for path in ROOT.rglob("__pycache__"):
        _fail(failures, path, None, "generated Python artifact")
    for path in ROOT.rglob("*.pyc"):
        _fail(failures, path, None, "generated Python artifact")


def _main() -> int:
    """
    Run KB lint checks.

    :return: Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all",
        action="store_true",
        help="accepted for future expansion; current default already checks all",
    )
    parser.parse_args()

    failures: List[str] = []
    files = _iter_md_files()
    pages = _wiki_pages()

    _check_wikilinks(files, failures)
    graph = _check_links(files, failures)
    _check_router_detail(files, failures)
    _check_leaf_contract(pages, failures)
    _check_reachability(graph, failures)
    _check_catalog(failures)
    _check_content_dir_indexes(failures)
    _check_glossary_catalog_signal(failures)
    _check_log_format(failures)
    _check_root_and_generated_artifacts(failures)

    if failures:
        print("KB lint failed:")
        for item in sorted(failures):
            print(f"- {item}")
        return 1

    print("KB lint passed.")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
