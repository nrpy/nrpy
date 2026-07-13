#!/usr/bin/env python3
"""
Mechanical lint checks for the markdown knowledge base.

Author: NRPy contributors
"""

import argparse
import os

# Regexes are needed for Markdown link, heading, and table checks.
import re
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from datetime import date as calendar_date
from pathlib import Path
from typing import Deque, Dict, List, Match, Optional, Set, Tuple
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
WIKI = ROOT / "wiki"
RAW = ROOT / "raw"
AGENTS = ROOT / "AGENTS.md"
CATALOG = WIKI / "catalog.md"
GLOSSARY = WIKI / "glossary.md"
SOURCE_MAP = WIKI / "source-map.md"
CONTRADICTIONS = WIKI / "contradictions.md"
SOURCES = RAW / "SOURCES.md"
RETIRED_OPERATION_OUTPUT = WIKI / "log.md"

REQUIRED_LEAF_SECTIONS = ["Summary", "Detail", "Sources", "See Also"]
SOURCE_MAP_COLUMNS = [
    "Source or aggregate",
    "Authority tier",
    "Status",
    "Ingest",
    "Dependent pages",
    "Covered subpaths",
    "Known gaps",
    "Last checked",
    "Next action",
]
SOURCE_AUTHORITY_TIERS = {
    "primary-code",
    "primary-test",
    "generated-evidence",
    "primary-doc",
    "external-spec",
    "background",
}
ROUTER_NAMES = {"index.md"}
SUPPORT_PAGES = {
    WIKI / "SCHEMA.md",
    WIKI / "workflows.md",
    WIKI / "glossary.md",
    WIKI / "catalog.md",
    WIKI / "source-map.md",
    WIKI / "contradictions.md",
    WIKI / "lint" / "CHECKS.md",
}
EXEMPT_INDEX_DIRS = {WIKI / "lint"}
SKIP_LINK_SCHEMES = {"http", "https", "mailto", "tel"}
PAGE_STATUSES = {"confirmed", "provisional", "contested", "stale"}
SOURCE_STATUSES = {"frozen", "living"}
INGEST_STATES = {"registered", "partial", "ingested", "exact seed"}
GLOBAL_HUB_NAMES = {
    "SCHEMA.md",
    "workflows.md",
    "glossary.md",
    "catalog.md",
    "source-map.md",
    "contradictions.md",
}
COMMISSIONED_COORDINATION_RE = re.compile(
    r"^(?:plan[1-5][.]md|plan_synth[.]md|tasks[1-6][.]md)$"
)
# Removed source-tracking metadata bans: no hash digest values of any
# algorithm, no Mtime/Hash manifest columns, no Mtime values, and no
# YYYY-MM-DD date literals in governed KB files. Prohibition/supersession
# statements may still name the removed metadata.
HASH_DIGEST_VALUE_RE = re.compile(
    r"`?\b(?:sha-?\d+|md-?5|blake-?\d\w*|xxh\d*|crc-?\d+)\b`?"
    r"(?:\s*[:=]\s*|\s+)`?[0-9a-f]{8,}\b`?",
    re.IGNORECASE,
)
MTIME_HASH_COLUMN_RE = re.compile(r"\|\s*(Mtime|Hash)\s*\|")
MTIME_WORD_RE = re.compile(r"\bmtimes?\b", re.IGNORECASE)
# A mention is only allowed when a prohibition/supersession keyword precedes
# it with no sentence boundary in between (e.g. "no `mtime` columns",
# "supersedes ... mtime source-tracking instructions"). Positive instructions
# such as "record mtime values" must fail.
METADATA_PROHIBITION_RE = re.compile(
    r"\b(?:no|not|never|without|removed?|supersed\w*|prohibit\w*|banned|bans?)\b"
    r"[^.!?]*\bmtimes?\b",
    re.IGNORECASE,
)
ISO_DATE_RE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")
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


def _is_relative_to(path: Path, base: Path) -> bool:
    """
    Return whether ``path`` is under ``base``.

    :param path: Candidate path.
    :param base: Candidate parent path.
    :return: Whether ``path`` is equal to or below ``base``.
    """
    try:
        path.resolve().relative_to(base.resolve())
    except ValueError:
        return False
    return True


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


def _mask_nonprose(text: str) -> str:
    """
    Mask fenced code and HTML comments while preserving text offsets.

    :param text: Markdown text.
    :return: Text with non-prose characters replaced by spaces.
    """

    def mask_value(value: str) -> str:
        return "".join("\n" if ch == "\n" else " " for ch in value)

    def repl(match: Match[str]) -> str:
        return mask_value(match.group(0))

    pattern = (
        r"^[ ]{0,3}`{3,}[^\n]*\n.*?^[ ]{0,3}`{3,}[ \t]*$"
        r"|^[ ]{0,3}~{3,}[^\n]*\n.*?^[ ]{0,3}~{3,}[ \t]*$"
        r"|<!--.*?-->"
    )
    masked = re.sub(pattern, repl, text, flags=re.DOTALL | re.MULTILINE)
    result: List[str] = []
    list_contexts: List[Tuple[int, int]] = []
    list_marker_re = re.compile(r"^([ \t]*)(?:[-+*]|\d{1,9}[.)])([ \t]+)")
    for line in masked.splitlines(keepends=True):
        content = line.rstrip("\r\n")
        if not content.strip():
            result.append(line)
            continue
        leading = re.match(r"^[ \t]*", content)
        indent = len(leading.group(0).expandtabs(4)) if leading else 0
        list_marker = list_marker_re.match(content)
        if list_marker:
            marker_indent = len(list_marker.group(1).expandtabs(4))
            content_indent = len(list_marker.group(0).expandtabs(4))
            parent_list = next(
                (
                    context
                    for context in reversed(list_contexts)
                    if context[1] <= marker_indent < context[1] + 4
                ),
                None,
            )
            if marker_indent > 3 and parent_list is None:
                result.append(mask_value(line))
                continue
            while list_contexts and marker_indent < list_contexts[-1][1]:
                list_contexts.pop()
            list_contexts.append((marker_indent, content_indent))
            result.append(line)
            continue
        if indent < 4:
            while list_contexts and indent <= list_contexts[-1][0]:
                list_contexts.pop()
            result.append(line)
            continue
        containing_list = next(
            (context for context in reversed(list_contexts) if context[1] <= indent),
            None,
        )
        if containing_list is not None and indent < containing_list[1] + 4:
            result.append(line)
        else:
            result.append(mask_value(line))
    return "".join(result)


def _mask_inline_code(text: str) -> str:
    """
    Mask Markdown inline-code spans while preserving text offsets.

    :param text: Markdown text with block-level non-prose already masked.
    :return: Text with inline-code characters replaced by spaces.
    """

    def repl(match: Match[str]) -> str:
        return "".join("\n" if ch == "\n" else " " for ch in match.group(0))

    pattern = r"(?<!`)(?P<ticks>`+)(?!`)(?P<body>.*?)(?<!`)(?P=ticks)(?!`)"
    return re.sub(pattern, repl, text, flags=re.DOTALL)


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
    for line in _mask_nonprose(_read(path)).splitlines():
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
    Extract Markdown links outside fenced code blocks and HTML comments.

    :param path: Markdown file to scan.
    :return: Links with local resolution metadata.
    """
    text = _read(path)
    search_text = _mask_inline_code(_mask_nonprose(text))
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
        _is_relative_to(path, WIKI)
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
            if not _is_relative_to(target, ROOT):
                _fail(
                    failures,
                    path,
                    link.line,
                    f"link target escapes repository: {link.target}",
                )
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
        if path == AGENTS or _is_relative_to(path, WIKI):
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
        elif not any(line.startswith("> Up:") for line in header_lines):
            _fail(failures, path, 2, "leaf header lacks Up link")

        prose = _mask_nonprose(text)
        found = [title for level, title in _headings(prose) if level == 2]
        if found != REQUIRED_LEAF_SECTIONS:
            _fail(
                failures,
                path,
                None,
                "leaf H2 sections must be exactly Summary, Detail, Sources, See Also",
            )
        for section in ("Sources", "See Also"):
            if not _section_body(prose, section):
                _fail(failures, path, None, f"leaf has empty {section}")

        registered = _manifest_literals()
        source_links = []
        for link in _extract_links(path):
            if not link.is_source_section:
                continue
            parsed = urlparse(link.target)
            external_literal = link.target.split("#", 1)[0]
            external = parsed.scheme in {"http", "https"} and _source_registered(
                external_literal, registered
            )
            local_source = (
                link.target_path is not None
                and _is_relative_to(link.target_path, ROOT)
                and not _is_relative_to(link.target_path, WIKI)
                and _source_registered(_rel(link.target_path), registered)
            )
            if external or local_source:
                source_links.append(link)
        if not source_links:
            _fail(failures, path, None, "leaf Sources has no registered source link")

        see_also = _mask_inline_code(_section_body(prose, "See Also"))
        see_also_targets: Set[Path] = set()
        for match in re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", see_also):
            target, _ = _path_from_link(path, match.group(2))
            if target is not None and _is_relative_to(target, WIKI):
                see_also_targets.add(target)
        parent = (path.parent / "index.md").resolve()
        if parent not in see_also_targets:
            _fail(failures, path, None, "leaf See Also lacks parent router")
        if not (see_also_targets - {parent, path.resolve()}):
            _fail(failures, path, None, "leaf See Also lacks neighboring wiki page")


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
                target == AGENTS.resolve() or _is_relative_to(target, WIKI)
            ):
                queue.append(target)

    for page in _wiki_pages():
        if page not in reachable:
            _fail(failures, page, None, "wiki page not reachable from AGENTS.md")


def _legal_router_target(router: Path, target: Path) -> bool:
    """
    Return whether a local router edge obeys router scope.

    :param router: Root or ``index.md`` router.
    :param target: Resolved local target.
    :return: Whether the edge is legal.
    """
    if target == AGENTS or target in SUPPORT_PAGES or target == SOURCES:
        return True
    if router == AGENTS:
        return target.parent == WIKI or (
            target.name == "index.md" and target.parent.parent == WIKI
        )
    if target.parent == router.parent:
        return True
    if target.name == "index.md" and target.parent.parent == router.parent:
        return True
    parent_router = (
        AGENTS if router.parent.parent == WIKI else router.parent.parent / "index.md"
    )
    return target == parent_router


def _check_router_edges(files: List[Path], failures: List[str]) -> None:
    """
    Check routers link only to children, parent/root, and global hubs.

    :param files: Governed Markdown files.
    :param failures: Mutable failure list.
    """
    for router in files:
        if not _is_router(router):
            continue
        for link in _extract_links(router):
            target = link.target_path
            if target is None or not _is_relative_to(target, ROOT):
                continue
            if not _legal_router_target(router, target):
                _fail(
                    failures,
                    router,
                    link.line,
                    f"illegal router edge: {link.target}",
                )


def _check_source_registration(pages: List[Path], failures: List[str]) -> None:
    """
    Check literal local source links are registered in the source manifest.

    :param pages: Wiki pages to inspect.
    :param failures: Mutable failure list.
    """
    registered = _manifest_literals()
    for path in pages:
        if not _is_normal_leaf(path):
            continue
        for link in _extract_links(path):
            target = link.target_path
            if not link.is_source_section or target is None:
                continue
            if _is_relative_to(target, WIKI) or not _is_relative_to(target, ROOT):
                continue
            literal = _rel(target)
            if not _source_registered(literal, registered):
                _fail(
                    failures,
                    path,
                    link.line,
                    f"literal source is not registered: {literal}",
                )


def _check_catalog(failures: List[str]) -> None:
    """
    Check catalog target equality and page status/date agreement.

    :param failures: Mutable failure list.
    """
    if not CATALOG.exists():
        return
    header, located_rows = _table_with_header(CATALOG, "Page")
    try:
        page_i = header.index("Page")
        type_i = header.index("Type")
        status_i = header.index("Status")
        date_i = header.index("Last reconciled")
    except ValueError:
        _fail(
            failures,
            CATALOG,
            None,
            "catalog lacks Page/Type/Status/Last reconciled columns",
        )
        return
    links: List[Tuple[Path, int, List[str]]] = []
    for line, row in located_rows:
        if len(row) <= max(page_i, type_i, status_i, date_i):
            _fail(failures, CATALOG, line, "catalog row has too few columns")
            continue
        page_links = list(
            re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", row[page_i])
        )
        if len(page_links) != 1:
            _fail(failures, CATALOG, line, "catalog Page cell must contain one link")
            continue
        for link in page_links:
            target, _ = _path_from_link(CATALOG, link.group(2))
            if target is None or not _is_relative_to(target, WIKI):
                _fail(failures, CATALOG, line, "catalog Page target is outside wiki")
                continue
            links.append((target, line, row))
    counts = Counter(target for target, _, _ in links)
    for page in _wiki_pages():
        count = counts.get(page, 0)
        if count == 0:
            _fail(failures, CATALOG, None, f"missing wiki page: {_rel(page)}")
        elif count > 1:
            _fail(failures, CATALOG, None, f"duplicate wiki page: {_rel(page)}")
    live_pages = set(_wiki_pages())
    for target, line, row in links:
        if target not in live_pages:
            _fail(
                failures,
                CATALOG,
                line,
                f"catalog target is not a live wiki page: {_rel(target)}",
            )
            continue
        expected_status: Optional[str]
        expected_date: Optional[str]
        if _is_router(target):
            expected_status, expected_date = "router", "n/a"
            if row[type_i].strip().lower() != "router":
                _fail(failures, CATALOG, line, "router catalog Type must be router")
        else:
            expected_status, expected_date = _page_header_metadata(target)
        catalog_status = row[status_i].strip().lower()
        catalog_date = row[date_i].strip()
        if expected_status is None:
            _fail(failures, target, 2, "page header lacks allowed Status")
        elif catalog_status != expected_status:
            _fail(
                failures,
                CATALOG,
                line,
                f"catalog status {catalog_status!r} disagrees with "
                f"{_rel(target)} status {expected_status!r}",
            )
        if expected_date is None:
            _fail(failures, target, 2, "page header lacks Last reconciled/checked date")
        elif catalog_date.lower() != expected_date.lower():
            _fail(
                failures,
                CATALOG,
                line,
                f"catalog date {catalog_date!r} disagrees with "
                f"{_rel(target)} date {expected_date!r}",
            )


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


def _parse_tables_with_lines(
    path: Path,
) -> List[Tuple[List[str], List[Tuple[int, List[str]]]]]:
    """
    Parse simple Markdown tables and retain one-based row locations.

    :param path: Markdown file containing zero or more tables.
    :return: Tables as header plus ``(line, cells)`` rows.
    """
    if not path.exists():
        return []
    tables: List[Tuple[List[str], List[Tuple[int, List[str]]]]] = []
    header: Optional[List[str]] = None
    rows: List[Tuple[int, List[str]]] = []
    prose = _mask_nonprose(_read(path))
    for line_number, line in enumerate(prose.splitlines(), start=1):
        if not line.startswith("|"):
            if header is not None:
                tables.append((header, rows))
                header = None
                rows = []
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells):
            continue
        if header is None:
            header = cells
        else:
            rows.append((line_number, cells))
    if header is not None:
        tables.append((header, rows))
    return tables


def _table_with_header(
    path: Path, required_column: str
) -> Tuple[List[str], List[Tuple[int, List[str]]]]:
    """
    Return first table containing a required column.

    :param path: Markdown file to inspect.
    :param required_column: Header cell identifying the desired table.
    :return: Header and located rows, or empty values.
    """
    for header, rows in _parse_tables_with_lines(path):
        if required_column in header:
            return header, rows
    return [], []


def _literal_code_values(text: str) -> Set[str]:
    """
    Return inline-code values from Markdown text.

    :param text: Markdown text.
    :return: Unique values between single backticks.
    """
    return {match.group(1).strip() for match in re.finditer(r"`([^`\n]+)`", text)}


def _manifest_literals() -> Set[str]:
    """
    Return literals appearing in manifest Source cells only.

    :return: Registered source identifiers and paths.
    """
    registered: Set[str] = set()
    for header, rows in _parse_tables_with_lines(SOURCES):
        if "Source" not in header:
            continue
        source_i = header.index("Source")
        for _, row in rows:
            if len(row) > source_i:
                registered.update(_literal_code_values(row[source_i]))
    return registered


def _source_registered(literal: str, registered: Set[str]) -> bool:
    """
    Return whether a source path is exact-registered or aggregate-covered.

    Aggregate coverage follows the finite path contracts named in the source
    manifest. It does not infer coverage from prose or arbitrary prefixes.

    :param literal: Repository-relative source path.
    :param registered: Manifest literals.
    :return: Whether registration is unambiguous.
    """
    if literal in registered or literal == "raw/SOURCES.md":
        return True
    path = Path(literal)
    suffix = path.suffix.lower()
    rules = [
        (
            "core-top-level-package-modules",
            path.parent == Path("nrpy")
            and (suffix in {".py", ".txt"} or path.name == "py.typed"),
        ),
        (
            "helpers-package-modules",
            _is_relative_to(ROOT / path, ROOT / "nrpy/helpers")
            and suffix in {".py", ".h"},
        ),
        (
            "helpers-validation-and-reference-metric-tests",
            any(
                _is_relative_to(ROOT / path, ROOT / prefix)
                for prefix in (
                    "nrpy/helpers",
                    "nrpy/validate_expressions",
                    "nrpy/tests",
                )
            ),
        ),
        (
            "equation-modules-and-trusted-values",
            _is_relative_to(ROOT / path, ROOT / "nrpy/equations"),
        ),
        (
            "infrastructure-modules-and-embedded-headers",
            _is_relative_to(ROOT / path, ROOT / "nrpy/infrastructures")
            and (suffix in {".py", ".h"} or (ROOT / path).is_dir()),
        ),
        (
            "carpetx-package-inventory",
            _is_relative_to(ROOT / path, ROOT / "nrpy/infrastructures/CarpetX")
            and suffix == ".py",
        ),
        (
            "example-generators-and-companion-scripts",
            _is_relative_to(ROOT / path, ROOT / "nrpy/examples"),
        ),
        (
            "ci-and-local-automation",
            _is_relative_to(ROOT / path, ROOT / ".github"),
        ),
    ]
    return any(name in registered and matches for name, matches in rules)


def _valid_date(value: str, allow_na: bool = False) -> bool:
    """
    Validate retained KB date syntax.

    Calendar-range validation is deterministic and does not consult the clock.

    :param value: Candidate date.
    :param allow_na: Whether ``n/a`` or ``-`` is accepted.
    :return: Whether the value is syntactically valid.
    """
    value = value.strip()
    if allow_na and value.lower() in {"n/a", "-"}:
        return True
    match = re.fullmatch(r"(\d{2})-(\d{2})-(\d{4})", value)
    if not match:
        return False
    month, day, year = (int(part) for part in match.groups())
    try:
        calendar_date(year, month, day)
    except ValueError:
        return False
    return True


def _page_header_metadata(path: Path) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract status and reconciliation/check date from page header.

    :param path: Wiki page.
    :return: Lowercase status and date, when present.
    """
    header = "\n".join(_read(path).splitlines()[:5])
    status_match = re.search(
        r"\bStatus:\s*(confirmed|provisional|contested|stale)\b",
        header,
        re.IGNORECASE,
    )
    date_match = re.search(
        r"\bLast (?:reconciled|checked|audited):\s*(\d{2}-\d{2}-\d{4}|n/a|-)\b",
        header,
        re.IGNORECASE,
    )
    status = status_match.group(1).lower() if status_match else None
    date = date_match.group(1) if date_match else None
    return status, date


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


def _check_retired_operation_output(failures: List[str]) -> None:
    """
    Check that the retired KB maintenance-output page has not been recreated.

    :param failures: Mutable failure list.
    """
    if RETIRED_OPERATION_OUTPUT.exists():
        _fail(
            failures,
            RETIRED_OPERATION_OUTPUT,
            None,
            "retired KB maintenance-output page must not exist; use git "
            "history for durable operations",
        )


def _governed_kb_files() -> List[Path]:
    """
    Return KB files governed by source-tracking metadata and date checks.

    :return: Governed markdown files, including ``raw/source-docs/**/*.md``.
    """
    files = set(_iter_md_files())
    source_docs = RAW / "source-docs"
    if source_docs.exists():
        files.update(p.resolve() for p in source_docs.rglob("*.md"))
    return sorted(files)


def _metadata_mention_allowed(lines: List[str], line_index: int) -> bool:
    """
    Return whether a removed-metadata mention is prohibition/supersession text.

    :param lines: File lines being checked.
    :param line_index: Zero-based index of the line with the mention.
    :return: Whether nearby context permits the mention.
    """
    context = " ".join(lines[max(0, line_index - 2) : line_index + 1])
    return METADATA_PROHIBITION_RE.search(context) is not None


def _check_source_tracking_metadata(failures: List[str]) -> None:
    """
    Check governed KB files for removed source-tracking metadata and dates.

    :param failures: Mutable failure list.
    """
    for path in _governed_kb_files():
        lines = _read(path).splitlines()
        for idx, line in enumerate(lines, start=1):
            if HASH_DIGEST_VALUE_RE.search(line):
                _fail(failures, path, idx, "hash digest value found")
            if MTIME_HASH_COLUMN_RE.search(line):
                _fail(
                    failures, path, idx, "source-tracking Mtime/Hash table column found"
                )
            if MTIME_WORD_RE.search(line) and not _metadata_mention_allowed(
                lines, idx - 1
            ):
                _fail(failures, path, idx, "source-tracking Mtime metadata found")
            if ISO_DATE_RE.search(line):
                _fail(
                    failures,
                    path,
                    idx,
                    "YYYY-MM-DD date literal found; KB dates use MM-DD-YYYY",
                )


def _check_status_vocabularies(failures: List[str]) -> None:
    """
    Check page, manifest, and source-map status/date vocabularies.

    :param failures: Mutable failure list.
    """
    for page in _wiki_pages():
        if _is_router(page):
            continue
        status, date = _page_header_metadata(page)
        if status not in PAGE_STATUSES:
            _fail(failures, page, 2, "page header has missing or invalid Status")
        if date is None or not _valid_date(date):
            _fail(
                failures,
                page,
                2,
                "page header has missing or invalid MM-DD-YYYY date",
            )

    for path in (SOURCES, SOURCE_MAP):
        for header, rows in _parse_tables_with_lines(path):
            if "Status" not in header:
                continue
            status_i = header.index("Status")
            ingest_i = header.index("Ingest") if "Ingest" in header else None
            date_header = "Accessed" if "Accessed" in header else "Last checked"
            date_i = header.index(date_header) if date_header in header else None
            for line, row in rows:
                if len(row) <= status_i:
                    _fail(failures, path, line, "source row lacks Status value")
                    continue
                if row[status_i].strip().lower() not in SOURCE_STATUSES:
                    _fail(
                        failures,
                        path,
                        line,
                        f"invalid source Status: {row[status_i]!r}",
                    )
                if ingest_i is not None:
                    if (
                        len(row) <= ingest_i
                        or row[ingest_i].strip().lower() not in INGEST_STATES
                    ):
                        value = row[ingest_i] if len(row) > ingest_i else ""
                        _fail(failures, path, line, f"invalid Ingest state: {value!r}")
                if date_i is not None:
                    if len(row) <= date_i or not _valid_date(row[date_i]):
                        value = row[date_i] if len(row) > date_i else ""
                        _fail(
                            failures,
                            path,
                            line,
                            f"invalid {date_header} date: {value!r}",
                        )


def _check_source_map_targets(failures: List[str]) -> None:
    """
    Check source-map source literals are registered and local targets exist.

    :param failures: Mutable failure list.
    """
    source_tables = [
        (header, rows)
        for header, rows in _parse_tables_with_lines(SOURCE_MAP)
        if "Source or aggregate" in header
    ]
    if not source_tables:
        _fail(failures, SOURCE_MAP, None, "source map lacks source table")
        return
    registered = _manifest_literals()
    for header, rows in source_tables:
        if header != SOURCE_MAP_COLUMNS:
            _fail(failures, SOURCE_MAP, None, "source-map columns do not match schema")
            continue
        source_i = header.index("Source or aggregate")
        authority_i = header.index("Authority tier")
        for line, row in rows:
            if len(row) != len(header):
                _fail(
                    failures, SOURCE_MAP, line, "source-map row has wrong field count"
                )
                continue
            tiers = {tier.strip() for tier in row[authority_i].split(";")}
            if not tiers or not tiers.issubset(SOURCE_AUTHORITY_TIERS):
                _fail(failures, SOURCE_MAP, line, "invalid source authority tier")
            literals = _literal_code_values(row[source_i])
            if len(literals) != 1:
                _fail(
                    failures,
                    SOURCE_MAP,
                    line,
                    "source-map source cell must contain one literal",
                )
                continue
            literal = next(iter(literals))
            if not _source_registered(literal, registered):
                _fail(
                    failures,
                    SOURCE_MAP,
                    line,
                    f"source-map literal is not registered: {literal}",
                )
            parsed = urlparse(literal)
            looks_local = not parsed.scheme and (
                "/" in literal or literal.startswith(".") or "." in Path(literal).name
            )
            if looks_local and not any(char in literal for char in "*?[]"):
                target = (ROOT / literal).resolve()
                if not _is_relative_to(target, ROOT) or not target.exists():
                    _fail(
                        failures,
                        SOURCE_MAP,
                        line,
                        f"source-map local target missing: {literal}",
                    )


def _affected_page_links(cell: str) -> List[str]:
    """
    Return repository-relative affected-page paths from a table cell.

    :param cell: Contradiction Affected pages cell.
    :return: Resolved wiki paths as strings.
    """
    paths: List[str] = []
    for match in re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", cell):
        target, _ = _path_from_link(CONTRADICTIONS, match.group(2))
        if target is not None and _is_relative_to(target, WIKI):
            paths.append(_rel(target))
    return paths


def _check_contradictions(failures: List[str]) -> None:
    """
    Check structured contradiction rows and active marker bidirectionality.

    :param failures: Mutable failure list.
    """
    required = [
        "ID",
        "Claim",
        "Claim status",
        "Source A",
        "Source B",
        "Authority decision",
        "Affected pages",
        "Page-status rationale",
        "Owner/trigger",
        "Resolution test",
        "Opened",
        "Resolved",
        "Notes",
    ]
    header, rows = _table_with_header(CONTRADICTIONS, "ID")
    if header != required:
        _fail(
            failures,
            CONTRADICTIONS,
            None,
            "contradiction table columns do not match schema",
        )
        return
    indexes = {name: header.index(name) for name in required}
    seen: Set[str] = set()
    active: Dict[str, Tuple[str, Set[str]]] = {}
    for line, row in rows:
        if len(row) != len(header):
            _fail(
                failures,
                CONTRADICTIONS,
                line,
                "contradiction row has wrong field count",
            )
            continue
        identifier = row[indexes["ID"]].strip("` ")
        if not re.fullmatch(r"CONTR-\d{4}", identifier):
            _fail(
                failures,
                CONTRADICTIONS,
                line,
                f"malformed contradiction ID: {identifier!r}",
            )
            continue
        if identifier in seen:
            _fail(
                failures,
                CONTRADICTIONS,
                line,
                f"duplicate contradiction ID: {identifier}",
            )
        seen.add(identifier)
        for field in required:
            if field == "Resolved":
                continue
            if not row[indexes[field]].strip():
                _fail(
                    failures,
                    CONTRADICTIONS,
                    line,
                    f"empty contradiction field: {field}",
                )
        status = row[indexes["Claim status"]].strip().lower()
        if status not in {"contested", "stale", "resolved"}:
            _fail(failures, CONTRADICTIONS, line, f"invalid Claim status: {status!r}")
        opened = row[indexes["Opened"]].strip()
        resolved = row[indexes["Resolved"]].strip()
        if not _valid_date(opened):
            _fail(failures, CONTRADICTIONS, line, f"invalid Opened date: {opened!r}")
        if status == "resolved":
            if not _valid_date(resolved):
                _fail(
                    failures, CONTRADICTIONS, line, "resolved row needs resolution date"
                )
        elif resolved not in {"-", "n/a"}:
            _fail(
                failures,
                CONTRADICTIONS,
                line,
                "active row Resolved must be '-' or 'n/a'",
            )
        affected = set(_affected_page_links(row[indexes["Affected pages"]]))
        if not affected:
            _fail(failures, CONTRADICTIONS, line, "Affected pages has no wiki links")
        if status in {"contested", "stale"}:
            active[identifier] = (status, affected)

    marker_re = re.compile(
        r"Claim status:\s*(contested|stale);\s*contradiction:\s*(CONTR-\d{4})[.]",
        re.IGNORECASE,
    )
    found_markers: Dict[str, Set[str]] = defaultdict(set)
    for page in _wiki_pages():
        if not _is_normal_leaf(page):
            continue
        text = _mask_inline_code(_mask_nonprose(_read(page)))
        for match in marker_re.finditer(text):
            status, identifier = match.group(1).lower(), match.group(2)
            found_markers[identifier].add(_rel(page))
            backlink = any(
                link.target_path == CONTRADICTIONS.resolve()
                and _slugify_heading(link.anchor) == identifier.lower()
                for link in _extract_links(page)
            )
            if not backlink:
                _fail(
                    failures,
                    page,
                    _line_for_offset(text, match.start()),
                    f"marker lacks backlink for {identifier}",
                )
            if identifier not in active:
                _fail(
                    failures,
                    page,
                    _line_for_offset(text, match.start()),
                    f"marker references missing or resolved contradiction: {identifier}",
                )
            elif active[identifier][0] != status:
                _fail(
                    failures,
                    page,
                    _line_for_offset(text, match.start()),
                    f"marker status disagrees with {identifier}",
                )
    for identifier, (_, affected) in active.items():
        markers = found_markers.get(identifier, set())
        for page_path in sorted(affected - markers):
            _fail(
                failures,
                ROOT / page_path,
                None,
                f"missing active marker for {identifier}",
            )
        for page_path in sorted(markers - affected):
            _fail(
                failures,
                ROOT / page_path,
                None,
                f"marker page absent from Affected pages for {identifier}",
            )


def _check_excluded_artifacts(failures: List[str]) -> None:
    """
    Check KB trees contain Markdown prose only and no retired artifacts.

    :param failures: Mutable failure list.
    """
    for base in (WIKI, RAW):
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix.lower() != ".md":
                _fail(failures, path, None, "non-Markdown artifact in KB tree")
            if path.is_file() and re.search(
                r"(?:^|[-_])(audit|scan|report|latest|token|plan|rank|log)(?:[-_.]|$)",
                path.name,
                re.IGNORECASE,
            ):
                _fail(failures, path, None, "excluded maintenance/planning artifact")


def _check_root_artifacts(failures: List[str]) -> None:
    """
    Check root maintenance artifacts.

    :param failures: Mutable failure list.
    """
    allowed_root_docs = {"AGENTS.md", "README.md", "CITATION.md", "coding_style.md"}
    for path in ROOT.glob("*.md"):
        if path.name in allowed_root_docs or COMMISSIONED_COORDINATION_RE.fullmatch(
            path.name
        ):
            continue
        if (
            re.fullmatch(r".*kb.*instructions.*[.]md", path.name, re.IGNORECASE)
            or ROOT_KB_ARTIFACT_RE.match(path.name)
            or re.search(
                r"(?:^|[-_])(audit|scan|report|latest|token|plan|rank|tasks?)(?:[-_.]|$)",
                path.name,
                re.IGNORECASE,
            )
        ):
            _fail(failures, path, None, "root-level KB maintenance artifact")


def _main() -> int:
    """
    Run KB lint checks.

    :return: Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all",
        action="store_true",
        help="compatibility alias; identical to default full deterministic coverage",
    )
    parser.parse_args()

    failures: List[str] = []
    files = _iter_md_files()
    pages = _wiki_pages()

    _check_wikilinks(files, failures)
    graph = _check_links(files, failures)
    _check_router_detail(files, failures)
    _check_router_edges(files, failures)
    _check_leaf_contract(pages, failures)
    _check_source_registration(pages, failures)
    _check_reachability(graph, failures)
    _check_catalog(failures)
    _check_content_dir_indexes(failures)
    _check_glossary_catalog_signal(failures)
    _check_status_vocabularies(failures)
    _check_source_map_targets(failures)
    _check_contradictions(failures)
    _check_excluded_artifacts(failures)
    _check_retired_operation_output(failures)
    _check_source_tracking_metadata(failures)
    _check_root_artifacts(failures)

    if failures:
        print("KB lint failed:")
        for item in sorted(failures):
            print(f"- {item}")
        return 1

    print("KB lint passed.")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
