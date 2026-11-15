#!/usr/bin/env python3
"""
Mechanical lint checks for the markdown knowledge base.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import argparse
import os

# Regexes are needed for Markdown link, heading, and table checks.
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
EXAMPLE_GENERATOR_CATALOG = WIKI / "examples" / "example-generator-catalog.md"
SOURCE_MAP = WIKI / "source-map.md"
CONTRADICTIONS = WIKI / "contradictions.md"
SOURCES = RAW / "SOURCES.md"
RETIRED_OPERATION_OUTPUT = WIKI / "log.md"

REQUIRED_LEAF_SECTIONS = ["Summary", "Detail", "Sources", "See Also"]
SOURCE_MAP_COLUMNS = [
    "Source or aggregate",
    "Glob pattern",
    "Authority tier",
    "Status",
    "Ingest",
    "Uncovered scope",
]
AGGREGATE_SOURCE_NAMES = {
    "core-top-level-package-modules",
    "helpers-validation-and-reference-metric-tests",
    "equation-modules-and-trusted-values",
    "infrastructure-modules-and-embedded-headers",
    "example-generators-and-companion-scripts",
    "ci-and-local-automation",
}
SOURCE_AUTHORITY_TIERS = {
    "primary-code",
    "primary-test",
    "generated-evidence",
    "primary-doc",
    "external-spec",
    "background",
}
REGISTERED_URL_PREFIXES = (
    "https://charm.readthedocs.io/en/v8.0.0/",
    "https://github.com/charmplusplus/charm/blob/v8.0.0/",
    "https://github.com/paralab/Dendro-5.01",
    "https://github.com/paralab/Dendro-GR",
    "https://github.com/zachetienne/nrpytutorial",
)
FILENAME_LABEL_RE = re.compile(
    r"(?:[A-Za-z0-9_./*-]+\.(?:md|py|c|h|cpp|hpp|cu|cuh|cxx|sh|yml|yaml|toml|ccl|ci|par|asc|txt|cfg|ini)"
    r"|(?:[A-Za-z0-9_./*-]+/)?(?:Makefile|Dockerfile|ThornList|py\.typed|\.pylintrc[\w.-]*|\.mypy\.ini|\.pydocstyle|\.darglint))"
)
ROUTER_NAMES = {"index.md"}
SUPPORT_PAGES = {
    WIKI / "SCHEMA.md",
    WIKI / "computational-physics-language.md",
    WIKI / "workflows.md",
    WIKI / "glossary.md",
    WIKI / "catalog.md",
    WIKI / "source-map.md",
    WIKI / "contradictions.md",
}
EXEMPT_INDEX_DIRS: Set[Path] = set()
SKIP_LINK_SCHEMES = {"http", "https", "mailto", "tel"}
PAGE_STATUSES = {"confirmed", "provisional", "contested", "stale"}
SOURCE_STATUSES = {"frozen", "frozen-historical", "living"}
INGEST_STATES = {"registered", "partial", "ingested"}

MAX_AGENTS_BYTES = 4_800
MAX_LEAF_BYTES = 20_000
MAX_UNSECTIONED_DETAIL_BYTES = 4_000
MAX_SUBSECTION_BYTES = 6_000
MIN_ROUTER_CHILDREN = 2
MAX_ROUTER_CHILDREN = 14
MAX_EVIDENCE_LINE_BYTES = 300
MAX_LINE_BYTES = 1_000
MAX_SEE_ALSO_LINKS = 6
MAX_TABLE_ROW_BYTES = 450
MAX_CATALOG_ROW_BYTES = 240
MAX_CONTRADICTION_ROW_BYTES = 220
MAX_EXTERNAL_NOTES_BYTES = 80
MAX_GLOSSARY_MEANING_BYTES = 180
RETIRED_AGGREGATE_NAMES = (
    "helpers-package-modules",
    "carpetx-package-inventory",
)
EVIDENCE_FILE_SUFFIX_RE = re.compile(
    r"(?:\.(?:py|c|h|cpp|hpp|cu|cuh|cxx|sh|yml|yaml|toml|ccl|ci|par|asc|md|txt|cfg|ini)$"
    r"|(?:^|/)(?:Makefile|Dockerfile|ThornList|py\.typed|\.pylintrc[\w.-]*|\.mypy\.ini|\.pydocstyle|\.darglint)$)"
)
ALLOWED_SEE_ALSO_PREFIXES = (
    "Depends on:",
    "Implements:",
    "Validated by:",
    "Example:",
    "Contrasts with:",
    "See also:",
)
ALLOWED_EVIDENCE_ROLES = {
    "descriptive behavior",
    "normative rule",
    "scientific/interface requirement",
    "CI behavior",
    "generated evidence",
}
BANNED_METAPHOR_RE = re.compile(
    r"\b(?:artifacts?|contracts?|schemas?|surfaces?|payloads?|plumbing|"
    r"pipelines?|bridges?|bundles?|gates?|harnesses?|fixtures?|consumers?|"
    r"owners?|ownership|owning|owned|owns|lifecycles?)\b",
    re.IGNORECASE,
)
MATH_IDIOM_EXEMPT_RE = re.compile(
    r"\b(?:apparent-horizon\s+surfaces?|horizon\s+surfaces?|adm\s+surfaces?|"
    r"two-surfaces?|2-surfaces?|coordinate\s+surfaces?|null\s+surfaces?|"
    r"caustic\s+surfaces?|spatial\s+hypersurfaces?|normal\s+bundles?|"
    r"tangent\s+bundles?|fiber\s+bundles?|tensor\s+contractions?|split-pipelines?|"
    r"memory\s+owners?(?:hip)?|posix\s+file\s+owners?(?:hip)?|file\s+owners?(?:hip)?|"
    r"distributed\s+subdomain\s+owners?(?:hip)?|subdomain\s+owners?(?:hip)?|"
    r"patch\s+owners?(?:hip)?)\b",
    re.IGNORECASE,
)

ROOT_COORDINATION_RE = re.compile(
    r"^(?:plan|tasks?)(?:[0-9]+)?(?:[-_][a-z0-9]+)*[.]md$",
    re.IGNORECASE,
)
ROOT_PROHIBITED_ARTIFACT_SEGMENT_RE = re.compile(
    r"(?:^|[._-])(?:audit|scan|report|log|latest|token|rank|ranking)s?[0-9]*"
    r"(?:$|[._-])",
    re.IGNORECASE,
)
# Volatile facts do not belong in the authored KB. Stable scientific
# cardinalities and identifier names remain valid; inventory counts and
# source-tracking values do not. Prohibition statements may name a removed
# maintenance-field category.
HASH_DIGEST_VALUE_RE = re.compile(
    r"`?\b(?:sha(?:-?3)?-?\d+|md-?5|blake-?\d\w*|xxh\d*|crc-?\d+)\b`?"
    r"(?:\s*[:=]\s*|\s+)`?[0-9a-f]{8,}\b`?"
    r"|\b(?:checksum|hash|digest)(?:\s+value)?\s*[:=]\s*`?\S+",
    re.IGNORECASE,
)
STORED_HEX_IDENTIFIER_RE = re.compile(
    r"(?<![0-9a-f])(?:[0-9a-f]{32}|[0-9a-f]{40}|[0-9a-f]{56}|"
    r"[0-9a-f]{64}|[0-9a-f]{96}|[0-9a-f]{128})(?![0-9a-f])",
    re.IGNORECASE,
)
STABLE_HEX_IDENTIFIER_PREFIX_RE = re.compile(
    r"\b(?:(?:api|interface|publication|release|source|format|schema)\s+)?"
    r"(?:identifier|id|version|label)\s*[:=]?\s*$",
    re.IGNORECASE,
)
VCS_HASH_VALUE_RE = re.compile(
    r"\b(?:commit|revision)\s+`?[0-9a-f]{7,40}\b`?|"
    r"https?://[^\s)]+/(?:commit|blob)/[0-9a-f]{7,40}(?:/|\b)",
    re.IGNORECASE,
)
URL_RE = re.compile(r"https?://[^\s)>`]+")
VOLATILE_DATE_RE = re.compile(
    r"(?<![\w./-])(?:\d{4}(?P<ymd_sep>[-/.])\d{2}(?P=ymd_sep)\d{2}|"
    r"\d{2}(?P<mdy_sep>[-/.])\d{2}(?P=mdy_sep)\d{4})(?![\w./-])|"
    r"(?<![\w-])(?:Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|"
    r"Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|"
    r"Nov(?:ember)?|Dec(?:ember)?)\.?\s+\d{1,2}(?:st|nd|rd|th)?[,]?\s+"
    r"\d{4}(?!\w)|"
    r"(?<![\w-])\d{1,2}(?:st|nd|rd|th)?\s+(?:Jan(?:uary)?|Feb(?:ruary)?|"
    r"Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|"
    r"Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)\.?[,]?\s+"
    r"\d{4}(?!\w)",
    re.IGNORECASE,
)
STABLE_DATE_IDENTIFIER_RE = re.compile(
    r"\b(?:api|interface|publication|release|source|format|schema)\s+"
    r"(?:identifier|id|version|label)\s*[:=]?\s*$",
    re.IGNORECASE,
)
VOLATILE_TIME_RE = re.compile(
    r"\b(?:[01]?\d|2[0-3]):[0-5]\d"
    + r"(?::[0-5]\d(?:[.]\d+)?)?(?:Z|[+-]\d{2}:?\d{2})?\b",
    re.IGNORECASE,
)
COUNT_VALUE_PATTERN = (
    r"(?:\d[\d,]*|zero|one|two|three|four|five|six|seven|eight|nine|ten|"
    r"eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|"
    r"nineteen|twenty)"
)
VOLATILE_INVENTORY_COUNT_RE = re.compile(
    r"(?<![\w.-])\d[\d,]*\s+(?:source\s+)?" + r"(?:files?|pages?|sources?|jobs?)\b",
    re.IGNORECASE,
)
VOLATILE_EXISTENTIAL_INVENTORY_RE = re.compile(
    rf"\bthere\s+(?:are|is)\s+{COUNT_VALUE_PATTERN}\s+(?:source\s+)?"
    r"(?:files?|pages?|sources?|jobs?)\b",
    re.IGNORECASE,
)
VOLATILE_CATALOG_COUNT_RE = re.compile(
    r"\b(?:inventory|catalog|aggregate)\b[^.!?\n|]{0,80}"
    rf"\b{COUNT_VALUE_PATTERN}(?:\s+[A-Za-z][\w-]*){{0,3}}\s+"
    r"(?:files?|pages?|sources?|generators?|jobs?|cases?|"
    r"projects?|applications?|builders?|dictionaries?|baselines?|comparisons?|"
    r"scripts?|modules?)\b",
    re.IGNORECASE,
)
VOLATILE_FIELD_RE = re.compile(
    r"(?:(?i:\bLast\s+(?:reconciled|checked|audited|reviewed|verified|"
    r"validated|updated|modified)\s*[:=])|"
    r"(?i:\b(?:Timestamp|Datestamp|Date stamp|Datetime|Date|Accessed|Created|"
    r"Updated|Modified|Opened|Resolved|Resolution|Reviewed|Checked|Audited|"
    r"Audit|Reconciliation|Run date|Validation date|Mtime|Hash|Checksum|Digest|"
    r"Source count|File count|Page count|Job count|Environment|Run result)"
    r"(?:\s+(?:at|on|date))?\s*[:=])|"
    r"(?i:\|\s*(?:Timestamp|Datestamp|Date stamp|Datetime|Date|Mtime|Hash|"
    r"Checksum|Digest|Accessed|Created|Updated|Modified|Reviewed|Checked|"
    r"Audited|Audit|Reconciliation|Environment|"
    r"Run result|Run date|Validation date|Source count|File count|Page count|"
    r"Job count|Last reconciled|Last checked|Last audited|Last reviewed|"
    r"Last verified|Last validated|Last updated|Last modified)\s*\|)|"
    r"\|\s*(?:Opened|Resolved|Resolution)\s*\|)",
)
VOLATILE_RESULT_TUPLE_RE = re.compile(
    r"\b(?:inspected|generated|built|run|result_checked)\s*=\s*"
    r"(?:pass|fail|not-run|not-applicable)\b",
    re.IGNORECASE,
)
ENVIRONMENT_ASSIGNMENT_RE = re.compile(
    r"\b(?:environment|platform|compiler|runtime|python|os|host)\s*[:=]\s*"
    + r"[^;|,]+",
    re.IGNORECASE,
)
VOLATILE_RUN_PROSE_RE = re.compile(
    r"(?:\b(?:was|were|remained)\s+(?:not\s+)?"
    r"(?:done|run|rerun|re-run|executed|exercised|inspected|built|generated|"
    r"validated)\b|"
    r"\b(?:has|have)\s+(?:not\s+)?been\s+"
    r"(?:audited|checked|reviewed|re-reviewed|run|rerun|re-run|executed|"
    r"exercised|inspected|built|generated|validated)\b|"
    r"\bnot\s+(?:yet\s+)?(?:audited|checked|reviewed|re-reviewed|reconciled|"
    r"validated)\b|\bunexecuted\b|"
    r"\bisolated\s+local\s+(?:[\w/-]+\s+){0,3}(?:builds?|runs?|results?)\b|"
    r"\bvalidation\s+artifacts?\s+(?:was|were)\s+temporary\b|"
    r"^\s*(?:[-*]\s+)?(?:the\s+)?(?:tests?|checks?|audits?|builds?|runs?|"
    r"validations?)\s+(?:passed|failed|succeeded|completed)"
    r"(?:\s+(?:on\s+(?:Linux|macOS|Ubuntu|Windows)"
    r"(?:\s+[0-9][\w.+-]*)*|in\s+(?:CI|continuous\s+integration)))?"
    r"\s*[.!?]?\s*$)",
    re.IGNORECASE,
)
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
ROOT_KB_ARTIFACT_RE = re.compile(r"^kb_audit_.*[.]md$")


@dataclass(frozen=True)
class _Link:
    """Markdown link with its source location, local path, anchor, and section flag."""

    source: Path
    line: int
    label: str
    target: str
    target_path: Optional[Path]
    anchor: str
    is_source_section: bool


_READ_CACHE: Dict[Path, str] = {}
_PROSE_CACHE: Dict[Path, str] = {}
_LINKS_CACHE: Dict[Path, List[_Link]] = {}
_REPO_FILES_CACHE: List[str] = []


def _clear_file_caches() -> None:
    """Reset per-run Markdown text, masked-prose, and extracted-link caches."""
    _READ_CACHE.clear()
    _PROSE_CACHE.clear()
    _LINKS_CACHE.clear()
    _REPO_FILES_CACHE.clear()


def _rel(path: Path) -> str:
    """
    Return a repository-relative POSIX path.

    :param path: Absolute or relative path inside the repository.
    :return: POSIX path relative to the repository root.
    """
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return Path(os.path.abspath(path)).relative_to(ROOT).as_posix()


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
    Read and cache a UTF-8 text file.

    :param path: File to read.
    :return: File contents.
    """
    resolved = path.resolve()
    if resolved not in _READ_CACHE:
        _READ_CACHE[resolved] = resolved.read_text(encoding="utf-8")
    return _READ_CACHE[resolved]


def _read_prose(path: Path) -> str:
    """
    Return cached Markdown text with fenced code and HTML comments masked.

    :param path: File to read.
    :return: Masked prose text preserving character offsets.
    """
    resolved = path.resolve()
    if resolved not in _PROSE_CACHE:
        _PROSE_CACHE[resolved] = _mask_nonprose(_read(resolved))
    return _PROSE_CACHE[resolved]


def _line_for_offset(text: str, offset: int) -> int:
    """
    Return one-based line number for character offset.

    :param text: Text containing the offset.
    :param offset: Zero-based character offset.
    :return: One-based line number.
    """
    return text.count("\n", 0, offset) + 1


def _has_prose_double_blank_line(text: str) -> bool:
    r"""
    Return whether prose outside fenced code blocks and HTML comments has two blank lines.

    Doctests:
    >>> _has_prose_double_blank_line("Alpha.\n\n```python\nx = 1\n\n\ny = 2\n```\n\nBeta.\n")
    False
    >>> _has_prose_double_blank_line("Alpha.\n\n\nBeta.\n")
    True

    :param text: Raw Markdown text.
    :return: True when prose contains ``\n\n\n``.
    """
    pattern = (
        r"^[ ]{0,3}`{3,}[^\n]*\n.*?^[ ]{0,3}`{3,}[ \t]*$"
        r"|^[ ]{0,3}~{3,}[^\n]*\n.*?^[ ]{0,3}~{3,}[ \t]*$"
        r"|<!--.*?-->"
    )
    without_fences = re.sub(pattern, "X", text, flags=re.DOTALL | re.MULTILINE)
    normalized = re.sub(r"[ \t]+\n", "\n", without_fences)
    return "\n\n\n" in normalized


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
    for line in _read_prose(path).splitlines():
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


def _mask_html_comments(text: str) -> str:
    """
    Mask HTML comments with spaces while preserving line breaks and offsets.

    :param text: Raw Markdown text.
    :return: Text with HTML comments replaced by spaces and newlines.
    """

    def repl(match: Match[str]) -> str:
        return "".join("\n" if ch == "\n" else " " for ch in match.group(0))

    return re.sub(r"<!--.*?-->", repl, text, flags=re.DOTALL)


def _raw_detail_and_chunks(raw: str, prose: str) -> Tuple[str, List[str]]:
    """
    Return unmasked Detail text (excluding HTML comments) and H3 chunks.

    :param raw: Unmasked file text.
    :param prose: Masked prose text with identical character offsets.
    :return: Tuple of Detail body and H3-partitioned chunks.
    """
    match = re.search(r"^## Detail\s*$", prose, re.MULTILINE)
    if not match:
        return "", []
    start = match.end()
    next_match = re.search(r"^##\s+", prose[start:], re.MULTILINE)
    end = start + next_match.start() if next_match else len(prose)
    uncommented = _mask_html_comments(raw)[start:end]
    prose_slice = prose[start:end]
    cuts = [0] + [m.start() for m in re.finditer(r"^###\s+", prose_slice, re.MULTILINE)]
    cuts.append(len(prose_slice))
    chunks = [
        uncommented[cuts[i] : cuts[i + 1]]
        for i in range(len(cuts) - 1)
        if uncommented[cuts[i] : cuts[i + 1]].strip()
    ]
    return uncommented.strip("\n"), chunks


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
    Extract and cache Markdown links outside fenced code blocks and HTML comments.

    :param path: Markdown file to scan.
    :return: Links with source locations, resolved paths, anchors, and section flags.
    """
    resolved = path.resolve()
    if resolved in _LINKS_CACHE:
        return _LINKS_CACHE[resolved]
    search_text = _mask_inline_code(_read_prose(resolved))
    pattern = re.compile(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)")
    links: List[_Link] = []
    for match in pattern.finditer(search_text):
        label, raw_target = match.groups()
        target = raw_target.strip()
        path_part, anchor = _path_from_link(resolved, target)
        links.append(
            _Link(
                source=resolved,
                line=_line_for_offset(search_text, match.start()),
                label=label.strip(),
                target=target,
                target_path=path_part,
                anchor=anchor,
                is_source_section=_in_sources_section(search_text, match.start()),
            )
        )
    _LINKS_CACHE[resolved] = links
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
    Validate Markdown links, reject multiline links, and return local link graph.

    :param files: Markdown files to scan.
    :param failures: Mutable failure list.
    :return: Directed graph of local links.
    """
    graph: Dict[Path, List[Path]] = defaultdict(list)
    anchor_cache: Dict[Path, Set[str]] = {}
    multiline_re = re.compile(r"(?<!!)\[([^\]\[]+)\]\(([^)]+)\)")
    for path in files:
        search_text = _mask_inline_code(_read_prose(path))
        for match in multiline_re.finditer(search_text):
            if "\n" in match.group(1) or "\n" in match.group(2):
                _fail(
                    failures,
                    path,
                    _line_for_offset(search_text, match.start()),
                    "multiline Markdown link must stay on one line",
                )
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
    Check router pages forbid Detail/Sources sections and carry a routing table.

    :param files: Markdown files to scan.
    :param failures: Mutable failure list.
    """
    for path in files:
        if (path == AGENTS or _is_relative_to(path, WIKI)) and _is_router(path):
            text = _read_prose(path)
            if re.search(r"^## Detail\s*$", text, re.MULTILINE):
                _fail(failures, path, None, "router has Detail section")
            if re.search(r"^## Sources\s*$", text, re.MULTILINE):
                _fail(failures, path, None, "router has Sources section")
            expected_table = (
                r"^\|\s*Task\s*\|\s*Read first\s*\|"
                if path == AGENTS
                else r"^\|\s*Page\s*\|\s*Go here when(?:\.\.\.|\u2026)?\s*\|"
            )
            if not re.search(expected_table, text, re.MULTILINE):
                _fail(
                    failures, path, None, "router lacks required routing table header"
                )


def _check_leaf_contract(
    pages: List[Path],
    registered: Set[str],
    source_statuses: Dict[str, str],
    failures: List[str],
) -> None:
    """
    Check normal leaves satisfy the required page layout and source rules.

    :param pages: Wiki pages to scan.
    :param registered: Registered source literals from ``raw/SOURCES.md``.
    :param source_statuses: Mapping from registered literal to source status.
    :param failures: Mutable failure list.
    """
    for path in pages:
        text = _read(path)
        lines = text.splitlines()
        page_status = _page_header_status(path)
        if page_status == "confirmed":
            for link in _extract_links(path):
                cite_literal = (
                    _rel(link.target_path)
                    if link.target_path is not None
                    and _is_relative_to(link.target_path, ROOT)
                    and not _is_relative_to(link.target_path, WIKI)
                    else (
                        link.target.split("#", 1)[0]
                        if urlparse(link.target).scheme in {"http", "https"}
                        else None
                    )
                )
                if (
                    cite_literal is not None
                    and source_statuses.get(cite_literal) == "frozen-historical"
                ):
                    source_line = (
                        lines[link.line - 1] if link.line <= len(lines) else ""
                    )
                    visible_line = re.sub(
                        r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)",
                        lambda m: (
                            " "
                            if FILENAME_LABEL_RE.fullmatch(m.group(1).strip())
                            else m.group(1)
                        ),
                        source_line,
                    )
                    if "historical" not in visible_line.lower():
                        _fail(
                            failures,
                            path,
                            link.line,
                            f"confirmed page cites frozen-historical source without "
                            f"explicit historical note: {cite_literal}",
                        )
        if not _is_normal_leaf(path):
            continue
        if not lines or not lines[0].startswith("# "):
            _fail(failures, path, 1, "leaf missing H1 title")
        expected_parent = (
            AGENTS.resolve()
            if path.parent == WIKI
            else (path.parent / "index.md").resolve()
        )
        header_lines = [line for line in lines[1:4] if line.startswith("> ")]
        if len(header_lines) < 2:
            _fail(failures, path, 2, "leaf missing two-line blockquote header")
        else:
            up_lines = [line for line in header_lines if line.startswith("> Up:")]
            if not up_lines:
                _fail(failures, path, 2, "leaf header lacks Up link")
            else:
                up_match = re.search(
                    r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)",
                    _mask_inline_code(up_lines[0]),
                )
                up_target = (
                    _path_from_link(path, up_match.group(2))[0] if up_match else None
                )
                if up_target != expected_parent:
                    _fail(
                        failures,
                        path,
                        3,
                        f"leaf Up link must resolve to {_rel(expected_parent)}",
                    )

        prose = _read_prose(path)
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

        source_links = []
        for link in _extract_links(path):
            if not link.is_source_section:
                continue
            parsed = urlparse(link.target)
            external_literal = link.target.split("#", 1)[0]
            external = parsed.scheme in {"http", "https"} and _source_registered(
                external_literal, registered
            )
            local_literal = (
                _rel(link.target_path)
                if link.target_path is not None
                and _is_relative_to(link.target_path, ROOT)
                and not _is_relative_to(link.target_path, WIKI)
                else None
            )
            if local_literal == "raw/SOURCES.md":
                _fail(
                    failures,
                    path,
                    link.line,
                    "normal leaf Sources must cite concrete sources, not raw/SOURCES.md",
                )
            local_source = local_literal is not None and _source_registered(
                local_literal, registered
            )
            if external or local_source:
                source_links.append(link)
        if not source_links:
            _fail(failures, path, None, "leaf Sources has no registered source link")
        for src_ln in _section_body(prose, "Sources").splitlines():
            if re.search(r"^\s*-\s*[,;:]|-\s*,\s*\[", src_ln):
                _fail(
                    failures,
                    path,
                    None,
                    f"malformed Sources bullet: {src_ln.strip()[:50]}",
                )

        see_also = _mask_inline_code(_section_body(prose, "See Also"))
        see_also_targets: Set[Path] = set()
        for match in re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", see_also):
            target, _ = _path_from_link(path, match.group(2))
            if target is not None and _is_relative_to(target, WIKI):
                see_also_targets.add(target)
        if not (see_also_targets - {expected_parent, path.resolve()}):
            _fail(failures, path, None, "leaf See Also lacks neighboring wiki page")


def _repo_source_files() -> List[str]:
    """
    Return cached repository-relative file paths outside KB, build, and scratch directories.

    :return: Sorted POSIX paths for repository source resolution.
    """
    if not _REPO_FILES_CACHE:
        pruned_dirs = {
            ".git",
            "wiki",
            "raw",
            "tmp",
            "project",
            "build",
            "dist",
            "__pycache__",
            ".mypy_cache",
            ".pytest_cache",
        }
        for dirpath, dirnames, filenames in os.walk(ROOT):
            dirnames[:] = [d for d in dirnames if d not in pruned_dirs]
            base = Path(dirpath)
            for fname in filenames:
                candidate = base / fname
                if candidate.is_file():
                    _REPO_FILES_CACHE.append(_rel(candidate))
        _REPO_FILES_CACHE.sort()
    return _REPO_FILES_CACHE


def _extract_evidence_tokens(
    span: str, is_evidence_line: bool
) -> List[Tuple[str, bool]]:
    """
    Extract ``(token, inside_parentheses)`` pairs from an Evidence line or ``(Ev: ...)`` span.

    Doctests:
    >>> _extract_evidence_tokens(
    ...     "Evidence (descriptive behavior): `solver.py::run` (`Ctx::out`); `bad_py::fn`.",
    ...     True,
    ... )
    [('solver.py::run', False), ('Ctx::out', True), ('bad_py::fn', False)]

    :param span: Stripped Evidence line or ``(Ev: ...)`` substring.
    :param is_evidence_line: Whether ``span`` is a full ``Evidence (<role>):`` line.
    :return: Ordered list of ``(token, in_parens)`` tuples.
    """
    start = 0
    if is_evidence_line:
        role_hdr = re.match(r"^Evidence \([^)]+\):", span)
        if role_hdr is not None:
            start = role_hdr.end()
    elif span.startswith("(Ev:"):
        start = 4
    body = span[start:]
    tokens: List[Tuple[str, bool]] = []
    depth = 0
    idx = 0
    while idx < len(body):
        ch = body[idx]
        if ch == "`":
            end = body.find("`", idx + 1)
            if end == -1:
                break
            tokens.append((body[idx + 1 : end], depth > 0))
            idx = end + 1
            continue
        if ch == "(":
            depth += 1
        elif ch == ")" and depth > 0:
            depth -= 1
        idx += 1
    return tokens


def _verify_code_symbol(
    path: Path,
    line_no: int,
    sym_raw: str,
    code_rel: str,
    failures: List[str],
) -> bool:
    """
    Verify ``Class.method`` / ``Class::method`` identifier segments in a resolved code file.

    :param path: Wiki page path.
    :param line_no: One-based Evidence line number.
    :param sym_raw: Symbol expression after ``::`` (or full parenthesized ``Class::member``).
    :param code_rel: Repository-relative path to the resolved code source file.
    :param failures: Mutable failure list.
    :return: True when syntax is valid and all identifier segments exist in ``code_rel``.
    """
    sym_base = sym_raw.split("(", 1)[0].strip()
    if not re.fullmatch(
        r"[A-Za-z_][A-Za-z0-9_]*(?:(?:\.|::)[A-Za-z_][A-Za-z0-9_]*)*", sym_base
    ):
        _fail(
            failures,
            path,
            line_no,
            f"invalid Evidence symbol syntax {sym_raw!r} in {code_rel}",
        )
        return False
    source_text = _read(ROOT / code_rel)
    missing = [seg for seg in re.split(r"::|\.", sym_base) if seg not in source_text]
    if missing:
        _fail(
            failures,
            path,
            line_no,
            f"Evidence symbol {missing[0]!r} not found in {code_rel}",
        )
        return False
    return True


def _verify_evidence_token(
    path: Path,
    line_no: int,
    tok_info: Tuple[str, bool],
    sources: Tuple[List[str], Set[str], Optional[str]],
    failures: List[str],
) -> Tuple[bool, Optional[str]]:
    """
    Verify one backticked Evidence token against source paths, URLs, and symbol identifiers.

    Doctests:
    >>> errs = []
    >>> ok, code_src = _verify_evidence_token(
    ...     WIKI / "core" / "c-codegen.md",
    ...     1,
    ...     ("c_codegen.py::c_codegen", False),
    ...     (["nrpy/c_codegen.py"], set(), None),
    ...     errs,
    ... )
    >>> (ok, code_src, errs)
    (True, 'nrpy/c_codegen.py', [])
    >>> errs = []
    >>> _verify_evidence_token(
    ...     WIKI / "core" / "c-codegen.md",
    ...     1,
    ...     ("c_codegen.py::MissingClass.c_codegen", False),
    ...     (["nrpy/c_codegen.py"], set(), None),
    ...     errs,
    ... )[0], bool(errs)
    (False, True)
    >>> errs = []
    >>> _verify_evidence_token(
    ...     WIKI / "infrastructures" / "dendro" / "constraints-and-diagnostic-norms.md",
    ...     1,
    ...     ("Ctx::write_grid_summary_data", True),
    ...     ([], set(), "nrpy/infrastructures/Dendro/solver_context.py"),
    ...     errs,
    ... ), errs
    ((False, 'nrpy/infrastructures/Dendro/solver_context.py'), [])
    >>> _verify_evidence_token(
    ...     WIKI / "infrastructures" / "dendro" / "constraints-and-diagnostic-norms.md",
    ...     1,
    ...     ("Ctx::write_grid_info", True),
    ...     ([], set(), "nrpy/infrastructures/Dendro/solver_context.py"),
    ...     errs,
    ... )[0], bool(errs)
    (False, True)
    >>> errs = []
    >>> _verify_evidence_token(
    ...     WIKI / "core" / "c-codegen.md",
    ...     1,
    ...     ("arXiv:0902.3652v2", False),
    ...     ([], {"https://arxiv.org/pdf/0902.3652v2"}, None),
    ...     errs,
    ... ), errs
    ((True, None), [])
    >>> _verify_evidence_token(
    ...     WIKI / "core" / "c-codegen.md",
    ...     1,
    ...     ("arXiv:9999.99999v1", False),
    ...     ([], {"https://arxiv.org/pdf/0902.3652v2"}, None),
    ...     errs,
    ... )[0], bool(errs)
    (False, True)

    :param path: Wiki page path.
    :param line_no: One-based Evidence line number.
    :param tok_info: Tuple of ``(token, in_parens)``.
    :param sources: Tuple of ``(local_sources, registered_sources, prev_code_file)``.
    :param failures: Mutable failure list.
    :return: Tuple of ``(resolved_citation, next_prev_code_file)``.
    """
    tok, in_parens = tok_info
    local_sources, registered, prev_code_file = sources
    parts = tok.split("::", 1)
    head = parts[0].strip()
    if head.startswith(("http://", "https://")):
        ext_url = head.split("#", 1)[0]
        ok = _source_registered(ext_url, registered)
        if not ok:
            _fail(
                failures,
                path,
                line_no,
                f"external URL is not registered: {ext_url}",
            )
        return ok, None
    if head.startswith("arXiv:"):
        aid = head.split(":", 1)[1]
        ok = re.fullmatch(
            r"(?:\d{4}\.\d{4,5}|[a-z-]+/\d{7})(?:v\d+)?", aid
        ) is not None and (
            f"https://arxiv.org/pdf/{aid}" in registered
            or f"https://arxiv.org/abs/{aid}" in registered
        )
        if not ok:
            _fail(
                failures,
                path,
                line_no,
                f"arXiv source is not registered: {head}",
            )
        return ok, None
    is_file_ref = not any(ch in head for ch in " *<{") and bool(
        EVIDENCE_FILE_SUFFIX_RE.search(head)
    )
    if not is_file_ref:
        if len(parts) == 2:
            if not in_parens or prev_code_file is None:
                _fail(
                    failures,
                    path,
                    line_no,
                    f"Evidence :: token must start with a source file path: {tok}",
                )
            else:
                _verify_code_symbol(path, line_no, tok, prev_code_file, failures)
        return False, prev_code_file
    matches = [src for src in local_sources if src == head or src.endswith("/" + head)]
    if len(matches) != 1:
        msg = (
            f"Evidence file reference missing from Sources: {head}"
            if not matches
            else f"ambiguous Evidence file reference {head!r} matches multiple Sources"
        )
        _fail(failures, path, line_no, msg)
        return False, None
    matched = matches[0]
    is_code = (
        matched.endswith((".py", ".c", ".h", ".cpp", ".hpp", ".cu", ".cuh", ".cxx"))
        and (ROOT / matched).is_file()
    )
    next_code = matched if is_code else None
    resolved = (
        _verify_code_symbol(path, line_no, parts[1].strip(), matched, failures)
        if (len(parts) == 2 and is_code)
        else True
    )
    return resolved, next_code


def _check_evidence_lines_and_sources(
    path: Path,
    prose: str,
    check_source_resolution: bool,
    registered: Set[str],
    failures: List[str],
) -> None:
    """
    Validate Evidence line roles/lengths, citation counts, and file/symbol resolution.

    :param path: Wiki page path.
    :param prose: Non-prose-masked page text.
    :param check_source_resolution: Whether to verify backticked file tokens against Sources.
    :param registered: Registered source literals from ``raw/SOURCES.md``.
    :param failures: Mutable failure list.
    """
    is_contradictions = path.resolve() == CONTRADICTIONS.resolve()
    has_sources = check_source_resolution or bool(_section_body(prose, "Sources"))
    local_sources = (
        _repo_source_files()
        if is_contradictions
        else (
            sorted(
                {
                    _rel(lk.target_path)
                    for lk in _extract_links(path)
                    if lk.is_source_section
                    and lk.target_path is not None
                    and _is_relative_to(lk.target_path, ROOT)
                    and not _is_relative_to(lk.target_path, WIKI)
                }
            )
            if has_sources
            else []
        )
    )
    prev_evidence_line: Optional[int] = None
    for line_no, ln in enumerate(prose.splitlines(), start=1):
        stripped = ln.strip()
        if not stripped:
            continue
        is_evidence_line = bool(re.match(r"^Evidence(?:\s*\(|:)", stripped))
        if stripped == "Claim evidence:":
            _fail(
                failures,
                path,
                line_no,
                "legacy multi-line Claim evidence block banned",
            )
        elif is_evidence_line:
            if prev_evidence_line is not None:
                _fail(
                    failures,
                    path,
                    line_no,
                    f"consecutive Evidence line follows line {prev_evidence_line}; merge tags",
                )
            prev_evidence_line = line_no
            if len(ln.encode("utf-8")) > MAX_EVIDENCE_LINE_BYTES:
                _fail(
                    failures,
                    path,
                    line_no,
                    f"Evidence line exceeds {MAX_EVIDENCE_LINE_BYTES} B",
                )
            role_match = re.match(r"^Evidence \(([^)]+)\):", stripped)
            if role_match is None or role_match.group(1) not in ALLOWED_EVIDENCE_ROLES:
                _fail(
                    failures,
                    path,
                    line_no,
                    f"Evidence line lacks allowed role: {stripped[:60]}",
                )
            if any(
                not urlparse(m.group(2).strip()).scheme
                for m in re.finditer(
                    r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", _mask_inline_code(stripped)
                )
            ):
                _fail(
                    failures,
                    path,
                    line_no,
                    "Evidence line must use backticked tags, not local Markdown links",
                )
        else:
            prev_evidence_line = None
        if has_sources or is_contradictions:
            spans = (
                [stripped]
                if is_evidence_line
                else [
                    stripped[m.start() : m.end()]
                    for m in re.finditer(r"\(Ev:[^)\n]+\)", _mask_inline_code(stripped))
                ]
            )
            for span in spans:
                resolved_count = sum(
                    1
                    for m in re.finditer(
                        r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)",
                        _mask_inline_code(span),
                    )
                    if urlparse(m.group(2).strip()).scheme in {"http", "https"}
                    and _source_registered(
                        m.group(2).strip().split("#", 1)[0], registered
                    )
                )
                prev_code_file: Optional[str] = None
                for tok_info in _extract_evidence_tokens(span, is_evidence_line):
                    resolved, prev_code_file = _verify_evidence_token(
                        path,
                        line_no,
                        tok_info,
                        (local_sources, registered, prev_code_file),
                        failures,
                    )
                    if resolved:
                        resolved_count += 1
                if resolved_count == 0:
                    _fail(
                        failures,
                        path,
                        line_no,
                        "Evidence span resolves zero source citations",
                    )


def _check_token_efficiency_invariants(
    pages: List[Path], registered: Set[str], failures: List[str]
) -> None:
    """
    Enforce byte caps, router fan-out, H3 chunking, Evidence/Sources bounds, and row widths.

    :param pages: Wiki pages to scan.
    :param registered: Registered source literals from ``raw/SOURCES.md``.
    :param failures: Mutable failure list.
    """
    if AGENTS.exists():
        ag_b = len(_read(AGENTS).encode("utf-8"))
        if ag_b > MAX_AGENTS_BYTES:
            _fail(
                failures,
                AGENTS,
                None,
                f"AGENTS.md is {ag_b} B (> {MAX_AGENTS_BYTES} B)",
            )

    for gov_path in dict.fromkeys(
        [AGENTS.resolve(), SOURCES.resolve()] + [p.resolve() for p in pages]
    ):
        if not gov_path.exists():
            continue
        for line_no, ln in enumerate(_read(gov_path).splitlines(), start=1):
            for retired_agg in RETIRED_AGGREGATE_NAMES:
                if retired_agg in ln:
                    _fail(
                        failures,
                        gov_path,
                        line_no,
                        f"retired aggregate reference {retired_agg!r}",
                    )

    for gov_doc in (AGENTS, SOURCES):
        if gov_doc.exists() and _has_prose_double_blank_line(_read(gov_doc)):
            _fail(failures, gov_doc, None, "contains double blank line (\\n\\n\\n)")

    for path in pages:
        raw = _read(path)
        prose = _read_prose(path)
        byte_len = len(raw.encode("utf-8"))
        if _has_prose_double_blank_line(raw):
            _fail(failures, path, None, "contains double blank line (\\n\\n\\n)")
        for line_no, raw_ln in enumerate(raw.splitlines(), start=1):
            if len(raw_ln.encode("utf-8")) > MAX_LINE_BYTES:
                _fail(
                    failures,
                    path,
                    line_no,
                    f"line exceeds {MAX_LINE_BYTES} B",
                )
        if _is_router(path) and path != AGENTS:
            children = {
                lk.target_path
                for lk in _extract_links(path)
                if lk.target_path is not None
                and lk.target_path != path
                and (
                    lk.target_path.parent == path.parent
                    or (
                        lk.target_path.name == "index.md"
                        and lk.target_path.parent.parent == path.parent
                    )
                )
            }
            if len(children) < MIN_ROUTER_CHILDREN:
                _fail(
                    failures,
                    path,
                    None,
                    f"router has {len(children)} child (< {MIN_ROUTER_CHILDREN})",
                )
            elif len(children) > MAX_ROUTER_CHILDREN:
                _fail(
                    failures,
                    path,
                    None,
                    f"router has {len(children)} children (> {MAX_ROUTER_CHILDREN})",
                )
            continue

        is_leaf = _is_normal_leaf(path)
        _check_evidence_lines_and_sources(path, prose, is_leaf, registered, failures)
        if not is_leaf:
            continue
        if byte_len > MAX_LEAF_BYTES:
            _fail(
                failures,
                path,
                None,
                f"leaf is {byte_len} B (> {MAX_LEAF_BYTES} B); split topic",
            )

        detail, chunks = _raw_detail_and_chunks(raw, prose)
        detail_b = len(detail.encode("utf-8"))
        h3_count = sum(1 for chunk in chunks if chunk.lstrip().startswith("### "))
        if detail_b > MAX_UNSECTIONED_DETAIL_BYTES and h3_count == 0:
            _fail(
                failures,
                path,
                None,
                f"Detail is {detail_b} B (> {MAX_UNSECTIONED_DETAIL_BYTES} B) with 0 H3s",
            )
        for idx, chunk in enumerate(chunks):
            chunk_b = len(chunk.encode("utf-8"))
            if chunk_b > MAX_SUBSECTION_BYTES:
                _fail(
                    failures,
                    path,
                    None,
                    f"Detail/H3 chunk {idx} is {chunk_b} B (> {MAX_SUBSECTION_BYTES} B)",
                )

        src_lines = [
            ln.strip()
            for ln in _section_body(prose, "Sources").splitlines()
            if ln.strip()
        ]
        for src_line in src_lines:
            if not src_line.startswith("- "):
                _fail(
                    failures,
                    path,
                    None,
                    f"Sources entry must use '- ' bullet: {src_line[:50]}",
                )
                continue
            src_links = re.findall(
                r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", _mask_inline_code(src_line[2:])
            )
            if len(src_links) != 1:
                _fail(
                    failures,
                    path,
                    None,
                    f"Sources bullet must contain exactly 1 link: {src_line[:50]}",
                )

        seen_sources: Set[str] = set()
        for lk in _extract_links(path):
            if not lk.is_source_section:
                continue
            key = (
                str(lk.target_path)
                if lk.target_path is not None
                else lk.target.split("#", 1)[0]
            )
            if key in seen_sources:
                _fail(failures, path, lk.line, f"duplicate Sources target: {lk.target}")
            seen_sources.add(key)

        expected_parent = (
            AGENTS.resolve()
            if path.parent == WIKI
            else (path.parent / "index.md").resolve()
        )
        sa_lines = [
            ln.strip()
            for ln in _section_body(prose, "See Also").splitlines()
            if ln.strip()
        ]
        for sa_line in sa_lines:
            if not sa_line.startswith("- "):
                _fail(
                    failures,
                    path,
                    None,
                    f"See Also entry must use '- ' bullet: {sa_line[:50]}",
                )
        sa_bullets = [ln[2:].strip() for ln in sa_lines if ln.startswith("- ")]
        if not 1 <= len(sa_bullets) <= MAX_SEE_ALSO_LINKS:
            _fail(
                failures,
                path,
                None,
                f"See Also has {len(sa_bullets)} bullets (expected 1..{MAX_SEE_ALSO_LINKS})",
            )
        seen_sa_targets: Set[Path] = set()
        for bullet in sa_bullets:
            if not bullet.startswith(ALLOWED_SEE_ALSO_PREFIXES):
                _fail(
                    failures,
                    path,
                    None,
                    f"See Also bullet lacks allowed prefix: {bullet[:50]}",
                )
            bullet_links = list(
                re.finditer(
                    r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", _mask_inline_code(bullet)
                )
            )
            if len(bullet_links) != 1:
                _fail(
                    failures,
                    path,
                    None,
                    f"See Also bullet must contain exactly 1 link: {bullet[:50]}",
                )
            for match in bullet_links:
                target_path, _ = _path_from_link(path, match.group(2))
                if target_path == expected_parent:
                    _fail(
                        failures,
                        path,
                        None,
                        "See Also duplicates header Up: parent router link",
                    )
                elif target_path == path.resolve():
                    _fail(failures, path, None, "See Also links leaf to itself")
                elif target_path is not None:
                    if target_path in seen_sa_targets:
                        _fail(
                            failures,
                            path,
                            None,
                            f"duplicate See Also target: {_rel(target_path)}",
                        )
                    seen_sa_targets.add(target_path)

    for table_file in (
        CATALOG,
        EXAMPLE_GENERATOR_CATALOG,
        GLOSSARY,
        SOURCE_MAP,
        CONTRADICTIONS,
        SOURCES,
    ):
        if not table_file.exists():
            continue
        row_cap = MAX_TABLE_ROW_BYTES
        if table_file in (CATALOG, EXAMPLE_GENERATOR_CATALOG):
            row_cap = MAX_CATALOG_ROW_BYTES
        elif table_file == CONTRADICTIONS:
            row_cap = MAX_CONTRADICTION_ROW_BYTES
        for line_no, ln in enumerate(_read_prose(table_file).splitlines(), start=1):
            if ln.startswith("|") and not re.match(r"^\|\s*---", ln):
                if len(ln.encode("utf-8")) > row_cap:
                    _fail(
                        failures,
                        table_file,
                        line_no,
                        f"table row exceeds {row_cap} B",
                    )
    if SOURCES.exists():
        for header, rows in _parse_tables_with_lines(SOURCES):
            if "Notes" not in header:
                continue
            notes_i = header.index("Notes")
            for line_no, row in rows:
                if (
                    len(row) > notes_i
                    and len(row[notes_i].encode("utf-8")) > MAX_EXTERNAL_NOTES_BYTES
                ):
                    _fail(
                        failures,
                        SOURCES,
                        line_no,
                        f"External Sources Notes cell exceeds {MAX_EXTERNAL_NOTES_BYTES} B",
                    )


def _check_computational_physics_language(
    pages: List[Path], failures: List[str]
) -> None:
    """
    Flag unquoted banned software metaphors in governed prose outside mathematical idioms.

    :param pages: Wiki pages to scan.
    :param failures: Mutable failure list.
    """
    for path in dict.fromkeys(
        [AGENTS.resolve(), SOURCES.resolve()] + [p.resolve() for p in pages]
    ):
        if path.name == "computational-physics-language.md" or not path.exists():
            continue

        def keep_prose_label(match: Match[str]) -> str:
            label = match.group(1).strip()
            return " " if FILENAME_LABEL_RE.fullmatch(label) else label

        masked = _mask_inline_code(_read_prose(path))
        for line_no, ln in enumerate(masked.splitlines(), start=1):
            cleaned = re.sub(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", keep_prose_label, ln)
            cleaned = re.sub(r'(?:"[^"\n]+"|\u201c[^\u201d\n]+\u201d)', " ", cleaned)
            cleaned = MATH_IDIOM_EXEMPT_RE.sub(" ", cleaned)
            for match in BANNED_METAPHOR_RE.finditer(cleaned):
                _fail(
                    failures,
                    path,
                    line_no,
                    f"banned prose metaphor {match.group(0)!r}",
                )


def _check_reachability(
    pages: List[Path], graph: Dict[Path, List[Path]], failures: List[str]
) -> None:
    """
    Check all wiki pages are reachable from ``AGENTS.md`` and linked from immediate parents.

    :param pages: Wiki pages to verify.
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

    for page in pages:
        if page not in reachable:
            _fail(failures, page, None, "wiki page not reachable from AGENTS.md")
        if _is_normal_leaf(page):
            parent = (
                AGENTS.resolve()
                if page.parent == WIKI
                else (page.parent / "index.md").resolve()
            )
            if page not in graph.get(parent, []):
                _fail(
                    failures,
                    page,
                    None,
                    f"leaf not linked from immediate parent router {_rel(parent)}",
                )
        elif _is_router(page) and page != AGENTS.resolve():
            parent = (
                AGENTS.resolve()
                if page.parent.parent == WIKI
                else (page.parent.parent / "index.md").resolve()
            )
            if page not in graph.get(parent, []):
                _fail(
                    failures,
                    page,
                    None,
                    f"router not linked from immediate parent router {_rel(parent)}",
                )


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
            target.name == "index.md" and _is_relative_to(target, WIKI)
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


def _unregistered_prose_urls(
    prose_masked: str, registered: Set[str]
) -> List[Tuple[int, str]]:
    """
    Find bare ``http``/``https`` URLs outside Markdown links and inline code that lack registration.

    Doctests:
    >>> _unregistered_prose_urls(
    ...     "See https://arxiv.org/pdf/0902.3652v2 and [link](https://unregistered.invalid/x).",
    ...     {"https://arxiv.org/pdf/0902.3652v2"},
    ... )
    []
    >>> _unregistered_prose_urls(
    ...     "Bare https://unregistered.invalid/doc.",
    ...     {"https://arxiv.org/pdf/"},
    ... )
    [(1, 'https://unregistered.invalid/doc')]

    :param prose_masked: Non-prose and inline-code-masked page text.
    :param registered: Registered source literals from ``raw/SOURCES.md``.
    :return: ``(line_no, url)`` pairs for unregistered bare URLs.
    """
    missing: List[Tuple[int, str]] = []
    for line_no, ln in enumerate(prose_masked.splitlines(), start=1):
        without_md_links = re.sub(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", " ", ln)
        for match in URL_RE.finditer(without_md_links):
            bare_url = match.group(0).rstrip(".,;:!?'\"\\]").split("#", 1)[0]
            if bare_url and not _source_registered(bare_url, registered):
                missing.append((line_no, bare_url))
    return missing


def _check_source_registration(
    pages: List[Path], registered: Set[str], failures: List[str]
) -> None:
    """
    Check external URLs across governed pages and literal local Sources links are registered.

    :param pages: Governed Markdown pages to inspect.
    :param registered: Registered source literals from ``raw/SOURCES.md``.
    :param failures: Mutable failure list.
    """
    for path in pages:
        if not path.exists():
            continue
        is_leaf = _is_normal_leaf(path)
        for link in _extract_links(path):
            parsed = urlparse(link.target)
            if parsed.scheme in {"http", "https"}:
                ext_url = link.target.split("#", 1)[0]
                if not _source_registered(ext_url, registered):
                    _fail(
                        failures,
                        path,
                        link.line,
                        f"external URL is not registered: {ext_url}",
                    )
                continue
            target = link.target_path
            if not is_leaf or not link.is_source_section or target is None:
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
        for line_no, bare_url in _unregistered_prose_urls(
            _mask_inline_code(_read_prose(path)), registered
        ):
            _fail(
                failures,
                path,
                line_no,
                f"external URL is not registered: {bare_url}",
            )


def _check_catalog(pages: List[Path], failures: List[str]) -> None:
    """
    Check catalog rows use ``Page | Scope | Query terms`` and cover each normal leaf once.

    :param pages: Wiki pages to verify.
    :param failures: Mutable failure list.
    """
    if not CATALOG.exists():
        return
    required = ["Page", "Scope", "Query terms"]
    header, located_rows = _table_with_header(CATALOG, "Page")
    if header != required:
        _fail(
            failures,
            CATALOG,
            None,
            "catalog columns must be exactly Page | Scope | Query terms",
        )
        return
    links: List[Tuple[Path, int]] = []
    for line, row in located_rows:
        if len(row) != len(required):
            _fail(failures, CATALOG, line, "catalog row has wrong column count")
            continue
        if not row[1].strip() or not row[2].strip():
            _fail(
                failures,
                CATALOG,
                line,
                "catalog Scope and Query terms must be nonempty",
            )
        page_links = list(re.finditer(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)", row[0]))
        if len(page_links) != 1:
            _fail(failures, CATALOG, line, "catalog Page cell must contain one link")
            continue
        target, _ = _path_from_link(CATALOG, page_links[0].group(2))
        if target is None or not _is_relative_to(target, WIKI):
            _fail(failures, CATALOG, line, "catalog Page target is outside wiki")
            continue
        ps_words = {
            word.lower()
            for word in re.findall(
                r"[A-Za-z0-9_]{3,}", f"{page_links[0].group(1)} {row[1]}"
            )
        }
        dup_words = sorted(
            {word.lower() for word in re.findall(r"[A-Za-z0-9_]{3,}", row[2])}
            & ps_words
        )
        if dup_words:
            _fail(
                failures,
                CATALOG,
                line,
                f"catalog Query terms repeat Page/Scope words: {', '.join(dup_words)}",
            )
        links.append((target, line))

    counts = Counter(target for target, _ in links)
    for page in pages:
        count = counts.get(page, 0)
        if count == 0 and _is_normal_leaf(page):
            _fail(failures, CATALOG, None, f"missing wiki page: {_rel(page)}")
        elif count > 1:
            _fail(failures, CATALOG, None, f"duplicate wiki page: {_rel(page)}")
    live_leaves = {page for page in pages if _is_normal_leaf(page)}
    for target, line in links:
        if target not in live_leaves:
            _fail(
                failures,
                CATALOG,
                line,
                f"catalog target is not a live normal leaf: {_rel(target)}",
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


def _split_table_cells(line: str) -> List[str]:
    r"""
    Split a Markdown table row on unescaped pipes outside inline code.

    Doctests:
    >>> _split_table_cells(r"| `a | b` | c \| d | e |")
    ['`a | b`', 'c \\| d', 'e']

    :param line: Raw Markdown table row starting and ending with ``|``.
    :return: Stripped cell strings.
    """
    inner = line.strip()
    if inner.startswith("|"):
        inner = inner[1:]
    if inner.endswith("|") and not inner.endswith("\\|"):
        inner = inner[:-1]
    masked = _mask_inline_code(inner)
    cuts = [-1] + [m.start() for m in re.finditer(r"(?<!\\)\|", masked)] + [len(inner)]
    return [inner[cuts[i] + 1 : cuts[i + 1]].strip() for i in range(len(cuts) - 1)]


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
    prose = _read_prose(path)
    for line_number, line in enumerate(prose.splitlines(), start=1):
        if not line.startswith("|"):
            if header is not None:
                tables.append((header, rows))
                header = None
                rows = []
            continue
        cells = _split_table_cells(line)
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


def _manifest_sources() -> Tuple[Set[str], Dict[str, str]]:
    """
    Return registered manifest literals and their declared statuses.

    :return: Tuple of registered source literals and literal-to-status mapping.
    """
    registered: Set[str] = set()
    statuses: Dict[str, str] = {}
    for header, rows in _parse_tables_with_lines(SOURCES):
        if "Source" not in header:
            continue
        source_i = header.index("Source")
        status_i = header.index("Status") if "Status" in header else None
        for _, row in rows:
            if len(row) <= source_i:
                continue
            values = _literal_code_values(row[source_i])
            registered.update(values)
            if status_i is not None and len(row) > status_i:
                status_val = row[status_i].strip().lower()
                for val in values:
                    statuses[val] = status_val
    return registered, statuses


def _source_registered(literal: str, registered: Set[str]) -> bool:
    """
    Return whether a source path or URL is exact-registered or aggregate-covered.

    Aggregate coverage follows the explicit path and URL-prefix rules named in
    the source manifest. It does not infer coverage from prose.

    :param literal: Repository-relative source path or external URL.
    :param registered: Manifest literals.
    :return: Whether registration is unambiguous.
    """
    if literal in registered or literal == "raw/SOURCES.md":
        return True
    if any(
        prefix in registered
        and (
            literal == prefix
            or literal.startswith(prefix if prefix.endswith("/") else prefix + "/")
        )
        for prefix in REGISTERED_URL_PREFIXES
    ):
        return True
    if urlparse(literal).scheme:
        return False
    path = Path(literal)
    suffix = path.suffix.lower()
    rules = [
        (
            "core-top-level-package-modules",
            path.parent == Path("nrpy")
            and (suffix in {".py", ".txt"} or path.name == "py.typed"),
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
            _is_relative_to(ROOT / path, ROOT / "nrpy/infrastructures"),
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


def _page_header_status(path: Path) -> Optional[str]:
    """
    Extract status from a page header.

    :param path: Wiki page.
    :return: Lowercase status, when present.
    """
    header = "\n".join(_read(path).splitlines()[:5])
    status_match = re.search(
        r"\bStatus:\s*(confirmed|provisional|contested|stale)\b",
        header,
        re.IGNORECASE,
    )
    return status_match.group(1).lower() if status_match else None


def _check_glossary_links(failures: List[str]) -> None:
    """
    Check glossary forbids ``## Sources`` and links a defining page or background note.

    :param failures: Mutable failure list.
    """
    if not GLOSSARY.exists():
        return
    if re.search(r"^## Sources\s*$", _read_prose(GLOSSARY), re.MULTILINE):
        _fail(failures, GLOSSARY, None, "glossary must not carry a Sources section")
    glossary_header, glossary_rows = _table_with_header(GLOSSARY, "Term")
    try:
        term_i = glossary_header.index("Term")
        meaning_i = glossary_header.index("Meaning")
    except ValueError:
        _fail(failures, GLOSSARY, None, "glossary lacks Term/Meaning table")
        return

    for line_no, row in glossary_rows:
        if len(row) <= max(term_i, meaning_i):
            continue
        term = row[term_i].strip()
        meaning = row[meaning_i].strip()
        if len(meaning.encode("utf-8")) > MAX_GLOSSARY_MEANING_BYTES:
            _fail(
                failures,
                GLOSSARY,
                line_no,
                f"glossary Meaning cell exceeds {MAX_GLOSSARY_MEANING_BYTES} B: {term}",
            )
        has_link = bool(re.search(r"\[[^\]]+\]\([^)]+\)", meaning))
        external = "external" in meaning.lower() or "background" in meaning.lower()
        if not (has_link or external):
            _fail(
                failures,
                GLOSSARY,
                line_no,
                f"glossary term lacks a defining-page link or external/background designation: {term}",
            )
    for sa_ln in _section_body(_read_prose(GLOSSARY), "See Also").splitlines():
        if sa_ln.strip().startswith("- ") and not sa_ln.strip()[2:].strip().startswith(
            ALLOWED_SEE_ALSO_PREFIXES
        ):
            _fail(
                failures,
                GLOSSARY,
                None,
                f"glossary See Also bullet lacks allowed prefix: {sa_ln.strip()[:50]}",
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
    Return authored KB files governed by volatile-data checks.

    Frozen source snapshots preserve imported evidence verbatim and are not
    authored KB claims. Their registration in ``raw/SOURCES.md`` is governed.

    :return: Authored Markdown files in the compiled KB and source manifest.
    """
    return _iter_md_files()


def _metadata_mention_allowed(lines: List[str], line_index: int) -> bool:
    """
    Return whether a removed maintenance-field mention is prohibition/supersession text.

    Doctests:
    >>> lines = ["Never store mtime values.", "Record mtime values."]
    >>> [_metadata_mention_allowed(lines, index) for index in range(2)]
    [True, False]

    :param lines: File lines being checked.
    :param line_index: Zero-based index of the line with the mention.
    :return: Whether same-line context permits the mention.
    """
    return METADATA_PROHIBITION_RE.search(lines[line_index]) is not None


def _volatile_metadata_issues(line: str) -> List[str]:
    """
    Return volatile-data policy violations found in one line.

    Doctests:
    >>> rejected = (
    ...     "Last reconciled: 09-12-2026",
    ...     "checked at 2026-09-12T14:03:22Z",
    ...     "Checked at 12:34.",
    ...     "commit 246043709e806021fcfc011fe657b8bf964cae4c",
    ...     "revision deadbee",
    ...     "sha256: 0123456789abcdef",
    ...     "inventory contains 390 files",
    ...     "The inventory covers 29 generators",
    ...     "The inventory configures 74 trusted comparisons.",
    ...     "| Source count | Last checked |",
    ...     "inspected=pass; generated=not-run",
    ...     "platform=linux; compiler=gcc-12",
    ...     "| Audit | Resolution |",
    ...     "No generated thorn build or restart validation was done.",
    ...     "Isolated local OpenMP builds and one-step startup were performed.",
    ...     "Validation artifacts were temporary and were not registered.",
    ... )
    >>> all(_volatile_metadata_issues(value) for value in rejected)
    True
    >>> allowed = (
    ...     "The state vector has 24 evolved components.",
    ...     "Finite-difference order 4 needs two centered halo points.",
    ...     "Compare CoordSystem_hash before restoring a checkpoint.",
    ...     "Brown, arXiv:0902.3652v2.",
    ...     (
    ...         "https://gist.githubusercontent.com/karpathy/"
    ...         "442a6bf555914893e9891c11519de94f/raw/"
    ...         "ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md"
    ...     ),
    ...     "https://example.org/archive/2026-09-12/report",
    ...     "Never store dates, timestamps, hashes, or file counts in KB prose.",
    ... )
    >>> [_volatile_metadata_issues(value) for value in allowed]
    [[], [], [], [], [], [], []]

    :param line: Authored KB line to inspect.
    :return: Violation messages for the line.
    """
    issues: List[str] = []
    # Complete locators are source identifiers. Opaque path components are not
    # maintenance snapshots even when they resemble dates or digests.
    snapshot_line = URL_RE.sub("", line)
    checks = (
        (HASH_DIGEST_VALUE_RE, snapshot_line, "hash digest value found"),
        (VCS_HASH_VALUE_RE, snapshot_line, "VCS revision value found"),
        (VOLATILE_TIME_RE, snapshot_line, "timestamp literal found"),
        (VOLATILE_INVENTORY_COUNT_RE, snapshot_line, "inventory count found"),
        (
            VOLATILE_EXISTENTIAL_INVENTORY_RE,
            snapshot_line,
            "inventory count found",
        ),
        (VOLATILE_CATALOG_COUNT_RE, snapshot_line, "catalog count found"),
        (VOLATILE_FIELD_RE, snapshot_line, "volatile maintenance field found"),
        (VOLATILE_RESULT_TUPLE_RE, snapshot_line, "recorded result tuple found"),
        (VOLATILE_RUN_PROSE_RE, snapshot_line, "recorded run/audit prose found"),
    )
    for pattern, subject, message in checks:
        if pattern.search(subject):
            issues.append(message)

    hex_match = STORED_HEX_IDENTIFIER_RE.search(snapshot_line)
    if hex_match and not STABLE_HEX_IDENTIFIER_PREFIX_RE.search(
        snapshot_line[: hex_match.start()]
    ):
        issues.append("hash-like value found")

    date_match = VOLATILE_DATE_RE.search(snapshot_line)
    if date_match and not STABLE_DATE_IDENTIFIER_RE.search(
        snapshot_line[: date_match.start()]
    ):
        issues.append("date literal found")
    if len(ENVIRONMENT_ASSIGNMENT_RE.findall(snapshot_line)) >= 2:
        issues.append("recorded environment tuple found")
    return issues


def _check_source_tracking_metadata(failures: List[str]) -> None:
    """
    Check governed KB files for forbidden maintenance fields and snapshot values.

    :param failures: Mutable failure list.
    """
    for path in _governed_kb_files():
        lines = _read(path).splitlines()
        for idx, line in enumerate(lines, start=1):
            for issue in _volatile_metadata_issues(line):
                _fail(failures, path, idx, issue)
            if MTIME_WORD_RE.search(line) and not _metadata_mention_allowed(
                lines, idx - 1
            ):
                _fail(failures, path, idx, "source-tracking Mtime value found")


def _check_status_vocabularies(pages: List[Path], failures: List[str]) -> None:
    """
    Check page, manifest, and source-map status vocabularies.

    :param pages: Wiki pages to inspect.
    :param failures: Mutable failure list.
    """
    for page in pages:
        if _is_router(page):
            continue
        status = _page_header_status(page)
        if status not in PAGE_STATUSES:
            _fail(failures, page, 2, "page header has missing or invalid Status")

    for path in (SOURCES, SOURCE_MAP):
        for header, rows in _parse_tables_with_lines(path):
            if "Status" not in header:
                continue
            status_i = header.index("Status")
            ingest_i = header.index("Ingest") if "Ingest" in header else None
            tier_i = header.index("Tier") if "Tier" in header else None
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
                if tier_i is not None:
                    if (
                        len(row) <= tier_i
                        or row[tier_i].strip().lower() not in SOURCE_AUTHORITY_TIERS
                    ):
                        value = row[tier_i] if len(row) > tier_i else ""
                        _fail(failures, path, line, f"invalid Tier value: {value!r}")


def _check_source_map_targets(registered: Set[str], failures: List[str]) -> None:
    """
    Check source-map rows contain registered multi-file aggregates only.

    :param registered: Registered source literals from ``raw/SOURCES.md``.
    :param failures: Mutable failure list.
    """
    source_tables = [
        (header, rows)
        for header, rows in _parse_tables_with_lines(SOURCE_MAP)
        if "Source or aggregate" in header
    ]
    if len(source_tables) != 1:
        _fail(failures, SOURCE_MAP, None, "source map must carry one aggregate table")
        return
    header, rows = source_tables[0]
    if header != SOURCE_MAP_COLUMNS:
        _fail(
            failures,
            SOURCE_MAP,
            None,
            "source-map column headings differ from SOURCE_MAP_COLUMNS",
        )
        return
    source_i = header.index("Source or aggregate")
    glob_i = header.index("Glob pattern")
    authority_i = header.index("Authority tier")
    scope_i = header.index("Uncovered scope")
    seen_aggregates: Set[str] = set()
    for line, row in rows:
        if len(row) != len(header):
            _fail(failures, SOURCE_MAP, line, "source-map row has wrong field count")
            continue
        tiers = {tier.strip() for tier in row[authority_i].split(";")}
        if not tiers or not tiers.issubset(SOURCE_AUTHORITY_TIERS):
            _fail(failures, SOURCE_MAP, line, "invalid source authority tier")
        if not row[glob_i].strip() or not row[scope_i].strip():
            _fail(
                failures,
                SOURCE_MAP,
                line,
                "source-map Glob pattern and Uncovered scope must be nonempty",
            )
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
        if literal not in AGGREGATE_SOURCE_NAMES:
            _fail(
                failures,
                SOURCE_MAP,
                line,
                f"source-map allows multi-file aggregate rows only: {literal}",
            )
        elif literal in seen_aggregates:
            _fail(failures, SOURCE_MAP, line, f"duplicate aggregate row: {literal}")
        seen_aggregates.add(literal)
        if literal not in registered:
            _fail(
                failures,
                SOURCE_MAP,
                line,
                f"source-map literal is not registered: {literal}",
            )
    for missing_agg in sorted(AGGREGATE_SOURCE_NAMES - seen_aggregates):
        _fail(failures, SOURCE_MAP, None, f"missing aggregate row: {missing_agg}")


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


def _check_contradictions(pages: List[Path], failures: List[str]) -> None:
    """
    Check 4-column active contradiction rows, H3 subsections, and leaf markers.

    :param pages: Wiki pages to inspect.
    :param failures: Mutable failure list.
    """
    required = [
        "ID",
        "Claim status",
        "Affected pages",
        "Summary",
    ]
    header, rows = _table_with_header(CONTRADICTIONS, "ID")
    if header != required:
        _fail(
            failures,
            CONTRADICTIONS,
            None,
            "contradiction table headings must be ID | Claim status | Affected pages | Summary",
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
            if not row[indexes[field]].strip():
                _fail(
                    failures,
                    CONTRADICTIONS,
                    line,
                    f"empty contradiction field: {field}",
                )
        status = row[indexes["Claim status"]].strip().lower()
        if status not in {"contested", "stale"}:
            _fail(
                failures,
                CONTRADICTIONS,
                line,
                f"invalid or unretired Claim status: {status!r}",
            )
        affected = set(_affected_page_links(row[indexes["Affected pages"]]))
        if not affected:
            _fail(failures, CONTRADICTIONS, line, "Affected pages has no wiki links")
        if status in {"contested", "stale"}:
            active[identifier] = (status, affected)

    subsection_ids = {
        title.strip()
        for level, title in _headings(_read_prose(CONTRADICTIONS))
        if level == 3 and re.fullmatch(r"CONTR-\d{4}", title.strip())
    }
    for identifier in sorted(seen - subsection_ids):
        _fail(
            failures,
            CONTRADICTIONS,
            None,
            f"missing ### {identifier} subsection for active contradiction",
        )
    for identifier in sorted(subsection_ids - seen):
        _fail(
            failures,
            CONTRADICTIONS,
            None,
            f"orphan ### {identifier} subsection without Register row",
        )
    contr_prose = _read_prose(CONTRADICTIONS)
    retired_ids: Set[str] = set()
    for ret_match in re.finditer(r"\bRetired IDs:\s*(.+)", contr_prose):
        tail = ret_match.group(1)
        retired_ids.update(re.findall(r"CONTR-\d{4}", tail))
        for start_s, end_s in re.findall(r"CONTR-(\d{4})\s+to\s+CONTR-(\d{4})", tail):
            for num in range(int(start_s), int(end_s) + 1):
                retired_ids.add(f"CONTR-{num:04d}")
    for reused in sorted(retired_ids & (seen | subsection_ids)):
        _fail(
            failures, CONTRADICTIONS, None, f"retired contradiction ID reused: {reused}"
        )
    for identifier in sorted(seen & subsection_ids):
        sub_m = re.search(
            rf"^### {re.escape(identifier)}\s*$(.*?)(?=^##|^###|\Z)",
            contr_prose,
            re.MULTILINE | re.DOTALL,
        )
        sub_body = sub_m.group(1) if sub_m else ""
        for req_field in (
            "Competing sources",
            "Recheck trigger:",
            "Resolution check:",
            "Evidence (",
        ):
            if req_field not in sub_body:
                _fail(
                    failures,
                    CONTRADICTIONS,
                    None,
                    f"### {identifier} lacks {req_field!r}",
                )

    marker_re = re.compile(
        r"Claim status:\s*(contested|stale);\s*contradiction:\s*(CONTR-\d{4})[.]",
        re.IGNORECASE,
    )
    found_markers: Dict[str, Set[str]] = defaultdict(set)
    for page in pages:
        if not _is_normal_leaf(page):
            continue
        text = _mask_inline_code(_read_prose(page))
        page_ids: Set[str] = set()
        for match in marker_re.finditer(text):
            status, identifier = match.group(1).lower(), match.group(2)
            page_ids.add(identifier)
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
                    f"marker references missing or retired contradiction: {identifier}",
                )
            elif active[identifier][0] != status:
                _fail(
                    failures,
                    page,
                    _line_for_offset(text, match.start()),
                    f"marker status disagrees with {identifier}",
                )
        if _page_header_status(page) == "confirmed" and len(page_ids) > 1:
            _fail(
                failures,
                page,
                None,
                f"confirmed leaf exceeds 1 contradiction marker ({len(page_ids)})",
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
    Check KB trees contain Markdown prose only and no retired maintenance files.

    :param failures: Mutable failure list.
    """
    for base in (WIKI, RAW):
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix.lower() != ".md":
                _fail(failures, path, None, "non-Markdown file in KB tree")
            if path.is_file() and re.search(
                r"(?:^|[-_])(audit|scan|report|latest|token|plan|rank|log)(?:[-_.]|$)",
                path.name,
                re.IGNORECASE,
            ):
                _fail(failures, path, None, "excluded maintenance/planning file")


def _check_root_artifacts(failures: List[str]) -> None:
    """
    Check for prohibited root-level maintenance files.

    :param failures: Mutable failure list.
    """
    root_tmp = ROOT / "tmp"
    if root_tmp.exists() or root_tmp.is_symlink():
        try:
            is_system_tmp = root_tmp.samefile(Path("/tmp"))
        except OSError:
            is_system_tmp = False
        if not is_system_tmp:
            _fail(
                failures,
                root_tmp,
                None,
                "root-level scratch directory tmp/ prohibited",
            )
    allowed_root_docs = {"AGENTS.md", "README.md", "CITATION.md", "coding_style.md"}
    for path in ROOT.glob("*.md"):
        if path.name in allowed_root_docs:
            continue
        if ROOT_PROHIBITED_ARTIFACT_SEGMENT_RE.search(path.name):
            _fail(failures, path, None, "root-level prohibited maintenance file")
            continue
        if ROOT_COORDINATION_RE.fullmatch(path.name):
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
            _fail(failures, path, None, "root-level prohibited maintenance file")


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

    _clear_file_caches()
    failures: List[str] = []
    files = _iter_md_files()
    pages = _wiki_pages()
    registered, source_statuses = _manifest_sources()

    _check_wikilinks(files, failures)
    graph = _check_links(files, failures)
    _check_router_detail(files, failures)
    _check_router_edges(files, failures)
    _check_leaf_contract(pages, registered, source_statuses, failures)
    _check_token_efficiency_invariants(pages, registered, failures)
    _check_computational_physics_language(pages, failures)
    _check_source_registration(
        ([AGENTS] if AGENTS.exists() else []) + pages, registered, failures
    )
    _check_reachability(pages, graph, failures)
    _check_catalog(pages, failures)
    _check_content_dir_indexes(failures)
    _check_glossary_links(failures)
    _check_status_vocabularies(pages, failures)
    _check_source_map_targets(registered, failures)
    _check_contradictions(pages, failures)
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
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)
    sys.exit(_main())
