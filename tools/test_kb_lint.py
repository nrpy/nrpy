#!/usr/bin/env python3
"""Focused regression tests for KB volatility policy."""

import unittest

import kb_lint


class VolatileMetadataPolicyTests(unittest.TestCase):
    """Keep transient inventory and reconciliation facts out of KB prose."""

    def test_rejects_volatile_literals_and_fields(self) -> None:
        rejected = (
            "Last reconciled: 09-12-2026",
            "checked at 2026-09-12T14:03:22Z",
            "Checked at 12:34.",
            "commit 246043709e806021fcfc011fe657b8bf964cae4c",
            "revision deadbee",
            "sha256: 0123456789abcdef",
            "inventory contains 390 files",
            "The inventory covers 29 generators",
            "The inventory configures 74 trusted comparisons.",
            "| Source count | Last checked |",
            "inspected=pass; generated=not-run",
            "platform=linux; compiler=gcc-12",
            "| Audit | Resolution |",
            "No generated thorn build or restart validation was done.",
            "Isolated local OpenMP builds and one-step startup were performed.",
            "Validation artifacts were temporary and were not registered.",
        )
        for line in rejected:
            with self.subTest(line=line):
                self.assertTrue(kb_lint._volatile_metadata_issues(line))

    def test_allows_stable_technical_content(self) -> None:
        allowed = (
            "The state vector has 24 evolved components.",
            "Finite-difference order 4 needs two centered halo points.",
            "Compare CoordSystem_hash before restoring a checkpoint.",
            "Brown, arXiv:0902.3652v2.",
            "https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md",
            "https://example.org/archive/2026-09-12/report",
            "Never store dates, timestamps, hashes, or file counts in KB prose.",
        )
        for line in allowed:
            with self.subTest(line=line):
                self.assertEqual(kb_lint._volatile_metadata_issues(line), [])

    def test_mtime_prohibition_does_not_whitelist_a_later_claim(self) -> None:
        lines = ["Never store mtime values.", "Record mtime values."]
        self.assertTrue(kb_lint._metadata_mention_allowed(lines, 0))
        self.assertFalse(kb_lint._metadata_mention_allowed(lines, 1))


if __name__ == "__main__":
    unittest.main()
