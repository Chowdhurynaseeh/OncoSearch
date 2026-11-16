"""Tests for dataset/query CRUD; no Streamlit server required."""
import importlib.util
import sqlite3
import tempfile
import unittest
from pathlib import Path

# crud_ui only requires Streamlit at runtime when rendering UI.
if importlib.util.find_spec("streamlit") is None:
    import sys
    import types
    sys.modules.setdefault("streamlit", types.ModuleType("streamlit"))

from crud_ui import (
    update_dataset, delete_dataset, update_query, delete_query,
    rename_panel, delete_panel,
)


class CRUDTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.datasets = self.root / "datasets"
        self.datasets.mkdir()
        self.db = self.root / "oncosearch.db"
        with sqlite3.connect(self.db) as conn:
            conn.execute("""
                CREATE TABLE datasets(
                    id TEXT PRIMARY KEY, name TEXT, original_filename TEXT,
                    stored_filename TEXT, source TEXT, size_bytes INTEGER,
                    created_at TEXT)
            """)
            conn.execute("""
                CREATE TABLE queries(
                    id TEXT PRIMARY KEY, panel TEXT, name TEXT,
                    sequence TEXT, created_at TEXT)
            """)
            conn.execute("INSERT INTO queries VALUES(?,?,?,?,?)", (
                "q1", "panelA", "motif1", "ACAC", "now",
            ))
            conn.execute("INSERT INTO queries VALUES(?,?,?,?,?)", (
                "q2", "panelA", "motif2", "AAAA", "now",
            ))
            conn.execute("INSERT INTO datasets VALUES(?,?,?,?,?,?,?)", (
                "d1", "Existing", "manual_input", "abc.fasta", "manual", 18, "now",
            ))
        (self.datasets / "abc.fasta").write_text(">manual_sequence\nACACACAC\n")

    def test_rename_dataset(self):
        update_dataset(self.db, "d1", "New name")
        with sqlite3.connect(self.db) as conn:
            self.assertEqual(conn.execute(
                "SELECT name FROM datasets WHERE id='d1'").fetchone()[0], "New name")
        self.assertEqual((self.datasets / "abc.fasta").read_text(), ">manual_sequence\nACACACAC\n")

    def test_update_manual_dna_with_backup(self):
        update_dataset(self.db, "d1", "Changed", self.datasets, "ACGTAC")
        self.assertIn("ACGTAC", (self.datasets / "abc.fasta").read_text())
        revisions = list((self.root / "revisions").glob("*abc.fasta"))
        self.assertEqual(len(revisions), 1)
        self.assertIn("ACACACAC", revisions[0].read_text())

    def test_query_update_duplicate_guard(self):
        update_query(self.db, "q1", "panelA", "renamed", "ACGT")
        with sqlite3.connect(self.db) as conn:
            self.assertEqual(conn.execute(
                "SELECT sequence FROM queries WHERE id='q1'").fetchone()[0], "ACGT")
        with self.assertRaises(ValueError):
            update_query(self.db, "q1", "panelA", "motif2", "CCCC")

    def test_query_delete_and_panel_rename_delete(self):
        delete_query(self.db, "q1")
        rename_panel(self.db, "panelA", "newPanel")
        with sqlite3.connect(self.db) as conn:
            self.assertEqual(conn.execute("SELECT panel FROM queries").fetchone()[0], "newPanel")
        delete_panel(self.db, "newPanel")
        with sqlite3.connect(self.db) as conn:
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM queries").fetchone()[0], 0)

    def test_delete_metadata_only_preserves_file(self):
        removed = delete_dataset(self.db, self.datasets, "d1", False)
        self.assertFalse(removed)
        self.assertTrue((self.datasets / "abc.fasta").exists())

    def test_delete_metadata_and_file(self):
        removed = delete_dataset(self.db, self.datasets, "d1", True)
        self.assertTrue(removed)
        self.assertFalse((self.datasets / "abc.fasta").exists())

    def test_deleting_shared_file_ref_preserves_file(self):
        with sqlite3.connect(self.db) as conn:
            conn.execute("INSERT INTO datasets VALUES(?,?,?,?,?,?,?)", (
                "d2", "Another", "manual_input", "abc.fasta", "manual", 18, "now",
            ))
        removed = delete_dataset(self.db, self.datasets, "d1", True)
        self.assertFalse(removed)
        self.assertTrue((self.datasets / "abc.fasta").is_file())

    def test_reject_path_traversal(self):
        with sqlite3.connect(self.db) as conn:
            conn.execute("UPDATE datasets SET stored_filename='../outside.txt' WHERE id='d1'")
        with self.assertRaises(ValueError):
            delete_dataset(self.db, self.datasets, "d1", True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
