"""CRUD forms for OncoSearch's existing SQLite schema.

No schema changes needed. UI uses parameterized SQLite queries, stable row IDs,
confirmation before deletion, and prevents deleting paths outside dataset_dir.
"""
from pathlib import Path
from uuid import uuid4
import os
import shutil
import sqlite3

import streamlit as st

DNA_ALPHABET = set("ACGTRYSWKMBDHVN")
QUERY_ALPHABET = set("ACGT")


def _rows(db_path, sql, params=()):
    with sqlite3.connect(db_path) as conn:
        conn.row_factory = sqlite3.Row
        return [dict(x) for x in conn.execute(sql, params).fetchall()]


def _safe_path(dataset_dir, stored_name):
    """Never act on outside-of-library paths or symlink targets."""
    base = Path(dataset_dir).resolve()
    if (not stored_name or Path(stored_name).name != stored_name
            or stored_name in (".", "..")):
        raise ValueError("Unsafe stored filename")
    candidate = base / stored_name
    if candidate.is_symlink():
        raise ValueError("Refusing to modify a symlink")
    resolved = candidate.resolve()
    if resolved.parent != base:
        raise ValueError("Unsafe file path")
    return resolved


def update_dataset(db_path, dataset_id, new_name, dataset_dir=None, new_sequence=None):
    """Rename a dataset; optionally replace manually entered DNA with backup."""
    new_name = new_name.strip()
    if not new_name:
        raise ValueError("Dataset name cannot be blank")
    with sqlite3.connect(db_path) as conn:
        row = conn.execute(
            "SELECT stored_filename, source FROM datasets WHERE id=?", (dataset_id,)
        ).fetchone()
        if not row:
            raise ValueError("Dataset no longer exists")
        stored_name, source = row
        if new_sequence is None:
            conn.execute("UPDATE datasets SET name=? WHERE id=?", (new_name, dataset_id))
            return

        if source != "manual" or dataset_dir is None:
            raise ValueError("Only manually entered DNA can be updated in place")
        sequence = "".join(new_sequence.upper().split())
        if not sequence or set(sequence) - DNA_ALPHABET:
            raise ValueError("Invalid DNA sequence; use DNA/IUPAC nucleotide letters")
        path = _safe_path(dataset_dir, stored_name)
        if not path.is_file():
            raise FileNotFoundError(path)
        current_seq = "".join(
            line.strip() for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith(">")
        ).upper()
        if current_seq == sequence:
            conn.execute("UPDATE datasets SET name=? WHERE id=?", (new_name, dataset_id))
            return
        other_refs = conn.execute(
            "SELECT COUNT(*) FROM datasets WHERE stored_filename=? AND id != ?",
            (stored_name, dataset_id),
        ).fetchone()[0]
        if other_refs:
            raise ValueError("File is shared by another dataset; copy it before editing")

        # Keep old version for reproducibility; app saves old DNA to data/revisions.
        revision_dir = Path(dataset_dir).resolve().parent / "revisions"
        revision_dir.mkdir(parents=True, exist_ok=True)
        backup_path = revision_dir / f"{uuid4().hex}_{stored_name}"
        shutil.copy2(path, backup_path)
        tmp_path = path.with_name(f".{uuid4().hex}.tmp")
        try:
            tmp_path.write_text(f">manual_sequence\n{sequence}\n", encoding="utf-8")
            os.replace(tmp_path, path)
            conn.execute("UPDATE datasets SET name=?, size_bytes=? WHERE id=?",
                         (new_name, path.stat().st_size, dataset_id))
        except Exception:
            # DB mutation failure should not silently lose the old sequence.
            shutil.copy2(backup_path, path)
            raise
        finally:
            tmp_path.unlink(missing_ok=True)


def delete_dataset(db_path, dataset_dir, dataset_id, remove_file=False):
    """Delete library entry; optionally remove physical file if unreferenced."""
    with sqlite3.connect(db_path) as conn:
        row = conn.execute("SELECT stored_filename FROM datasets WHERE id=?",
                           (dataset_id,)).fetchone()
        if row is None:
            raise ValueError("Dataset no longer exists")
        stored_name = row[0]
        path = _safe_path(dataset_dir, stored_name)
        conn.execute("DELETE FROM datasets WHERE id=?", (dataset_id,))
        references = conn.execute(
            "SELECT COUNT(*) FROM datasets WHERE stored_filename=?", (stored_name,)
        ).fetchone()[0]
    if remove_file and references == 0:
        path.unlink(missing_ok=True)
    return bool(remove_file and references == 0)


def update_query(db_path, query_id, panel, name, sequence):
    panel, name, sequence = panel.strip(), name.strip(), sequence.strip().upper()
    if not panel or not name or not sequence or set(sequence) - QUERY_ALPHABET:
        raise ValueError("Panel/name required; exact query may contain only A/C/G/T")
    with sqlite3.connect(db_path) as conn:
        existing = conn.execute("SELECT 1 FROM queries WHERE id=?", (query_id,)).fetchone()
        if not existing:
            raise ValueError("Query no longer exists")
        duplicate = conn.execute(
            "SELECT 1 FROM queries WHERE panel=? AND name=? AND id != ?",
            (panel, name, query_id),
        ).fetchone()
        if duplicate:
            raise ValueError("This query name already exists in the selected panel")
        conn.execute("UPDATE queries SET panel=?, name=?, sequence=? WHERE id=?",
                     (panel, name, sequence, query_id))


def delete_query(db_path, query_id):
    with sqlite3.connect(db_path) as conn:
        conn.execute("DELETE FROM queries WHERE id=?", (query_id,))


def rename_panel(db_path, old_panel, new_panel):
    new_panel = new_panel.strip()
    if not new_panel:
        raise ValueError("New panel name cannot be empty")
    with sqlite3.connect(db_path) as conn:
        if conn.execute("SELECT 1 FROM queries WHERE panel=?", (new_panel,)).fetchone():
            raise ValueError("Panel name already exists")
        conn.execute("UPDATE queries SET panel=? WHERE panel=?", (new_panel, old_panel))


def delete_panel(db_path, panel):
    with sqlite3.connect(db_path) as conn:
        conn.execute("DELETE FROM queries WHERE panel=?", (panel,))


def _manual_sequence(dataset_dir, stored_filename):
    path = _safe_path(dataset_dir, stored_filename)
    if not path.exists():
        return ""
    # Manual sequences were saved as a single FASTA record.
    return "".join(
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )


def render_dataset_manager(db_path, dataset_dir):
    rows = _rows(db_path, "SELECT * FROM datasets ORDER BY created_at DESC")
    if not rows:
        st.info("No datasets saved yet. Upload a FASTQ or add a manual sequence.")
        return

    st.caption("Update or delete any dataset. For FASTQ/FASTA, update changes the display name, not the original sequencing file.")
    headings = st.columns([3.5, 2.9, 1.0, 0.9, 0.9])
    for col, name in zip(headings, ["Dataset name", "Original file", "Size", "Update", "Delete"]):
        col.markdown(f"**{name}**")
    for row in rows:
        cells = st.columns([3.5, 2.9, 1.0, 0.9, 0.9], vertical_alignment="center")
        cells[0].write(row["name"])
        cells[1].caption(f"{row['original_filename']} · {row['source']}")
        cells[2].caption(f"{row['size_bytes'] / 1048576:.3f} MB")
        if cells[3].button("✏️", key=f"dataset_edit_{row['id']}", help="Update this dataset"):
            st.session_state["editing_dataset"] = row["id"]
            st.session_state.pop("deleting_dataset", None)
            st.rerun()
        if cells[4].button("🗑️", key=f"dataset_del_{row['id']}", help="Delete this dataset"):
            st.session_state["deleting_dataset"] = row["id"]
            st.session_state.pop("editing_dataset", None)
            st.rerun()

    editing = next((x for x in rows if x["id"] == st.session_state.get("editing_dataset")), None)
    if editing:
        st.divider()
        st.markdown(f"#### ✏️ Update dataset — {editing['name']}")
        with st.form(f"edit_dataset_form_{editing['id']}"):
            new_name = st.text_input("Dataset display name", value=editing["name"])
            new_seq = None
            if editing["source"] == "manual":
                st.caption("Manual DNA can also be edited. A backup of its old sequence is saved in data/revisions/.")
                new_seq = st.text_area(
                    "DNA sequence", value=_manual_sequence(dataset_dir, editing["stored_filename"]),
                    height=110,
                )
            save = st.form_submit_button("Save updates", type="primary")
        if save:
            try:
                update_dataset(db_path, editing["id"], new_name, dataset_dir, new_seq)
            except (ValueError, FileNotFoundError, OSError, sqlite3.Error) as exc:
                st.error(str(exc))
            else:
                st.session_state.pop("editing_dataset", None)
                st.success("Dataset updated.")
                st.rerun()
        if st.button("Cancel update", key="cancel_dataset_update"):
            st.session_state.pop("editing_dataset", None)
            st.rerun()

    deleting = next((x for x in rows if x["id"] == st.session_state.get("deleting_dataset")), None)
    if deleting:
        st.divider()
        st.warning(f"Delete dataset '{deleting['name']}'? This cannot be undone.")
        remove_file = st.checkbox(
            "Also permanently delete its FASTQ/FASTA file from server disk",
            value=False, key=f"remove_dataset_file_{deleting['id']}",
        )
        a, b = st.columns(2)
        if a.button("Confirm deletion", type="primary", key="confirm_dataset_delete"):
            try:
                delete_dataset(db_path, dataset_dir, deleting["id"], remove_file)
            except (ValueError, OSError, sqlite3.Error) as exc:
                st.error(str(exc))
            else:
                st.session_state.pop("deleting_dataset", None)
                st.rerun()
        if b.button("Cancel", key="cancel_dataset_delete"):
            st.session_state.pop("deleting_dataset", None)
            st.rerun()


def render_query_manager(db_path):
    st.subheader("Saved Search Queries")
    rows = _rows(db_path, "SELECT * FROM queries ORDER BY created_at DESC, rowid DESC")
    if not rows:
        st.info("No saved queries. Create a search panel above.")
        return

    headings = st.columns([2.5, 2.5, 3.4, 0.9, 0.9])
    for col, name in zip(headings, ["Panel", "Query name", "DNA sequence", "Update", "Delete"]):
        col.markdown(f"**{name}**")
    for row in rows:
        cells = st.columns([2.5, 2.5, 3.4, 0.9, 0.9], vertical_alignment="center")
        cells[0].caption(row["panel"])
        cells[1].write(row["name"])
        cells[2].code(row["sequence"], language=None, wrap_lines=True)
        if cells[3].button("✏️", key=f"query_edit_{row['id']}", help="Update this query"):
            st.session_state["editing_query"] = row["id"]
            st.session_state.pop("deleting_query", None)
            st.rerun()
        if cells[4].button("🗑️", key=f"query_delete_{row['id']}", help="Delete this query"):
            st.session_state["deleting_query"] = row["id"]
            st.session_state.pop("editing_query", None)
            st.rerun()

    editing = next((x for x in rows if x["id"] == st.session_state.get("editing_query")), None)
    if editing:
        st.divider()
        st.markdown(f"#### ✏️ Update query — {editing['name']}")
        with st.form(f"edit_query_form_{editing['id']}"):
            panel = st.text_input("Panel", value=editing["panel"])
            name = st.text_input("Query name", value=editing["name"])
            sequence = st.text_input("DNA search pattern", value=editing["sequence"])
            save = st.form_submit_button("Save query", type="primary")
        if save:
            try:
                update_query(db_path, editing["id"], panel, name, sequence)
            except (ValueError, sqlite3.Error) as exc:
                st.error(str(exc))
            else:
                st.session_state.pop("editing_query", None)
                st.rerun()
        if st.button("Cancel update", key="cancel_query_update"):
            st.session_state.pop("editing_query", None)
            st.rerun()

    deleting = next((x for x in rows if x["id"] == st.session_state.get("deleting_query")), None)
    if deleting:
        st.divider()
        st.warning(f"Delete query '{deleting['name']}' from panel '{deleting['panel']}'?")
        a, b = st.columns(2)
        if a.button("Confirm query deletion", type="primary"):
            delete_query(db_path, deleting["id"])
            st.session_state.pop("deleting_query", None)
            st.rerun()
        if b.button("Cancel deletion", key="cancel_query_deletion"):
            st.session_state.pop("deleting_query", None)
            st.rerun()

    st.divider()
    with st.expander("Manage entire query panel", expanded=False):
        panels = sorted({row["panel"] for row in rows})
        selected = st.selectbox("Select panel", panels, key="manage_panel")
        new_name = st.text_input("Rename panel to", key="new_panel_name")
        if st.button("Rename entire panel"):
            try:
                rename_panel(db_path, selected, new_name)
            except (ValueError, sqlite3.Error) as exc:
                st.error(str(exc))
            else:
                st.rerun()
        st.warning(f"Deleting panel '{selected}' removes all of its queries.")
        confirmed = st.checkbox("I understand — delete all queries in this panel",
                                key="confirm_delete_panel")
        if st.button("Delete entire panel", disabled=not confirmed):
            delete_panel(db_path, selected)
            st.rerun()
