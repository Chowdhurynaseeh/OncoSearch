"""OncoSearch: persistent FASTQ/FASTA library and pattern matching dashboard.

Drop-in replacement for the user's existing app.py. Retains the original
SQLite datasets and queries schemas and five existing tabs.
"""
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
import gzip
import sqlite3

import streamlit as st

from crud_ui import render_dataset_manager, render_query_manager
from search_page import render_search_page


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
DATASET_DIR = DATA_DIR / "datasets"
DB_PATH = DATA_DIR / "oncosearch.db"
DATASET_DIR.mkdir(parents=True, exist_ok=True)

DNA_ALPHABET = set("ACGTRYSWKMBDHVN")
QUERY_ALPHABET = set("ACGT")
ALLOWED_EXTENSIONS = (
    ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fna.gz",
    ".fastq", ".fq", ".fasta", ".fa", ".fna",
)


def initialize_database():
    """Create tables if missing. Never delete or migrate existing user data."""
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS datasets (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                original_filename TEXT NOT NULL,
                stored_filename TEXT NOT NULL,
                source TEXT NOT NULL,
                size_bytes INTEGER NOT NULL,
                created_at TEXT NOT NULL
            )
        """)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS queries (
                id TEXT PRIMARY KEY,
                panel TEXT NOT NULL,
                name TEXT NOT NULL,
                sequence TEXT NOT NULL,
                created_at TEXT NOT NULL
            )
        """)


def fetch_all(sql, params=()):
    with sqlite3.connect(DB_PATH) as conn:
        conn.row_factory = sqlite3.Row
        return [dict(row) for row in conn.execute(sql, params).fetchall()]


def now_utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def identify_extension(filename):
    low = filename.lower()
    return next((ext for ext in ALLOWED_EXTENSIONS if low.endswith(ext)), None)


def open_sequence_file(path):
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def validate_file(path, extension):
    """Quick format check of first record; full parsing occurs at search time."""
    try:
        with open_sequence_file(path) as handle:
            if extension.startswith(".fastq") or extension.startswith(".fq"):
                header = handle.readline().rstrip("\r\n")
                seq = handle.readline().strip().upper()
                sep = handle.readline().rstrip("\r\n")
                qual = handle.readline().rstrip("\r\n")
                return (header.startswith("@") and bool(seq)
                        and set(seq) <= DNA_ALPHABET and sep.startswith("+")
                        and len(qual) == len(seq))
            header = handle.readline().rstrip("\r\n")
            seq = handle.readline().strip().upper()
            return (header.startswith(">") and bool(seq)
                    and set(seq) <= DNA_ALPHABET)
    except (OSError, UnicodeError, EOFError):
        return False


def save_metadata(name, original, stored, source, size):
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("""
            INSERT INTO datasets
            (id, name, original_filename, stored_filename, source, size_bytes, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (str(uuid4()), name, original, stored, source, size, now_utc()))


st.set_page_config(page_title="OncoSearch", page_icon="🧬", layout="wide")
initialize_database()

st.title("🧬 OncoSearch")
st.markdown(
    "### **Group 7 — Algorithm Term Project**"
)

st.markdown(
    "**Group Members:** Abdullah Al Naseeh Chowdhury "
    "and Proloy Kumar Mondol"
)

st.caption(
    "CSCE 5150 — Analysis of Computer Algorithms | "
    "University of North Texas"
)
st.info(
    "OncoSearch provides an interactive environment for real FASTQ/FASTA "
    "sequence analysis, custom DNA query management, exact and multi-pattern "
    "searching, comparative algorithm benchmarking, and step-by-step "
    "visualization of DNA pattern matching techniques."
)

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📤 Upload FASTQ / FASTA",
    "🧬 Manual DNA Input",
    "🎯 Search Query Panel",
    "📂 Saved Dataset Library",
    "🔎 Search & Benchmark",
])

# TAB 1 — UPLOAD REAL FILES
with tab1:
    st.subheader("Upload Sequencing Data")
    with st.form("upload_form", clear_on_submit=False):
        dataset_name = st.text_input("Dataset Name", placeholder="e.g. NCBI_Ecoli_SRR1770413_R1")
        uploaded_files = st.file_uploader(
            "Choose FASTQ / FASTA files",
            type=["fastq", "fq", "gz", "fasta", "fa", "fna"],
            accept_multiple_files=True,
            help="R1/R2 can be uploaded separately; .gz files are accepted.",
        )
        submit_upload = st.form_submit_button("Save Dataset to Server", type="primary")

    if submit_upload:
        if not uploaded_files:
            st.warning("Choose at least one FASTQ / FASTA file.")
        else:
            for uploaded in uploaded_files:
                extension = identify_extension(uploaded.name)
                if not extension:
                    st.error(f"Unsupported file extension: {uploaded.name}")
                    continue
                stored = f"{uuid4()}{extension}"
                destination = DATASET_DIR / stored
                try:
                    uploaded.seek(0)
                    with destination.open("wb") as out:
                        while chunk := uploaded.read(1024 * 1024):
                            out.write(chunk)
                    if not validate_file(destination, extension):
                        raise ValueError("Invalid first FASTQ/FASTA record")
                    save_metadata(
                        dataset_name.strip() or uploaded.name,
                        uploaded.name,
                        stored,
                        "uploaded",
                        destination.stat().st_size,
                    )
                    st.success(f"Saved: {uploaded.name}")
                except (OSError, ValueError, sqlite3.Error) as exc:
                    destination.unlink(missing_ok=True)
                    st.error(f"Could not save {uploaded.name}: {exc}")

# TAB 2 — MANUAL SEQUENCES
with tab2:
    st.subheader("Enter DNA Sequence Manually")
    with st.form("manual_form", clear_on_submit=False):
        manual_name = st.text_input("Sequence Name", placeholder="e.g. KMP_repetitive_demo")
        manual_input = st.text_area(
            "Paste DNA sequence (or one FASTA record)",
            placeholder="ACACACAC",
            height=160,
        )
        submit_manual = st.form_submit_button("Save DNA Sequence", type="primary")

    if submit_manual:
        lines = [s.strip() for s in manual_input.splitlines() if s.strip()]
        headers = [s for s in lines if s.startswith(">")]
        seq = "".join(s for s in lines if not s.startswith(">")).upper()
        if len(headers) > 1:
            st.error("Please enter only one FASTA record.")
        elif not seq:
            st.warning("Enter a DNA sequence before saving.")
        elif set(seq) - DNA_ALPHABET:
            st.error("Sequence must contain only DNA/IUPAC nucleotide letters.")
        else:
            stored = f"{uuid4()}.fasta"
            destination = DATASET_DIR / stored
            try:
                destination.write_text(f">manual_sequence\n{seq}\n", encoding="utf-8")
                save_metadata(
                    manual_name.strip() or "Manual DNA Sequence", "manual_input",
                    stored, "manual", destination.stat().st_size,
                )
                st.success(f"Saved {len(seq):,}-base DNA sequence.")
            except (OSError, sqlite3.Error) as exc:
                destination.unlink(missing_ok=True)
                st.error(f"Save failed: {exc}")

# TAB 3 — QUERY CREATION AND CRUD
with tab3:
    st.subheader("Create Search Query Panel")
    st.caption("Enter one query per line: query_name,DNA_sequence (A/C/G/T only)")
    with st.form("query_form"):
        panel_name = st.text_input("Panel Name", value="", placeholder="Enter a new panel name")
        query_lines = st.text_area(
            "Search Patterns",
            value="",
            placeholder="short_pattern,CACCCG\nmedium_pattern,CACCCGGCATCA",
            height=140,
        )
        submit_query = st.form_submit_button("Save Query Panel", type="primary")

    if submit_query:
        parsed = []
        errors = []
        seen_names = set()
        for line_no, raw in enumerate(query_lines.splitlines(), start=1):
            if not raw.strip():
                continue
            if "," not in raw:
                errors.append(f"Line {line_no}: expected name,sequence")
                continue
            name, dna = (s.strip() for s in raw.split(",", 1))
            dna = dna.upper()
            if not name or not dna or set(dna) - QUERY_ALPHABET:
                errors.append(f"Line {line_no}: invalid name or DNA sequence")
            elif name in seen_names:
                errors.append(f"Line {line_no}: duplicate query name '{name}'")
            else:
                seen_names.add(name)
                parsed.append((name, dna))
        if not panel_name.strip():
            errors.append("Panel name cannot be empty.")
        if errors:
            for error in errors:
                st.error(error)
        elif not parsed:
            st.warning("No valid queries supplied.")
        else:
            with sqlite3.connect(DB_PATH) as conn:
                existing_names = {
                    row[0] for row in conn.execute(
                        "SELECT name FROM queries WHERE panel=?", (panel_name.strip(),)
                    )
                }
                overlap = sorted(existing_names & seen_names)
                if overlap:
                    st.error(f"Query names already in this panel: {', '.join(overlap)}. Update existing queries below instead.")
                else:
                    conn.executemany(
                        "INSERT INTO queries (id,panel,name,sequence,created_at) VALUES(?,?,?,?,?)",
                        [(str(uuid4()), panel_name.strip(), name, dna, now_utc())
                         for name, dna in parsed],
                    )
                    st.success(f"✅ Saved {len(parsed)} queries to '{panel_name.strip()}'. Check Saved Search Queries below.")

    st.divider()
    render_query_manager(DB_PATH)

# TAB 4 — DATASET CRUD + PREVIEW
with tab4:
    st.subheader("Saved Dataset Library")
    render_dataset_manager(DB_PATH, DATASET_DIR)
    st.divider()
    datasets = fetch_all("SELECT * FROM datasets ORDER BY created_at DESC")
    if datasets:
        picked = st.selectbox(
            "Select dataset to preview",
            range(len(datasets)),
            format_func=lambda i: f"{datasets[i]['name']} — {datasets[i]['original_filename']}",
            key="dataset_preview_select",
        )
        selected = datasets[picked]
        base = DATASET_DIR.resolve()
        path = (base / selected["stored_filename"]).resolve()
        if path.parent != base or not path.is_file():
            st.warning("Stored file is missing or outside data/datasets.")
        else:
            try:
                with open_sequence_file(path) as handle:
                    preview = "".join(handle.readline() for _ in range(12))
                st.code(preview, language="text")
            except (OSError, UnicodeError, EOFError) as exc:
                st.error(f"Could not preview file: {exc}")

# TAB 5 — ALGORITHM ENGINE (unchanged)
with tab5:
    render_search_page(fetch_all, DATASET_DIR, DATA_DIR)

st.divider()
st.caption("OncoSearch | Exact DNA searching for education and reproducible benchmarking; not clinical variant calling.")
