"""Streamlit search and algorithm-demo panel for the existing OncoSearch app."""
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
import csv
import io
import json
import time
import html
import random
import statistics

import pandas as pd
import streamlit as st

from teaching_studio import render_teaching_studio

from search_engine import (
    METHODS, COMPLEXITY, Matcher, Target, lps_table, prepare_targets,
    load_sample, run_benchmark, trace_steps
)


def _csv_preview(rows):
    if not rows:
        return ""
    f = io.StringIO()
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    return f.getvalue()


def _draw_sequence(seq, step, pattern):
    """Small bounded text strip. User-supplied DNA is HTML-escaped."""
    highlight = step.get("pos", -1)
    align = step.get("align")
    parts = []
    for i, base in enumerate(seq):
        focused = i == highlight
        in_window = align is not None and align <= i < align + len(pattern)
        completed_hit = any(h["start"] <= i < h["end"]
                            for h in step.get("hits", []))
        bg = "#15803d" if completed_hit else (
            ("#166534" if step.get("match") else "#991b1b") if focused else (
                "#334155" if in_window else "#1e293b"))
        border = "2px solid #e2e8f0" if focused else "1px solid #475569"
        parts.append(
            f'<span title="position {i}" style="display:inline-block;padding:5px 6px;'
            f'font-family:monospace;background:{bg};border:{border};'
            f'color:#fff;border-radius:4px;margin:2px">{html.escape(base)}</span>'
        )
    st.markdown('<div style="line-height:2.9">' + ''.join(parts) + '</div>',
                unsafe_allow_html=True)
    st.caption(f"Highlighted position: {highlight} (0-based). Gray = current pattern window.")
    st.code(step["note"], language="text")


def _show_match_window(seq, hit):
    """Display an exact-match window with context and 0-based coordinates."""
    start, end = int(hit["start"]), int(hit["end"])
    if not (0 <= start < end <= len(seq)):
        st.warning("Match coordinate is outside the displayed sequence.")
        return
    left = max(0, start - 12)
    right = min(len(seq), end + 12)
    cells = []
    for i in range(left, right):
        in_match = start <= i < end
        background = "#15803d" if in_match else "#334155"
        border = "2px solid #bbf7d0" if in_match else "1px solid #64748b"
        cells.append(
            '<span title="0-based index ' + str(i) + '" style="'
            'display:inline-block;min-width:22px;text-align:center;'
            'margin:2px;padding:6px 4px;border-radius:5px;'
            f'font-family:monospace;font-weight:700;color:white;'
            f'background:{background};border:{border}">'
            + html.escape(seq[i]) + '</span>'
        )
    st.success(f"🎯 EXACT MATCH — {hit['name']} ({hit['strand']}) "
               f"at position {start}")
    st.markdown('<div style="line-height:3.0">' + ''.join(cells) + '</div>',
                unsafe_allow_html=True)
    st.caption(f"Pattern: {hit['sequence']} | Match: [{start}, {end}) "
               f"(0-based, end exclusive) | Displayed region: [{left}, {right})")


@st.dialog("🎯 Exact DNA Match")
def _open_match_popup(seq, hit):
    _show_match_window(seq, hit)


def _show_history(result_dir):
    with st.expander("📚 Saved experiment history", expanded=False):
        history = sorted(result_dir.glob("search_*.json"), reverse=True)
        if not history:
            st.caption("No saved searches yet.")
            return
        chosen = st.selectbox("Previously saved run", history,
                              format_func=lambda p: p.name, key="saved_search_select")
        data = json.loads(chosen.read_text(encoding="utf-8"))
        st.write("**Created:**", data.get("created_at"))
        st.write("**Input datasets:**", ', '.join(data.get("files", [])))
        st.dataframe(pd.DataFrame(data["metrics"])[[
            "algorithm", "median_search_ms", "preprocessing_ms",
            "total_matches", "correct_vs_first"]], hide_index=True)
        st.download_button("Download saved JSON", data=chosen.read_bytes(),
                           file_name=chosen.name, mime="application/json",
                           key="download_saved_run")

        # Saved search reports are immutable results; allow deletion after confirmation.
        if st.button("🗑️ Delete selected experiment", key="delete_selected_run"):
            st.session_state["pending_run_delete"] = chosen.name
            st.rerun()
        if st.session_state.get("pending_run_delete") == chosen.name:
            st.warning(f"Permanently delete saved report '{chosen.name}'?")
            yes, no = st.columns(2)
            if yes.button("Confirm delete report", type="primary", key="confirm_run_delete"):
                try:
                    # chosen is selected from result_dir.glob, not user-provided path.
                    chosen.unlink()
                except OSError as exc:
                    st.error(f"Cannot delete report: {exc}")
                else:
                    st.session_state.pop("pending_run_delete", None)
                    st.rerun()
            if no.button("Cancel", key="cancel_run_delete"):
                st.session_state.pop("pending_run_delete", None)
                st.rerun()


def render_search_page(fetch_all, dataset_dir, data_dir):
    st.subheader("🔎 Multi-Pattern Search & Algorithm Laboratory")
    st.caption("Search saved FASTQ/FASTA or pasted DNA with five techniques. "
               "Results and benchmark evidence can be stored on the BioAI server.")
    result_dir = Path(data_dir) / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    datasets = fetch_all("SELECT * FROM datasets ORDER BY created_at DESC")
    panels = fetch_all("SELECT DISTINCT panel FROM queries ORDER BY panel")
    if not datasets or not panels:
        st.warning("Save at least one dataset and one query panel in the earlier tabs.")
        _show_history(result_dir)
        return

    dataset_choices = {f"{i+1}. {d['name']} — {d['original_filename']}": d
                       for i, d in enumerate(datasets)}
    selected_labels = st.multiselect("Select one or more datasets (R1/R2 allowed)",
                                     list(dataset_choices),
                                     default=list(dataset_choices)[:1])
    panel_name = st.selectbox("Saved query panel", [p["panel"] for p in panels])
    queries = fetch_all("SELECT name, sequence FROM queries WHERE panel = ? ORDER BY rowid",
                        (panel_name,))
    st.dataframe(pd.DataFrame(queries), hide_index=True, use_container_width=True)

    c1, c2, c3 = st.columns(3)
    with c1:
        max_reads = st.number_input("Max records per file (0 = all, subject to memory safeguards)",
                                    min_value=0, max_value=250000, value=1000, step=100)
    with c2:
        min_quality = st.slider("Minimum mean Phred (FASTQ)", 0, 40, 0)
    with c3:
        both = st.checkbox("Search reverse complement", value=True)

    methods = st.multiselect("Compare these algorithms", METHODS, default=METHODS)
    repeats = st.slider("Timed runs for median", 1, 5, 3)
    if max_reads == 0:
        st.warning("Full-file mode loads accepted reads into memory on the BioAI server. "
                   "A 40-million-base safety limit prevents excessive memory use. "
                   "For large sequencing studies, use bounded subsets or a future streaming mode.")
    st.caption("Search is exact (0 mismatches). N is literal, not a wildcard. "
               "Each FASTQ read is scanned separately. R1/R2 are not merged.")

    if st.button("▶ Run DNA Search & Benchmark", type="primary"):
        if not selected_labels or not methods:
            st.error("Select at least one dataset and one algorithm.")
        else:
            try:
                targets = prepare_targets([(q["name"], q["sequence"]) for q in queries], both)
                paths = []
                base = Path(dataset_dir).resolve()
                for label in selected_labels:
                    d = dataset_choices[label]
                    p = (base / d["stored_filename"]).resolve()
                    if not p.is_relative_to(base) or not p.is_file():
                        raise ValueError(f"Stored file is missing or invalid: {label}")
                    paths.append((d["original_filename"], p))
                with st.spinner("Loading read samples and benchmarking exact matches..."):
                    reads, stats = load_sample(paths, int(max_reads), min_quality)
                    if not reads:
                        raise ValueError("No reads passed filters; lower minimum Phred threshold.")
                    metrics, previews = run_benchmark(reads, targets, methods,
                                                      repetitions=repeats)
                st.session_state["last_search_reads_v3"] = reads
                report = {
                    "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "files": [x[0] for x in paths],
                    "panel": panel_name,
                    "targets": [asdict(t) for t in targets],
                    "sampling": stats,
                    "methods": methods,
                    "metrics": metrics,
                    "preview_hits": previews,
                    "example_read": next((r["example_read"] for r in metrics
                                          if r["example_read"]), reads[0]["seq"]),
                }
                result_path = result_dir / f"search_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid4().hex[:7]}.json"
                result_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
                st.session_state["last_search_report"] = report
                st.success(f"Search complete. Result saved on server: {result_path.name}")
            except Exception as ex:
                st.error(f"Search failed: {ex}")

    report = st.session_state.get("last_search_report")
    if report:
        st.divider()
        st.markdown("### 📊 Search Results & Performance")
        stats = report["sampling"]
        a, b, c, d = st.columns(4)
        a.metric("Reads analyzed", f"{stats['used_reads']:,}")
        b.metric("Bases scanned", f"{stats['bases']:,}")
        c.metric("Filtered reads", f"{stats['filtered_reads']:,}")
        d.metric("Target variants (incl. RC)", len(report["targets"]))

        df = pd.DataFrame(report["metrics"])
        st.dataframe(df[["algorithm", "strategy", "preprocessing_ms", "median_search_ms",
                         "total_with_setup_ms", "total_matches", "reads_with_hits",
                         "operations", "operation_type", "correct_vs_first"]],
                     hide_index=True, use_container_width=True)
        st.bar_chart(df, x="algorithm", y="median_search_ms", horizontal=True)
        if not df["correct_vs_first"].all():
            st.error("Algorithm disagreement detected. Inspect implementations!")
        else:
            st.success("All selected algorithms returned identical read-level match coordinates.")

        st.caption("Adaptive (ours) is a transparent course-level heuristic: ≥4 queries → Aho; "
                   "otherwise a strongly self-overlapping pattern → KMP; else Horspool. "
                   "No theoretical novelty or guaranteed speedup is claimed.")
        st.caption("Runtime = median of sequential search-only runs on already-parsed reads; "
                   "preprocessing measured separately, no FASTQ disk I/O included. "
                   "Python built-in operations are not instrumented. "
                   "Aho counts transition checks, not character equality comparisons.")

        st.markdown("### 📚 Theory vs Measured Work")
        theory = [{"Algorithm": m, "Preprocessing": COMPLEXITY[m][0],
                   "Search (worst-case unless noted)": COMPLEXITY[m][1],
                   "Counted operation": COMPLEXITY[m][2]}
                  for m in report["methods"]]
        st.dataframe(pd.DataFrame(theory), hide_index=True, use_container_width=True)
        st.caption("N = total read bases; P = searched patterns including reverse complements; "
                   "M = maximum pattern length; L = total pattern lengths; Z = matches. "
                   "Python str.find is CPython-optimized; we do not assign it a simplistic worst-case bound.")

        summary = report["metrics"][0]
        target_df = pd.DataFrame([
            {"Query": tgt["name"], "Strand": tgt["strand"],
             "Pattern": tgt["sequence"], "Match count": summary["counts"][i]}
            for i, tgt in enumerate(report["targets"])
        ])
        st.markdown("### 🎯 Matches by Target")
        st.dataframe(target_df, hide_index=True, use_container_width=True)
        selected_method = st.selectbox("Inspect hits for algorithm", report["methods"],
                                       key="inspection_method")
        preview = report["preview_hits"][selected_method]
        st.dataframe(pd.DataFrame(preview), hide_index=True, use_container_width=True)
        if preview:
            st.download_button("⬇ Download preview hits CSV", _csv_preview(preview),
                               "oncosearch_preview_hits.csv", mime="text/csv")
        st.info("Read-level exact motifs are not validated clinical variant calls. "
                "Reference alignment, quality at the matched bases and biological context "
                "are required for variant interpretation.")

        st.divider()
        render_teaching_studio(report, st.session_state.get("last_search_reads_v3"))

    st.divider()
    with st.expander("🧪 Random vs Repetitive vs Worst-Case Stress Test", expanded=False):
        st.caption("Controlled synthetic DNA complements real FASTQ results. "
                   "A fixed random seed makes results reproducible.")
        workload = st.selectbox("Sequence type", ["Random A/C/G/T", "Repetitive AC", "Adversarial AAA…AC"])
        if st.button("Run synthetic benchmark"):
            samples = []
            for n in [100, 1000, 5000, 20000]:
                if workload == "Random A/C/G/T":
                    rng = random.Random(17 + n)
                    seq = ''.join(rng.choices("ACGT", k=n))
                    pattern = "ACGTACGT"
                elif workload == "Repetitive AC":
                    seq = ("AC" * (n // 2 + 1))[:n]
                    pattern = "ACACACAC"
                else:
                    seq = "A" * (n - 1) + "C"
                    pattern = "AAAAAAC"
                target = [Target("synthetic", pattern)]
                for method in METHODS:
                    engine = Matcher(method, target)
                    times = []
                    for _ in range(3):
                        start = time.perf_counter_ns()
                        hits, ops = engine.scan(seq)
                        times.append((time.perf_counter_ns() - start) / 1e6)
                    samples.append({"DNA length": n, "Algorithm": method,
                                    "Time (ms)": round(statistics.median(times), 4),
                                    "Operations": ops, "Matches": len(hits)})
            st.session_state["synthetic_results"] = samples
        if "synthetic_results" in st.session_state:
            sample_df = pd.DataFrame(st.session_state["synthetic_results"])
            st.line_chart(sample_df.pivot(index="DNA length", columns="Algorithm",
                                          values="Time (ms)"))
            st.dataframe(sample_df, hide_index=True, use_container_width=True)
            st.caption("Comparison counts are only directly comparable for Naive/KMP/Horspool. "
                       "Aho counts transition checks; str.find does not expose internal counts.")

    _show_history(result_dir)
