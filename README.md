# 🧬 OncoSearch v2

**Interactive DNA Pattern Matching and Algorithm Visualization**

**CSCE 5150 — Analysis of Computer Algorithms | University of North Texas**

OncoSearch is a Streamlit-based web application for searching DNA patterns and primer sequences in real FASTQ/FASTA datasets.

It combines classical string-matching algorithms, performance benchmarking, and interactive step-by-step visualization to demonstrate how DNA searching techniques work.

## ✨ Features

- **FASTQ/FASTA Support:** Upload real sequencing datasets, including compressed files.
- **Manual DNA Input:** Enter custom DNA sequences for testing.
- **Query Management:** Create, update, and delete DNA patterns and query panels.
- **Multi-Pattern Search:** Search multiple sequences simultaneously.
- **Reverse Complement:** Search both DNA strand orientations.
- **Quality Filtering:** Filter sequencing reads using mean Phred scores.
- **Algorithm Benchmarking:** Compare runtime, preprocessing, operations, and exact matches.
- **Interactive Teaching Studio:** Visualize comparisons, pattern movement, mismatches, and algorithm-specific shifts.
- **Match Detection:** Real-time match notifications, highlighted sliding windows, and detailed popups.
- **Experiment History:** Save benchmark summaries and results.

## 🧠 Implemented Algorithms

| Algorithm | Technique | Time Complexity |
|---|---|---|
| Naïve Search | Character-by-character matching | O(N × M) |
| KMP | Longest Prefix Suffix (LPS) | O(N + M) |
| Boyer–Moore–Horspool | Bad-character shifting | O(N × M) worst case |
| Aho–Corasick | Trie and failure links | O(N + Z) search |
| Python `str.find()` | Optimized built-in baseline | Implementation-dependent |
| Adaptive (Ours) | Pattern-dependent algorithm selection | Depends on selected algorithm |

Here, N represents sequence length, M represents pattern length, and Z represents reported matches. The listed single-pattern complexities exclude additional costs for searching multiple patterns separately.

**Adaptive (Ours)** is a rule-based selection heuristic, not a new theoretical string-matching algorithm.

## 🎬 Interactive Algorithm Visualization

The Algorithm Teaching Studio demonstrates how each algorithm searches DNA sequences.

Features include:

- Moving pattern alignment beneath the DNA sequence
- Character-by-character comparison
- Match and mismatch highlighting
- KMP LPS fallback visualization
- Horspool shift visualization
- Aho–Corasick trie-state transitions
- Play, Pause, Next, Previous, and speed controls
- Live operation and match counters
- Exact-match popup with surrounding DNA context
- Navigation between multiple matching positions

The visualization supports full selected-read traces, subject to memory safeguards.

## 🛠️ Installation

### 1. Clone the Repository

```bash
git clone https://github.com/YOUR-USERNAME/OncoSearch.git
cd OncoSearch
```

### 2. Create a Python Environment

Python 3.11 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
python -m pip install -r requirements.txt
```

### 4. Run the Application

```bash
streamlit run app.py
```

Open the local URL shown by Streamlit, usually:

http://localhost:8501

## 🔬 How to Use

1. Upload a FASTQ/FASTA dataset or enter a DNA sequence manually.
2. Create a search query panel containing one or more DNA patterns.
3. Select the dataset, query panel, and searching algorithms.
4. Configure quality filtering and reverse-complement searching.
5. Click **Run DNA Search & Benchmark**.
6. Explore match counts, matching positions, runtime, and operation counts.
7. Open the **Algorithm Teaching Studio** to visualize the searching process.

### Example

**DNA sequence:**

```text
CACCCGCACCCG
```

**Search query:**

```text
example_pattern,CACCCG
```

**Expected matches (0-based):**

```text
Position 0
Position 6
```

The Teaching Studio highlights both matches and demonstrates how different algorithms locate them.

## 📊 Experimental Evaluation

OncoSearch supports evaluation using:

- Real FASTQ sequencing reads
- Random DNA sequences
- Repetitive DNA sequences
- Adversarial DNA patterns
- Different sequence lengths and query counts

Evaluation metrics include:

- Exact match count and positions
- Execution time
- Preprocessing time
- Character comparisons or state-transition checks
- Algorithm agreement
- Input-size scalability

## 📁 Project Structure

```text
OncoSearch/
├── app.py
├── crud_ui.py
├── search_engine.py
├── search_page.py
├── teaching_trace.py
├── teaching_studio.py
├── test_search_engine.py
├── test_match_popup_events.py
├── test_crud.py
├── test_teaching_studio.py
├── requirements.txt
├── .gitignore
└── README.md
```

## 🧪 Run Tests

```bash
python -m unittest -q \
  test_search_engine.py \
  test_match_popup_events.py \
  test_crud.py \
  test_teaching_studio.py
```

## 🚀 Deployment

This application can be deployed using Streamlit Community Cloud.

1. Push the project source code to GitHub.
2. Visit https://share.streamlit.io/
3. Select your GitHub repository and branch.
4. Set the main file path to `app.py`.
5. Deploy the application.

**Important:** Do not upload private FASTQ datasets, SQLite databases, credentials, or other sensitive files to GitHub.

The current application supports uploads and deletions without per-user authentication. Public deployment should use safe demo data and restrict destructive actions. Streamlit Community Cloud's local filesystem is not durable storage; persistent multi-user deployment requires an external database and file storage.

## ⚠️ Limitations

- Only exact DNA pattern matching is supported; mismatches are not currently allowed.
- Primer-sequence matches do not establish experimentally validated primer binding.
- Read-level matching positions are not genomic reference coordinates.
- Sequence matches are not evidence of a confirmed cancer mutation or diagnosis.
- Quality filtering uses the average quality of the entire sequencing read.
- Large datasets and visual traces are subject to memory safeguards.
- Python's built-in `str.find()` does not expose its internal comparison steps.

## 🎓 Academic Information

**Course:** CSCE 5150 — Analysis of Computer Algorithms  
**Institution:** University of North Texas  
**Project:** OncoSearch — Interactive DNA Pattern Matching and Algorithm Visualization

---

**Disclaimer:** OncoSearch is developed for educational and algorithmic research purposes. It is not intended for clinical diagnosis.
