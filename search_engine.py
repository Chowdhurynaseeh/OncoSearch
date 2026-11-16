"""Exact nucleotide pattern-matching engines (no external matching libraries)."""
from collections import deque
from dataclasses import dataclass
from pathlib import Path
import gzip
import hashlib
import statistics
import time

METHODS = ["Python str.find", "Naive", "KMP", "Horspool", "Aho-Corasick", "Adaptive (ours)"]

COMPLEXITY = {
    "Python str.find": ("Internal optimized routine", "Implementation-dependent", "N/A"),
    "Naive": ("O(1)", "O(P × N × M) worst", "Character comparisons"),
    "KMP": ("O(L)", "O(P × N + Z)", "Character comparisons"),
    "Horspool": ("O(L)", "O(P × N × M) worst", "Character comparisons"),
    "Aho-Corasick": ("O(L) typical", "O(N + Z) after construction", "State-transition checks"),
    "Adaptive (ours)": ("O(L) typical", "Depends on chosen underlying algorithm", "Underlying engine operations"),
}
# N = total bases scanned (sum of read lengths); P = number of patterns;
# M = max pattern length; L = sum of pattern lengths; Z = reported matches.


def reverse_complement(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def lps_table(p):
    lps = [0] * len(p)
    j = 0
    for i in range(1, len(p)):
        while j and p[i] != p[j]:
            j = lps[j - 1]
        if p[i] == p[j]:
            j += 1
            lps[i] = j
    return lps


@dataclass(frozen=True)
class Target:
    name: str
    sequence: str
    strand: str = "+"


def prepare_targets(rows, both_strands=False):
    result = []
    for name, seq in rows:
        s = seq.strip().upper()
        if not s or set(s) - set("ACGT"):
            raise ValueError(f"Query {name!r} must contain only A/C/G/T")
        result.append(Target(name, s, "+"))
        rc = reverse_complement(s)
        if both_strands and rc != s:
            result.append(Target(name, rc, "-"))
    if not result:
        raise ValueError("Choose at least one saved query.")
    return result


class Matcher:
    """Preprocess patterns once and re-use them across all FASTQ reads.

    scan(read) -> ([(target_index, zero_based_position), ...], operation_count)
    The operation_count is None for native str.find.
    """
    def __init__(self, method, targets):
        if method not in METHODS:
            raise ValueError(method)
        self.method = method
        self.targets = targets
        self.patterns = [t.sequence for t in targets]
        self.tables = []
        self.trie = None
        self.strategy = method
        self.delegate = None
        if method == "Adaptive (ours)":
            # Course-level heuristic, not a novel string-matching theorem.
            if len(self.patterns) >= 4:
                choice = "Aho-Corasick"
            elif any(max(lps_table(p), default=0) >= len(p) // 2 and len(p) >= 4
                     for p in self.patterns):
                choice = "KMP"
            else:
                choice = "Horspool"
            self.delegate = Matcher(choice, targets)
            self.strategy = f"Selected {choice}"
        elif method == "KMP":
            self.tables = [lps_table(p) for p in self.patterns]
        elif method == "Horspool":
            self.tables = [self._shift_table(p) for p in self.patterns]
        elif method == "Aho-Corasick":
            self.trie = self._build_trie(self.patterns)

    @staticmethod
    def _shift_table(p):
        return {c: len(p) - 1 - i for i, c in enumerate(p[:-1])}

    @staticmethod
    def _build_trie(patterns):
        nodes = [{"next": {}, "fail": 0, "out": []}]
        for idx, pattern in enumerate(patterns):
            state = 0
            for ch in pattern:
                if ch not in nodes[state]["next"]:
                    nodes[state]["next"][ch] = len(nodes)
                    nodes.append({"next": {}, "fail": 0, "out": []})
                state = nodes[state]["next"][ch]
            nodes[state]["out"].append(idx)
        queue = deque()
        for state in nodes[0]["next"].values():
            queue.append(state)
        while queue:
            parent = queue.popleft()
            for ch, child in nodes[parent]["next"].items():
                queue.append(child)
                f = nodes[parent]["fail"]
                while f and ch not in nodes[f]["next"]:
                    f = nodes[f]["fail"]
                nodes[child]["fail"] = nodes[f]["next"].get(ch, 0)
                nodes[child]["out"].extend(nodes[nodes[child]["fail"]]["out"])
        return nodes

    def scan(self, text):
        if self.delegate is not None:
            return self.delegate.scan(text)
        if self.method == "Aho-Corasick":
            return self._aho(text)
        results = []
        operations = 0
        for idx, pattern in enumerate(self.patterns):
            if self.method == "Python str.find":
                positions, ops = self._builtin(text, pattern)
            elif self.method == "Naive":
                positions, ops = self._naive(text, pattern)
            elif self.method == "KMP":
                positions, ops = self._kmp(text, pattern, self.tables[idx])
            else:
                positions, ops = self._horspool(text, pattern, self.tables[idx])
            results.extend((idx, pos) for pos in positions)
            if ops is not None:
                operations += ops
        results.sort()
        return results, None if self.method == "Python str.find" else operations

    @staticmethod
    def _builtin(t, p):
        out = []
        start = 0
        while True:
            pos = t.find(p, start)
            if pos < 0:
                return out, None
            out.append(pos)
            start = pos + 1  # retain overlapping occurrences

    @staticmethod
    def _naive(t, p):
        n, m, result, ops = len(t), len(p), [], 0
        for i in range(n - m + 1):
            for j in range(m):
                ops += 1
                if t[i + j] != p[j]:
                    break
            else:
                result.append(i)
        return result, ops

    @staticmethod
    def _kmp(t, p, lps):
        result, ops = [], 0
        i = j = 0
        while i < len(t):
            ops += 1
            if t[i] == p[j]:
                i += 1
                j += 1
                if j == len(p):
                    result.append(i - j)
                    j = lps[j - 1]
            elif j:
                j = lps[j - 1]
            else:
                i += 1
        return result, ops

    @staticmethod
    def _horspool(t, p, shifts):
        n, m, result, ops = len(t), len(p), [], 0
        i = 0
        while i <= n - m:
            j = m - 1
            while j >= 0:
                ops += 1
                if t[i + j] != p[j]:
                    break
                j -= 1
            if j < 0:
                result.append(i)
            i += shifts.get(t[i + m - 1], m)
        return result, ops

    def _aho(self, t):
        nodes = self.trie
        state = checks = 0
        result = []
        for pos, ch in enumerate(t):
            while state and ch not in nodes[state]["next"]:
                checks += 1
                state = nodes[state]["fail"]
            checks += 1
            state = nodes[state]["next"].get(ch, 0)
            for idx in nodes[state]["out"]:
                result.append((idx, pos - len(self.patterns[idx]) + 1))
        result.sort()
        return result, checks


def iter_records(path):
    """Yield (read_id, sequence, avg_phred_or_None) without loading entire file."""
    path = Path(path)
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    name = path.name.lower()
    with opener(path, "rt", encoding="utf-8") as fh:
        if any(x in name for x in (".fastq", ".fq")):
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip().upper()
                plus = fh.readline()
                quality = fh.readline().rstrip("\r\n")
                if not (header.startswith("@") and plus.startswith("+") and seq
                        and len(seq) == len(quality)):
                    raise ValueError(f"Malformed FASTQ record near {header[:50]!r}")
                avg = sum(ord(c) - 33 for c in quality) / len(quality)
                yield header[1:].strip().split()[0], seq, avg
        elif any(x in name for x in (".fasta", ".fa", ".fna")):
            rid, fragments, record_len = None, [], 0
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if rid is not None:
                        yield rid, "".join(fragments).upper(), None
                    rid, fragments, record_len = line[1:].split()[0], [], 0
                elif line:
                    if rid is None:
                        raise ValueError("FASTA must start with >header")
                    record_len += len(line)
                    if record_len > 5_000_000:
                        raise ValueError("FASTA record exceeds 5 Mb; split into smaller regions for this demo")
                    fragments.append(line)
            if rid is not None:
                yield rid, "".join(fragments).upper(), None
        else:
            raise ValueError(f"Unsupported sequence format: {path.name}")


def load_sample(dataset_paths, max_per_file=1000, min_mean_q=0):
    reads, seen, filtered, total_bases = [], 0, 0, 0
    for label, path in dataset_paths:
        for index, (rid, seq, quality) in enumerate(iter_records(path)):
            if max_per_file > 0 and index >= max_per_file:
                break
            seen += 1
            if quality is not None and quality < min_mean_q:
                filtered += 1
                continue
            if total_bases + len(seq) > 40_000_000:
                raise ValueError("This in-memory benchmark has a 40-million-base safety limit. Set a smaller Max records value.")
            reads.append({"source": label, "id": rid, "seq": seq,
                          "mean_q": None if quality is None else round(quality, 2)})
            total_bases += len(seq)
    return reads, {"considered_reads": seen, "filtered_reads": filtered,
                   "used_reads": len(reads), "bases": total_bases}


def run_benchmark(reads, targets, methods=METHODS, repetitions=3, hit_preview_limit=150):
    results = []
    previews = {}
    for method in methods:
        start = time.perf_counter_ns()
        engine = Matcher(method, targets)
        preprocess_ms = (time.perf_counter_ns() - start) / 1e6
        # Warm up on identical reads; do not include warm-up in the timed samples.
        for read in reads[:min(10, len(reads))]:
            engine.scan(read["seq"])
        trials = []
        for _ in range(repetitions):
            t0 = time.perf_counter_ns()
            total_hits = 0
            for read in reads:
                hits, _ = engine.scan(read["seq"])
                total_hits += len(hits)
            trials.append((time.perf_counter_ns() - t0) / 1e6)
        # Validation and hit preview are outside the timed benchmark.
        digest = hashlib.sha256()
        ops_total = 0
        counts = [0] * len(targets)
        preview = []
        matching_read_count = 0
        example_read = None
        for r_idx, read in enumerate(reads):
            hits, operations = engine.scan(read["seq"])
            ops_total += operations or 0
            if hits:
                matching_read_count += 1
                if example_read is None:
                    example_read = read["seq"]
            for target_idx, position in hits:
                counts[target_idx] += 1
                digest.update(f"{r_idx}:{target_idx}:{position};".encode("ascii"))
                if len(preview) < hit_preview_limit:
                    preview.append({
                        "read_id": read["id"], "file": read["source"],
                        "query": targets[target_idx].name,
                        "strand": targets[target_idx].strand,
                        "position_0based": position,
                        "matched_dna": targets[target_idx].sequence,
                        "mean_read_phred": read["mean_q"],
                    })
        results.append({
            "algorithm": method,
            "strategy": engine.strategy,
            "median_search_ms": round(statistics.median(trials), 4),
            "preprocessing_ms": round(preprocess_ms, 4),
            "total_with_setup_ms": round(statistics.median(trials) + preprocess_ms, 4),
            "total_matches": sum(counts),
            "reads_with_hits": matching_read_count,
            "operations": None if method == "Python str.find" else ops_total,
            "operation_type": COMPLEXITY[method][2],
            "signature": digest.hexdigest(),
            "counts": counts,
            "example_read": example_read,
        })
        previews[method] = preview
    reference_signature = results[0]["signature"] if results else ""
    for row in results:
        row["correct_vs_first"] = row["signature"] == reference_signature
    return results, previews


def trace_steps(text, targets, method, max_events=140):
    """Bounded per-character trace; `hits` lists *complete* exact matches."""
    text = text.upper()
    engine = Matcher(method, targets)
    if method == "Adaptive (ours)":
        return trace_steps(text, targets, engine.delegate.method,
                           max_events=max_events)

    steps = []

    def hit(target_idx, start):
        target = targets[target_idx]
        return {"name": target.name, "sequence": target.sequence,
                "strand": target.strand, "start": start,
                "end": start + len(target.sequence)}

    if method == "Aho-Corasick":
        node, state = engine.trie, 0
        for i, char in enumerate(text):
            failures = 0
            while state and char not in node[state]["next"]:
                state = node[state]["fail"]
                failures += 1
            state = node[state]["next"].get(char, 0)
            found = [hit(k, i - len(targets[k].sequence) + 1)
                     for k in node[state]["out"]]
            steps.append({"pos": i, "align": None, "match": bool(found),
                          "hits": found,
                          "note": f"Read {char}; state={state}; "
                                  f"failure jumps={failures}; "
                                  f"complete matches={len(found)}"})
            if len(steps) >= max_events:
                break
        return steps

    if method == "Python str.find":
        matches, _ = engine.scan(text)
        for target_idx, pos in matches[:max_events]:
            steps.append({"pos": pos, "align": pos, "match": True,
                          "hits": [hit(target_idx, pos)],
                          "note": f"str.find found complete match at {pos}; "
                                  "native comparisons are not exposed."})
        return steps

    p = targets[0].sequence
    if method == "Naive":
        for i in range(len(text) - len(p) + 1):
            for j in range(len(p)):
                ok = text[i + j] == p[j]
                found = [hit(0, i)] if ok and j == len(p) - 1 else []
                note = (f"Compare text[{i+j}]={text[i+j]} "
                        f"to pattern[{j}]={p[j]}: "
                        f"{'match' if ok else 'mismatch'}")
                if found:
                    note += f" — COMPLETE MATCH at {i}"
                steps.append({"pos": i + j, "pattern_pos": j,
                              "align": i, "match": ok,
                              "hits": found, "note": note})
                if len(steps) >= max_events:
                    return steps
                if not ok:
                    break

    elif method == "KMP":
        table = lps_table(p)
        i = j = 0
        while i < len(text) and len(steps) < max_events:
            compared_i, compared_j = i, j
            ok = text[i] == p[j]
            align = i - j
            found = []
            note = f"Compare text[{i}]={text[i]} vs pattern[{j}]={p[j]}"
            if ok:
                i += 1
                j += 1
                if j == len(p):
                    found = [hit(0, i - j)]
                    note += (f" — COMPLETE MATCH at {i-j}; "
                             f"set j=LPS[{j-1}]={table[j-1]}")
                    j = table[j-1]
                else:
                    note += f" — advance to (i={i}, j={j})"
            elif j:
                note += (f" — mismatch, jump j={j} → "
                         f"LPS[{j-1}]={table[j-1]} (no text restart)")
                j = table[j-1]
            else:
                i += 1
                note += f" — mismatch at j=0, advance text to i={i}"
            steps.append({"pos": compared_i, "pattern_pos": compared_j,
                          "align": align, "match": ok, "hits": found,
                          "note": note})

    elif method == "Horspool":
        i, m = 0, len(p)
        shifts = engine.tables[0]
        while i <= len(text) - m:
            j = m - 1
            while j >= 0:
                ok = text[i + j] == p[j]
                found = [hit(0, i)] if ok and j == 0 else []
                note = (f"Compare text[{i+j}]={text[i+j]} "
                        f"vs pattern[{j}]={p[j]}")
                if found:
                    note += f" — COMPLETE MATCH at {i}"
                steps.append({"pos": i + j, "pattern_pos": j,
                              "align": i, "match": ok,
                              "hits": found, "note": note})
                if len(steps) >= max_events:
                    return steps
                if not ok:
                    break
                j -= 1
            i += shifts.get(text[i + m - 1], m)
            steps[-1]["note"] += f"; shift alignment to {i}"
    return steps
