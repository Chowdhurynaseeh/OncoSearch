"""Deterministic teaching events for exact DNA string matching.

No searching dependency beyond our own search_engine. This module is pure Python
and can be unit-tested without Streamlit. The visualizer never invents matches:
all final hit coordinates must agree with Matcher.scan().
"""
from collections import deque

from search_engine import Matcher, Target, lps_table


DEFAULT_MAX_EVENTS = 100_000  # Safety cap, not a 100-step teaching limit.


def _hit(target, start):
    return {
        "name": target.name, "strand": target.strand,
        "sequence": target.sequence,
        "start": start, "end": start + len(target.sequence),
    }


def build_teaching_trace(text, targets, method, max_events=DEFAULT_MAX_EVENTS):
    """Return fully specified animation events for a single selected read.

    For Naive/KMP/Horspool traces, pass one target. Aho-Corasick accepts all
    targets and displays actual automaton state transitions and failure links.
    `max_events` is an explicit browser-memory safeguard. `complete` indicates
    whether we reached the end without truncating. Event counters track *search*
    operations, not preprocessing operations.
    """
    text = text.strip().upper()
    if set(text) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("DNA contains unsupported characters")
    if not targets:
        raise ValueError("Select at least one target")
    if not all(set(t.sequence) <= set("ACGT") and t.sequence for t in targets):
        raise ValueError("Target sequences must contain only A/C/G/T")
    requested_method = method
    matcher = Matcher(method, targets)
    effective_method = matcher.delegate.method if matcher.delegate else method
    if effective_method not in ("Aho-Corasick", "Python str.find") and len(targets) != 1:
        raise ValueError("Choose one target for this algorithm's character trace")

    steps = []
    found = []
    ops = 0
    truncated = False

    def add(kind, note, *, pos=None, align=None, p_index=None, ok=None,
            shift=None, state=None, prefix="", fail_from=None, fail_to=None,
            hits=None):
        nonlocal truncated
        if len(steps) >= max_events:
            truncated = True
            return False
        matches = hits or []
        found.extend(matches)
        steps.append({
            "kind": kind, "message": note, "text_index": pos,
            "alignment": align, "pattern_index": p_index,
            "matched": ok, "shift": shift, "state": state,
            "active_prefix": prefix, "fail_from": fail_from,
            "fail_to": fail_to, "hits": matches,
            "operations": ops, "found_count": len(found),
        })
        return True

    def finish(extra=None):
        return {
            "algorithm": requested_method,
            "effective_algorithm": effective_method,
            "strategy": matcher.strategy,
            "target_count": len(targets),
            "sequence_length": len(text),
            "patterns": [t.sequence for t in targets],
            "lps": (lps_table(targets[0].sequence) if effective_method == "KMP" else []),
            "shift_table": (matcher.tables[0] if effective_method == "Horspool" else {}),
            "events": steps, "hits": found, "complete": not truncated,
            "operations": ops, "metadata": extra or {},
        }

    if not text:
        return finish()
    if effective_method in ("Naive", "Horspool") and len(targets[0].sequence) > len(text):
        add("info", "The search pattern is longer than this read, so no alignment is possible.", pos=0)
        return finish()

    if effective_method == "Python str.find":
        hits, _ = matcher.scan(text)
        for target_idx, start in sorted(hits, key=lambda x: (x[1], x[0])):
            t = targets[target_idx]
            h = _hit(t, start)
            if not add("exact_match", "Python str.find returned a complete match. Its internal comparisons are not exposed.",
                       pos=start, align=start, hits=[h]):
                break
        if not steps:
            add("info", "No exact matches. Python's native character-comparison trace is unavailable.", pos=0)
        return finish()

    if effective_method == "Naive":
        pattern = targets[0].sequence
        m = len(pattern)
        for a in range(max(0, len(text) - m + 1)):
            for j in range(m):
                ops += 1
                ok = text[a + j] == pattern[j]
                if not add("compare", f"DNA[{a+j}]={text[a+j]} vs Pattern[{j}]={pattern[j]} → {'MATCH' if ok else 'MISMATCH'}",
                           pos=a+j, align=a, p_index=j, ok=ok):
                    return finish()
                if not ok:
                    break
            else:
                h = _hit(targets[0], a)
                if not add("exact_match", f"Every base matches. Exact pattern found at 0-based position {a}.",
                           pos=a+m-1, align=a, p_index=m-1, ok=True, hits=[h]):
                    return finish()
            if a < len(text) - m:
                if not add("shift", f"Naïve shifts the whole pattern from start {a} → {a+1}. It does not reuse earlier comparisons.",
                           pos=a+1, align=a+1, shift=1):
                    return finish()
        return finish()

    if effective_method == "KMP":
        pattern = targets[0].sequence
        lps = matcher.tables[0]
        i = j = 0
        while i < len(text):
            a = i-j
            ops += 1
            ok = text[i] == pattern[j]
            if not add("compare", f"DNA[{i}]={text[i]} vs Pattern[{j}]={pattern[j]} → {'MATCH' if ok else 'MISMATCH'}",
                       pos=i, align=a, p_index=j, ok=ok):
                return finish()
            if ok:
                i += 1
                j += 1
                if j == len(pattern):
                    start = i - j
                    h = _hit(targets[0], start)
                    if not add("exact_match", f"Full KMP match at {start}! The next j uses LPS[{j-1}]={lps[j-1]} so overlapping matches remain possible.",
                               pos=i-1, align=start, p_index=j-1, ok=True, hits=[h]):
                        return finish()
                    new_j = lps[j-1]
                    if new_j:
                        if not add("fallback", f"After full match, reuse prefix: j {j} → {new_j}, next alignment = {i-new_j}.",
                                   pos=min(i, len(text)-1), align=i-new_j, p_index=new_j,
                                   shift=i-new_j-start):
                            return finish()
                    j = new_j
            elif j:
                new_j = lps[j-1]
                if not add("fallback", f"Mismatch: LPS[{j-1}]={new_j}, so j {j} → {new_j}. DNA index i={i} stays fixed; the pattern slides {j-new_j} positions.",
                           pos=i, align=i-new_j, p_index=new_j, shift=j-new_j):
                    return finish()
                j = new_j
            else:
                i += 1
                if not add("advance", f"Mismatch at j=0. Move DNA index i to {i}; pattern starts at {i}.",
                           pos=min(i, len(text)-1), align=i, p_index=0, shift=1):
                    return finish()
        return finish()

    if effective_method == "Horspool":
        pattern = targets[0].sequence
        m = len(pattern)
        shifts = matcher.tables[0]
        a = 0
        while a <= len(text) - m:
            j = m-1
            while j >= 0:
                ops += 1
                ok = text[a+j] == pattern[j]
                if not add("compare", f"Right-to-left: DNA[{a+j}]={text[a+j]} vs Pattern[{j}]={pattern[j]} → {'MATCH' if ok else 'MISMATCH'}",
                           pos=a+j, align=a, p_index=j, ok=ok):
                    return finish()
                if not ok:
                    break
                j -= 1
            if j < 0:
                h = _hit(targets[0], a)
                if not add("exact_match", f"Full Horspool match at position {a}.",
                           pos=a, align=a, p_index=0, ok=True, hits=[h]):
                    return finish()
            last_ch = text[a+m-1]
            jump = shifts.get(last_ch, m)
            next_a = a+jump
            if not add("shift", f"Bad-character shift: current window's final DNA base '{last_ch}' → shift {jump}. Window start {a} → {next_a}.",
                       pos=min(next_a, len(text)-1), align=next_a, shift=jump):
                return finish()
            a = next_a
        return finish()

    if effective_method == "Aho-Corasick":
        trie = matcher.trie
        # Compute the literal trie prefix associated with each state.
        prefixes = {0: ""}
        q = deque([0])
        while q:
            parent = q.popleft()
            for char, child in trie[parent]["next"].items():
                prefixes[child] = prefixes[parent] + char
                q.append(child)
        state = 0
        for i, ch in enumerate(text):
            while state and ch not in trie[state]["next"]:
                old = state
                state = trie[old]["fail"]
                ops += 1
                if not add("fallback", f"No '{ch}' edge from trie state {old}; follow failure link {old} → {state} without moving the DNA pointer.",
                           pos=i, state=state, prefix=prefixes[state], fail_from=old, fail_to=state):
                    return finish({"states": len(trie)})
            ops += 1
            old = state
            state = trie[state]["next"].get(ch, 0)
            if not add("transition", f"Consume DNA[{i}]='{ch}'. Trie state {old} → {state}; active prefix='{prefixes[state] or '∅'}'.",
                       pos=i, state=state, prefix=prefixes[state], fail_from=old):
                return finish({"states": len(trie)})
            for target_idx in trie[state]["out"]:
                target = targets[target_idx]
                start = i-len(target.sequence)+1
                h = _hit(target, start)
                if not add("exact_match", f"Aho output state {state}: '{target.name}' found at [{start}, {i+1}).",
                           pos=i, align=start, state=state, prefix=prefixes[state], hits=[h]):
                    return finish({"states": len(trie)})
        return finish({"states": len(trie)})
    raise ValueError(f"Unhandled algorithm: {effective_method}")
