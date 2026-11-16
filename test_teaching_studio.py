import random
import unittest
from search_engine import Matcher, Target
from teaching_trace import build_teaching_trace


class TeachingTraceTests(unittest.TestCase):
    def compare_single(self, dna, pat):
        t = [Target("primer", pat)]
        for method in ["Python str.find", "Naive", "KMP", "Horspool", "Aho-Corasick"]:
            with self.subTest(dna=dna[:20], pattern=pat, method=method):
                trace = build_teaching_trace(dna, t, method)
                self.assertTrue(trace["complete"])
                coords = sorted((0, h["start"]) for h in trace["hits"])
                engine_hits, ops = Matcher(method, t).scan(dna)
                self.assertEqual(coords, engine_hits)
                if ops is not None:
                    self.assertEqual(trace["operations"], ops)
                self.assertTrue(all(dna[h['start']:h['end']] == pat for h in trace['hits']))
                for ev in trace["events"]:
                    self.assertIn(ev["kind"], {"compare", "shift", "fallback", "advance", "transition", "exact_match", "info"})
                    if ev["hits"]:
                        self.assertEqual(ev["kind"], "exact_match")
                self.assertEqual(trace["events"][-1]["found_count"], len(trace["hits"]))

    def test_known_patterns(self):
        for dna, p in [("ACACACAC", "ACAC"), ("CACCCGCACCCG", "CACCCG"),
                       ("A"*80+"C", "AAAAAC"), ("CGTACGTA", "ACGT"),
                       ("TTTTT", "AAA"), ("AC", "ACGT"), ("A", "A")]:
            self.compare_single(dna, p)

    def test_random_exact_coordinate_and_operation_parity(self):
        rng = random.Random(141)
        for _ in range(150):
            n = rng.randint(1, 120)
            m = rng.randint(1, 16)
            dna = ''.join(rng.choices('ACGTN', k=n))
            pat = ''.join(rng.choices('ACGT', k=m))
            self.compare_single(dna, pat)

    def test_aho_multiple_targets(self):
        targets = [Target("p", "ACAC"), Target("q", "CAC"), Target("r", "CG"),
                   Target("rc", "GTAC", "-")]
        for dna in ["ACACACAC", "CGTACGTACG", "A"*300+"CGTACAC"]:
            result = build_teaching_trace(dna, targets, "Aho-Corasick")
            expected, operations = Matcher("Aho-Corasick", targets).scan(dna)
            self.assertEqual(sorted((next(i for i, t in enumerate(targets) if t.name==h['name']), h['start']) for h in result["hits"]), expected)
            self.assertEqual(result["operations"], operations)

    def test_kmp_shows_real_fallback(self):
        trace = build_teaching_trace("A" * 30 + "C", [Target("x", "AAAAC")], "KMP")
        fallbacks = [x for x in trace["events"] if x["kind"] == "fallback"]
        self.assertGreater(len(fallbacks), 3)
        self.assertTrue(any("DNA index i=" in x["message"] for x in fallbacks))

    def test_horspool_shifts_and_exact_overlap(self):
        trace = build_teaching_trace("AAAAAAAA", [Target("x", "AAA")], "Horspool")
        self.assertEqual([h["start"] for h in trace["hits"]], list(range(6)))
        self.assertTrue(any(x['kind'] == 'shift' for x in trace['events']))

    def test_100_step_limit_removed(self):
        trace = build_teaching_trace("A"*300+"C", [Target("x", "AAAAC")], "Naive")
        self.assertTrue(trace["complete"])
        self.assertGreater(len(trace["events"]), 100)
        self.assertEqual([h["start"] for h in trace["hits"]], [296])

    def test_explicit_safety_limit(self):
        trace = build_teaching_trace("A"*300+"C", [Target("x", "AAAAC")], "Naive", max_events=10)
        self.assertFalse(trace["complete"])
        self.assertEqual(len(trace["events"]), 10)


if __name__ == "__main__":
    unittest.main()
