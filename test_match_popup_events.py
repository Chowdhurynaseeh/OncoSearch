import unittest
from search_engine import METHODS, Target, trace_steps


class MatchPopupEventTests(unittest.TestCase):
    def test_exact_hits_include_coordinates_all_methods(self):
        text = "CACCCGCACCCG"
        targets = [Target("short_motif", "CACCCG")]
        for method in METHODS:
            with self.subTest(method=method):
                steps = trace_steps(text, targets, method, max_events=100)
                hits = [hit for event in steps for hit in event.get("hits", [])]
                self.assertEqual([h["start"] for h in hits], [0, 6])
                self.assertEqual([h["end"] for h in hits], [6, 12])
                self.assertTrue(all(text[h["start"]:h["end"]] == h["sequence"]
                                    for h in hits))

    def test_partial_char_match_does_not_trigger_popup(self):
        text = "CAAAAA"
        targets = [Target("short_motif", "CACCCG")]
        for method in METHODS:
            with self.subTest(method=method):
                steps = trace_steps(text, targets, method, max_events=100)
                self.assertFalse(any(event.get("hits") for event in steps))

    def test_multi_target_aho_hits_have_correct_query_labels(self):
        targets = [Target("first", "ACAC"), Target("second", "CAC")]
        steps = trace_steps("ACACACAC", targets, "Aho-Corasick")
        hits = [h for s in steps for h in s.get("hits", [])]
        self.assertEqual([(h["name"], h["start"]) for h in hits],
                         [("first", 0), ("second", 1), ("first", 2),
                          ("second", 3), ("first", 4), ("second", 5)])


if __name__ == "__main__":
    unittest.main()
