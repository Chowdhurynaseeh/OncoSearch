import gzip
import random
import tempfile
import unittest
from pathlib import Path

from search_engine import (
    METHODS, Matcher, Target, lps_table, prepare_targets,
    reverse_complement, iter_records, load_sample, run_benchmark, trace_steps,
)


class SearchAlgorithmTests(unittest.TestCase):
    def verify_all(self, text, patterns, both=False):
        targets = prepare_targets([(f"q{i}", p) for i, p in enumerate(patterns)], both)
        expected = Matcher("Python str.find", targets).scan(text)[0]
        for algorithm in METHODS:
            got, _ = Matcher(algorithm, targets).scan(text)
            self.assertEqual(got, expected, f"Mismatch in {algorithm}")

    def test_lps(self):
        self.assertEqual(lps_table("ACAC"), [0, 0, 1, 2])
        self.assertEqual(lps_table("AAAA"), [0, 1, 2, 3])

    def test_overlapping(self):
        self.verify_all("ACACACAC", ["ACAC", "AC"])

    def test_longer_than_text(self):
        self.verify_all("AC", ["ACG", "AAAAAC"])

    def test_no_matches_and_n(self):
        self.verify_all("NNNNNN", ["ACGT", "A", "N" if False else "CG"])

    def test_adversarial(self):
        self.verify_all("A" * 1000 + "C", ["AAAAAC", "AAAC", "AAC"])

    def test_both_strands(self):
        self.assertEqual(reverse_complement("ACGA"), "TCGT")
        self.verify_all("TCGTACGA", ["ACGA"], both=True)

    def test_random(self):
        rng = random.Random(17)
        for _ in range(100):
            seq = "".join(rng.choices("ACGTN", k=rng.randrange(0, 130)))
            pats = ["".join(rng.choices("ACGT", k=rng.randrange(1, 15)))
                    for _ in range(rng.randrange(1, 9))]
            self.verify_all(seq, pats, both=True)

    def test_fastq_and_wrapped_fasta(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "real.fastq.gz"
            with gzip.open(path, "wt") as f:
                f.write("@r1\nACACAC\n+\nIIIIII\n@r2\nAAAAAC\n+\n!!!!!!\n")
            fasta = Path(td) / "test.fasta"
            fasta.write_text(">chr1\nACAC\nACAC\n")
            self.assertEqual(len(list(iter_records(path))), 2)
            self.assertEqual(list(iter_records(fasta))[0][1], "ACACACAC")
            reads, stats = load_sample([("fastq", path), ("fa", fasta)], 10, 25)
            self.assertEqual(stats["used_reads"], 2)
            self.assertEqual(stats["filtered_reads"], 1)
            targets = [Target("p", "ACAC")]
            metrics, previews = run_benchmark(reads, targets, repetitions=1)
            self.assertTrue(all(m["correct_vs_first"] for m in metrics))
            self.assertEqual(metrics[0]["total_matches"], 5)

    def test_animation(self):
        tar = [Target("motif", "ACAC")]
        for method in METHODS:
            steps = trace_steps("ACACACAC", tar, method)
            self.assertTrue(steps, method)
            self.assertLessEqual(len(steps), 140)


if __name__ == "__main__":
    unittest.main(verbosity=2)
