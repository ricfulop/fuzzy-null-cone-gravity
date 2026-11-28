"""Synthetic stdlib tests; never import or run either scientific checker."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import ast
import copy
import importlib.util
import json
import sys
import unittest

import leading_certificate as certificate


def synthetic_rows():
    rows = []
    for L, m in certificate.SECTORS:
        labels = certificate.shift_labels(L, m)
        columns = len(labels)
        rank = 3 if m == L else 4
        rows.append({
            "L": L, "N": L + 1, "m": m, "labels": labels,
            "normal_row_labels": certificate.normal_labels(L),
            "leading_DeWitt_Gram": [["0"] * columns for _ in range(columns)],
            "normal_coefficients": [["0"] * columns for _ in range(12)],
            "exact_leading_rank": rank, "exact_leading_normal_map_rank": rank,
            "nonzero_principal_minor_indices": list(range(rank)),
            "nonzero_principal_minor_determinant": "1",
            "positive_Sylvester_minors": ["1"] * rank,
        })
    return rows


class CertificateTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        for name in certificate.SOURCE_FILES:
            (self.directory / name).write_text("# synthetic source " + name)
        self.sources = certificate.source_fingerprints(self.directory)
        self.path = self.directory / "exact-leading.json"
        self.rows = synthetic_rows()
        certificate.write_certificate(self.path, self.rows, self.sources)

    def rewrite_payload(self, rows):
        data = (json.dumps(rows) + "\n").encode()
        self.path.write_bytes(data)
        sidecar = certificate.manifest_path(self.path)
        manifest = json.loads(sidecar.read_bytes())
        manifest["payload_sha256"] = certificate.sha256_bytes(data)
        sidecar.write_text(json.dumps(manifest))

    def test_roundtrip_keeps_data_list_and_records_exact_input_digests(self):
        rows, evidence = certificate.load_certificate(self.path, self.sources)
        self.assertEqual(rows, self.rows)
        self.assertIsInstance(json.loads(self.path.read_bytes()), list)
        self.assertEqual(evidence["payload_sha256"],
                         certificate.sha256_bytes(self.path.read_bytes()))
        self.assertEqual(evidence["manifest_sha256"], certificate.sha256_bytes(
            certificate.manifest_path(self.path).read_bytes()))
        self.assertEqual(evidence["source_sha256"], self.sources)

    def test_sibling_loader_works_without_its_directory_on_sys_path(self):
        # Extract only the stdlib loader, never import/execute scientific code.
        source_path = Path(__file__).with_name("all_modes.py")
        tree = ast.parse(source_path.read_text())
        loader = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                      and node.name == "_load_leading_certificate")
        isolated = ast.fix_missing_locations(ast.Module(body=[loader], type_ignores=[]))
        helper_path = self.directory / "leading_certificate.py"
        helper_path.write_text("MARKER = 'loaded adjacent helper'\n")
        namespace = {"Path": Path, "importlib": importlib,
                     "__file__": str(self.directory / "all_modes.py")}
        before = list(sys.path)
        self.assertNotIn(str(self.directory), before)
        exec(compile(isolated, str(source_path), "exec"), namespace)
        module = namespace["_load_leading_certificate"]()
        self.assertEqual(module.MARKER, "loaded adjacent helper")
        self.assertEqual(Path(module.__file__).resolve(), helper_path.resolve())
        self.assertEqual(sys.path, before)

    def test_missing_legacy_manifest_requires_regeneration(self):
        certificate.manifest_path(self.path).unlink()
        with self.assertRaisesRegex(ValueError, "legacy exports.*regenerate"):
            certificate.load_certificate(self.path, self.sources)

    def test_stale_source_is_rejected_before_data_read(self):
        (self.directory / "all_modes.py").write_text("# changed synthetic source")
        expected = certificate.source_fingerprints(self.directory)
        original = Path.read_bytes
        reads = []

        def read(path):
            reads.append(path)
            return original(path)

        with patch.object(Path, "read_bytes", read):
            with self.assertRaisesRegex(ValueError, "source hashes differ"):
                certificate.load_certificate(self.path, expected)
        self.assertEqual(reads, [certificate.manifest_path(self.path)])

    def test_metadata_layout_change_is_rejected(self):
        sidecar = certificate.manifest_path(self.path)
        manifest = json.loads(sidecar.read_bytes())
        manifest["basis_layout"][0]["normal_row_labels"].reverse()
        sidecar.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "basis layout changed"):
            certificate.load_certificate(self.path, self.sources)

    def test_hash_mismatch_precedes_data_json_parsing(self):
        self.path.write_bytes(b"not even JSON")
        with patch.object(certificate.json, "loads", wraps=json.loads) as parse:
            with self.assertRaisesRegex(ValueError, "data hash does not match"):
                certificate.load_certificate(self.path, self.sources)
        self.assertEqual(parse.call_count, 1)  # Only the sidecar was parsed.

    def test_replacement_after_read_cannot_change_parsed_data_or_provenance(self):
        sidecar = certificate.manifest_path(self.path)
        original_data = self.path.read_bytes()
        original_manifest = sidecar.read_bytes()
        original = Path.read_bytes
        reads = []

        def replace_after_read(path):
            data = original(path)
            reads.append(path)
            path.write_bytes(b"replacement must not be consumed")
            return data

        with patch.object(Path, "read_bytes", replace_after_read), \
             patch.object(certificate.json, "loads", wraps=json.loads) as parse:
            rows, evidence = certificate.load_certificate(self.path, self.sources)
        self.assertEqual(reads, [sidecar, self.path])
        self.assertEqual([call.args[0] for call in parse.call_args_list],
                         [original_manifest, original_data])
        self.assertEqual(rows, self.rows)
        self.assertEqual(evidence["payload_sha256"], certificate.sha256_bytes(original_data))
        self.assertEqual(evidence["manifest_sha256"], certificate.sha256_bytes(original_manifest))

    def test_missing_duplicate_or_reordered_sectors_are_rejected(self):
        for rows in (self.rows[:-1], [self.rows[0]] * 7, list(reversed(self.rows))):
            with self.subTest(rows=[(row["L"], row["m"]) for row in rows]):
                self.rewrite_payload(rows)
                with self.assertRaisesRegex(ValueError, "seven ordered sectors|sector order"):
                    certificate.load_certificate(self.path, self.sources)

    def test_authenticated_payload_still_needs_ordered_labels_and_matrix_shapes(self):
        cases = []
        for field in ("labels", "normal_row_labels"):
            rows = copy.deepcopy(self.rows)
            rows[0][field].reverse()
            cases.append((field, rows))
        for field in ("normal_coefficients", "leading_DeWitt_Gram"):
            for mutation in ("row", "column", "entry"):
                rows = copy.deepcopy(self.rows)
                if mutation == "row":
                    rows[0][field].pop()
                elif mutation == "column":
                    rows[0][field][0].pop()
                else:
                    rows[0][field][0][0] = 0
                cases.append((field + mutation, rows))
        for name, rows in cases:
            with self.subTest(name=name):
                self.rewrite_payload(rows)
                with self.assertRaises(ValueError):
                    certificate.load_certificate(self.path, self.sources)

    def test_integer_basis_labels_cannot_be_replaced_by_booleans_or_floats(self):
        for value in (False, 0.0):
            with self.subTest(value=value):
                rows = copy.deepcopy(self.rows)
                rows[0]["labels"][0][1] = value
                self.rewrite_payload(rows)
                with self.assertRaisesRegex(ValueError, "shift label order changed"):
                    certificate.load_certificate(self.path, self.sources)

    def test_principal_minor_indices_must_be_unique_in_bounds_integers(self):
        for indices in ([0, 1, 1, 3], [0, 1, 2, 99], [0, 1, 2, -1],
                        [False, 1, 2, 3], [0, 1, 2]):
            with self.subTest(indices=indices):
                rows = copy.deepcopy(self.rows)
                rows[0]["nonzero_principal_minor_indices"] = indices
                self.rewrite_payload(rows)
                with self.assertRaisesRegex(ValueError, "minor indices are invalid"):
                    certificate.load_certificate(self.path, self.sources)

    def test_matching_but_wrong_sector_ranks_are_rejected_before_publication(self):
        for index in (0, 2):  # A rank-four sector and a rank-three sector.
            with self.subTest(index=index):
                rows = copy.deepcopy(self.rows)
                rows[index]["exact_leading_rank"] = 2
                rows[index]["exact_leading_normal_map_rank"] = 2
                output = self.directory / "invalid.json"
                with self.assertRaisesRegex(ValueError, "unexpected exact-leading sector rank"):
                    certificate.write_certificate(output, rows, self.sources)
                self.assertFalse(output.exists())
                self.rewrite_payload(rows)
                with self.assertRaisesRegex(ValueError, "unexpected exact-leading sector rank"):
                    certificate.load_certificate(self.path, self.sources)


if __name__ == "__main__":
    unittest.main(verbosity=2)
