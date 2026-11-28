"""Mocked supervisor regressions. All child, ps, and signal operations are replaced."""

from contextlib import ExitStack
from datetime import datetime, timedelta, timezone
from itertools import chain, count
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import io
import json
import unittest
from unittest.mock import Mock, call, patch

import run_bounded as runner


def manufacture_bundle(path):
    """Synthetic wire data with all 78 member files. No scientific modules are loaded."""
    path.parent.mkdir()
    folder = path.parent / "atlas-columns"
    folder.mkdir()
    rows, cols = runner.labels(3, 4), runner.labels(0, 2)
    hashes = {"config_sha256": "a"*64, "source_sha256": {n: "b"*64 for n in runner.SCIENCE_NAMES}}
    ident = {**hashes, "runtime": {"python": "3.12.11", "sympy": "1.14.0"},
             "archived_regression_reference": {"path": "paper1/"+runner.MAP_PATH,
                                                "sha256": runner.ARCHIVE_PINS[runner.MAP_PATH]}}
    catalog = {**ident, "atlas_complete": False, "labels": ["field","ell","m"],
               "columns": cols, "rows": rows, "axes": list(runner.AXES),
               "weights": "Cartesian generator, unit w, exact real/imaginary rational strings"}
    def save(name, obj):
        member = folder/name
        member.write_text(json.dumps(obj))
        return {"file": "atlas-columns/"+name, "sha256": runner.sha(member)}

    cat_ref = save("catalog.json", catalog)
    files = []
    z_columns = [i for i,c in enumerate(cols) if c[2] == 0][:9]
    z_rows = [i for i,c in enumerate(rows) if c[2] == 0][:8]
    for i, label in enumerate(cols):
        maps = {axis: [] for axis in runner.AXES}
        if i in z_columns:
            maps["z"].append({"row": z_rows[z_columns.index(i)%8], "real": "1", "imag": "0"})
        if i == 0:
            row = next(j for j,r in enumerate(rows) if r[2] == 1)
            maps["x"] = [{"row": row, "real": "1/2", "imag": "0"}]
            maps["y"] = [{"row": row, "real": "0", "imag": "-1/2"}]
        files.append(save(f"column-{i:03d}.json", {
            "column": i, "label": label, "atlas_complete": False, "catalog_sha256": cat_ref["sha256"],
            "maps": maps, "map_seconds": dict.fromkeys(runner.TIMERS, 0.1),
            "rotation_zero": True, "round_reconstruction_checked": True, "metric_chain_rule_checked": True}))
    evidence = [
        ("explicit_low_harmonic_basis", {"fixture_count": 6, "signed_modes": [-2,-1,0,1,2]}),
        ("catalog_and_fourier_selection", {"input_columns": 78, "discarded_labels": 160,
                                          "generator_columns": 234, "nonzero_entries": 11}),
        ("archived_axisymmetric_eight_rows", {"input_columns": 22, "nonzero_output_rows": 8,
            "nonzero_entries": 9, "all_other_entries_zero": True, "temporal_W_rows_included": True,
            "round_trace_X_from_C_axial_included": True}),
        ("all_column_complex_conjugacy", {"column_pairs_checked": 234, "coefficient_comparisons": 11,
              "negative_m_convention": "conjugate positive mode; no factorial reweighting"}),
        ("all_column_cartesian_transverse_phases", {"input_columns": 78, "coefficient_comparisons": 1,
                                                    "raising_factor": "-I", "lowering_factor": "I"}),
        ("all_column_rotation_zero", {"rotation_columns_checked": 234, "all_discarded_zero": True}),
        ("deliberate_coefficient_defects", {"case_count": 3,
                 "cases": [{"case": n, "passed": True} for n in runner.MUTATIONS]}),
    ]
    independent = {"case_count": 7, "status": "independent atlas coefficient controls passed",
                   "cases": [{"case": n, "passed": True, "evidence": v} for n,v in evidence],
                   "scope": "exact coefficient maps and conventions; no Euler pairing or seed",
                   "archived_source": {"path": "paper1/"+runner.REFERENCE_PATH,
                       "function": "expected_discarded_generator",
                       "convention": "unit w; eight output rows, nine nonzero matrix entries"}}
    path.write_text(json.dumps({**ident, "atlas_complete": True, "state": "complete N=3 discarded coefficient maps",
       "generator_count": 3, "coefficient_column_count": 78, "discarded_row_count": 160,
       "columns": files, "catalog": cat_ref, "checks": dict.fromkeys(runner.CHECKS,True),
       "scope": dict.fromkeys(runner.SCOPE,False), "independent_controls": independent}))
    return hashes


class GroupProbeTests(unittest.TestCase):
    def test_permission_denied_requires_a_readable_process_table(self):
        cases = [
            ("1 1 50\n2 2 80\n", "empty"),
            ("1 1 50\n111 77 80\n", "live"),
            ("", "unknown"),
            ("1 1 50\ntruncated\n", "unknown"),
            (runner.subprocess.TimeoutExpired("ps", 1), "unknown"),
        ]
        for table, expected in cases:
            with self.subTest(table=table), \
                 patch.object(runner.os, "killpg", side_effect=PermissionError()), \
                 patch.object(runner.subprocess, "check_output") as ps:
                if isinstance(table, str):
                    ps.return_value = table
                else:
                    ps.side_effect = table
                state = runner.group_state(77)
                self.assertEqual(state["status"], expected)
                self.assertEqual(state["probe"], "EPERM")
                if expected == "live":
                    self.assertEqual(state["members"], [{"pid": 111, "rss_bytes": 81920}])
                self.assertGreater(ps.call_args.kwargs["timeout"], 0)

    def test_stubborn_descendant_escalated_after_launcher_exit(self):
        process = Mock(pid=77)
        process.poll.return_value = 0
        clock = count()
        events = []
        with patch.object(runner.os, "killpg", side_effect=PermissionError()) as signals, \
             patch.object(runner, "group_state", return_value={"status": "live"}), \
             patch.object(runner.time, "monotonic", side_effect=lambda: next(clock)), \
             patch.object(runner.time, "sleep"):
            runner.stop_group(process, events)
        self.assertEqual(signals.call_args_list, [call(77, runner.signal.SIGTERM),
                                                 call(77, runner.signal.SIGKILL)])
        self.assertEqual(events[-1]["status"], "live")


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.patches = ExitStack()
        self.addCleanup(self.patches.close)
        self.out = Path(self.patches.enter_context(TemporaryDirectory()))
        self.handlers = {}
        self.received_during_install = False

        def install(sig, handler):
            self.handlers[sig] = handler
            if self.received_during_install and callable(handler):
                handler(sig, None)
            return runner.signal.SIG_DFL

        self.replace(runner.signal, "signal", side_effect=install)
        self.killpg = self.replace(runner.os, "killpg", side_effect=ProcessLookupError())
        self.ps = self.replace(runner.subprocess, "check_output", return_value="55 55 50\n")
        self.child = Mock(pid=77, returncode=0, stdout=io.BytesIO())
        self.child.poll.return_value = 0
        self.launch = self.replace(runner.subprocess, "Popen", return_value=self.child)
        self.stop_implementation = runner.stop_group
        self.stop = self.replace(runner, "stop_group")
        self.console = self.replace(runner.os, "write", side_effect=lambda fd, data: len(data))
        console_stream = self.replace(runner.sys, "stdout", new=Mock())
        console_stream.fileno.return_value = 1
        self.replace(runner.os, "get_blocking", return_value=True)
        self.replace(runner.os, "set_blocking")
        self.record = {"wall_seconds": 10, "memory_bytes": 2048,
                       "scientific_invocations": 0, "command": ["dummy"],
                       "source_sha256_before": {}, "checker_input_sha256": {}}
        self.certificate = self.out / runner.CERTIFICATE
        self.record["checker_input_sha256"] = manufacture_bundle(self.certificate)

    def replace(self, obj, name, **kwargs):
        return self.patches.enter_context(patch.object(obj, name, **kwargs))

    def supervise(self):
        rc = runner.supervise(self.out, self.record, {},
                              datetime.now(timezone.utc) + timedelta(minutes=10))
        execution = json.loads((self.out / "execution.json").read_bytes())
        samples = json.loads((self.out / "resource-samples.json").read_bytes())
        self.assertEqual(execution["metrics"], json.loads((self.out / "metrics.json").read_bytes()))
        self.assertEqual(execution["rss_sample_count"], len(samples))
        self.assertEqual(rc == 0, execution["status"] == "passed")
        self.assertEqual(set(self.handlers.values()), {runner.signal.SIG_DFL})
        if self.launch.called and execution["scientific_invocations"]:
            self.assertGreater(self.child.wait.call_args.kwargs["timeout"], 0)
        return execution, samples

    def test_normal_exit_eperm_empty_and_exact_raw_log(self):
        self.killpg.side_effect = PermissionError()
        self.child.stdout = io.BytesIO(b"stage\nlast")
        self.console.side_effect = lambda fd, data: 1   # partial writes
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["group_disposition"]["status"], "empty")
        self.assertEqual((self.out / "verification.log").read_bytes(), b"stage\nlast")
        self.assertGreaterEqual(self.console.call_count, len(b"stage\nlast"))
        self.stop.assert_not_called()

    def test_timeout_preserves_measured_sample_when_cleanup_raises(self):
        self.record["wall_seconds"] = 1
        clock = chain([0, 0.5], count(2))
        self.replace(runner.time, "monotonic", side_effect=lambda: next(clock))
        self.child.poll.return_value = None
        self.child.returncode = None
        self.ps.return_value = "77 77 12\n"
        self.killpg.side_effect = PermissionError()
        self.stop.side_effect = PermissionError("manufactured termination error")
        self.child.wait.side_effect = runner.subprocess.TimeoutExpired("dummy", 0.5)
        evidence, samples = self.supervise()
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["group_disposition"]["status"], "live")
        self.assertEqual([e["action"] for e in evidence["cleanup_errors"]],
                         ["terminate_group", "wait_child"])
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["rss_bytes"], 12 * 1024)

    def test_timeout_with_eperm_after_term_keeps_both_records(self):
        self.record["wall_seconds"] = 1
        clock = chain([0, 0.5], count(2))
        self.replace(runner.time, "monotonic", side_effect=lambda: next(clock))
        self.child.poll.side_effect = [None, -15]
        self.child.returncode = -15
        self.stop.side_effect = self.stop_implementation
        self.killpg.side_effect = [None, None, PermissionError(), PermissionError()]
        self.ps.side_effect = ["77 77 12\n", "1 1 50\n", "1 1 50\n"]
        evidence, samples = self.supervise()
        self.assertEqual(evidence["group_disposition"]["status"], "empty")
        self.assertEqual(evidence["group_disposition"]["probe"], "EPERM")
        self.assertEqual(evidence["cleanup_errors"], [])
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["sampled_max_rss_bytes"], samples[0]["rss_bytes"])
        self.assertEqual(self.killpg.call_args_list,
                         [call(77, 0), call(77, runner.signal.SIGTERM), call(77, 0), call(77, 0)])

    def test_zero_child_exit_with_unknown_group_is_failure(self):
        self.killpg.side_effect = PermissionError()
        self.ps.side_effect = runner.subprocess.TimeoutExpired("ps", 1)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["exit_code"], 0)
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["group_disposition"]["status"], "unknown")
        self.stop.assert_called_once()

    def test_signal_before_spawn(self):
        self.received_during_install = True
        evidence, _ = self.supervise()
        self.launch.assert_not_called()
        self.assertTrue(evidence["received_signals"])
        self.assertEqual(evidence["scientific_invocations"], 0)
        self.assertEqual(evidence["group_disposition"]["status"], "not_spawned")

    def test_signal_during_spawn_still_cleans_child_and_drains_log(self):
        def launch(*args, **kwargs):
            self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return self.child

        self.launch.side_effect = launch
        self.killpg.side_effect = [None, ProcessLookupError()]
        self.child.stdout = io.BytesIO(b"before stop\n")
        evidence, samples = self.supervise()
        self.assertEqual(evidence["received_signals"], [int(runner.signal.SIGTERM)])
        self.assertEqual(evidence["scientific_invocations"], 1)
        self.assertEqual(evidence["stop_reason"], "supervisor received signal")
        self.assertEqual(evidence["group_disposition"]["status"], "empty")
        self.assertFalse(evidence["relay_active_after_join"])
        self.assertEqual((self.out / "verification.log").read_bytes(), b"before stop\n")
        self.assertEqual(samples, [])
        self.stop.assert_called_once()

    def test_stalled_relay_join_is_bounded_and_fails_closed(self):
        relay = self.replace(runner.threading, "Thread").return_value
        relay.is_alive.return_value = True
        evidence, _ = self.supervise()
        self.assertGreater(relay.join.call_args.kwargs["timeout"], 0)
        self.assertTrue(evidence["relay_active_after_join"])
        self.assertEqual(evidence["status"], "failed")
        # The mock did not start the worker that ordinarily owns these streams.
        relay_target = runner.threading.Thread.call_args.kwargs["target"]
        relay_target()


    def test_zero_exit_without_certificate_fails(self):
        self.certificate.unlink()
        evidence, _ = self.supervise()
        self.assertEqual(evidence["stop_reason"], "terminal certificate validation failed")
        self.assertEqual(evidence["metrics"]["atlas_bundle_passed"], 0)
        self.assertNotIn("checkpoint_column_count", evidence["metrics"])

    def test_zero_exit_observed_after_deadline_fails(self):
        with patch.object(runner.time, "monotonic", side_effect=[0, 11, 11]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["exit_code"], 0)
        self.assertEqual(evidence["stop_reason"], "wall budget reached")
        self.assertEqual(evidence["metrics"]["atlas_bundle_passed"], 0)

    def test_final_schema_rejections(self):
        baseline = self.certificate.read_bytes()
        for keys, value in [
            (["config_sha256"], "wrong"), (["source_sha256"], {}), (["atlas_complete"], 1),
            (["runtime", "python"], "3.13.0"), (["scope", runner.SCOPE[0]], True),
            (["checks", runner.CHECKS[0]], 1), (["archived_regression_reference", "sha256"], ""),
            (["columns", 0, "file"], "../elsewhere"), (["columns", 1, "file"], "atlas-columns/column-000.json"),
            (["columns", 0, "sha256"], "0"*64), (["columns"], []),
            (["independent_controls", "case_count"], 6),
            (["independent_controls", "cases", 6, "evidence", "cases", 1, "passed"], False),
            (["independent_controls", "cases", 1, "evidence", "nonzero_entries"], 99),
        ]:
            with self.subTest(keys=keys):
                obj = json.loads(baseline)
                target = obj
                for k in keys[:-1]:
                    target = target[k]
                target[keys[-1]] = value
                self.certificate.write_text(json.dumps(obj))
                with self.assertRaises(ValueError):
                    runner.read_certificate(self.certificate, self.record["checker_input_sha256"])

    def test_hash_consistent_malformed_member_rejected(self):
        final = json.loads(self.certificate.read_bytes())
        for member_key, changes in [
            ("catalog", [(["rows",0,1], False), (["columns",0,0], "invalid"), (["axes"], ["z","x","y"])]),
            ("column", [(["column"], False), (["label",2], True), (["rotation_zero"], 1),
                        (["atlas_complete"], True), (["catalog_sha256"], "wrong"), (["map_seconds","rotation_x"], float("inf")),
                        (["maps","x",0,"imag"], "0.0"), (["maps","x",0,"real"], "2/4"),
                        (["maps","x",0,"real"], "0"), (["maps","x",0,"row"], True),
                        (["maps","x",0,"row"], 160), (["maps","x",0,"row"], 0),
                        (["maps","x"], [{"row":60,"real":"1","imag":"0"}]*2)])]:
            ref = final["catalog"] if member_key == "catalog" else final["columns"][0]
            member = self.certificate.parent/ref["file"]
            baseline = member.read_bytes()
            for keys, value in changes:
                with self.subTest(member=member_key, keys=keys, value=value):
                    obj = json.loads(baseline)
                    target = obj
                    for k in keys[:-1]:
                        target = target[k]
                    target[keys[-1]] = value
                    member.write_text(json.dumps(obj))
                    ref["sha256"] = runner.sha(member)
                    self.certificate.write_text(json.dumps(final))
                    with self.assertRaises(ValueError):
                        runner.read_certificate(self.certificate, self.record["checker_input_sha256"])
            member.write_bytes(baseline)
            ref["sha256"] = runner.sha(member)

    def test_unlisted_missing_or_symlink_column_rejected(self):
        folder = self.certificate.parent/"atlas-columns"
        member = folder/"column-077.json"
        extra = folder/"column-078.json"
        extra.write_bytes(member.read_bytes())
        def reject():
            with self.assertRaises(ValueError):
                runner.read_certificate(self.certificate, self.record["checker_input_sha256"])
        reject()
        extra.unlink()
        outside = self.out/"other.json"
        member.rename(outside)
        reject()
        member.symlink_to(outside)
        reject()

    def test_pass_metrics_follow_authenticated_terminal_records(self):
        def inspect(metrics):
            self.assertEqual(metrics["atlas_bundle_passed"], 1)
            self.assertEqual(metrics["checkpoint_column_count"], 78)
            self.assertEqual(metrics["independent_control_count"], 7)
            marker = json.loads((self.out / "terminal-evidence.json").read_bytes())
            self.assertEqual(marker["status"], "passed")
            self.assertEqual(marker["execution_sha256"], runner.sha(self.out / "execution.json"))
            execution = json.loads((self.out / "execution.json").read_bytes())
            self.assertEqual(execution["certificate_sha256"], runner.sha(self.certificate))
        emit = self.replace(runner, "publish_metrics", side_effect=inspect)
        self.supervise()
        emit.assert_called_once()

    def test_resource_stop_keeps_sample_and_fails(self):
        self.child.poll.side_effect = [None, 0]
        self.ps.return_value = "77 77 9\n"
        evidence, samples = self.supervise()
        self.assertEqual(evidence["stop_reason"], "sampled process-group RSS exceeded threshold")
        self.assertEqual(samples[0]["rss_bytes"], 9 * 1024)
        self.assertEqual(evidence["metrics"]["atlas_bundle_passed"], 0)

    def test_existing_terminal_json_is_preserved(self):
        path = self.out / "standalone.json"
        runner.save_json(path, {"value": 1})
        with self.assertRaises(FileExistsError):
            runner.save_json(path, {"value": 2})
        self.assertEqual(json.loads(path.read_text()), {"value": 1})
        self.assertFalse(path.with_suffix(".json.tmp").exists())

    def test_console_backpressure_preserves_archived_log(self):
        self.console.side_effect = BlockingIOError("console reader stalled")
        self.child.stdout = io.BytesIO(b"archived regardless of console\n" * 1024)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["relay_errors"], [])
        self.assertEqual((self.out / "verification.log").read_bytes(),
                         b"archived regardless of console\n" * 1024)

    def test_empty_rss_sample_set_is_unknown(self):
        evidence, samples = self.supervise()
        self.assertEqual(samples, [])
        self.assertIsNone(evidence["sampled_max_rss_bytes"])

    def test_sampling_timeout_is_capped_by_remaining_budget(self):
        self.record["wall_seconds"] = 1
        self.child.poll.side_effect = [None, 0]
        self.replace(runner.time, "sleep")
        with patch.object(runner.time, "monotonic", side_effect=[0, 0.75, 0.8, 0.9, 1]):
            evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(self.ps.call_args.kwargs["timeout"], 0.25)

    def test_signal_after_terminal_boundary_is_outside_the_decision(self):
        original_save = runner.save_json
        def inject(path, payload):
            if path.name == "metrics.json":
                self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return original_save(path, payload)
        self.replace(runner, "save_json", side_effect=inject)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "passed")
        self.assertEqual(evidence["received_signals"], [])
        self.assertEqual(evidence["cancellation_boundary"], "before terminal evidence publication")

    def test_terminal_marker_authenticates_both_records(self):
        evidence, _ = self.supervise()
        marker = json.loads((self.out / "terminal-evidence.json").read_bytes())
        self.assertEqual(marker["status"], evidence["status"])
        for key, filename in [("execution_sha256", "execution.json"), ("metrics_sha256", "metrics.json")]:
            self.assertEqual(marker[key], hashlib.sha256((self.out / filename).read_bytes()).hexdigest())

    def test_failed_execution_write_cannot_commit_pass_metrics(self):
        original_save = runner.save_json
        def fail_execution(path, payload):
            if path.name == "execution.json":
                raise OSError("synthetic execution-write failure")
            return original_save(path, payload)
        self.replace(runner, "save_json", side_effect=fail_execution)
        with self.assertRaisesRegex(OSError, "synthetic execution-write"):
            self.supervise()
        self.assertFalse((self.out / "terminal-evidence.json").exists())

    def test_terminal_marker_write_failure_never_emits_pass(self):
        save = runner.save_json
        def fail_marker(path, obj):
            if path.name == "terminal-evidence.json":
                raise OSError("mock terminal marker failure")
            return save(path, obj)
        self.replace(runner, "save_json", side_effect=fail_marker)
        emit = self.replace(runner, "publish_metrics")
        with self.assertRaisesRegex(OSError, "marker failure"):
            self.supervise()
        emit.assert_not_called()
        self.assertFalse((self.out/"terminal-evidence.json").exists())

    def test_duplicate_json_keys_are_rejected(self):
        original = self.certificate.read_text()
        self.certificate.write_text('{"atlas_complete":false,'+original[1:])
        with self.assertRaisesRegex(ValueError, "Duplicate JSON field"):
            runner.read_certificate(self.certificate, self.record["checker_input_sha256"])

    def test_signal_during_certificate_validation_still_fails(self):
        original_read = runner.read_certificate
        def cancel(*args, **kwargs):
            self.handlers[runner.signal.SIGTERM](runner.signal.SIGTERM, None)
            return original_read(*args, **kwargs)
        self.replace(runner, "read_certificate", side_effect=cancel)
        evidence, _ = self.supervise()
        self.assertEqual(evidence["status"], "failed")
        self.assertEqual(evidence["received_signals"], [int(runner.signal.SIGTERM)])


class DispatchTests(unittest.TestCase):
    """Only synthetic files enter the supervisor. Git, process and clock are mocked."""

    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(TemporaryDirectory()))
        self.here = self.root / "verification/atlas"
        self.here.mkdir(parents=True)
        self.commit = "a" * 40
        self.out = self.root / "run-output/mock-run"
        self.config = json.loads((runner.HERE / "run-config.json").read_bytes())
        pins = {}
        for name in runner.ARCHIVE_PINS:
            file = self.root/name
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_bytes(("mock archive "+name).encode())
            pins[name] = runner.sha(file)
        self.stack.enter_context(patch.object(runner, "ARCHIVE_PINS", pins))
        self.config["regression_archive"]["paper1_paths"] = pins
        for name in runner.SOURCE_NAMES:
            (self.here / name).write_text("mock source " + name)
        (self.here / "run-config.json").write_text(json.dumps(self.config))
        self.committed = {path.relative_to(self.root).as_posix(): path.read_bytes()
                          for path in self.root.rglob("*") if path.is_file()}

        def git(args, **kwargs):
            self.assertEqual(args[0], "git")
            if args[3] == "rev-parse":
                return self.commit
            self.assertEqual(args[3], "show")
            full_id, path = args[4].split(":paper1/")
            self.assertIn(full_id, (self.commit, runner.ARCHIVE_COMMIT))
            return self.committed[path]

        self.stack.enter_context(patch.object(runner, "ROOT", self.root))
        self.stack.enter_context(patch.object(runner, "HERE", self.here))
        self.stack.enter_context(patch.object(runner.subprocess, "check_output", side_effect=git))
        self.popen = self.stack.enter_context(patch.object(runner.subprocess, "Popen", side_effect=AssertionError("unexpected child")))
        self.dispatch = self.stack.enter_context(patch.object(runner, "supervise", return_value=0))
        self.stack.enter_context(patch.object(runner, "report"))
        self.stack.enter_context(patch.object(runner.platform, "platform", return_value="mock host"))
        self.now = self.stack.enter_context(patch.object(runner, "datetime"))
        self.now.now.return_value = datetime(2026, 11, 28, tzinfo=timezone.utc)
        self.stack.enter_context(patch.dict(runner.os.environ, {"ORX_RUN_ID": "mock-run",
                            "ORX_COMMIT": self.commit, "PYTHONOPTIMIZE": "0"}))

    def test_only_snapshot_bytes_dispatched_and_fresh_directory_required(self):
        self.assertEqual(runner.main(), 0)
        out, record, env, cutoff = self.dispatch.call_args.args
        self.assertEqual((out, cutoff), (self.out, runner.CUTOFF))
        self.assertEqual((record["wall_seconds"], record["memory_bytes"]), (900, 2 * 1024**3))
        argv = record["command"]
        self.assertEqual(argv[:8], ["uv", "run", "--isolated", "--no-project", "--python", "3.12",
                                    "--with", "sympy==1.14.0"])
        frozen = out / "sources/verification/atlas"
        self.assertEqual(argv[argv.index("python") + 1:],
                         [str(frozen / "check_atlas.py"), "--config",
                          str(frozen / "run-config.json"), "--output", str(out / runner.CERTIFICATE)])
        self.assertEqual(record["checker_input_sha256"],
                         {"source_sha256": {name: runner.sha(frozen/name) for name in runner.SCIENCE_NAMES},
                          "config_sha256": runner.sha(frozen / "run-config.json")})
        self.assertEqual(env["PYTHONOPTIMIZE"], "0")
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "MKL_NUM_THREADS"):
            self.assertEqual(env[key], "1")
        errors = []
        runner.check_hashes_after(out, record, errors)
        self.assertEqual(errors, [])
        self.assertEqual(record["source_sha256_before"], record["source_sha256_after"])
        self.assertEqual(record["source_sha256_before"], record["source_snapshot_sha256_after"])
        (out / "sources" / runner.MAP_PATH).write_bytes(b"changed copy")
        (self.here / "check_atlas.py").write_bytes(b"changed original")
        runner.check_hashes_after(out, record, errors)
        self.assertEqual(len(errors), 2)
        with self.assertRaises(FileExistsError):
            runner.main()
        self.assertEqual(self.dispatch.call_count, 1)

    def check_preflight_failure(self):
        self.assertEqual(runner.main(), 1)
        self.dispatch.assert_not_called()
        self.popen.assert_not_called()
        evidence = json.loads((self.out / "execution.json").read_bytes())
        self.assertEqual(evidence["scientific_invocations"], 0)
        self.assertEqual(evidence["group_disposition"]["status"], "not_spawned")
        self.assertNotIn("checkpoint_column_count", evidence["metrics"])
        return evidence["supervisor_error"]

    def test_report_mismatch_cannot_dispatch(self):
        (self.root / runner.MAP_PATH).write_bytes(b"wrong report")
        self.assertIn("pin mismatch", self.check_preflight_failure())

    def test_modified_source_cannot_dispatch(self):
        (self.here / "check_atlas.py").write_bytes(b"wrong source")
        self.assertIn("pin mismatch", self.check_preflight_failure())

    def test_expired_window_cannot_dispatch(self):
        self.now.now.return_value = runner.CUTOFF
        self.assertIn("cutoff", self.check_preflight_failure())

    def test_missing_identity_cannot_dispatch(self):
        for key in ("ORX_RUN_ID", "ORX_COMMIT"):
            with self.subTest(key=key), patch.dict(runner.os.environ, {key: ""}):
                with self.assertRaises(ValueError):
                    runner.main()
        self.dispatch.assert_not_called()
        self.assertFalse(self.out.exists())


if __name__ == "__main__":
    unittest.main()
