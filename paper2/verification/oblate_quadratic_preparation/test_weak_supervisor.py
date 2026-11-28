import errno
from pathlib import Path
import signal
import subprocess
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

import supervise_weak_validation as supervisor


class FakeChild:
    def __init__(self, alive_polls=0):
        self.alive_polls = alive_polls
        self.polls = 0

    def poll(self):
        self.polls += 1
        if self.polls <= self.alive_polls:
            return None
        return -signal.SIGTERM


class FakeClock:
    def __init__(self):
        self.value = 0.0

    def monotonic(self):
        self.value += 1.0
        return self.value

    def sleep(self, seconds):
        self.value += seconds


class SupervisorCleanupTests(unittest.TestCase):
    def test_empty_group_permission_probe_is_treated_as_empty_race(self):
        events = []
        with mock.patch.object(supervisor, "process_group_members", return_value=[]), \
             mock.patch.object(
                 supervisor.os, "killpg",
                 side_effect=PermissionError(errno.EPERM, "Operation not permitted"),
             ):
            exists, members = supervisor.probe_group(1234, events)

        self.assertFalse(exists)
        self.assertEqual(members, [])
        self.assertEqual(events[-1]["status"], "permission_error_empty")

    def test_live_group_permission_probe_remains_live(self):
        events = []
        members = [{"pid": 7, "rss_bytes": 4096}]
        with mock.patch.object(supervisor, "process_group_members", return_value=members), \
             mock.patch.object(
                 supervisor.os, "killpg",
                 side_effect=PermissionError(errno.EPERM, "Operation not permitted"),
             ):
            exists, observed = supervisor.probe_group(1234, events)

        self.assertTrue(exists)
        self.assertEqual(observed, members)
        self.assertEqual(events[-1]["status"], "permission_error_live")

    def test_empty_group_permission_signal_records_empty_race(self):
        events = []
        with mock.patch.object(supervisor, "process_group_members", return_value=[]), \
             mock.patch.object(
                 supervisor.os, "killpg",
                 side_effect=PermissionError(errno.EPERM, "Operation not permitted"),
             ):
            sent = supervisor.signal_group(1234, signal.SIGKILL, events)

        self.assertFalse(sent)
        self.assertEqual(events[-1]["status"], "permission_error_empty")
        self.assertEqual(events[-1]["members"], [])

    def test_permission_signal_records_unknown_when_member_lookup_fails(self):
        events = []
        with mock.patch.object(
            supervisor, "process_group_members",
            side_effect=subprocess.TimeoutExpired(["ps"], 5),
        ), mock.patch.object(
            supervisor.os, "killpg",
            side_effect=PermissionError(errno.EPERM, "Operation not permitted"),
        ):
            sent = supervisor.signal_group(1234, signal.SIGKILL, events)

        self.assertFalse(sent)
        self.assertEqual(events[-2]["status"], "failed")
        self.assertEqual(events[-1]["status"], "permission_error_unknown")
        self.assertIsNone(events[-1]["members"])

    def test_permission_probe_records_unknown_when_member_lookup_fails(self):
        events = []
        with mock.patch.object(
            supervisor, "process_group_members",
            side_effect=subprocess.TimeoutExpired(["ps"], 5),
        ):
            exists, members = supervisor.probe_group(1234, events)

        self.assertTrue(exists)
        self.assertIsNone(members)
        self.assertEqual(events[-1]["status"], "unknown")

    def test_normal_termination_stops_after_empty_group(self):
        events = []
        fake_clock = FakeClock()
        members = [
            [{"pid": 7, "rss_bytes": 4096}],
            [],
        ]

        with mock.patch.object(supervisor, "process_group_members", side_effect=members), \
             mock.patch.object(
                 supervisor.os, "killpg",
                 side_effect=[None, None, ProcessLookupError],
             ) as killpg, \
             mock.patch.object(supervisor.time, "monotonic", fake_clock.monotonic), \
             mock.patch.object(supervisor.time, "sleep", fake_clock.sleep):
            live, final_members = supervisor.terminate_group(
                1234, FakeChild(), events, grace_seconds=5
            )

        self.assertFalse(live)
        self.assertEqual(final_members, [])
        self.assertEqual(killpg.call_args_list[0].args, (1234, signal.SIGTERM))
        self.assertNotIn(
            (1234, signal.SIGKILL),
            [call.args for call in killpg.call_args_list],
        )

    def test_escalation_records_sigkill_when_group_stays_live(self):
        events = []
        fake_clock = FakeClock()
        live_members = [{"pid": 7, "rss_bytes": 4096}]

        with mock.patch.object(
            supervisor, "process_group_members", return_value=live_members
        ), mock.patch.object(supervisor.os, "killpg") as killpg, \
             mock.patch.object(supervisor.time, "monotonic", fake_clock.monotonic), \
             mock.patch.object(supervisor.time, "sleep", fake_clock.sleep):
            live, final_members = supervisor.terminate_group(
                1234, FakeChild(alive_polls=10), events, grace_seconds=3
            )

        self.assertTrue(live)
        self.assertEqual(final_members, live_members)
        self.assertIn(
            (1234, signal.SIGTERM),
            [call.args for call in killpg.call_args_list],
        )
        self.assertIn(
            (1234, signal.SIGKILL),
            [call.args for call in killpg.call_args_list],
        )

    def test_final_disposition_records_persistent_live_denial(self):
        events = []
        fake_clock = FakeClock()
        live_members = [{"pid": 7, "rss_bytes": 4096}]

        with mock.patch.object(
            supervisor, "process_group_members", return_value=live_members,
        ), mock.patch.object(
            supervisor.os, "killpg",
            side_effect=PermissionError(errno.EPERM, "Operation not permitted"),
        ), mock.patch.object(supervisor.time, "sleep", fake_clock.sleep):
            live, members, status = supervisor.final_disposition(
                1234, FakeChild(alive_polls=10), events
            )

        self.assertTrue(live)
        self.assertEqual(members, live_members)
        self.assertEqual(status, "live")
        self.assertEqual(events[-1]["status"], "live")
        self.assertTrue(
            any(event["status"] == "permission_error_live" for event in events)
        )


if __name__ == "__main__":
    unittest.main()
