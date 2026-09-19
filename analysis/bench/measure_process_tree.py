#!/usr/bin/env python3
"""Run a command and sample aggregate memory for its Linux process session."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--interval",
        type=float,
        default=0.1,
        help="Sampling interval in seconds (default: 0.1).",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ["--"]:
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if args.interval <= 0:
        parser.error("--interval must be positive")
    return args


def read_stat(pid: int) -> tuple[str, int, int] | None:
    """Return (state, process group, session) from /proc/PID/stat."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        fields = stat[stat.rfind(")") + 2 :].split()
        return fields[0], int(fields[2]), int(fields[3])
    except (OSError, IndexError, ValueError):
        return None


def read_memory_kib(pid: int) -> tuple[int, int]:
    """Return (RSS, PSS) in KiB; PSS is zero when unavailable."""
    rss_kib = 0
    pss_kib = 0
    try:
        for line in Path(f"/proc/{pid}/smaps_rollup").read_text().splitlines():
            if line.startswith("Rss:"):
                rss_kib = int(line.split()[1])
            elif line.startswith("Pss:"):
                pss_kib = int(line.split()[1])
    except (OSError, ValueError):
        try:
            for line in Path(f"/proc/{pid}/status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    rss_kib = int(line.split()[1])
                    break
        except (OSError, ValueError):
            pass
    return rss_kib, pss_kib


def live_session_processes(session_id: int) -> list[tuple[int, int]]:
    """Return (PID, process group) for live processes in the session."""
    processes: list[tuple[int, int]] = []
    try:
        proc_entries = os.scandir("/proc")
    except OSError:
        return processes

    with proc_entries:
        for entry in proc_entries:
            if not entry.name.isdigit():
                continue
            pid = int(entry.name)
            stat = read_stat(pid)
            if stat is None or stat[2] != session_id:
                continue
            if stat[0] in {"X", "x", "Z"}:
                continue
            processes.append((pid, stat[1]))
    return processes


def sample_session(session_id: int) -> tuple[int, int, int]:
    rss_kib = 0
    pss_kib = 0
    processes = live_session_processes(session_id)
    for pid, _process_group in processes:
        process_rss, process_pss = read_memory_kib(pid)
        rss_kib += process_rss
        pss_kib += process_pss
    return rss_kib, pss_kib, len(processes)


def main() -> int:
    args = parse_args()
    started = time.monotonic()
    process = subprocess.Popen(args.command, start_new_session=True)
    session_id = process.pid

    def forward_signal(signum: int, _frame: object) -> None:
        process_groups = {
            process_group
            for _pid, process_group in live_session_processes(session_id)
        }
        for process_group in process_groups:
            try:
                os.killpg(process_group, signum)
            except (PermissionError, ProcessLookupError):
                pass

    for signum in (signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, forward_signal)

    peak_rss_kib = 0
    peak_pss_kib = 0
    max_processes = 0
    samples = 0

    while True:
        leader_status = process.poll()
        rss_kib, pss_kib, process_count = sample_session(session_id)
        peak_rss_kib = max(peak_rss_kib, rss_kib)
        peak_pss_kib = max(peak_pss_kib, pss_kib)
        max_processes = max(max_processes, process_count)
        samples += 1
        if leader_status is not None and process_count == 0:
            break
        time.sleep(args.interval)

    elapsed = time.monotonic() - started
    print(
        "process_tree_memory"
        f" peak_rss_kib={peak_rss_kib}"
        f" peak_pss_kib={peak_pss_kib}"
        f" max_processes={max_processes}"
        f" samples={samples}"
        f" interval_seconds={args.interval:.3f}"
        f" elapsed_seconds={elapsed:.3f}"
        f" exit_status={process.returncode}",
        flush=True,
    )
    return int(process.returncode or 0)


if __name__ == "__main__":
    sys.exit(main())
