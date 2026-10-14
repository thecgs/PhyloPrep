"""Dependency-free, throttled progress reporting for CLI and browser callbacks."""
from contextlib import contextmanager
from contextvars import ContextVar
import json
import sys
import time

_CURRENT = ContextVar('phyloprep_progress', default=None)


class ProgressReporter:
    def __init__(self, scope, total, unit='tasks', callback=None, interval=5.0, clock=None, stream=None):
        self.scope = scope
        self.total = total
        self.unit = unit
        self.completed = 0
        self.counts = {}
        self.label = 'Starting'
        self.callback = callback
        self.interval = interval
        self.clock = clock or time.monotonic
        self.stream = stream
        self.started = self.clock()
        self.last_update = float('-inf')

    def update(self, label=None, completed=None, status=None, force=True):
        now = self.clock()
        if not force and now - self.last_update < self.interval:
            return
        if label is not None:
            self.label = label
        if completed is not None:
            self.completed = completed
        if status:
            self.counts[status] = self.counts.get(status, 0) + 1
        elapsed = max(0.0, now - self.started)
        # ETA is only an estimate based on completed units; stage runtimes differ.
        eta = elapsed * (self.total - self.completed) / self.completed if self.completed else None
        event = dict(scope=self.scope, label=self.label, completed=self.completed,
                     total=self.total, unit=self.unit, elapsed_seconds=round(elapsed, 1),
                     eta_seconds=round(eta, 1) if eta is not None and self.unit == 'loci' else None,
                     counts=dict(self.counts))
        self.last_update = now
        if self.callback is not None:
            self.callback(json.dumps(event))
        counts = ', '.join(f'{key}={value}' for key, value in sorted(self.counts.items()))
        text = f'[{self.scope}] {self.completed}/{self.total} {self.unit} | {self.label} | elapsed {elapsed:.0f}s'
        if self.unit == 'loci' and self.completed:
            text += f' | estimated remaining {eta:.0f}s'
        if counts:
            text += f' | {counts}'
        print(text, file=self.stream if self.stream is not None else sys.stderr, flush=True)

    def tick(self):
        self.update(force=False)


@contextmanager
def progress_context(reporter):
    token = _CURRENT.set(reporter)
    try:
        yield reporter
    finally:
        _CURRENT.reset(token)


def heartbeat():
    reporter = _CURRENT.get()
    if reporter is not None:
        reporter.tick()


def report_phase(args, label):
    reporter = getattr(args, '_progress', None)
    if reporter is not None:
        reporter.update(label)
