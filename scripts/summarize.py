#!/usr/bin/env python3
"""Riassume esclusivamente le ripetizioni della campagna specificata."""
import csv
import statistics
import sys
from pathlib import Path
campaign = Path(sys.argv[1])
with (campaign / 'manifest.csv').open() as f:
    runs = list(csv.DictReader(f))
times = {}
for run in runs:
    p = int(run['processes'])
    path = campaign / run['directory'] / 'results' / f'M{p}.txt'
    lines = path.read_text().splitlines()
    if len(lines) != 1:
        raise SystemExit(f'Attesa una sola misura valida: {path}')
    times.setdefault(p, []).append(float(lines[0].split('\t')[5]))
base = statistics.median(times[1])
with (campaign / 'summary.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['processes', 'repeats', 'median_ms', 'min_ms', 'max_ms', 'speedup', 'efficiency'])
    for p, values in sorted(times.items()):
        median = statistics.median(values)
        speedup = base / median
        writer.writerow([p, len(values), median, min(values), max(values), speedup, speedup / p])
print(f"Riepilogo: {campaign / 'summary.csv'}")
