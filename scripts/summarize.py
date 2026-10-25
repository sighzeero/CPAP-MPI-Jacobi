#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Riepilogo delle ripetizioni: compatibile con Python 2.7 e Python 3."""
from __future__ import print_function
import csv
import os
import sys

def median(values):
    ordered = sorted(values)
    half = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[half]
    return (ordered[half - 1] + ordered[half]) / 2.0

campaign = sys.argv[1]
with open(os.path.join(campaign, 'manifest.csv'), 'r') as f:
    runs = list(csv.DictReader(f))
times = {}
for run in runs:
    p = int(run['processes'])
    path = os.path.join(campaign, run['directory'], 'results', 'M{0}.txt'.format(p))
    with open(path, 'r') as f:
        lines = f.read().splitlines()
    if len(lines) != 1:
        raise SystemExit('Attesa una sola misura valida: {0}'.format(path))
    times.setdefault(p, []).append(float(lines[0].split('\t')[5]))
base = median(times[1])
path = os.path.join(campaign, 'summary.csv')
# csv su Python 2 richiede un file binario; Python 3 gestisce newline.
output = open(path, 'wb') if sys.version_info[0] == 2 else open(path, 'w', newline='')
with output as f:
    writer = csv.writer(f)
    writer.writerow(['processes', 'repeats', 'median_ms', 'min_ms', 'max_ms', 'speedup', 'efficiency'])
    for p, values in sorted(times.items()):
        elapsed = median(values)
        speedup = base / elapsed
        writer.writerow([p, len(values), elapsed, min(values), max(values), speedup, speedup / p])
print('Riepilogo: {0}'.format(path))
