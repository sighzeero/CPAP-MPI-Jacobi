#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Riepilogo delle ripetizioni: compatibile con Python 2.7 e Python 3."""
from __future__ import print_function
import csv
import os
import re
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
iterations = {}
converged = {}
configurations = {}
for run in runs:
    p = int(run['processes'])
    config = (run.get('nodes', ''), run.get('map_by', ''), run.get('bind_to', ''))
    if p in configurations and configurations[p] != config:
        raise SystemExit('Configurazioni diverse per {0} processi: non aggregare le misure.'.format(p))
    configurations[p] = config
    path = os.path.join(campaign, run['directory'], 'results', 'M{0}.txt'.format(p))
    with open(path, 'r') as f:
        lines = f.read().splitlines()
    if len(lines) != 1:
        raise SystemExit('Attesa una sola misura valida: {0}'.format(path))
    times.setdefault(p, []).append(float(lines[0].split('\t')[5]))
    with open(os.path.join(campaign, run['directory'], 'run.log'), 'r') as f:
        match = re.search(r'processes=\d+ iterations=(\d+) relative_update=\S+ converged=([01])', f.read())
    if not match:
        raise SystemExit('Stato del solver assente nel log: {0}'.format(run['directory']))
    iterations.setdefault(p, []).append(int(match.group(1)))
    converged.setdefault(p, []).append(int(match.group(2)))
base = median(times[1])
path = os.path.join(campaign, 'summary.csv')
# csv su Python 2 richiede un file binario; Python 3 gestisce newline.
output = open(path, 'wb') if sys.version_info[0] == 2 else open(path, 'w', newline='')
with output as f:
    writer = csv.writer(f)
    writer.writerow(['processes', 'repeats', 'median_ms', 'min_ms', 'max_ms', 'speedup', 'efficiency', 'min_iterations', 'max_iterations', 'converged_runs', 'nodes', 'map_by', 'bind_to'])
    for p, values in sorted(times.items()):
        elapsed = median(values)
        speedup = base / elapsed
        writer.writerow([p, len(values), elapsed, min(values), max(values), speedup, speedup / p, min(iterations[p]), max(iterations[p]), sum(converged[p])] + list(configurations[p]))
print('Riepilogo: {0}'.format(path))

if len(set(k for values in iterations.values() for k in values)) > 1:
    print("Attenzione: iterazioni diverse fra le run; il confronto non e a lavoro fisso.")
