#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
"""Soluzioni note e residui; --direct evita mpirun per il caso a un rank."""
import os
import shutil
import threading
import shlex
import subprocess
import sys
import tempfile
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
binary = os.path.join(root, 'src', 'jacobi_mpi_bin')
direct = '--direct' in sys.argv

def verify(matrix, expected, coefficients, rhs, cwd, ranks):
    command = [str(binary), str(matrix), 'M', str(ranks), '1']
    if not direct:
        command = ['mpirun'] + shlex.split(os.environ.get('MPI_TEST_FLAGS', '--bind-to none')) + ['-np', str(ranks)] + command
    run = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True)
    expired = threading.Event()
    def terminate():
        expired.set()
        try:
            run.kill()
        except OSError:
            pass
    timer = threading.Timer(60, terminate)
    timer.start()
    try:
        stdout, stderr = run.communicate()
    finally:
        timer.cancel()
    if expired.is_set() or run.returncode:
        raise AssertionError('{0}\n{1}\n{2}'.format(command, stdout, stderr))
    line = next(line for line in stdout.splitlines() if line.startswith('Results:'))
    x = [float(v) for v in line.split(':', 1)[1].split()]
    assert len(x) == len(expected)
    error = max(abs(a-b) for a,b in zip(x, expected))
    residual = max(abs(sum(a*b for a,b in zip(row,x))-b) for row,b in zip(coefficients,rhs))
    assert error < 1e-4 and residual < 1e-3, (ranks,error,residual,stdout)
    print('OK {0}: rank={1}, errore={2:.3g}, residuo={3:.3g}'.format(
        os.path.basename(matrix), ranks, error, residual))

tmp = tempfile.mkdtemp(prefix='jacobi-test-')
try:
    cwd = tmp
    os.mkdir(os.path.join(cwd, 'results'))
    n=10
    a=[[10.0 if i==j else 1.0 if abs(i-j)==1 else 0.0 for j in range(n)] for i in range(n)]
    for p in ([1] if direct else [1,2,4]):
        verify(os.path.join(root, 'big', 'test10.mtx'),[1.0]*n,a,[sum(row) for row in a],cwd,p)
    # Dipendenze verso rank non adiacenti, matrice generale e b esplicito.
    n=12
    a=[[0.0]*n for _ in range(n)]
    for i in range(n):
        a[i][i]=10.0
        a[i][(i+7)%n]=0.5
        a[i][(i+1)%n]=-1.0
    expected=[float(i+1) for i in range(n)]
    rhs=[sum(v*x for v,x in zip(row,expected)) for row in a]
    matrix = os.path.join(cwd, 'general.mtx')
    entries=['{0} {1} {2:.10f}\n'.format(i+1, j+1, a[i][j]) for i in range(n) for j in range(n) if a[i][j]]
    with open(matrix, 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n' +
                '{0} {0} {1}\n'.format(n, len(entries)) + ''.join(entries))
    with open(os.path.join(cwd, 'general_b.mtx'), 'w') as f:
        f.write('{0} 1\n'.format(n) + '\n'.join(map(str, rhs)) + '\n')
    for p in ([1] if direct else [1,2,4]):
        verify(matrix,expected,a,rhs,cwd,p)

finally:
    shutil.rmtree(tmp)
