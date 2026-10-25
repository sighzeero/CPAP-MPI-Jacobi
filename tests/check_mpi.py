#!/usr/bin/env python3
"""Soluzioni note e residui; --direct evita mpirun per il caso a un rank."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
root = Path(__file__).resolve().parents[1]
binary = root / 'src/jacobi_mpi_bin'
direct = '--direct' in sys.argv

def verify(matrix, expected, coefficients, rhs, cwd, ranks):
    command = [str(binary), str(matrix), 'M', str(ranks), '1']
    if not direct:
        command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-np', str(ranks)] + command
    run = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=60)
    if run.returncode:
        raise AssertionError(f'{command}\n{run.stdout}\n{run.stderr}')
    line = next(line for line in run.stdout.splitlines() if line.startswith('Results:'))
    x = [float(v) for v in line.split(':', 1)[1].split()]
    assert len(x) == len(expected)
    error = max(abs(a-b) for a,b in zip(x, expected))
    residual = max(abs(sum(a*b for a,b in zip(row,x))-b) for row,b in zip(coefficients,rhs))
    assert error < 1e-4 and residual < 1e-3, (ranks,error,residual,run.stdout)
    print(f'OK {matrix.name}: rank={ranks}, errore={error:.3g}, residuo={residual:.3g}')

with tempfile.TemporaryDirectory(prefix='jacobi-test-') as tmp:
    cwd=Path(tmp); (cwd/'results').mkdir()
    n=10
    a=[[10.0 if i==j else 1.0 if abs(i-j)==1 else 0.0 for j in range(n)] for i in range(n)]
    for p in ([1] if direct else [1,2,4]):
        verify(root/'big/test10.mtx',[1.0]*n,a,[sum(row) for row in a],cwd,p)
    # Dipendenze verso rank non adiacenti, matrice generale e b esplicito.
    n=12
    a=[[0.0]*n for _ in range(n)]
    for i in range(n):
        a[i][i]=10.0
        a[i][(i+7)%n]=0.5
        a[i][(i+1)%n]=-1.0
    expected=[float(i+1) for i in range(n)]
    rhs=[sum(v*x for v,x in zip(row,expected)) for row in a]
    matrix=cwd/'general.mtx'
    entries=[f'{i+1} {j+1} {a[i][j]:.10f}\n' for i in range(n) for j in range(n) if a[i][j]]
    matrix.write_text('%%MatrixMarket matrix coordinate real general\n'+f'{n} {n} {len(entries)}\n'+''.join(entries))
    (cwd/'general_b.mtx').write_text(f'{n} 1\n'+'\n'.join(map(str,rhs))+'\n')
    for p in ([1] if direct else [1,2,4]):
        verify(matrix,expected,a,rhs,cwd,p)
