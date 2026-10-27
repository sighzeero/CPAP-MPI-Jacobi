#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Genera A con diagonale N e tutti gli elementi fuori diagonale -1."""
from __future__ import print_function
import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=5000, help='Ordine della matrice (default: 5000)')
    parser.add_argument('--output', help='File .mtx (default: big/dense_N.mtx)')
    args = parser.parse_args()
    if args.size < 2:
        parser.error('--size deve essere almeno 2')
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output = os.path.abspath(args.output or os.path.join(root, 'big', 'dense_{0}.mtx'.format(args.size)))
    if not output.endswith('.mtx'):
        parser.error('--output deve terminare con .mtx')
    if os.path.exists(output):
        parser.error('File gia esistente: {0}'.format(output))
    # Un b preesistente cambierebbe il problema b=A*1 previsto per il test.
    if os.path.exists(output[:-4] + '_b.mtx'):
        parser.error('Esiste un file _b.mtx associato: scegliere un altro nome')
    parent = os.path.dirname(output)
    if not os.path.isdir(parent):
        os.makedirs(parent)
    n = args.size
    try:
        with open(output, 'w') as f:
            f.write('%%MatrixMarket matrix coordinate real symmetric\n')
            f.write('% Dense symmetric positive definite: diagonal N, off-diagonal -1\n')
            f.write('{0} {0} {1}\n'.format(n, n * (n + 1) // 2))
            # Si salva il triangolo inferiore: il loader ricostruisce il superiore.
            # Memoria del generatore costante, indipendente dal numero di elementi.
            for i in range(1, n + 1):
                for j in range(1, i):
                    f.write('{0} {1} -1\n'.format(i, j))
                f.write('{0} {0} {1}\n'.format(i, n))
    except BaseException:
        if os.path.exists(output):
            os.remove(output)
        raise
    print('Matrice creata: {0}'.format(output))
    print('Ordine: {0}; nonzeri nella matrice completa: {1}'.format(n, n * n))
    print('Dimensione del file: {0:.1f} MB'.format(os.path.getsize(output) / 1e6))
    print('Senza file b associato, il solver genera b=A*1: soluzione attesa x=1.')


if __name__ == '__main__':
    main()
