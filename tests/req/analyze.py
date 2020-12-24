#!/bin/python
from mpmath import mp
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'h2.log'

mp.dps = 400

dD, dU, dA, dR = [], [], [], []

n = 0
U, A, R = 0, 0, 0
for line in open(filename, 'r').readlines():
    if '<<' in line: continue

    n += 1
    data = line.split()

    if n > 0:
        Uold = U
        Aold = A
        Rold = R

    D = int(data[2])
    U = mp.mpf(data[3])
    A = mp.mpf(data[4])
    R = mp.mpf(data[5])

    if n > 0:
        dD.append(D)
        dU.append(-mp.log10(abs(Uold - U)))
        dA.append(-mp.log10(abs(Aold - A)))
        dR.append(-mp.log10(abs(Rold - R)))

plt.plot(dD, dU, 'x', label = 'U')
plt.plot(dD, dA, '.', label = 'A')
plt.plot(dD, dR, 'o', label = 'R')

plt.legend()
plt.show()
