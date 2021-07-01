#/bin/python
import os
import numpy as np
import numpy.polynomial as poly
import mpmath as mp
import argparse
import subprocess
import sys

# This is where the executable should be located
exe = '/home/jgarcia/git/h2p-ricpad/build/h2p-ricpad'

# Number of digits in the computations
ndigits = 150

# Tolerance for the Newton-Raphson method
tol = "1E-30"

# Value of h in the numerical differentiation. Should be at least sqrt(tol) or
# smaller
h = "1E-60"

# Other parameters
d = 0
Dmax0 = 10
Dmin0 = 2

# List of values of R that will be computed
Rlist = list(np.arange(0.1, 5, 0.1)) + list(np.arange(5, 10, 0.2)) + \
    list(np.arange(10, 20, 1)) + list(np.arange(20, 50, 2)) + \
    list(np.arange(50, 101, 5))

mp.mp.dps = ndigits

Dmax = Dmax0
Dmin = Dmin0

filename = sys.argv[1]

def get_points(U, A, R, m, s, Dmin, d, of):
    '''
    command = exe + \
            (' fixed  --Dmin {} --Dmax {} --d {} --U0 {} --A0 {} --R0 {}' +
            ' --m {} --s {} --ndigits {} --tol {} --h {}' +
            ' --output-file "{}"' + 
            ' --use-E --nr-max-iter 10000').format(
                Dmin, Dmax, d, U, A, R, m, s, ndigits, tol, h, of
                )
                '''
    args = [
            exe, 
            'fixed',
            f'--Dmin', f'{Dmin}',
            f'--Dmax', f'{Dmax}',
            f'--d', f'{d}',
            f'--U0', f'{U}',
            f'--A0', f'{A}',
            f'--R0', f'{R}',
            f'--m', f'{m}',
            f'--s', f'{s}',
            f'--ndigits', f'{ndigits}',
            f'--tol', f'{tol}',
            f'--h', f'{h}',
            f'--output-file', f'{of}',
            '--use-E',
            '--nr-max-iter', f'1000'
            ]

    print(' '.join(args))

    subprocess.run(args, check = True)


def analyze(filename):
    lines = open(filename, 'r').readlines()
    first = lines[:2]

    Us, As = [], []
    for line in first:
        data = line.strip().split()
        Us.append(mp.mpf(data[3]))
        As.append(mp.mpf(data[4]))

    ldu = abs((Us[1] - Us[0])/Us[1])
    lda = abs((As[1] - As[0])/As[1])

    # Raise a ValueError if the relative difference between the first two steps is too large
    if ldu > 0.5 or lda > 0.5:
        print(f' >> ldu = {mp.nstr(ldu,4)}, lda = {mp.nstr(lda,4)}')
        raise ValueError

    last = lines[-2:]

    Us, As = [], []
    for line in last:
        data = line.strip().split()
        Us.append(mp.mpf(data[3]))
        As.append(mp.mpf(data[4]))

    ldu = int(mp.floor(-mp.log10(abs(Us[1] - Us[0]))))
    lda = int(mp.floor(-mp.log10(abs(As[1] - As[0]))))

    return mp.nstr(Us[1], ldu), mp.nstr(As[1], ldu)

def sweep_state(l, m, I):
    global Dmin, Dmax 
    if m > l:
        raise ValueError('l should be greater or equal than m')

    nu = I + l
    ns = nu - int(np.ceil((l-m)/2))
    Nm = l - m
    Nl = I - 1
    s = Nm % 2

    U = -4/nu**2/2
    A = -l*(l+1)

    print(120*'-')
    print('Running with:')
    print('m = {}, l = {}, I = {}'.format(m, l, I))
    print('Nl = {}, Nm = {}, nu = {}, ns = {}'.format(Nl, Nm, nu, ns))
    print('E0 = {}, A0 = {}'.format(U, A))
    print(120*'-')
    print()

    of = '{}_{}_{}.tmp'.format(l,m,I)
    logfile = '{}_{}_{}.log'.format(l,m,I)

    Uvals = [U]
    Avals = [A]
    Rvals = [0]

    for R in Rlist:
        if R == 0:
            with open(logfile, 'a') as write_file:
                write_file.write('{:7.3f}    {}    {}\n'.format(R, U, A))

            continue

        ntries = 0
        # Try with the initial values; if the first two steps are too far apart from each other, try again with Dmin += 1.

        while ntries < 10:
            ntries += 1
            try:
                try:
                    os.remove(of)
                except OSError:
                    pass

                deg = len(Avals) - 1
                print('Uvals: ', Uvals)
                print('Avals: ', Avals)
                print('Rvals: ', Rvals)

                if deg >= 1 and deg < 3:
                    Apoly = poly.Polynomial.fit(Rvals, Avals, deg=deg)
                    Astart = Apoly(R)

                    Upoly = poly.Polynomial.fit(Rvals, Uvals, deg=deg)
                    Ustart = Upoly(R)
                else:
                    Astart = A
                    Ustart = U

                # This may fail with subprocess.CalledProcessError due to NR not converging
                try: 
                    get_points(Ustart, Astart, R, m, s, Dmin, d, of)
                except subprocess.CalledProcessError:
                    pass

                # This may fail with ValueError due to NR converging to wrong roots, or IndexError due to the file being empty
                Ustr, Astr = analyze(of)

                U = float(Ustr)
                A = float(Astr)

                ldu = max(abs((U - Ustart)/U), abs((U - Ustart)/Ustart))
                if Astart != 0 and deg == 2:
                    lda = max(abs((A - Astart)/A), abs((A - Astart)/Astart))
                else: 
                    lda = 0

                print(f'ldu = {ldu}, lda = {lda}')
                print('A = ', A, 'Astart = ', Astart)
                if ldu > 0.2 or lda > 0.6:
                    U = Uvals[-1]
                    A = Avals[-1]
                    raise ValueError

                Rvals.append(R)
                Uvals.append(U)
                Avals.append(A)

                if len(Rvals) > 3:
                    Rvals.pop(0)
                    Uvals.pop(0)
                    Avals.pop(0)

                with open(logfile, 'a') as write_file:
                    write_file.write('{:6.2f}    {}    {}\n'.format(R, Ustr, Astr))

                break
            except ValueError as e:
                print(e)
                print('Failed. Increasing Dmin.')
                Dmin += 1
                Dmax += 1

            except IndexError as e:
                print(e)
                print('IndexError caught. Increasing Dmin.')
                Dmin += 1
                Dmax += 1

        if ntries == 10:
            print(120*'-')
            print()
            print('At l = {}, m = {}, I = {}.'.format(l, m, I))
            print(
                 '10 successive tries have failed. Please reduce step size.'
                 )
            print()
            print(120*'-')

    try:
        os.remove(of)
    except OSError:
        pass

if __name__ == '__main__':
    for line in open(filename, 'r').readlines():
        Dmax = Dmax0
        Dmin = Dmin0

        if line[0] == '#': continue
        data = line.split()

        l = int(data[0])
        m = int(data[1])
        I = int(data[2])

        sweep_state(l, m, I)

