#/bin/python
import os

exe = '/home/jgarcia/git/h2p/build/h2p-ricpad'

ndigits = 600
Dmax = 40
tol = "1E-90"
h = "1E-200"

points = [
        # U, R, A, m, s, Dmin
        [-1.1, 2, 0.8, 0, 0, 2],
        [-0.08, 8, -10.7, 1, 0, 4], 
        [-0.07, 8, -11, 2, 1, 4],
        [-0.06, 10, -29, 4, 1, 2],
        [-0.04, 10, -41, 3, 1, 5]
        ]

def get_points(U, A, R, m, s, Dmin, d, num, of):
    os.makedirs(num, exist_ok=True)

    os.chdir(num)
    command = exe + \
            (' fixed --Dmin {} --Dmax {} --d {} --U0 {} --R0 {} --A0 {}' +
            ' --m {} --s {} --ndigits {} --tol {} --h {}' +
            ' --output-file "{}"' + 
            ' --use-E').format(
                Dmin, Dmax, d, U, A, R, m, s, ndigits, tol, h, of
                )
    os.system(command)
    os.chdir('..')

for n, data in enumerate(points):
    for d in range(2):
        of = '{}-{}.log'.format(n, d)
        get_points(*data, d, str(n), of)
