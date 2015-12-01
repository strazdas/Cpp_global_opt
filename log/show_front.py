# coding: utf-8
# should construct simplexes and call show potential function
import sys
from numpy import array as a, matrix as m, arange, sqrt, isnan, pi, cos, sin, mean
from itertools import permutations
from mpl_toolkits.mplot3d import axes3d


def show_front(filename='front.txt'):
    from matplotlib import pyplot as plt
    f = open(filename)
    reading_front = True
    front = []
    tols = []
    for line in f:
        if 'Tolerance' in line:
            reading_front = False
            continue
        if reading_front:
            vals = line.split('->')[-1]
            p = []
            for v in vals.split():
                if v.strip():
                    p.append(float(v.strip()))
            front.append(p)
        else:
            p = []
            for v in line.split(' '):
                if v.strip():
                    p.append(float(v.strip()))
            tols.append(p)

    for p in front:
        plt.plot(p[0], p[1], 'bo')
    for p in tols:
        plt.plot(p[0], p[1], 'ro')
    plt.show()


def l2norm(a1, a2):
    '''Euclidean norm, which converts arguments to arrays automatically.'''
    if isinstance(a1, (int, float)):  # len(X) < 2:
        return abs(a(a1)-a(a2))
    return sqrt(sum([e**2 for e in (a(a1)-a(a2))]))


if __name__ == '__main__':
    if len(sys.argv) == 2:
        show_front(sys.argv[1])
    else:
        show_front()
