# coding: utf-8
# A. Žilinsko apatinės ribos įvertinimo strategija.
#
# Simpleksams rasti apatinę ribą rasti geriausią žinomą reikšmę ir paklaidą
# Patikrink ar sprendinio tikslumas pakankamas. Paklaida: skirtumas tarp reikšmės ir apatinės ribos.
# Parink simpleksą tolimesniam dalinimui
# Padalink simpleksą
# Išmesk dominuojamus simpleksus (apatinė riba didesnė, nei žinoma reikšmė).
#
# Panaudoti A. Žilinsko simplekso apatinės ribos skaičiavimo strategiją.
# Dalinimas mažiausios apatinės ribos simplekso - optimaliai minimizuoja apatinę ribą.


from itertools import permutations
from datetime import datetime
from math import sqrt
from numpy import mean, array as a
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import numpy as np


def l2norm(a1, a2):
    '''Euclidean norm, which converts arguments to arrays automatically.'''
    if isinstance(a1, (int, float)):  # len(X) < 2:
        return abs(a(a1)-a(a2))
    return sqrt(sum([e**2 for e in (a(a1)-a(a2))]))


class Simplex:
    def __init__(self, verts=[], values=[], L=None):
        self.D = len(verts) - 1
        self.verts = verts
        self.values = values
        self.L = L
        self.min_value = min(values)
        self.diameter, self.le_vert1, self.le_vert2 = self.find_longest_edge(verts)
        self.min_lb = self.find_lower_bound_minimum()
        self.tolerance = abs(self.min_lb - self.min_value)
        self.min_value_id = values.index(self.min_value)
        self.max_value_id = values.index(max(values))

    def find_longest_edge(self, verts):
        max_len = None
        le_vert1, le_vert2 = None, None
        for i, j in permutations(range(len(verts)), 2):
            length = l2norm(verts[i], verts[j])
            if max_len == None or max_len < length:
                max_len = length
                le_vert1 = i
                le_vert2 = j
        return max_len, le_vert1, le_vert2

    def find_lower_bound_minimum(self):
        # Some how min lower bound is increasing
        # Draw this case and explain why
        '''A. Žilinsko apatinės ribos įvertinimo strategija'''
        i = self.le_vert1
        j = self.le_vert2
        return (self.values[i] + self.values[j]) / 2. - self.L * l2norm(self.verts[i], self.verts[j]) / sqrt(2)

    def divide(self, lb_func, min_value, tolerance):
        middle_point = (a(self.verts[self.le_vert1]) + a(self.verts[self.le_vert2])) / 2.
        middle_point_value = lb_func(middle_point)
        new_verts = self.verts[:self.le_vert1] + [middle_point] + self.verts[self.le_vert1 + 1:]
        new_values = self.values[:self.le_vert1] + [middle_point_value] + self.values[self.le_vert1 + 1:]
        new_simplex1 = Simplex(new_verts, new_values, self.L)
        new_verts = self.verts[:self.le_vert2] + [middle_point] + self.verts[self.le_vert2 + 1:]
        new_values = self.values[:self.le_vert2] + [middle_point_value] + self.values[self.le_vert2 + 1:]
        new_simplex2 = Simplex(new_verts, new_values, self.L)
        if min_value >= new_simplex1.min_value:
            min_value = new_simplex1.min_value
            tolerance = new_simplex1.tolerance
        if min_value >= new_simplex2.min_value:
            min_value = new_simplex2.min_value
            tolerance = new_simplex2.tolerance
        return [new_simplex1, new_simplex2], min_value, tolerance

    def __str__(self):
        verts = '\n'.join([str(self.verts[i][0]) + ' ' + str(self.verts[i][1]) + '  ->  ' + str(self.values[i]) for i in range(len(self.verts))])
        representation = 'min value: %s   min lower bound: %s    tolerance: %s\n' % (self.min_value, self.min_lb, self.tolerance)
        return representation + verts

def lower_bound_surface(func, X1, X2):
    LB = []
    for i, x1_row in enumerate(X1):
        row = []
        for j, x1 in enumerate(x1_row):
            row.append(func([x1, X2[i,j]]))
        LB.append(row)
    return a(LB)

class LBFunction:
    def __init__(self, simplex, L):
        self.simplex = simplex
        self.L = L

    def __call__(self, point):
        verts = self.simplex.verts
        values = self.simplex.values
        L = self.L
        return max([values[i] - L*l2norm(point, verts[i]) for i in range(simplex.D + 1)])


def draw_simplex_3d_euclidean_bounds(simplex, L=1., points=[]):
    '''2d->1d simplex lower bound surface.'''
    func = LBFunction(simplex, L)

    t = a([simplex.verts[0], simplex.verts[1], simplex.verts[2]])
    y = a([simplex.values[0], simplex.values[1], simplex.values[2]])   # (obj1, obj2) for A, B, C

    X1 = np.arange(min(t[:,0]), max(t[:,0]), (max(t[:,0]) - min(t[:,0]))/60.)
    X2 = np.arange(min(t[:,1]), max(t[:,1]), (max(t[:,1]) - min(t[:,1]))/60.)

    X1, X2 = np.meshgrid(X1, X2)
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    LB = lower_bound_surface(func, X1, X2)

    ax.plot_surface(X1, X2, LB[:,:], linewidth=0, rstride=1, cstride=1, cmap=cm.coolwarm)
    ax.plot_wireframe(np.hstack((t[:,0], t[0,0])), np.hstack(([t[:,1], t[0,1]])), np.hstack((y[:],y[0])), zorder=1000)  ## Line surface

    # points = [simplex[-1]['mins_ABC'][1]]
    for p in points:
        ax.plot([p[0]], [p[1]], [p[2]], 'go')

    # points = [[1.5892857095626787, 1.5892857194077434, 1.186929245703222]]
    # for p in points:
    #     ax.plot([p[0]], [p[1]], [p[2]], 'ro')
    plt.show()

def optimize(simplex, L, error=None):
    if error is None:
        error = 10**(-8 + 1.4 * simplex.D)
    # Simpleksams rasti apatinę ribą ir jos tikslumą (skirtumas tarp reikšmės ir apatinės ribos)
    min_value = simplex.min_value
    tolerance = simplex.tolerance   # Updates min_value and tolerance
    simplexes = [simplex]
    lb_func = LBFunction(simplex, L)

    min_L = ((simplex.values[simplex.max_value_id] - simplex.values[simplex.min_value_id]) /
                l2norm(simplex.verts[simplex.max_value_id], simplex.verts[simplex.min_value_id]))
    if L < min_L:
        raise ValueError("Given L=%f is too small for this simplex (min_L=%f)" % (L, min_L))

    # Apskaičiuoti Lipschitco konstantą iš Simplekso viršūnių ir jeigu ji yra
    # didesnė, nei nurodytoji, rodyti klaidą


    # Patikrinti ar sprendinio tikslumas pakankamas
    i = 0
    while tolerance > error:
        i += 1
        print(tolerance) #, min_value, simplexes)
        # print(tolerance, min_value, [simplex for simplex in simplexes if simplex.tolerance == tolerance])
        # Parinkti simpleksus tolimesniam dalinimui
        selected = min(simplexes, key=lambda x: x.min_lb)  # select_simplexes_for_partition(simplexes)
        # Padalinti simpleksus
        new_simplexes, min_value, tolerance = selected.divide(lb_func, min_value, tolerance)   # finds meta for newly created simplexes and updates min_value and tolerance

        simplexes.remove(selected)
        simplexes += new_simplexes

        # simplex_with_tolerance =  min(simplexes, key=lambda x: x.tolerance)
        # tolerance = simplex_with_tolerance.tolerance
        # print(simplex_with_tolerance.tolerance)

        # simplex_with_min_value = min(simplexes, key=lambda x: x.min_value)
        # print(simplex_with_min_value.min_value)

        # print('After adding', tolerance, min_value, simplexes)

        # Išmesti dominuojamus simpleksus (apatinė riba didesnė, nei žinoma reikšmė)
        # for simplex in simplexes:
        #     if simplex.min_lb > min_value:
        #         simplexes.remove(simplex)
        # print('At the loop end', tolerance, min_value, simplexes)
    min_value_simplex = min(simplexes, key=lambda x: x.min_value)
    vert_id = min_value_simplex.values.index(min_value)
    print(min_value_simplex.verts[vert_id], '  ', min_value_simplex.values[vert_id])
    print("Reikejo: %d" % (i,))
    return list(min_value_simplex.verts[vert_id]) + [min_value_simplex.values[vert_id]]
    # print 'Where the center bound is and its value'



# def lower_bound_function(simplex, L):
#     return function
#
# def find_simplex_min_lb(simplex, L):
#     return optimize(lb(simplex, L), [simplex], L)


if __name__ == '__main__':
    ## Simplex([(0, 0), (1, 1), (1, 0)], [2.60757, 3.26901, 5.02497], L)
    # 0 1  (0.218556); 0 0  (2.60757); 1 1  (3.26901);  (1.41421,0.0237495)
    ## Simplex([(0, 1), (0, 0), (1, 1)], [0.218556, 2.60757, 3.26901], L)
    # 0 1  (0.218556); 0.5 0.5  (0.938293); 0 0  (2.60757);  (1,-2.85539)
    #
    # 1 1 1  (1.61907); 1 0 0  (4.12837); 1 1 0  (5.181); 0 0 0  (7.69918);  (1.73205,-2.82922)
    # 1 0 1  (0.566444); 0.5 0.5 0.5  (1.65913); 1 0 0  (4.12837); 0 0 0  (7.69918);  (1.41421,-7.97019)
    # 1 0 1  (0.566444); 0.5 0.5 0.5  (1.65913); 0 0 1  (4.13725); 0 0 0  (7.69918);  (1.41421,-7.97019)

    L = 5
    start = datetime.now()
    # simplex = Simplex([(0, 0), (1, 1), (1, 0)], [2.60757, 3.26901, 5.02497], L)
    # simplex = Simplex([(0, 1), (0, 0), (1, 1)], [0.218556, 2.60757, 3.26901], L)
    # simplex = Simplex([(1, 1, 1), (1, 0, 0), (1, 1, 0), (0, 0, 0)], [1.61907, 4.12837, 5.181, 7.69918], L)
    simplex = Simplex([
        (1, 1, 1, 0, 0),
        (1, 0, 0, 0, 0),
        (1, 1, 0, 0, 1),
        (0, 0, 0, 0, 1),
        (0, 1, 0, 0, 0),
        (0, 1, 0, 0, 1),
    ],
    [1.61907, 1.81907, 4.12837, 5.181, 7.69918, 8.242], L)
    point = optimize(simplex, L)
    end = datetime.now()
    print(end - start)
    draw_simplex_3d_euclidean_bounds(simplex, L, [point])
