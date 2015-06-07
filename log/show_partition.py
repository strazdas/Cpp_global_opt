# should construct simplexes and call show potential function
from numpy import array as a, matrix as m, arange, sqrt, isnan, pi, cos, sin, mean
from itertools import permutations


def show_partition(filename='partition.txt'):
    from matplotlib import pyplot as plt
    draw_from_iteration = 24
    iteration = 0
    f = open(filename)
    simplexes = []
    selected_mode = False
    selected = []
    title = ''
    # ok = False
    for line in f:
        # if 'Iteration 56' in line:
        #     ok = True
        # if ok:
        if 'Iteration' in line:
            # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,6))
            # ax1.axis([-0.05, 1.05, -0.05, 1.05])
            title = line
            iteration = int(line.strip().strip('Iteration:').strip())
            continue
        if 'Selected' in line:
            selected_mode = True
            continue
            # Ignore till empty line found
        parts = line.split(';')
        simplex = []

        if line == '\n':
            selected_mode = False
            if simplexes and iteration >= draw_from_iteration:
                show_potential(simplexes, selected, title=title)
            simplexes = []
            selected = []
            continue
                    # [[0.75, 0.5], [0.875, 0.375], [1.0, 0.5]]
                    # plt.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'b-')
        else:
            if not selected_mode:
                add_to = simplexes
            else:
                add_to = selected
            for part in parts:
                if ',' not in part:
                    simplex.append([float(e) for e in part.split() if not ('(') in e])
                    if '(' in part:
                        simplex[-1].append({'obj': (float(part.split()[-1].strip().strip('()')),)})
                if '(' in part and ',' in part:
                    size, value = part.strip().strip('()').split(',')
                    simplex.append({'size': float(size), 'value': float(value)})
            if simplex:
                add_to.append(simplex)

    show_potential(simplexes, selected, title=title)


def l2norm(a1, a2):
    '''Euclidean norm, which converts arguments to arrays automatically.'''
    if isinstance(a1, (int, float)):  # len(X) < 2:
        return abs(a(a1)-a(a2))
    return sqrt(sum([e**2 for e in (a(a1)-a(a2))]))


def show_potential(simplexes, selected=[], show=True, title=''):
    from matplotlib import pyplot as plt
    # Draw two plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,6))
    fig.suptitle(title)

    for simplex in simplexes:
        ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'bo')
    for simplex in selected:
        ax2.plot([simplex[-1]['size']], [simplex[-1]['value']], 'ro')
    for i in range(len(selected[:-1])):
        ax2.plot([selected[i][-1]['size'], selected[i+1][-1]['size']],
                 [selected[i][-1]['value'], selected[i+1][-1]['value']], 'r-')
    ax2.set_ylabel('Function values on simplices vertices')
    ax2.set_xlabel('Simplices diameter')

    for simplex in simplexes:
        s = simplex[:-1]
        for j in range(3):
            ax1.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'b-')

    for simplex in selected:
        s = sort_vertexes_longest_edge_first(simplex)[:-1]
        for j in range(3):
            ax1.plot([s[j-1][0], s[j][0]], [s[j-1][1], s[j][1]], 'r-', linewidth=2)

        edge_lengths = []   # [(vertex_index, vertex_index, edge_length),]
        for i, j in permutations(range(len(s)), 2):
            if j > i:
                edge_lengths.append((i, j, l2norm(s[i][:-1], s[j][:-1])))
        le_i, le_j, le_length = max(edge_lengths, key=lambda x: x[-1])

        division_point = [(s[le_i][0]+s[le_j][0])/2., (s[le_i][1]+s[le_j][1])/2.]
        # ax1.plot([division_point[0], simplex[2][0]], [division_point[1], simplex[2][0]], 'r-')
        ax1.plot([division_point[0]], [division_point[1]], 'ro')
        ax1.plot([division_point[0], s[2][0]], [division_point[1], s[2][1]], 'r--')
        # Get longest edge, divide it
        # division line needed

    for simplex in simplexes:
        for j in range(3):
            ax1.plot([simplex[j][0]], [simplex[j][1]], 'bo')

    # ax2.axis([min([simplexes]) -0.05, 1.05, -0.05, 1.05])
    max_size = max([s[-1]['size'] for s in simplexes])
    ax2.set_xlim([-0.05, max_size + 0.05])
    ax1.axis([-0.05, 1.05, -0.05, 1.05])
    if show:
        plt.show()
    return ax1, ax2


def sort_vertexes_longest_edge_first(simplex):
    '''nD->nD Moves longest edge vertexes to the simplex vertex list beginning.'''
    # Find simplex edges lengths
    edge_lengths = []   # [(vertex_index, vertex_index, edge_length),]
    for i, j in permutations(range(len(simplex[:-1])), 2):
        if j > i:
            edge_lengths.append((i, j, l2norm(simplex[i][:-1], simplex[j][:-1])))


    # Get longest edge vertexes ids
    le_i, le_j, le_length = max(edge_lengths, key=lambda x: x[-1])

    # Move longest edge vertexes to simplex vertex list beginning
    vi = simplex[le_i]
    vj = simplex[le_j]
    simplex.remove(vi)
    simplex.remove(vj)
    simplex.insert(0, vj)
    simplex.insert(0, vi)
    return simplex

if __name__ == '__main__':
    show_partition()
