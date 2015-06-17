# should construct simplexes and call show potential function
from numpy import array as a, matrix as m, arange, sqrt, isnan, pi, cos, sin, mean
from itertools import permutations

from os import listdir

# List directories (exclude bin)
#   iterate through classes      filenames and summarise them

ignore = ['bin', '.ropeproject', 'show_stats.py']

def show_stats(root_path="."):
    for path in listdir(root_path):
        if path not in ignore:
            print('===== ' + path + ' =====')
            for cls in range(1, 9):

                stats = {'alg': '', 'cls': '', 'calls': [], 'subregions': [], 'duration': []}
                for fid in range(1, 101):
                    filename = str(cls) + "_" + str(fid)
                    try:
                        f = open(path + '/' + filename)
                    except:  # Temporarily
                        continue
                    file_content = f.read().strip()
                    for o in file_content.split(','):
                        if ':' in o:
                            key, value = [e.strip() for e in o.strip().split(':')]
                            if key == 'calls' or key == 'subregions':
                                stats[key].append(int(value))
                            if key == 'duration':
                                stats[key].append(float(value))
                print " ", cls,
                if stats['calls']:
                    l = len(stats['calls'])
                    stats['calls'] = sorted(stats['calls'])
                    stats['subregions'] = sorted(stats['subregions'])
                    print "fc50: %5d   fc100: %5d   calls: %9.3f   runs: %3d   parts50: %6d   parts100: %6d" % (stats['calls'][l/2], stats['calls'][-1], mean(stats['calls']), len(stats['calls']), stats['subregions'][l/2], stats['subregions'][-1])
                else:
                    print


if __name__ == '__main__':
    show_stats()
