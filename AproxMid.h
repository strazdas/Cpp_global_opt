#ifndef APROXMID_H
#define APROXMID_H
/* Approximate min_lb_point by averaging weighted verts (good for Euclidean space */
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include <sys/time.h>
#include <vector>
#include "utils.h"
#include "subsimplex.h"
#include <math.h>
// Depencencies:   main  <  Algorithm  <  Simplex  <  AproxMid  <  Subsimplex


using namespace std;

class AproxMid {
    AproxMid(const AproxMid& other) {};
    AproxMid& operator=(const AproxMid& other) {};
public:
    AproxMid(vector<Point*> verts, vector<double> Ls, int crit_id) {
        // Note: Should pass simplex here and use store dist between verts in it
        _verts = verts;
        _L = Ls[crit_id];
        _D = _verts.size() - 1;
        _V = _verts.size();
        _crit_id = crit_id;
    };
    vector<Point*> _verts;   // Simplex vertexes
    double _L;
    int _D;
    int _V;
    int _crit_id;

    double get_lb_value(double* p) {
        /* Returns lb value at a given point */
        double cone_value;
        double dist;
        double lb_value = -numeric_limits<double>::max();;

        for (int i=0; i < _verts.size(); i++) {
            // dist = l2norm(point, verts[i]);
            dist = 0;
            //// L2norm:
            for (int j=0; j < _verts[i]->size(); j++){
                dist += pow(_verts[i]->_X[j] - p[j], 2);
            };
            dist = sqrt(dist);
            //// L1norm:
            // for (int j=0; j < _verts[i]->size(); j++){
            //     dist += fabs(_verts[i]->_X[j] - p[j]);
            // };

            cone_value = _verts[i]->_values[_crit_id] - _L*dist;
            if (lb_value < cone_value) {
                lb_value = cone_value;
            };
        };
        return lb_value;
    };

    Point* minimize() {       // Returns estimated simplex_lb_min_value
        double X[_D];
        double value = numeric_limits<double>::max();

        sort(_verts.begin(), _verts.end(), Point::compare_by_value);

        // combinatorically select vert pairs
        double dist_ab, ca, cb;
        double pairs_count = 0;
        for (int a=0; a <= _D; a++)
            for (int b=a; b <= _D; b++)
                if (a != b) {
                    pairs_count += 1;
                    dist_ab = l2norm(_verts[a], _verts[b]);    // Distance in variable space

                    ca = 0.5 + (_verts[b]->_values[_crit_id] - _verts[a]->_values[_crit_id]) / (2 * _L * dist_ab);
                    cb = 0.5 - (_verts[b]->_values[_crit_id] - _verts[a]->_values[_crit_id]) / (2 * _L * dist_ab);

                    for (int k=0; k < _D; k++) {
                        X[k] += ca * _verts[a]->_X[k];
                        X[k] += cb * _verts[b]->_X[k];
                    };
                };

        // Divide X by pairs_count
        for (int k=0; k < _D; k++) {
            X[k] = X[k] / pairs_count;
        };
        Point* mid = new Point(X, _D);
        mid->add_value(_verts[0]->_values[_crit_id] - _L * l2norm(_verts[0], mid));  // dist to min_vert point  or  value in mid point
        return mid;
    };

    virtual ~AproxMid() {
    };
};

#endif
