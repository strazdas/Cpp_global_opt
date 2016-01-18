#ifndef CONTE_H
#define CONTE_H 
/* CONTE algorithm purpose is to solve ASIMPL inner optimization problem, which
 * is to optimize simplex's lower bound by dividing it into subsimplexes. */
/* CONTinuous lower bound coverage to find its minimum Estimate ~ CONTE */  
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
// Depencencies are:   main  <  Disimplv  <  Simplex  <  Conte  <  Subsimplex


using namespace std;

class Conte {   // For efficiency do not store links to parents
    // Conte uses CS (ComponentS) and CS_LEN (number of components) consts from utils.h
    Conte(const Conte& other) {};
    Conte& operator=(const Conte& other) {};
public:
    Conte(vector<Point*> verts, vector<double> Ls, int crit_id) {
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
            // for (int j=0; j < _verts[i]->size(); j++){
            //     dist += pow(_verts[i]->_X[j] - p[j], 2);
            // };
            // dist = sqrt(dist);
            //// L1norm:
            for (int j=0; j < _verts[i]->size(); j++){
                dist += fabs(_verts[i]->_X[j] - p[j]);
            };

            cone_value = _verts[i]->_values[_crit_id] - _L*dist;
            if (lb_value < cone_value) {
                lb_value = cone_value;
            };
        };
        return lb_value;
    };

    Point* minimize() {       // Returns estimated simplex_lb_min_value
        double best_value;
        double best_p[_D];
        double p[_D];
        double coef;
        double value;

        for (int i=0; i < CS_LEN[_D-2]; i++) {   // Iterates thgrough coeficients
            for (int j=0; j < _V; j++) { p[j] = 0; };
            for (int k=0; k < _V; k++) {  // Iterates through coeficient components
                coef = CS[_D-2][_V*i + k];
                // Iterates through point components
                for (int j=0; j < _D; j++) {
                    p[j] += coef * _verts[k]->_X[j];
                };
            };
            value = get_lb_value(p);

            if (value < best_value) {
                best_value = value;
                for (int j=0; j < _V; j++) {
                    best_p[j] = p[j];
                };
            };
        };

        // Construct Point object
        Point* lb_estimate = new Point(best_p, _V);
        lb_estimate->add_value(best_value);
        return lb_estimate;
    };

    virtual ~Conte() {
    };
};

#endif
