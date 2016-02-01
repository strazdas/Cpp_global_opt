#ifndef PROPCONTE_H
#define PROPCONTE_H
/* Proportional CONTE algorithm purpose is to solve ASIMPL inner optimization
 * problem, which is to optimize simplex's lower bound by dividing it into
 * subsimplexes. CONTinuous lower bound coverage to find its minimum Estimate */
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
// Depencencies are:   main  <  Disimplv  <  Simplex  <  PropConte  <  Subsimplex


using namespace std;

class PropConte {
    PropConte(const PropConte& other) {};
    PropConte& operator=(const PropConte& other) {};
public:
    PropConte(vector<Point*> verts, vector<double> Ls, int crit_id, double diameter) {
        _verts = verts;
        _L = Ls[crit_id];
        _D = _verts.size() - 1;
        _V = _verts.size();
        _crit_id = crit_id;
        _diameter = diameter;
        _min_step_size = 1. / 15;
        _steps = ceil(_diameter / (sqrt(_D) * _min_step_size));
        if (_steps < _V) {
            _steps = _V;
        };
        _step_size = 1. / _steps;
    };
    vector<Point*> _verts;   // Simplex vertexes
    double _L;
    int _D;
    int _V;
    int _crit_id;
    double _diameter;
    double _min_step_size;
    double _steps;
    double _step_size;

    double get_lb_value(double* p) {
        /* Returns lb value at a given point */
        double cone_value;
        double dist;
        double lb_value = -numeric_limits<double>::max();

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

    Point* minimize() {      // Returns estimated simplex_lb_min_value
        double best_value = numeric_limits<double>::max();  // min_lb_estimate
        double best_p[_D];   // min_lb_estimate point

        // Tmp vars
        double p[_D];
        double coef;
        double c[_V];
        double value;
        double left1, left2, left3, left4;

        if (_D == 2) {
            for (double v1=0; v1 <= 1; v1 += _step_size) {
                left1 = 1. - v1;
                for (double v2=0; v2 <= left1; v2 += _step_size) {
                    c[0] = v1;
                    c[1] = v2;
                    c[2] = left1 - v2;

                    // Point
                    for (int j=0; j < _V; j++) { p[j] = 0; };
                    for (int k=0; k < _V; k++) {  // Iterates through coeficient components
                        coef = c[k];
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
            };
        };
        if (_D == 3) {
            for (double v1=0; v1 <= 1; v1 += _step_size) {
                left1 = 1. - v1;
                for (double v2=0; v2 <= left1; v2 += _step_size) {
                    left2 = left1 - v2;
                    for (double v3=0; v3 <= left2; v3 += _step_size) {
                        c[0] = v1;
                        c[1] = v2;
                        c[2] = v3;
                        c[3] = left2 - v3;

                        // Point
                        for (int j=0; j < _V; j++) { p[j] = 0; };
                        for (int k=0; k < _V; k++) {  // Iterates through coeficient components
                            coef = c[k];
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
                };
            };
        };
        if (_D == 4) {
            for (double v1=0; v1 <= 1; v1 += _step_size) {
                left1 = 1. - v1;
                for (double v2=0; v2 <= left1; v2 += _step_size) {
                    left2 = left1 - v2;
                    for (double v3=0; v3 <= left2; v3 += _step_size) {
                        left3 = left2 - v3;
                        for (double v4=0; v4 <= left3; v4 += _step_size) {
                            c[0] = v1;
                            c[1] = v2;
                            c[2] = v3;
                            c[3] = v4;
                            c[4] = left3 - v4;

                            // Point
                            for (int j=0; j < _V; j++) { p[j] = 0; };
                            for (int k=0; k < _V; k++) {  // Iterates through coeficient components
                                coef = c[k];
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
                    };
                };
            };
        };
        if (_D == 5) {
            for (double v1=0; v1 <= 1; v1 += _step_size) {
                left1 = 1. - v1;
                for (double v2=0; v2 <= left1; v2 += _step_size) {
                    left2 = left1 - v2;
                    for (double v3=0; v3 <= left2; v3 += _step_size) {
                        left3 = left2 - v3;
                        for (double v4=0; v4 <= left3; v4 += _step_size) {
                            left4 = left3 - v4;
                            for (double v5=0; v5 <= left4; v5 += _step_size) {
                                c[0] = v1;
                                c[1] = v2;
                                c[2] = v3;
                                c[3] = v4;
                                c[4] = v5;
                                c[5] = left4 - v5;

                                // Point
                                for (int j=0; j < _V; j++) { p[j] = 0; };
                                for (int k=0; k < _V; k++) {  // Iterates through coeficient components
                                    coef = c[k];
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
                        };
                    };
                };
            };
        };

        // Construct Point object
        Point* lb_estimate = new Point(best_p, _V);
        lb_estimate->add_value(best_value);
        return lb_estimate;
    };

    virtual ~PropConte() {
    };
};

#endif
