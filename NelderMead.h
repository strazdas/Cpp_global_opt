#ifndef NELDERMEAD_H
#define NELDERMEAD_H
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
#include "Eigen/Dense"
#include "Eigen/LU"
// Depencencies are:   main  <  Disimplv  <  Simplex  <  NelderMead  <  Subsimplex


using namespace std;

class NelderMead {
    NelderMead(const NelderMead& other) {};
    NelderMead& operator=(const NelderMead& other) {};
public:
    NelderMead(vector<Point*> verts, vector<double> Ls, int crit_id, double diameter) {
        _verts = verts;
        _L = Ls[crit_id];
        _D = _verts.size() - 1;
        _V = _verts.size();
        _crit_id = crit_id;

        // Nelder-Mead parameters
        _step = 0.1 * diameter;
        _no_improve_thr = 10e-6;
        _no_improv_break = 10;
        _max_iters = 0;
        _iters = 0;
        _alpha = 1.;
        _gamma = 2.;
        _rho = -0.5;
        _sigma = 0.5;

        // is_in_simplex method cached variables
        _last_vert = _verts[_D];
        _p.resize(_D);
        _T.resize(_D, _D);
        for (int i=0; i < _D; i++) {
            for (int j=0; j < _D; j++) {
                _T(j, i) = _verts[i]->_X[j] - _last_vert->_X[j];
            };
        };
        _T = _T.inverse();
    };
    vector<Point*> _verts;     // Simplex vertexes
    Point* _last_vert;
    double _L;
    int _D;
    int _V;
    int _crit_id;

    double _step;
    double _no_improve_thr;
    double _no_improv_break;
    double _max_iters;       // Terminate NelderMead if reached this iteration (0 - ignore)
    double _iters;
    double _alpha;
    double _gamma;
    double _rho;
    double _sigma;

    Eigen::MatrixXd _T;      // Cache of Already transposed _T
    Eigen::VectorXd _p;      // Tmp var
    Eigen::VectorXd _lamdas; // Tmp var

    static double _max_iteration;   // Max iteration which was reached

    double get_lb_value(double* p) {
        /* Returns lb value at a given point */
        double cone_value;
        double dist;
        double lb_value = -numeric_limits<double>::max();

        if (!is_in_simplex(p)) {
            return numeric_limits<double>::max();
        };

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

    bool is_in_simplex(double* p) {
        // https://math.stackexchange.com/questions/1226707/how-to-check-if-point-x-in-mathbbrn-is-in-a-n-simplex
        for (int i=0; i < _D; i++) {
            _p(i) = p[i] - _last_vert->_X[i];
        };

        Eigen::VectorXd _lamdas = _T * _p;
        double lamdas_sum = 0;
        for (int i=0; i < _lamdas.size(); i++) {
            if (_lamdas[i] < 0) {
                return false;
            };
            lamdas_sum += _lamdas[i];
        };
        if (lamdas_sum > 1.0) {
            return false;
        };
        return true;
    };

    static double ascending_first_value(double* v1, double* v2) {
        return v1[0] < v2[0];
    };

    Point* to_point(double* x) {
        if (_iters > NelderMead::_max_iteration) {
            NelderMead::_max_iteration = _iters;
        };
        Point* best_point = new Point(&x[1], _D);
        best_point->add_value(x[0]);
        return best_point;
    };

    Point* minimize() {
        //// Starting point strategies:
        // {+} Imti suvidurkintas viršūnes kaip pradinį tašką
        // { } Viršūnės gali sutapti su pradinėmis viršūnėmis
        // { } Imti dvimačius sankirtų taškus kaip pradinį artinį
        double x_start[_V];
        for (int i=0; i < _V; i++) { // Verts
            for (int j=0; j < _D; j++) { // Coords
                x_start[j+1] += _verts[i]->_X[j] / _V;
            };
        };
        x_start[0] = get_lb_value(x_start+1);

        double prev_best = x_start[0];
        int no_improv = 0;
        vector<double*> res;    // list of [(func, x1, x2, .., xn), ..]
        res.push_back(x_start);

        //// Construct starting simplex from starting point
        double x[_V * _D];
        for (int i=0; i < _D; i++) {
            for (int j=0; j < _D; j++) {
                if (j == i) {
                    x[_V*i + j+1] = x_start[j+1] + _step;
                } else {
                    x[_V*i + j+1] = x_start[j+1];
                };
            };
            x[_V*i] = get_lb_value(&x[_V*i+1]);
            res.push_back(x + _V*i);
        };

        //// Nelder-mead downhill simplex method
        double best;
        double x0[_V];  // centroid
        double xr[_V];  // new vertex using alpha
        double xe[_V];  // new vertex using gamma
        double xc[_V];  // new vertex using rho
        double rscore, escore, cscore;
        while (true) {
            sort(res.begin(), res.end(), NelderMead::ascending_first_value);
            best = res[0][0];

            // cout.precision(16);
            // cout << "  Current simplex: " << endl;
            // for (int i=0; i < _V; i++) {  // verts - 1
            //     cout << "   ";
            //     for (int j=0; j < _D; j++) {  // coords
            //         cout << res[i][j+1] << ", ";
            //     };
            //     cout << "  -->  " << res[i][0] << endl;
            // };

            // break after max_iters
            if (_max_iters != 0 and _iters >= _max_iters) {
                return to_point(res[0]);
            };
            _iters += 1;
            // cout << _iters << ". Best: " << res[0][1] << ", " << res[0][2] << " --> " << res[0][0] << endl;

            // break after no_improv_break iterations with no improvement
            if (best < prev_best - _no_improve_thr) {
                no_improv = 0;
                prev_best = best;
            } else {
                no_improv += 1;
            };

            if (no_improv >= _no_improv_break) {
                return to_point(res[0]);
            };

            // centroid
            for (int i=0; i < _D; i++) { x0[i+1] = 0; };
            for (int i=0; i < _D; i++) {  // verts - 1
                for (int j=0; j < _D; j++) {  // coords
                    x0[j+1] += res[i][j+1] / _D;
                };
            };

            //// reflection
            for (int i=1; i < _V; i++) {
                xr[i] = x0[i] + _alpha*(x0[i] - res[_V-1][i]);
            };
            rscore = get_lb_value(&xr[1]);
            xr[0] = rscore;
            if ((res[0][0] <= rscore) and (rscore < res[_V-2][0])) {
                for (int i=0; i < _V; i++) {
                    res[_V-1][i] = xr[i];
                };
                // cout << "  Reflection: " << res[_V-1][1] << ", " << res[_V-1][2] << " --> " << res[_V-1][0] << endl;
                continue;
            };

            //// expansion
            if (rscore < res[0][0]) {
                for (int i=1; i < _V; i++) {
                    xe[i] = x0[i] + _gamma*(x0[i] - res[_V-1][i]);
                };
                escore = get_lb_value(&xe[1]);
                xe[0] = escore;
                if (escore < rscore) {
                    for (int i=0; i < _V; i++) {
                        res[_V-1][i] = xe[i];
                    };
                    // cout << "  Expansion: " << res[_V-1][1] << ", " << res[_V-1][2] << " --> " << res[_V-1][0] << endl;
                    continue;
                } else {
                    for (int i=0; i < _V; i++) {
                        res[_V-1][i] = xr[i];
                    };
                    // cout << "  Expansion: " << res[_V-1][1] << ", " << res[_V-1][2] << " --> " << res[_V-1][0] << endl;
                    continue;
                };
            };


            //// contraction
            for (int i=1; i < _V; i++) {
                xc[i] = x0[i] + _rho*(x0[i] - res[_V-1][i]);
            };
            cscore = get_lb_value(&xc[1]);
            xc[0] = cscore;
            if (cscore < res[_V-1][0]) {
                for (int i=0; i < _V; i++) {
                    res[_V-1][i] = xc[i];
                };
                // cout << "  Contraction: " << res[_V-1][1] << ", " << res[_V-1][2] << " --> " << res[_V-1][0] << endl;
                continue;
            };

            //// reduction
            // Cache first vert in xr
            for (int i=1; i < _V; i++) {
                xr[i] = res[0][i];
            };
            // update each vertex, each coordinate
            for (int i=0; i < _V; i++) { // verts
                for (int j=1; j < _V; j++) { // coords
                    res[i][j] = xr[j] + _sigma*(res[i][j] - xr[j]);
                };
                res[i][0] = get_lb_value(&res[i][1]);
            };
            // cout << "  Reduction: " << endl;
            // for (int i=0; i < _V; i++) {  // verts - 1
            //     cout << "   ";
            //     for (int j=0; j < _D; j++) {  // coords
            //         cout << res[i][j+1] << ", ";
            //     };
            //     cout << "  -->  " << res[i][0] << endl;
            // };
        };
    };

    virtual ~NelderMead() {
    };
};
double NelderMead::_max_iteration = 0;

#endif
