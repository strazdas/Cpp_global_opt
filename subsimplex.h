#ifndef SUBSIMPLEX_H
#define SUBSIMPLEX_H 
/* Subsimplex is a simplex used for solving inner optimization problem, which
 * is finding simplex lower bound minimum estimate */

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

using namespace std;


class SubsimplexTree;
class SubsimplexTreeNode;

class Subsimplex {  // Suitable for inner problems, e.g. Simplex lower bound min finding
    Subsimplex(const Subsimplex& other){}
    Subsimplex& operator=(const Subsimplex& other){}
public:
    Subsimplex(double L) {
        // _parent = 0;
        _D = 0;
        _L = L;
        _min_vert = 0;
        _diameter = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        _min_lb = 0;
        _min_lb_value = 0;
        _min_vert = 0;
        _max_vert = 0;
        _max_vert_value = -numeric_limits<double>::max();
        _min_vert_value = numeric_limits<double>::max();
        _should_be_divided = false;
    };

    Subsimplex(vector<Point*> verts, double L) {
        _L = L;
        _verts = verts;
        _D = 0;
        _min_vert = 0;
        _diameter = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        _min_lb = 0;
        _min_lb_value = 0;
        _min_vert = 0;
        _max_vert = 0;
        _max_vert_value = numeric_limits<double>::min();
        _min_vert_value = numeric_limits<double>::max();
        _should_be_divided = false;
        init_parameters();
    };

    int _D;                         // Variable space dimension

    double _L;                      // Predefined inner problem Lipschitz constant

    vector<Point*> _verts;          // Points with coordinates and values

    // Subsimplex* _parent;         // For efficiency do not store links to parents

    Point* _le_v1;                  // Longest edge vertexes
    Point* _le_v2; 
    double _diameter;               // Longest edge length

    Point* _min_vert;               // Pointer to vertex with lowest function value 
    double _min_vert_value;         // _min_vert function value 

    Point* _max_vert;
    double _max_vert_value;

    Point* _min_lb;                 // Subsimplexes estimate of accurate lower bound
    double _min_lb_value;
    double _tolerance;              // ||_min_lb - _min_vert||  or ||_min_lb_value - _min_vert_value||
    
    bool _should_be_divided;        // Should be divided in next iteration

    void init_parameters() {   // Called when all _verts have been added
        _D = _verts.size() - 1; 

        // Find _diameter and _le_v1, _le_v2
        double edge_length;         // Temporary variable
        for (int a=0; a < _verts.size(); a++) {
            for (int b=0; b < _verts.size(); b++){
                if (b > a) {
                    edge_length = l2norm(_verts[a], _verts[b]); 
                    if (edge_length > _diameter) {
                        _diameter = edge_length;
                        _le_v1 = _verts[a];
                        _le_v2 = _verts[b];
                    };
                };
            };
        }; 
 
        // Sort vertexes and set _min_vert, _max_vert
        sort(_verts.begin(), _verts.end(), Point::compare_by_value);
        _min_vert = _verts[0];
        _min_vert_value = _min_vert->_values[0];
        _max_vert = _verts[_verts.size()-1];
        _max_vert_value = _max_vert->_values[0];

        // Find longest edge lower bound
        _min_lb_value = find_simplex_lb_min_value_estimate(_le_v1, _le_v2, _L, _diameter);
        // cout << "_min_vert_value is " << _min_vert_value << endl;
        // if (_min_vert_value == 2.22507e-308) {
        //     cout << "_min_vert_value is: " << _min_vert_value << endl;
        //     exit(0);
        // };
        // Warrning (reasoning needed):  using accurate lb would work better, because function
        // Strictly changes by L But it wouldn't be efficient, because calculation would cost more?
        _tolerance = fabs(_min_vert_value - _min_lb_value); 
        // cout << "Calculated tolerance " << _tolerance << " _min_lb_value " << _min_lb_value << endl;
    };

    double find_simplex_lb_min_value_estimate(Point* _le_v1, Point* _le_v2, double _L, double _diameter) {
        /* This is A.Ž. simplex lower bound estimation algorithm Elleme, which do not
         * garantee tolerance strict decrease */
        /* Elleme: Efficient Lower Longest Edge bound Minimum Estimation */

        /* Error: was done here: all experiments using Elbme are not valid */ 
        // double x = (_le_v1->_values[0] + _le_v2->_values[0]) / 2 - _L * (l2norm(_le_v1, _le_v2)) / sqrt(2);
        
        // Use the bulk estimate
        double v1 = _le_v1->_values[0] - _L*_diameter;
        double v2 = _le_v2->_values[0] - _L*_diameter;
        if (v2 > v1){
            return v2;
        };
        return v1;
    };

    // double find_edge_lb_min_value(Point* _le_v1, Point* _le_v2, double _L) {
    // /* Finds accurate edge lower bound value */
    // };

    static bool wont_be_divided(Subsimplex* s) {
        return !s->_should_be_divided;
    };

    static double compare_by_min_lb_value(Subsimplex* s1, Subsimplex* s2) {  // sorts ascending
        return s1->_min_lb_value < s2->_min_lb_value; 
    };

    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
    };

    int size() {
        return _verts.size();
    };

    void print(){
        cout << " Subsimplex  (" << _diameter << ", " << _min_lb_value << "):" << endl;
        for (int i=0; i < _verts.size(); i++) {
            _verts[i]->print();
        };
    };

    static void print(vector<Subsimplex*> subsimplexes, string label="Printing subsimplexes:"){
        cout << label << endl;
        for (int i=0; i < subsimplexes.size(); i++){
            subsimplexes[i]->print();
        };
    };

    virtual ~Subsimplex(){    // Several subsimplexes can have same points - pakartotinis jų panaikinimas sukelia problemų, reikia struktūros, kurioje būtų saugomi unikalūs taškai, kuriuos būtų galima sunaikinti
        // for (int i=0; i < _verts.size(); i++) {
        //     delete _verts[i];
        // };
        _verts.clear();
    };  
};

#endif
