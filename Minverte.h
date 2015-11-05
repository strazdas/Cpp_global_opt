#ifndef MINVERTE_H
#define MINVERTE_H 
/* Minverte algorithm purpose is to solve ASIMPL inner optimization problem, which
 * is to optimize simplex's lower bound by dividing it into subsimplexes. */
/* MIN VERTex lower bound minimum Estimate ~ MINVERTE */  
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
// Depencencies are:   main  <  Disimplv  <  Simplex  <  Minverte  <  Subsimplex


using namespace std;

class Minverte {   // For efficiency do not store links to parents
    Minverte(const Minverte& other) {};
    Minverte& operator=(const Minverte& other) {};
public:
    Minverte(vector<Point*> verts, double L) {
        _verts = verts;
        _L = L;
        _D = _verts.size() - 1;
        _partition.push_back(new Subsimplex(_verts, _L));
        _accuracy = 1e-1;    // should be proportional to problem dimension
        _max_iter = 10;      // This means max number of function evalutaions
    };
    vector<Point*> _verts;   // Simplex vertexes used for initialization
    double _L;
    int _D;
    vector<Subsimplex*> _partition;
    vector<Point*> _points;  // Cashe of Point objects to be deleted
    double _accuracy;        // Stopping criteria (tolerance threshold)
    int _max_iter;           // This means max number of function evalutaions
    
    double get_lb_value(vector<Point*> verts, double L, Point* point) {
        /* Returns lb value at a given point */
        double lb_value = -numeric_limits<double>::max();  // Note check numeric_limits<double>min() usage, because the value is 2e-308 ~ 0, and not negative.
        double dist;
        int min_i = numeric_limits<int>::min();
        // 2e-308 value is calculated here, is it correct?
        for (int i=0; i < verts.size(); i++) {
            dist = l2norm(point, verts[i]);
            // cout << lb_value << " is greater thant " << verts[i]->_values[0] - L*dist << endl;
            if (lb_value < verts[i]->_values[0] - L*dist) {
                lb_value = verts[i]->_values[0] - L*dist;
                min_i = i;
            };
        };
        return lb_value;
    };

    vector<Subsimplex*> divide_subsimplex(Subsimplex* subsimplex, string strategy="longest_half") {
        vector<Subsimplex*> divided_subsimplexes;

        if (strategy== "longest_half") {
            // Find middle point
            double c[_D];
            for (int i=0; i < _D; i++) {
                c[i] = (subsimplex->_le_v1->_X[i] + subsimplex->_le_v2->_X[i]) / 2.;
            };

            Point* middle_point = new Point(c, _D);

            _points.push_back(middle_point);
            middle_point->add_value(get_lb_value(_verts, _L, middle_point));

            // Construct two new simplexes using this middle point.
            Subsimplex* left_subsimplex = new Subsimplex(subsimplex->_L);
            Subsimplex* right_subsimplex = new Subsimplex(subsimplex->_L);

            for (int i=0; i < _D + 1; i++){
                // Point* point = _func->get(new Point(triangle[i], n)); 
                if (subsimplex->_verts[i] != subsimplex->_le_v1){
                    right_subsimplex->add_vertex(subsimplex->_verts[i]);
                } else {
                    right_subsimplex->add_vertex(middle_point);
                };
                if (subsimplex->_verts[i] != subsimplex->_le_v2) {
                    left_subsimplex->add_vertex(subsimplex->_verts[i]);
                } else {
                    left_subsimplex->add_vertex(middle_point);
                };
            };
            left_subsimplex->init_parameters();
            right_subsimplex->init_parameters();

            divided_subsimplexes.push_back(left_subsimplex);
            divided_subsimplexes.push_back(right_subsimplex);
        };
        return divided_subsimplexes;
    };

    // int get_simplex_to_divide_id(vector<Subsimplex*> subsimplexes) {
    //     // Selects subsimplex with lowest estimated lower bound minimum
    //     int simplex_with_min_lb_value_id = 0;
    //     for (int i=0; i < subsimplexes.size(); i++) {
    //         if (subsimplexes[simplex_with_min_lb_value_id]->_min_lb_value < subsimplexes[i]->_min_lb_value) {
    //             simplex_with_min_lb_value_id = i;
    //         };
    //     };
    //     return simplex_with_min_lb_value_id;
    // };

    vector<Subsimplex*> divide_subsimplexes(vector<Subsimplex*> simpls) {
        vector<Subsimplex*> divided_subsimplexes;
        for (int i=0; i < simpls.size(); i++) {
            vector<Subsimplex*> two_subsimplexes = divide_subsimplex(simpls[i]);
            for (int j=0; j < two_subsimplexes.size(); j++) {
                divided_subsimplexes.push_back(two_subsimplexes[j]);
            };
        };
        return divided_subsimplexes;
    };


    Point* minimize() {     // Returns estimated simplex_lb_min_value
        // Warrning: unimplemented (warning: old idea, minverte now is simply min vert)
        // Minverte idea is: to have several subsimplexes groups by number of subsimplex divisions
        // In each iteration divide one best simplex from each group (if value match, divide them bouth)
        // Best subsimplexes from each group can be found by itearting through all simplexes and comparing values.

        // double tolerance = _partition[0]->_tolerance;
        // int iter = 0;
        // // while (tolerance > _accuracy) {
        // while (_max_iter > iter) {
        //     // Pop subsimplex whith lowest min lb value from partition
        //     // Subsimplex* to_divide = _partition[0];
        //     // _partition.erase(_partition.begin() + 0);
        //
        //     // Divide selected subsimplexes
        //     // vector<Subsimplex*> divided_subsimplexes = divide_subsimplex(to_divide);
        //     vector<Subsimplex*> divided_subsimplexes = divide_subsimplexes(_partition);
        //
        //     // Remove all _partition subsimplexes.    Note: Does this work?
        //     for (int i=0; i < _partition.size(); i++) {
        //         delete _partition[i];
        //     };
        //     _partition.clear();
        //
        //     // Erase to_divide
        //     // delete to_divide;
        //
        //     // insert new simplexes as bubble sort
        //     // Iterate through
        //
        //     for (int i=0; i < divided_subsimplexes.size(); i++) { // Note: structure which has instert 
        //         _partition.push_back(divided_subsimplexes[i]);    // command would allow to insert as bubble sort  
        //     };
        //
        //     // sort partition by min lb value 
        //     sort(_partition.begin(), _partition.end(), Subsimplex::compare_by_min_lb_value);
        //
        //     // tolerance = _partition[0]->_tolerance;
        //     iter += 1;
        //     // cout << it << ". tol " << tolerance << " acc " << _accuracy << " estimate " << _partition[0]->_min_vert_value << " min_lb " << _partition[0]->_min_lb_value << " diameter " << _partition[0]->_diameter << endl;
        // };

        sort(_verts.begin(), _verts.end(), Point::compare_by_value);
        return _verts[0];
    };

    virtual ~Minverte() {
        // for (int i=0; i < _points.size(); i++) {
        //     delete _points[i];
        // };
        // _points.clear();
        // for (int i=0; i < _partition.size(); i++) {
        //     delete _partition[i];
        // };
        _partition.clear();
    };
};

#endif
